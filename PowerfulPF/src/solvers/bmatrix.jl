"""
    bmatrix.jl

Efficient B-matrix construction using incidence matrix representation.

The incidence matrix A (L×N) captures branch connectivity with +1 at from-bus
and -1 at to-bus for each branch. This enables unified matrix construction:

    B = Aᵀ D A   where D = diag(susceptances)

This file provides:
- `gather_branch_susceptances`: Collect susceptance values matching incidence rows
- `build_B_from_incidence`: Build N×N B-matrix from incidence and susceptances
"""

using SparseArrays: SparseMatrixCSC, sparse, nnz, nzrange, rowvals, nonzeros
using PowerfulCore: build_incidence_matrix, extract_from_to_idx, get_foreign_indices, in_service

"""
    gather_branch_susceptances(sys::SystemModel, branch_indices::Vector{Tuple{Symbol, Int}}) -> Vector{Float64}

Gather susceptance values for each branch in incidence order.

# Arguments
- `sys::SystemModel`: System containing Branch and Transformer2W models
- `branch_indices`: Mapping from incidence row l to (model_name, instance_idx)

# Returns
Vector of susceptances where:
- For lines: b[l] = 1/x
- For transformers: b[l] = 1/(x * tap) where tap = windv1/windv2

# Notes
- Uses `x` field from Branch model (series reactance in per-unit)
- For transformers, accounts for tap ratio in DC power flow approximation
- The tap adjustment matches MATPOWER's DC power flow convention
"""
function gather_branch_susceptances(sys::SystemModel, branch_indices::Vector{Tuple{Symbol, Int}})
    b = Vector{Float64}(undef, length(branch_indices))

    for (l, (model_name, instance_idx)) in enumerate(branch_indices)
        _, x = series_impedance(sys, model_name, instance_idx)
        if model_name == :Transformer2W
            # Tap ratio = windv1/windv2 (effective turns ratio)
            model = sys.models[model_name]
            b[l] = dc_susceptance(x, model.windv1[instance_idx] / model.windv2[instance_idx])
        else
            b[l] = dc_susceptance(x)
        end
    end

    return b
end

"""
    dc_bus_injections(sys::SystemModel) -> Vector{Float64}

Net active power injection (pu) at each internal bus index in the DC model: PV
`pg`, minus PQ `pl·f(vm) + ip·vm + yp·vm²` (`dc_load_power`) and the conductance `g·vm²` of every fixed shunt model
(`dc_shunt_power`, `fixed_shunt_models`; a switched shunt has none) at the case bus voltage `vm`, over in-service devices. The slack generator is not included; the reference bus
absorbs the mismatch. Same injections as the `DCPowerFlow` model kernels.
"""
function dc_bus_injections(sys::SystemModel)
    P = zeros(Float64, length(sys.models.Bus))
    vm = sys.models.Bus.vm

    if haskey(sys.models, :PV)
        pv = sys.models.PV
        active = in_service(pv)
        for (k, b) in enumerate(get_foreign_indices(sys, :PV, :i))
            active[k] && (P[b] += pv.pg[k])
        end
    end
    if haskey(sys.models, :PQ)
        pq = sys.models.PQ
        active = in_service(pq)
        for (k, b) in enumerate(get_foreign_indices(sys, :PQ, :i))
            active[k] && (P[b] -= dc_load_power(pq.pl[k], pq.ip[k], pq.yp[k], vm[b],
                                                sys.settings.constant_power_vmin))
        end
    end
    for name in fixed_shunt_models(sys)
        shunt = sys.models[name]
        params = shunt_admittance_params(shunt)
        active = in_service(shunt)
        for (k, b) in enumerate(get_foreign_indices(sys, name, :i))
            active[k] && (P[b] -= dc_shunt_power(shunt_admittance(params, k)[1], vm[b]))
        end
    end

    # Every other model in the system: PV, PQ and the shunts above contribute injections; Branch
    # and Transformer2W contribute the angle-dependent flow through the B-matrix instead (not here);
    # Slack is deliberately excluded (the reference bus absorbs the mismatch); the topological
    # models (AreaInterchange, Zone, Owner) have none. Anything else is not accounted for by either
    # path, and its power would be silently missing from the DC injections.
    known = union!(Set((:Bus, :PV, :PQ, :Slack, :Branch, :Transformer2W, :AreaInterchange, :Zone, :Owner)),
                   fixed_shunt_models(sys))
    unhandled = [name for name in keys(sys.models) if name ∉ known]
    isempty(unhandled) || throw(ArgumentError(
        "dc_bus_injections does not know how to include $(join(unhandled, ", ")) in the DC power flow bus injections"))

    return P
end

"""
    series_impedance(sys, model_name, k) -> (r, x)

Series resistance and reactance (per-unit) of instance `k` of a 2-terminal
branch model: `r`/`x` for Branch, the primary-secondary `r1_2`/`x1_2` for
Transformer2W.
"""
function series_impedance(sys::SystemModel, model_name::Symbol, k::Int)
    model = sys.models[model_name]
    model_name == :Branch && return (model.r[k], model.x[k])
    model_name == :Transformer2W && return (model.r1_2[k], model.x1_2[k])
    error("Unknown branch model type: $model_name")
end

"""
    build_B_from_incidence(A::SparseMatrixCSC{Int8}, b::Vector{T}) where T

Build bus susceptance matrix B = Aᵀ D A using COO format.

# Arguments
- `A::SparseMatrixCSC{Int8}`: Branch-bus incidence matrix (L×N)
- `b::Vector{T}`: Branch susceptances (length L), where b[l] = 1/x_l

# Returns
`SparseMatrixCSC{T}`: Bus susceptance matrix (N×N)

# Algorithm
Uses COO format to avoid intermediate allocations:
```
For each branch l from bus i to bus j:
    B[i,i] += b_l
    B[j,j] += b_l
    B[i,j] -= b_l
    B[j,i] -= b_l
```

This is equivalent to B = Aᵀ diag(b) A but without allocating the L×N intermediate.
"""
function build_B_from_incidence(
    A::SparseMatrixCSC{Int8},
    b::Vector{T}
) where T
    L, N = size(A)
    @assert length(b) == L "Susceptance vector length ($length(b)) must match incidence rows ($L)"

    # Branch endpoints in O(nnz) instead of O(L×N) per-branch extraction
    from_indices, to_indices = extract_from_to_idx(A)

    # Pre-allocate COO arrays: 4 entries per branch
    I_row = Vector{Int}(undef, 4 * L)
    J_col = Vector{Int}(undef, 4 * L)
    V_val = Vector{T}(undef, 4 * L)

    # Build COO triplets
    @inbounds for l in 1:L
        from_idx = from_indices[l]
        to_idx = to_indices[l]

        b_l = b[l]
        base = 4 * (l - 1)

        # Diagonal: B[from, from] += b
        I_row[base + 1] = from_idx
        J_col[base + 1] = from_idx
        V_val[base + 1] = b_l

        # Diagonal: B[to, to] += b
        I_row[base + 2] = to_idx
        J_col[base + 2] = to_idx
        V_val[base + 2] = b_l

        # Off-diagonal: B[from, to] -= b
        I_row[base + 3] = from_idx
        J_col[base + 3] = to_idx
        V_val[base + 3] = -b_l

        # Off-diagonal: B[to, from] -= b
        I_row[base + 4] = to_idx
        J_col[base + 4] = from_idx
        V_val[base + 4] = -b_l
    end

    # Create sparse matrix (duplicate entries are summed)
    return sparse(I_row, J_col, V_val, N, N)
end

"""
    dc_network(sys, slack_idx) -> NamedTuple

DC network model shared by DC power flow and PTDF/LODF:
- `A`, `branch_indices`: incidence matrix (L×N) and the (model, index) of each row
- `b`: branch susceptances in the same order
- `from_indices`, `to_indices`: endpoints of each branch
- `non_slack`: bus indices other than `slack_idx`
- `F`: KLU factorization of B = Aᵀ diag(b) A with the slack row/column removed
"""
function dc_network(sys::SystemModel, slack_idx::Int)
    A, branch_indices = build_incidence_matrix(sys)
    b = gather_branch_susceptances(sys, branch_indices)
    from_indices, to_indices = extract_from_to_idx(A)
    non_slack = collect(setdiff(1:size(A, 2), slack_idx))
    F = klu(build_B_from_incidence(A, b)[non_slack, non_slack])
    return (; A, branch_indices, b, from_indices, to_indices, non_slack, F)
end
