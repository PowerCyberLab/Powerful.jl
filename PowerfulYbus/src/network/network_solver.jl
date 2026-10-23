# ============================================================================
# Network Solver
# ============================================================================
#
# Solves the network equations at every bus with a Newton iteration in
# rectangular coordinates, in sparse form. Constant-impedance loads and the
# generator Norton admittances are in Y; the constant-power and constant-current
# parts of loads and current sources are nonlinear injections
# (network_injections.jl); static PV/Slack generators are constraints
# (static_generators.jl).
#
# ============================================================================

using ForwardDiff: ForwardDiff
using KLU
using LinearAlgebra: ldiv!, mul!
using SparseArrays

"""
    NetworkJacobian

Admittance matrix of one topology with the factorizations of its network
equations: the complex matrix of the linear equations, or the real Newton Jacobian
of the nonlinear ones (each built on first use).

# Fields
- `Y::SparseMatrixCSC{ComplexF64,Int}`: Admittance matrix of the topology
- `Y_fact`: KLU factorization of `Y` with the Slack rows replaced by identity rows:
  the matrix of the linear network equations
- `J::SparseMatrixCSC{Float64,Int}`: Real Jacobian work matrix (2n × 2n): the 2×2
  block `[G -B; B G]` of each entry G + jB of `Y` (plus a structural diagonal
  block for every bus), with bus k at rows/columns `2k-1` (real part) and `2k`
  (imaginary part); its values are rewritten when it is refactored. Empty until
  first needed
- `Y_nzval::Vector{Float64}`: The values of `J` due to `Y` alone
- `diag::Vector{NTuple{4,Int}}`: Indices into `nonzeros(J)` of the diagonal 2×2
  block of each bus, ordered (rr, ri, ir, ii)
- `pv_rows`, `slack_rows::Vector{Vector{NTuple{2,Int}}}`: For each static generator
  bus (`StaticGeneratorBuses` order), the indices into `nonzeros(J)` of its two rows,
  as (real row, imaginary row) pairs per stored column; the constraints replace them
- `fact`: KLU factorization of `J`, reused across Newton steps and network solves
  while it keeps converging (see [`solve_network!`](@ref))
- `refresh::Bool`: Whether the next Newton step must rebuild and refactor `J`
"""
mutable struct NetworkJacobian
    Y::SparseMatrixCSC{ComplexF64, Int}
    Y_fact::Union{Nothing, KLU.KLUFactorization{ComplexF64, Int}}
    J::SparseMatrixCSC{Float64, Int}
    Y_nzval::Vector{Float64}
    diag::Vector{NTuple{4, Int}}
    pv_rows::Vector{Vector{NTuple{2, Int}}}
    slack_rows::Vector{Vector{NTuple{2, Int}}}
    fact::Union{Nothing, KLU.KLUFactorization{Float64, Int}}
    refresh::Bool
end

NetworkJacobian(Y::SparseMatrixCSC{ComplexF64, Int}) =
    NetworkJacobian(Y, nothing, spzeros(Float64, 0, 0), Float64[], NTuple{4, Int}[],
                    Vector{NTuple{2, Int}}[], Vector{NTuple{2, Int}}[], nothing, true)

"""
    _linear_factorization!(nj::NetworkJacobian, static) -> KLUFactorization

Factorization of the linear network equations (built on first use): `Y` with the
row of each Slack bus replaced by the identity row (V = vs∠va there).
"""
function _linear_factorization!(nj::NetworkJacobian, static::StaticGeneratorBuses)
    nj.Y_fact === nothing || return nj.Y_fact
    A = copy(nj.Y)
    slack = falses(size(A, 1))
    slack[static.slack_bus] .= true
    rows, vals = rowvals(A), nonzeros(A)
    for c in axes(A, 2), k in nzrange(A, c)
        slack[rows[k]] && (vals[k] = rows[k] == c ? one(ComplexF64) : zero(ComplexF64))
    end
    for b in static.slack_bus
        iszero(A[b, b]) && (A[b, b] = one(ComplexF64))  # a Slack bus with no stored diagonal
    end
    return nj.Y_fact = klu(A)
end

"""
    _real_jacobian!(nj::NetworkJacobian, static) -> nj

Build the real Jacobian pattern of `nj` (on first use), directly in CSC form from
the columns of `Y`: column c of `Y` (with its diagonal) gives columns 2c-1 and 2c,
rows 2r-1 and 2r for each stored row r.
"""
function _real_jacobian!(nj::NetworkJacobian, static::StaticGeneratorBuses)
    size(nj.J, 1) == 0 || return nj
    Y = nj.Y
    n = size(Y, 1)
    rows, vals = rowvals(Y), nonzeros(Y)
    has_diag(c) = any(k -> rows[k] == c, nzrange(Y, c))
    colptr = Vector{Int}(undef, 2n + 1)
    colptr[1] = 1
    for c in 1:n
        m = 2 * (length(nzrange(Y, c)) + !has_diag(c))
        colptr[2c] = colptr[2c - 1] + m
        colptr[2c + 1] = colptr[2c] + m
    end
    nz = colptr[end] - 1
    rowval, nzval = Vector{Int}(undef, nz), Vector{Float64}(undef, nz)
    diag = Vector{NTuple{4, Int}}(undef, n)
    for c in 1:n
        pa, pb = colptr[2c - 1], colptr[2c]  # next slots in columns 2c-1 and 2c
        range = nzrange(Y, c)
        k, diag_done = first(range), false
        while k <= last(range) || !diag_done
            # Merge the diagonal into the sorted rows of column c
            if !diag_done && (k > last(range) || rows[k] >= c)
                r = c
                G, B = (k <= last(range) && rows[k] == c) ? reim(vals[k]) : (0.0, 0.0)
                (k <= last(range) && rows[k] == c) && (k += 1)
                diag_done = true
                diag[c] = (pa, pb, pa + 1, pb + 1)
            else
                r = rows[k]
                G, B = reim(vals[k])
                k += 1
            end
            rowval[pa], nzval[pa], rowval[pa + 1], nzval[pa + 1] = 2r - 1, G, 2r, B
            rowval[pb], nzval[pb], rowval[pb + 1], nzval[pb + 1] = 2r - 1, -B, 2r, G
            pa += 2; pb += 2
        end
    end
    nj.J = SparseMatrixCSC(2n, 2n, colptr, rowval, nzval)
    nj.Y_nzval = copy(nzval)
    nj.diag = diag
    # Row k of Y is column k of its transpose: the columns where bus k's rows are stored
    Yt = sparse(transpose(Y))
    nj.pv_rows = [_bus_rows(nj.J, Yt, k) for k in static.pv_bus]
    nj.slack_rows = [_bus_rows(nj.J, Yt, k) for k in static.slack_bus]
    return nj
end

# (real row, imaginary row) indices into nonzeros(J) of bus k, per stored column: the
# columns of bus k's row of Y (column k of Yt) and its diagonal block; the 2×2 block
# expansion stores both rows of a bus in the same columns, adjacently
function _bus_rows(J::SparseMatrixCSC, Yt::SparseMatrixCSC, k::Int)
    cols = sort!(unique!(push!(rowvals(Yt)[nzrange(Yt, k)], k)))
    pairs = NTuple{2, Int}[]
    for c in cols, jc in (2c - 1, 2c)
        ia = _nz_index(J, 2k - 1, jc)
        rowvals(J)[ia + 1] == 2k || error("Rows of bus $k are not paired in column $jc")
        push!(pairs, (ia, ia + 1))
    end
    return pairs
end

# Index into nonzeros(A) of stored entry (r, c) (row indices are sorted per column)
function _nz_index(A::SparseMatrixCSC, r::Int, c::Int)
    range = nzrange(A, c)
    k = searchsortedfirst(view(rowvals(A), range), r)
    (k <= length(range) && rowvals(A)[range[k]] == r) ||
        error("Entry ($r, $c) is not stored in the sparse matrix")
    return range[k]
end

"""
    NetworkCache{T, TC, L}

Cache for the network solution.

# Type Parameters
- `T`: Numeric type for buffers (Float64 or Dual)
- `TC`: TopologyCache type (a Float64 cache is shared with Dual buffers)
- `L`: Tuple type of the load models (`network_load_models`)

# Fields
- `topology_cache::TC`: Y matrix per topology
- `E_buf::InternalVoltageBuffer{T}`: Generator internal voltages and impedances
- `cs_buf::CurrentSourceBuffer{T}`: Current source injections (REGCA1, etc.)
- `loads::L`: Load models with nonlinear current, as `(model, instances, buses)`
- `static::StaticGeneratorBuses`: Static PV/Slack generator constraints
- `gen_buses::Vector{Int}`: Generator buses (order of the generator buffers)
- `Y_gen::Vector{ComplexF64}`: Generator admittances 1/Z_gen
- `jacobians::Dict{TopologyKey, NetworkJacobian}`: Admittance and Newton Jacobian
  per topology (built on first use; shared with Dual caches)
- `V_bus::Vector{Complex{T}}`: Bus voltages: the Newton iterate, warm start and
  solution
- `V_prev::Vector{Complex{T}}`: Iterate before the current line search
- `I_inj`, `mismatch::Vector{Complex{T}}`: Nonlinear injections and KCL mismatch F
- `residual::Vector{Complex{T}}`: Equation residual per bus: F, or at a static
  generator bus its two constraints packed as `complex(r1, r2)`
- `blocks::Vector{NTuple{4,T}}`: Injection derivative 2×2 block per bus
- `step::Vector{T}`: Real Newton right-hand side / step (length 2n)
- `tangent::Matrix{ComplexF64}`: Voltage derivatives of a Dual cache, one column per
  partial (n × number of partials; empty for Float64 and nested Duals)
- `I_gen`, `V_terminal`, `P_e`, `Q_e`: Generator outputs
"""
struct NetworkCache{T, TC<:TopologyCache, L<:Tuple}
    topology_cache::TC
    E_buf::InternalVoltageBuffer{T}
    cs_buf::CurrentSourceBuffer{T}
    loads::L
    static::StaticGeneratorBuses
    gen_buses::Vector{Int}
    Y_gen::Vector{ComplexF64}
    jacobians::Dict{TopologyKey, NetworkJacobian}
    V_bus::Vector{Complex{T}}
    V_prev::Vector{Complex{T}}
    I_inj::Vector{Complex{T}}
    mismatch::Vector{Complex{T}}
    residual::Vector{Complex{T}}
    blocks::Vector{NTuple{4, T}}
    step::Vector{T}
    tangent::Matrix{ComplexF64}
    I_gen::Vector{Complex{T}}
    V_terminal::Vector{Complex{T}}
    P_e::Vector{T}
    Q_e::Vector{T}
end

"""
    NetworkCache(sys; T=Float64, V_pf=nothing)

Create a NetworkCache from system model.

Loads keep their constant-impedance part in Y; their constant-power and
constant-current parts, and current sources, are nonlinear injections at their
buses. Static PV/Slack generators are constraints at their buses. `V_pf` (power
flow bus voltages) is the initial Newton iterate; without it the iteration starts
flat (1 pu).
"""
function NetworkCache(sys; T::Type=Float64, V_pf=nothing)
    topo_cache = TopologyCache(sys)
    base = topo_cache.base
    n_bus, gen_buses = base.n_bus, base.gen_buses
    n_gen = length(gen_buses)

    E_buf = InternalVoltageBuffer(T, n_gen)
    setup_internal_voltage_buffer!(E_buf, gen_buses, sys)
    populate_generator_impedances!(E_buf, sys.models, sys.external)
    Y_gen = ComplexF64[1 / z for z in E_buf.Z_gen]

    cs_buf = CurrentSourceBuffer(T, 0)
    setup_current_source_buffer!(cs_buf, sys)

    V = isnothing(V_pf) ? ones(Complex{T}, n_bus) : Vector{Complex{T}}(V_pf)
    length(V) == n_bus || error("V_pf has $(length(V)) entries for $n_bus buses")

    return NetworkCache(topo_cache, E_buf, cs_buf, network_load_models(sys),
                        static_generator_buses(sys), gen_buses, Y_gen,
                        Dict{TopologyKey, NetworkJacobian}(),
                        V, copy(V), zeros(Complex{T}, n_bus), zeros(Complex{T}, n_bus),
                        zeros(Complex{T}, n_bus), fill(ntuple(_ -> zero(T), 4), n_bus),
                        zeros(T, 2n_bus), _tangent_buffer(T, n_bus),
                        zeros(Complex{T}, n_gen), zeros(Complex{T}, n_gen),
                        zeros(T, n_gen), zeros(T, n_gen))
end

"""
    similar_network_cache(cache::NetworkCache, ::Type{T}) -> NetworkCache{T}

A cache with `T`-typed work buffers (e.g. ForwardDiff Duals) that shares the
topology cache, device mappings and per-topology matrices of `cache`.
"""
function similar_network_cache(cache::NetworkCache, ::Type{T}) where T
    E, cs = cache.E_buf, cache.cs_buf
    n_bus, n_gen = length(cache.V_bus), length(cache.I_gen)
    E_buf = InternalVoltageBuffer{T}(zeros(Complex{T}, n_gen), E.gen_bus_indices, E.gen_model_map,
                                     Complex{T}.(E.Z_gen), Set{Symbol}())
    cs_buf = CurrentSourceBuffer{T}(zeros(T, length(cs.Ip)), zeros(T, length(cs.Iq)),
                                    cs.src_bus, cs.src_model_map)
    return NetworkCache(cache.topology_cache, E_buf, cs_buf, cache.loads,
                        cache.static, cache.gen_buses, cache.Y_gen, cache.jacobians,
                        Complex{T}.(cache.V_bus), Complex{T}.(cache.V_bus),
                        zeros(Complex{T}, n_bus), zeros(Complex{T}, n_bus),
                        zeros(Complex{T}, n_bus), fill(ntuple(_ -> zero(T), 4), n_bus),
                        zeros(T, 2n_bus), _tangent_buffer(T, n_bus),
                        zeros(Complex{T}, n_gen), zeros(Complex{T}, n_gen),
                        zeros(T, n_gen), zeros(T, n_gen))
end

"""
    network_jacobian!(cache::NetworkCache) -> NetworkJacobian

The admittance and Newton Jacobian of the current topology, built on first use.
"""
function network_jacobian!(cache::NetworkCache)
    key = cache.topology_cache.current[]
    nj = get(cache.jacobians, key, nothing)
    nj === nothing || return nj
    return cache.jacobians[key] = NetworkJacobian(current_Y(cache.topology_cache).Y)
end

"""
    NetworkSolution{T}

Result of network solution.

# Fields
- `converged::Bool`: Whether the network solve converged to tolerance
- `iterations::Int`: Newton steps taken (1 for a direct solve of linear equations)
- `residual::T`: Final max equation residual (pu): the KCL current mismatch, or the
  constraint residual at static generator buses (0 for a direct solve, which is
  exact up to round-off)
- `I_gen::Vector{Complex{T}}`: Generator stator currents (E - V_terminal)/Z_gen
- `V_bus::Vector{Complex{T}}`: Bus voltages
- `V_terminal::Vector{Complex{T}}`: Generator terminal voltages
- `P_e::Vector{T}`: Generator electrical power
- `Q_e::Vector{T}`: Generator reactive power
- `detail::String`: Where not `converged` and the solve names what failed (time, iterations,
  residual, worst rows), else empty
"""
struct NetworkSolution{T}
    converged::Bool
    iterations::Int
    residual::T
    I_gen::Vector{Complex{T}}
    V_bus::Vector{Complex{T}}
    V_terminal::Vector{Complex{T}}
    P_e::Vector{T}
    Q_e::Vector{T}
    detail::String
end

"""
    solve_network!(cache::NetworkCache{Float64}, sys; maxiters=50, abstol=1e-8)

Solve the network equations at every bus by Newton iteration.

This solve serves the `EliminatedNetwork` treatment of the DAE residual (`YBusResidual`),
the baseline the simultaneous treatment is compared with, and the tests of the network
solver. The fixed-step integration (`solve_trapezoidal`) and the DAE integration with the
network simultaneous solve the network together with the models' algebraic rows
(`solve_algebraic!`).

# Equations
With the bus voltages V, KCL at every bus is

    F(V) = Y V - I_E - I_inj(V) = 0

where `Y` is the sparse admittance matrix (lines, transformers, shunts,
constant-impedance loads and generator admittances), `I_E`
the generator Norton currents E/Z_gen at the generator buses, and `I_inj` the
nonlinear injections: minus the constant-power and constant-current load currents
(`calc_load_current`) and plus the current-source currents (`source_current`).
They depend on |V| and are not holomorphic, so the Newton iteration runs in
rectangular coordinates with exact 2×2 derivative blocks. At a static generator bus
F is the current the generator injects, and its two equations are replaced by the
generator's constraints: Re(V F*) = Σ pg and |V|² = vs² (PV), or V = vs∠va
(Slack).

Without nonlinear injections or PV buses the equations are linear and holomorphic:
they are solved directly, `Y V = I_E` with the Slack rows as identity rows, by the
complex sparse LU factorization kept per topology (`iterations = 1`, `residual = 0`:
exact up to round-off, not re-evaluated).

# Iteration (nonlinear equations)
The iteration starts from the previous solution. Each step solves `J Δ = r` for the
equation residual `r` (sparse KLU) and takes the largest step `V - α Δ`
(α = 1, 1/2, ...) that decreases `‖r‖²` (Armijo backtracking). The factorization of
`J` is kept across steps and network solves (factoring dominates the cost): it is
rebuilt at the current iterate when a step with it contracts the residual by less
than a factor of 10, or finds no decrease; a fresh `J` that finds no decrease ends
the iteration. It stops when the max equation residual (current mismatch, or
constraint residual at static generator buses) is below `abstol` (pu), or after
`maxiters` steps; non-convergence is reported in the returned `NetworkSolution`
(`converged = false`), which callers must check.

After the solve, generator currents are I_gen = (E - V_terminal)/Z_gen.

ForwardDiff derivatives of the solution come from the Dual method
`solve_network!(cache::NetworkCache{<:Dual}, primal, sys)`.
"""
function solve_network!(cache::NetworkCache{Float64}, sys;
                        maxiters::Int=50, abstol::Real=1e-8)
    rn = network_jacobian!(cache)
    vmin = sys.settings.constant_power_vmin

    if _is_nonlinear(cache)
        converged, iterations, residual = _newton!(cache, rn, cache.E_buf.E_internal, vmin,
                                                   maxiters, abstol)
    else
        _solve_linear!(cache, rn)
        converged, iterations, residual = true, 1, 0.0
    end

    @debug "Network Newton" iterations residual converged
    return _network_solution!(cache, converged, iterations, residual)
end

"""
    solve_network!(cache::NetworkCache{D}, primal::NetworkCache{Float64}, sys;
                   maxiters=50, abstol=1e-8) where {D<:ForwardDiff.Dual}

Network solution with ForwardDiff derivatives by the implicit function theorem.

`primal` (the Float64 cache sharing `cache`'s topology cache) solves the network at the
values of the Dual inputs: the generator internal voltages `cache.E_buf.E_internal` and
the current-source currents `cache.cs_buf`. The derivatives are then one linear solve
with the Jacobian `J` of the equation residual `r` at that solution:

    dV = -J⁻¹ (∂r/∂E dE + ∂r/∂I_src dI_src)

obtained as the derivative part of one Newton step from the solution with zero
voltage derivatives: the Dual residual there has `∂r/∂(E, I_src)·d(E, I_src)` as
partials. The voltage values are those of the `primal` solution. For linear network
equations `J` is the complex matrix factored per topology (no new factorization);
otherwise the real Jacobian at the solution is factored. The voltages of `cache` are
set from `primal` on every call, so no derivative carries over from an earlier
evaluation.
"""
function solve_network!(cache::NetworkCache{D}, primal::NetworkCache{Float64}, sys;
                        maxiters::Int=50, abstol::Real=1e-8) where {D<:ForwardDiff.Dual}
    cache.topology_cache === primal.topology_cache ||
        error("The Dual network cache and its Float64 primal cache have different topology caches")
    E, E0 = cache.E_buf.E_internal, primal.E_buf.E_internal
    @inbounds for i in eachindex(E)
        E0[i] = _value(E[i])
    end
    cs, cs0 = cache.cs_buf, primal.cs_buf
    @inbounds for s in eachindex(cs.Ip)
        cs0.Ip[s] = _value(cs.Ip[s])
        cs0.Iq[s] = _value(cs.Iq[s])
    end
    sol = solve_network!(primal, sys; maxiters, abstol)

    V = cache.V_bus
    copyto!(V, primal.V_bus)  # zero voltage derivatives
    if sol.converged
        rn = network_jacobian!(cache)
        vmin = sys.settings.constant_power_vmin
        if _is_nonlinear(cache)
            _real_jacobian!(rn, cache.static)
            _evaluate_residual!(cache, rn, E, vmin)
            _newton_direction!(cache, rn, vmin)
            @inbounds for i in eachindex(V)
                V[i] -= _derivative_part(complex(cache.step[2i - 1], cache.step[2i]))
            end
        else
            # Linear equations: J is the complex matrix factored per topology
            _evaluate_residual!(cache, rn, E, vmin)
            _subtract_tangent!(V, _linear_factorization!(rn, cache.static), cache.residual,
                               cache.tangent)
        end
    end
    return _network_solution!(cache, sol.converged, sol.iterations, D(sol.residual))
end

# The derivative part of a (possibly nested) Dual: its Float64 value set to zero
_derivative_part(z) = z - _value(z)

# Tangent buffer of a cache with number type T: one column per partial of a single-level
# Dual (solved in place), none otherwise
_tangent_buffer(::Type, n::Int) = zeros(ComplexF64, 0, 0)
_tangent_buffer(::Type{ForwardDiff.Dual{Tag, Float64, N}}, n::Int) where {Tag, N} =
    zeros(ComplexF64, n, N)

# V -= J⁻¹ ∂r (derivative parts only) with the complex factorization F of J: every
# partial of the complex Dual residual r solved in place in the columns of W
function _subtract_tangent!(V::AbstractVector{Complex{ForwardDiff.Dual{Tag, Float64, N}}}, F,
                            r::AbstractVector, W::Matrix{ComplexF64}) where {Tag, N}
    @inbounds for i in 1:N, k in eachindex(r)
        W[k, i] = complex(ForwardDiff.partials(real(r[k]), i), ForwardDiff.partials(imag(r[k]), i))
    end
    ldiv!(F, W)
    @inbounds for k in eachindex(V)
        V[k] -= complex(ForwardDiff.Dual{Tag}(0.0, ntuple(i -> real(W[k, i]), Val(N))...),
                        ForwardDiff.Dual{Tag}(0.0, ntuple(i -> imag(W[k, i]), Val(N))...))
    end
    return V
end
# Nested Duals: solve value and partials recursively
_subtract_tangent!(V, F, r, W) = (V .-= _derivative_part.(_solve_components(F, r)))

# Generator outputs from the bus voltages of `cache`: terminal voltages, stator currents
# I_gen = (E - V_terminal)/Z_gen and electrical powers
function _network_solution!(cache::NetworkCache{T}, converged::Bool, iterations::Int,
                            residual, detail::String = "") where T
    E_internal, Z_gen, V = cache.E_buf.E_internal, cache.E_buf.Z_gen, cache.V_bus
    @inbounds for (i, b) in enumerate(cache.gen_buses)
        cache.V_terminal[i] = V[b]
        cache.I_gen[i] = (E_internal[i] - V[b]) / Z_gen[i]
    end
    compute_electrical_power!(cache.P_e, cache.Q_e, E_internal, cache.V_terminal, cache.E_buf)

    # References to cache buffers (no copy needed - consumed immediately by caller)
    return NetworkSolution{T}(converged, iterations, residual, cache.I_gen, cache.V_bus,
                              cache.V_terminal, cache.P_e, cache.Q_e, detail)
end

"""
Direct solve of the linear network equations `Y V = I_E` (Slack rows: V = vs∠va)
with the complex KLU factorization of the topology, into `cache.V_bus`.
"""
function _solve_linear!(cache::NetworkCache{Float64}, nj::NetworkJacobian)
    F = _linear_factorization!(nj, cache.static)
    V, st, E = cache.V_bus, cache.static, cache.E_buf.E_internal
    fill!(V, zero(ComplexF64))
    @inbounds for (i, b) in enumerate(cache.gen_buses)
        V[b] += E[i] * cache.Y_gen[i]
    end
    @inbounds for (k, b) in enumerate(st.slack_bus)
        V[b] = st.slack_V[k]
    end
    ldiv!(F, V)
    return V
end

"""
Newton iteration on the network equations (see [`solve_network!`](@ref)); returns
`(converged, iterations, residual)`. It serves the `EliminatedNetwork` baseline and the tests of
the network solver only.
"""
function _newton!(cache::NetworkCache{Float64}, rn::NetworkJacobian, E_internal, vmin,
                  maxiters::Int, abstol::Real)
    _real_jacobian!(rn, cache.static)
    residual, merit = _evaluate_residual!(cache, rn, E_internal, vmin)
    iterations = 0
    while !(residual < abstol) && iterations < maxiters && isfinite(residual)
        fresh = _newton_direction!(cache, rn, vmin)
        accepted, residual_new, merit_new = _line_search!(cache, rn, E_internal, vmin, merit)
        iterations += 1
        if !accepted
            fresh && break  # no decrease even along the Newton direction of a fresh Jacobian
            rn.refresh = true
            residual, merit = _evaluate_residual!(cache, rn, E_internal, vmin)
            continue
        end
        # Slow contraction with a reused Jacobian: rebuild it for the next step
        !fresh && residual_new > residual / 10 && (rn.refresh = true)
        residual, merit = residual_new, merit_new
    end
    return residual < abstol, iterations, residual
end

"""
Evaluate the KCL mismatch and the equation residual at the current iterate; returns
the max absolute residual component and the merit `‖r‖²`.
"""
function _evaluate_residual!(cache::NetworkCache, rn::NetworkJacobian, E_internal, vmin)
    _network_mismatch!(cache, rn, E_internal, vmin)
    r = _equation_residual!(cache)
    merit = sum(z -> abs2(real(z)) + abs2(imag(z)), r; init=zero(real(eltype(r))))
    return _max_abs_component(r), merit
end

"""
Backtracking line search along the Newton step `cache.step` from `V`: tries
`V - α Δ` for α = 1, 1/2, ..., 1/512 and accepts the first with
`‖r‖² ≤ (1 - 2·10⁻⁴ α) merit0` (Armijo). Returns `(accepted, max residual, merit)`;
on rejection `V` is restored.
"""
function _line_search!(cache::NetworkCache, rn::NetworkJacobian, E_internal, vmin, merit0)
    V, V_prev, step = cache.V_bus, cache.V_prev, cache.step
    copyto!(V_prev, V)
    α = 1.0
    residual, merit = zero(merit0), merit0
    for _ in 1:10
        @inbounds for i in eachindex(V)
            V[i] = V_prev[i] - α * complex(step[2i - 1], step[2i])
        end
        residual, merit = _evaluate_residual!(cache, rn, E_internal, vmin)
        merit <= (1 - 2e-4 * α) * merit0 && return true, residual, merit
        α /= 2
    end
    copyto!(V, V_prev)
    return false, residual, merit
end

"""
KCL mismatch `Y V - I_E - I_inj(V)` at every bus, into `cache.mismatch`.
"""
function _network_mismatch!(cache::NetworkCache, rn::NetworkJacobian, E_internal, vmin)
    V, I_inj, mismatch = cache.V_bus, cache.I_inj, cache.mismatch
    mul!(mismatch, rn.Y, V)
    @inbounds for (i, b) in enumerate(cache.gen_buses)
        mismatch[b] -= E_internal[i] * cache.Y_gen[i]
    end
    fill!(I_inj, zero(eltype(I_inj)))
    _add_load_currents!(I_inj, cache.loads, V, vmin)
    cs = cache.cs_buf
    @inbounds for (s, pos) in enumerate(cs.src_bus)
        pos == 0 && continue
        v = V[pos]
        i = source_current(real(v), imag(v), cs.Ip[s], cs.Iq[s])
        I_inj[pos] += complex(i.i_r, i.i_i)
    end
    mismatch .-= I_inj
    return mismatch
end

"""
Equation residual per bus into `cache.residual`: the KCL mismatch F, or at
a static generator bus its constraints packed as `complex(r1, r2)`:
PV `r1 = Re(V F*) - Σ pg`, `r2 = |V|² - vs²`; Slack `r1 + j r2 = V - vs∠va`.
"""
function _equation_residual!(cache::NetworkCache)
    V, F, r, st = cache.V_bus, cache.mismatch, cache.residual, cache.static
    copyto!(r, F)
    @inbounds for (k, pos) in enumerate(st.pv_bus)
        v, f = V[pos], F[pos]
        r[pos] = complex(real(v) * real(f) + imag(v) * imag(f) - st.pv_p[k], abs2(v) - st.pv_vs[k]^2)
    end
    @inbounds for (k, pos) in enumerate(st.slack_bus)
        r[pos] = V[pos] - st.slack_V[k]
    end
    return r
end

_max_abs_component(x::AbstractVector{<:Complex}) =
    mapreduce(z -> max(abs(real(z)), abs(imag(z))), max, x; init=zero(real(eltype(x))))

"""
Derivative 2×2 blocks of the mismatch due to the nonlinear injections, per bus, into `cache.blocks` (loads add their current derivative, sources
subtract theirs).
"""
function _injection_blocks!(cache::NetworkCache{T}, vmin) where T
    V, blocks = cache.V_bus, cache.blocks
    fill!(blocks, ntuple(_ -> zero(T), 4))
    _add_load_jacobians!(blocks, cache.loads, V, vmin)
    cs = cache.cs_buf
    @inbounds for (s, pos) in enumerate(cs.src_bus)
        pos == 0 && continue
        v = V[pos]
        d = source_current_jacobian(real(v), imag(v), cs.Ip[s], cs.Iq[s])
        b = blocks[pos]
        blocks[pos] = (b[1] - d.rr, b[2] - d.ri, b[3] - d.ir, b[4] - d.ii)
    end
    nothing
end

"""
Real right-hand side of the Newton step, the equation residual interleaved
real/imaginary, into `cache.step`.
"""
function _residual_rhs!(cache::NetworkCache)
    @inbounds for (k, r) in enumerate(cache.residual)
        cache.step[2k - 1] = real(r)
        cache.step[2k] = imag(r)
    end
    nothing
end

"""
Whether the network equations are nonlinear (nonlinear injections or PV buses);
otherwise the Jacobian is constant per topology.
"""
_is_nonlinear(cache::NetworkCache) =
    !isempty(cache.loads) || any(>(0), cache.cs_buf.src_bus) || !isempty(cache.static.pv_bus)

# Float64 value of a (possibly nested) ForwardDiff Dual
_value(x::Float64) = x
_value(x::ForwardDiff.Dual) = _value(ForwardDiff.value(x))
_value(z::Complex) = complex(_value(real(z)), _value(imag(z)))

"""
Assemble the values `nz` of the real Jacobian at the current iterate: `Y`, the
injection derivative blocks and the static generator constraint rows (Float64
values of a Dual iterate).
"""
function _assemble_jacobian!(nz::Vector{Float64}, rn::NetworkJacobian, cache::NetworkCache, vmin)
    _injection_blocks!(cache, vmin)
    copyto!(nz, rn.Y_nzval)
    @inbounds for (k, idx) in enumerate(rn.diag)
        b = cache.blocks[k]
        for m in 1:4
            nz[idx[m]] += _value(b[m])
        end
    end
    _constraint_rows!(nz, rn, cache)
    return nz
end

"""
Newton direction `J \\ r` into `cache.step` with the sparse KLU factorization of the
current topology. `J` is rebuilt at the current iterate and refactored (with
pivoting) when the topology has none yet or `rn.refresh` is set; otherwise the kept
factorization is reused. Returns whether `J` is exact at the current iterate (just
built, or the equations are linear).
"""
function _newton_direction!(cache::NetworkCache{Float64}, rn::NetworkJacobian, vmin)
    fresh = !_is_nonlinear(cache)
    if rn.fact === nothing || rn.refresh
        _assemble_jacobian!(nonzeros(rn.J), rn, cache, vmin)
        rn.fact = klu(rn.J)
        rn.refresh = false
        @debug "Network Jacobian refactored" n=size(rn.J, 1) nnz=nnz(rn.J)
        fresh = true
    end
    _residual_rhs!(cache)
    ldiv!(rn.fact, cache.step)
    return fresh
end

"""
Newton direction for a Dual-typed (ForwardDiff) cache: the Float64 Jacobian at the
current voltages is factored with sparse KLU (a copy, so the Float64 cache's kept
factorization is untouched) and applied to the value and to each partial of the
Dual right-hand side. At a network solution with zero voltage derivatives this is
the implicit-function derivative step of `solve_network!(::NetworkCache{<:Dual}, ...)`.
"""
function _newton_direction!(cache::NetworkCache{T}, rn::NetworkJacobian, vmin) where T
    J = copy(rn.J)
    _assemble_jacobian!(nonzeros(J), rn, cache, vmin)
    _residual_rhs!(cache)
    cache.step .= _solve_components(klu(J), cache.step)
    return true
end

# J \ x for a Float64 (real or complex) factorization and a (possibly nested) Dual
# vector x: the solve is linear, so it applies to the value and to each partial
# separately
_solve_components(F, x::Vector{Float64}) = F \ x
_solve_components(F, x::Vector{ComplexF64}) = F \ x
function _solve_components(F, x::AbstractVector{ForwardDiff.Dual{Tag, V, N}}) where {Tag, V, N}
    v = _solve_components(F, ForwardDiff.value.(x))
    p = ntuple(i -> _solve_components(F, ForwardDiff.partials.(x, i)), Val(N))
    return [ForwardDiff.Dual{Tag}(v[k], ntuple(i -> p[i][k], Val(N))...) for k in eachindex(x)]
end
function _solve_components(F, x::AbstractVector{Complex{ForwardDiff.Dual{Tag, V, N}}}) where {Tag, V, N}
    v = _solve_components(F, [complex(ForwardDiff.value(real(z)), ForwardDiff.value(imag(z))) for z in x])
    p = ntuple(i -> _solve_components(F, [complex(ForwardDiff.partials(real(z), i),
                                                   ForwardDiff.partials(imag(z), i)) for z in x]), Val(N))
    return [complex(ForwardDiff.Dual{Tag}(real(v[k]), ntuple(i -> real(p[i][k]), Val(N))...),
                    ForwardDiff.Dual{Tag}(imag(v[k]), ntuple(i -> imag(p[i][k]), Val(N))...))
            for k in eachindex(x)]
end

"""
Replace the Jacobian rows (values `nz` of the sparse Jacobian) of the static
generator buses by the derivatives of their constraints. With rows a (real) and
b (imaginary) of bus k holding ∂F_r/∂x and ∂F_i/∂x:
PV `∂r1/∂x = V_r ∂F_r/∂x + V_i ∂F_i/∂x + F_r e_a + F_i e_b`, `∂r2/∂x = 2V_r e_a + 2V_i e_b`;
Slack rows become e_a and e_b.
"""
function _constraint_rows!(nz::Vector{Float64}, rn::NetworkJacobian, cache::NetworkCache)
    st = cache.static
    @inbounds for (k, pos) in enumerate(st.pv_bus)
        v, f = _value(cache.V_bus[pos]), _value(cache.mismatch[pos])
        vr, vi = real(v), imag(v)
        for (ia, ib) in rn.pv_rows[k]
            nz[ia] = vr * nz[ia] + vi * nz[ib]
            nz[ib] = 0.0
        end
        d = rn.diag[pos]
        nz[d[1]] += real(f); nz[d[2]] += imag(f)
        nz[d[3]] = 2vr; nz[d[4]] = 2vi
    end
    @inbounds for (k, pos) in enumerate(st.slack_bus)
        for (ia, ib) in rn.slack_rows[k]
            nz[ia] = 0.0; nz[ib] = 0.0
        end
        d = rn.diag[pos]
        nz[d[1]] = 1.0; nz[d[4]] = 1.0
    end
    nothing
end

"""
    network_failure_message(net_sol::NetworkSolution, t) -> String

Describe a non-converged network solve at time `t` (Newton steps used and the final
equation residual; the `detail` of the solution where it has one, naming the worst rows),
for warnings and errors.
"""
network_failure_message(net_sol::NetworkSolution, t) =
    isempty(net_sol.detail) ?
    "Y-bus network solve did not converge at t = $t: $(net_sol.iterations) Newton " *
    "steps, final max equation residual $(net_sol.residual) pu" :
    net_sol.detail

"""
    compute_electrical_power!(P_e, Q_e, E_internal, V_terminal, E_buf)

Compute electrical power using stator current from terminal voltages.

# Sign Convention (ANDES-compatible)
P_e is computed from the internal EMF perspective:
- S = E_internal * conj(I_stator)
- P_e = real(S)

Sign interpretation:
- P_e > 0: Generator producing power (E drives current into network)
- P_e < 0: Generator absorbing power (current flows back into E)

This matches the ANDES convention where:
- Generator injects power: appears as -P_e in bus power balance
- For swing equation: P_m - P_e - D*(ω-1) = 2H * dω/dt

# Physics
The stator current is:
- I_stator = (E_internal - V_terminal) / Z_gen  (positive = from E to V)
"""
function compute_electrical_power!(P_e::AbstractVector{T}, Q_e::AbstractVector{T},
                                   E_internal::AbstractVector{Complex{T}},
                                   V_terminal::AbstractVector{Complex{T}},
                                   E_buf::InternalVoltageBuffer{T}) where T
    @inbounds for i in eachindex(P_e)
        # Compute stator current: I_stator = (E - V_terminal) / Z_gen
        Z_gen = E_buf.Z_gen[i]
        I_stator = (E_internal[i] - V_terminal[i]) / Z_gen

        # Compute electrical power from internal EMF: S = E * conj(I)
        # P_e > 0 means power flows from E to network (generation)
        # P_e < 0 means power flows from network to E (absorption/motoring)
        S = E_internal[i] * conj(I_stator)
        P_e[i] = real(S)
        Q_e[i] = imag(S)

        @debug "compute_electrical_power! Gen $i" E=E_internal[i] V=V_terminal[i] Z=Z_gen I=I_stator S P_e=P_e[i] Q_e=Q_e[i]
    end
end

export NetworkCache, NetworkSolution, NetworkJacobian
export solve_network!, similar_network_cache, network_jacobian!, network_failure_message
export compute_electrical_power!
