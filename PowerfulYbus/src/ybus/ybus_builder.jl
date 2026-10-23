# ============================================================================
# Y-Bus Matrix Builder
# ============================================================================
#
# Builds the full n×n admittance matrix from system models.
# ============================================================================

using SparseArrays
using LinearAlgebra
using PowerfulPF: fixed_shunt_models, shunt_admittance_params, shunt_admittance

"""
    FullYMatrix{T}

Admittance matrix of the network and its generator buses.

# Fields
- `Y::SparseMatrixCSC{Complex{T}, Int}`: Full n×n admittance matrix
- `n_bus::Int`: Number of buses
- `gen_buses::Vector{Int}`: Internal indices of the buses of the in-service machines
  (`online_machines`), sorted
"""
struct FullYMatrix{T}
    Y::SparseMatrixCSC{Complex{T}, Int}
    n_bus::Int
    gen_buses::Vector{Int}
end

"""
    device_bus_indices(sys, model_name::Symbol) -> Vector{Int}

Internal Bus index of each instance of `model_name`, through its first terminal
(from `foreign_refs`).
"""
device_bus_indices(sys, model_name::Symbol) =
    terminal_vertex_indices(sys.order, model_name, first(get_terminals(eltype(sys.models[model_name]))))

"""
    online_machines(sys) -> Vector{Tuple{Symbol, Int, Int}}

`(model_name, instance_index, internal_bus_index)` of every in-service (per
`in_service`) `MachineWithInternalVoltage` instance. The single source of generators
for the Y matrix, the bus partition and the internal voltage buffer.
"""
function online_machines(sys)
    machines = Tuple{Symbol, Int, Int}[]
    for (model_name, model) in pairs(sys.models)
        network_partition_role(eltype(model)) isa MachineWithInternalVoltage || continue
        isempty(model) && continue
        active = in_service(model)
        for (idx, bus) in enumerate(device_bus_indices(sys, model_name))
            active[idx] && push!(machines, (model_name, idx, bus))
        end
    end
    return machines
end

"""
    generators_by_bus(sys) -> Dict{Int, Tuple{Symbol, Int}}

Internal bus index → `(model_name, instance_index)` of the in-service generator
(`online_machines`) at that bus. Errors if a bus has two.
"""
function generators_by_bus(sys)
    gens = Dict{Int, Tuple{Symbol, Int}}()
    for (model_name, idx, bus) in online_machines(sys)
        haskey(gens, bus) &&
            error("Two generators at bus index $bus: $(gens[bus]) and ($model_name, $idx)")
        gens[bus] = (model_name, idx)
    end
    return gens
end

"""
    build_full_ybus(sys::SystemModel{T}; topology_key=nothing) where T

Build the full n×n admittance matrix from system models.

# Arguments
- `sys`: SystemModel containing all network elements
- `topology_key`: Optional key specifying which topology (for fault/line-out scenarios)

# Returns
`FullYMatrix{T}` containing the admittance matrix and the generator buses.

# Algorithm
1. Initialize sparse Y matrix (n_bus × n_bus)
2. Add Branch contributions (series + shunt)
3. Add Transformer contributions (with tap ratio)
4. Add the fixed shunt contributions (every `AbstractShunt` model: `fixed_shunt_models`)
5. Add generator admittances and constant-impedance loads

Refuses a model whose Y-bus role is `FixedNetworkElement` that it does not stamp
(`check_stamped_network_elements`), naming it.

# Extended help
Each series branch of admittance `y = 1/(r + jx)` with line-charging susceptance
split as shunts `jb/2` at each end contributes the standard two-port stamp:
`Y[from,from] += y + jb/2`, `Y[to,to] += y + jb/2`, `Y[from,to] -= y`, `Y[to,from] -= y`
(a transformer's off-nominal tap scales the off-diagonal terms and its own-bus term
asymmetrically). Generator and constant-impedance load admittances are added to their
own bus's diagonal only.
"""
function build_full_ybus(sys; topology_key=nothing)
    T = Float64  # Default precision
    check_stamped_network_elements(sys)

    n_bus = length(sys.models[:Bus])

    # Initialize sparse Y matrix using COO format for efficient construction
    I_idx = Int[]
    J_idx = Int[]
    V_val = Complex{T}[]

    # Process branches (transmission lines)
    if haskey(sys.models, :Branch)
        _add_branch_contributions!(I_idx, J_idx, V_val, sys)
    end

    # Process transformers (two-winding, if present)
    if haskey(sys.models, :Transformer2W)
        _add_transformer_contributions!(I_idx, J_idx, V_val, sys)
    end

    # Process the fixed shunts (Shunt, SwitchedShunt)
    for name in fixed_shunt_models(sys)
        _add_shunt_contributions!(I_idx, J_idx, V_val, sys, name)
    end

    # Add generator admittances (their Norton sources are network solve injections)
    _add_generator_contributions!(I_idx, J_idx, V_val, sys)

    # Add PQ Z component (constant admittance yp, yq); the other load parts are
    # network solve injections
    @debug "build_full_ybus COO size before PQ Z" coo_size=length(I_idx)
    _add_pq_z_contributions!(I_idx, J_idx, V_val, sys)
    @debug "build_full_ybus COO size after PQ Z" coo_size=length(I_idx)

    # Build sparse matrix from COO
    Y = sparse(I_idx, J_idx, V_val, n_bus, n_bus)

    @debug "build_full_ybus Y[1,1]" Y_11=Y[1,1]

    return FullYMatrix{T}(Y, n_bus, generator_buses(sys))
end

"""
Add branch (transmission line) contributions to Y matrix.
"""
function _add_branch_contributions!(I_idx, J_idx, V_val, sys)
    jay = im
    branches = sys.models[:Branch]
    isempty(branches) && return
    from_bus = get_foreign_indices(sys, :Branch, :i)
    to_bus = get_foreign_indices(sys, :Branch, :j)

    for idx in eachindex(branches)
        from_int = from_bus[idx]
        to_int = to_bus[idx]

        # Check if line is in service
        if !branches.status[idx]
            continue
        end

        # Get line parameters (all branches must have r, x, b)
        r = branches.r[idx]
        x = branches.x[idx]
        b = branches.b[idx]

        # Series admittance
        z = complex(r, x)
        y_series = 1.0 / z

        # Shunt admittance (half at each end)
        y_shunt = jay * b / 2

        # Diagonal contributions
        # Y[from, from] += y_series + y_shunt
        push!(I_idx, from_int); push!(J_idx, from_int); push!(V_val, y_series + y_shunt)
        # Y[to, to] += y_series + y_shunt
        push!(I_idx, to_int); push!(J_idx, to_int); push!(V_val, y_series + y_shunt)

        # Off-diagonal contributions
        # Y[from, to] -= y_series
        push!(I_idx, from_int); push!(J_idx, to_int); push!(V_val, -y_series)
        # Y[to, from] -= y_series
        push!(I_idx, to_int); push!(J_idx, from_int); push!(V_val, -y_series)
    end
end

"""
Add transformer contributions to Y matrix (with tap ratio and phase shift).

Uses Transformer2W fields:
- `r1_2`, `x1_2`: Resistance and reactance between windings
- `windv1`, `windv2`: Winding voltage ratios (pu of the bus base voltages); tap t = windv1/windv2
- `ang1`: Phase shift angle in degrees
- `mag1`, `mag2`: Magnetizing conductance and susceptance (pu on system base), a shunt at
  the winding 1 bus on the bus side of the tap
- `stat`: Status (true = in service)

Same admittances as the power flow kernel (`calc_transformer_currents`):
Y_ii = y/t² + y_mag, Y_ij = -y/t·e^{jφ}, Y_ji = -y/t·e^{-jφ}, Y_jj = y.

"""
function _add_transformer_contributions!(I_idx, J_idx, V_val, sys)
    jay = im
    transformers = sys.models[:Transformer2W]
    isempty(transformers) && return
    from_bus = get_foreign_indices(sys, :Transformer2W, :i)
    to_bus = get_foreign_indices(sys, :Transformer2W, :j)

    for idx in eachindex(transformers)
        from_int = from_bus[idx]
        to_int = to_bus[idx]

        # Check if transformer is in service (stat is BitVector)
        if !transformers.status[idx]
            continue
        end

        # Get transformer parameters (Transformer2W field names)
        r = transformers.r1_2[idx]
        x = transformers.x1_2[idx]

        # Off-nominal tap ratio t = windv1/windv2 (both converted to pu at parse)
        tap_mag = transformers.windv1[idx] / transformers.windv2[idx]

        # Phase shift in degrees from ang1
        phase_deg = transformers.ang1[idx]

        # Complex tap ratio: a = tap * exp(-j*phase)
        # Negative sign in the exponential
        tap = (1.0 / tap_mag) * exp(-jay * phase_deg * π / 180)

        # Series admittance
        z = complex(r, x)
        y_series = 1.0 / z

        # Y matrix contributions with tap ratio
        # Y = C_line * y * C_line' where C_line = C_from - C_to
        # Expanded:
        # Y[from,from] += y_series * |tap|^2
        # Y[to,to] += y_series
        # Y[from,to] -= y_series * conj(tap)
        # Y[to,from] -= y_series * tap

        # Magnetizing admittance: shunt at the winding 1 bus, not scaled by the tap
        y_mag = complex(transformers.mag1[idx], transformers.mag2[idx])

        tap_sq = abs2(tap)
        push!(I_idx, from_int); push!(J_idx, from_int); push!(V_val, y_series * tap_sq + y_mag)
        push!(I_idx, to_int); push!(J_idx, to_int); push!(V_val, y_series)
        push!(I_idx, from_int); push!(J_idx, to_int); push!(V_val, -y_series * conj(tap))
        push!(I_idx, to_int); push!(J_idx, from_int); push!(V_val, -y_series * tap)
    end
end

"""
Add the shunt (capacitor/reactor) contributions of the shunt model `name` (`fixed_shunt_models`)
to the Y matrix diagonal: the admittance G + jB of each in-service instance (`shunt_admittance`).
"""
function _add_shunt_contributions!(I_idx, J_idx, V_val, sys, name::Symbol)
    shunts = sys.models[name]
    isempty(shunts) && return
    shunt_bus = get_foreign_indices(sys, name, :i)
    params = shunt_admittance_params(shunts)

    for idx in eachindex(shunts)
        bus_internal = shunt_bus[idx]

        # Check if shunt is in service
        shunts.status[idx] || continue

        g, b = shunt_admittance(params, idx)

        # Add to diagonal
        push!(I_idx, bus_internal); push!(J_idx, bus_internal); push!(V_val, complex(g, b))
    end
end

"""
    check_stamped_network_elements(sys)

Throw an `ArgumentError` naming each model of `sys` whose Y-bus role is
`FixedNetworkElement` and that `build_full_ybus` does not stamp: its admittance would be
missing from Y. Stamped are the buses (the nodes), the branches and two-winding
transformers, the fixed shunt models (`fixed_shunt_models`), and the faults, stamped in
the topologies where they are applied (`process_ybus_events`). A model without a declared
role is left to `check_ybus_support`, which refuses it by name.
"""
function check_stamped_network_elements(sys)
    stamped = Set{Symbol}((:Bus, :Branch, :Transformer2W, :Fault))
    union!(stamped, fixed_shunt_models(sys))
    unstamped = [name for (name, model) in pairs(sys.models)
                 if declares_network_partition_role(get_base_model_type(model)) &&
                    network_partition_role(get_base_model_type(model)) isa FixedNetworkElement && !(name in stamped)]
    isempty(unstamped) || throw(ArgumentError(
        "The Y-bus builder does not stamp these network elements (FixedNetworkElement), whose " *
        "admittance would be missing from Y: " * join(unstamped, ", ")))
    return nothing
end

"""
    _add_generator_contributions!(I_idx, J_idx, V_val, sys)

Add the admittance Y_gen = 1/(ra + j*xd') (per unit on the system base) of each in-service machine
(`online_machines`, e.g. GENCLS, GENROU) to the Y matrix diagonal (COO triplets):
its Norton source is Y_gen·E. Static PV/Slack generators have no admittance in Y:
they are constraints of the network solve (`StaticGenerator`).
"""
function _add_generator_contributions!(I_idx, J_idx, V_val, sys)
    for (model_name, idx, bus_internal) in online_machines(sys)
        y_gen = 1 / get_machine_impedance(sys.models[model_name], sys.external[model_name], idx)
        push!(I_idx, bus_internal); push!(J_idx, bus_internal); push!(V_val, y_gen)
    end
end

"""
    _add_pq_z_contributions!(I_idx, J_idx, V_val, sys)

Add the constant-impedance part of the in-service PQ loads, Y_load = yp - j*yq
(yq > 0 is inductive), to the Y matrix diagonal (COO triplets). The constant-power and
constant-current parts are current injections of the network solve
(`nc_load_current`).
"""
function _add_pq_z_contributions!(I_idx, J_idx, V_val, sys)
    haskey(sys.models, :PQ) || return
    pq = sys.models[:PQ]
    isempty(pq) && return
    load_bus = device_bus_indices(sys, :PQ)

    for idx in eachindex(pq)
        pq.status[idx] || continue
        y_load = complex(pq.yp[idx], -pq.yq[idx])
        iszero(y_load) && continue
        push!(I_idx, load_bus[idx]); push!(J_idx, load_bus[idx]); push!(V_val, y_load)
    end

    @debug "Finished adding PQ Z components"
end

"""
    generator_buses(sys) -> Vector{Int}

Sorted internal indices of the buses with an in-service machine
(`online_machines`), the buses of the generator Norton sources. Their order is the
order of the generator buffers and of `NetworkSolution` generator arrays.
"""
generator_buses(sys) = sort!([bus for (_, _, bus) in online_machines(sys)])

export FullYMatrix, build_full_ybus
