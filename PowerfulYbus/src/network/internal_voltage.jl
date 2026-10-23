# ============================================================================
# Internal Voltage Gathering
# ============================================================================
#
# Gathers internal voltages E_internal from generator states.
# E_internal = |E| * exp(j*delta) in system reference frame.
#
# ============================================================================

using PowerfulTDS: GENCLS, GENROU, GENSAL, GENTPJ1

# Machines behind their subtransient impedance ra + j*Xd2 whose stator current (Id, Iq) is a
# machine variable: the Norton source E = V + (ra + j*Xd2)*I of the network is exact whatever
# the machine's other reactances (GENTPJ1's Xq2 among them), since I is the machine's own current
const SubtransientMachine = Union{GENROU, GENSAL, GENTPJ1}

"""
    InternalVoltageBuffer{T}

Buffer for storing gathered internal voltages from generators.

# Fields
- `E_internal::Vector{Complex{T}}`: Internal voltages E∠δ
- `gen_bus_indices::Vector{Int}`: Internal bus indices for each generator
- `gen_model_map::Vector{Tuple{Symbol,Int}}`: Per-instance mapping (model_name, model_idx) for each E_internal index
- `Z_gen::Vector{Complex{T}}`: Stator impedances (r_a + j*xd_p) for each generator
- `gathered_models::Set{Symbol}`: Pre-allocated set for tracking gathered models (reused each call)
"""
struct InternalVoltageBuffer{T}
    E_internal::Vector{Complex{T}}
    gen_bus_indices::Vector{Int}
    gen_model_map::Vector{Tuple{Symbol,Int}}
    Z_gen::Vector{Complex{T}}
    gathered_models::Set{Symbol}
end

"""
    InternalVoltageBuffer(::Type{T}, n_gen::Int) where T

Create buffer for `n_gen` generators.
"""
function InternalVoltageBuffer(::Type{T}, n_gen::Int) where T
    return InternalVoltageBuffer{T}(
        zeros(Complex{T}, n_gen),
        zeros(Int, n_gen),
        Tuple{Symbol,Int}[],
        zeros(Complex{T}, n_gen),
        Set{Symbol}()
    )
end

"""
    setup_internal_voltage_buffer!(buf::InternalVoltageBuffer, gen_buses, sys)

Setup buffer with generator bus mapping ordered by `gen_buses` (internal bus indices).

The gen_model_map is built to match the order of gen_buses, ensuring
E_internal[i] corresponds to the generator at gen_buses[i].
"""
function setup_internal_voltage_buffer!(buf::InternalVoltageBuffer{T},
                                        gen_buses::AbstractVector{Int}, sys) where T
    n_gen = length(gen_buses)
    resize!(buf.E_internal, n_gen)
    resize!(buf.gen_bus_indices, n_gen)
    resize!(buf.Z_gen, n_gen)
    resize!(buf.gen_model_map, n_gen)
    copyto!(buf.gen_bus_indices, gen_buses)

    gens = generators_by_bus(sys)
    for (E_idx, bus) in enumerate(gen_buses)
        buf.gen_model_map[E_idx] = gens[bus]
    end
    return buf
end

"""
    populate_generator_impedances!(buf::InternalVoltageBuffer, models, external)

Populate Z_gen (stator impedances, per unit on the system base) of the generators in
`buf.gen_model_map` from `models` (anything indexed by model name, e.g. `sys.models` or
`sys_ws.models_by_name`) and their external parameters `external` (`sys.external`), which hold
the ratio of the machine base to the system base (`Mratio`).

Must be called after setup_internal_voltage_buffer!.
Uses the `get_machine_impedance` trait to get Z_gen = r_a + j*xd' for each generator type.
"""
function populate_generator_impedances!(buf::InternalVoltageBuffer{T}, models, external) where T
    for (E_idx, (model_name, model_idx)) in enumerate(buf.gen_model_map)
        buf.Z_gen[E_idx] = get_machine_impedance(models[model_name], external[model_name], model_idx)
    end
    return buf
end

# ============================================================================
# Model-specific get_machine_impedance implementations
# ============================================================================

"""
    get_machine_impedance(model::StructArray{<:GENCLS}, external, idx::Int) -> Complex

Get machine impedance Z = r_a + j*xd' for GENCLS generator, per unit on the system base: the
machine's own (on its MVA base) divided by `external.Mratio` = Sn/Sb.
"""
function get_machine_impedance(model::StructArray{<:GENCLS}, external, idx::Int)
    return complex(model.r_a[idx], model.xd_p[idx]) / external.Mratio[idx]
end

"""
    get_machine_impedance(model::StructArray{<:SubtransientMachine}, external, idx::Int) -> Complex

Get machine impedance Z = ra + j*Xd2 for a GENROU, GENSAL or GENTPJ1 generator, per unit on
the system base (the machine's own divided by `external.Mratio` = Sn/Sb).
Uses subtransient reactance Xd2 (xd'') for Y-bus formulation.
"""
function get_machine_impedance(model::StructArray{<:SubtransientMachine}, external, idx::Int)
    return complex(model.ra[idx], model.Xd2[idx]) / external.Mratio[idx]
end


"""
    gather_internal_voltages!(buf::InternalVoltageBuffer, sys_ws, dae)

Set `buf.E_internal` from the machine variables in `dae` (x, y), by the internal
voltage gatherer of each machine model in `sys_ws` (`get_internal_voltage` of its
in-service instances, in network generator order).

# Arguments
- `buf`: InternalVoltageBuffer to fill
- `sys_ws`: YBusSystemBuffers with the internal voltage gatherers
- `dae`: DAE with x, y vectors
"""
function gather_internal_voltages!(buf::InternalVoltageBuffer{T},
                                   sys_ws, dae) where T
    # Extract vectors (supports both full DAE and algebraic-only)
    x = hasfield(typeof(dae), :x) ? dae.x : T[]
    y = dae.y
    for gather! in sys_ws.internal_voltage_gatherers
        gather!(buf.E_internal, x, y)
    end
    nothing
end

"""
    gather_internal_voltages_from_pf!(buf::InternalVoltageBuffer, pf_sol, sys)

Initialize internal voltages from power flow solution.

This creates the initial E_internal values based on power flow results,
assuming generators are operating at their power flow setpoints.

For classical model approximation:
- E_internal ≈ V_bus + j*xd' * I_gen
"""
function gather_internal_voltages_from_pf!(buf::InternalVoltageBuffer{T},
                                           V_bus::AbstractVector{Complex{T}},
                                           sys) where T
    # For power flow initialization, use bus voltage as E_internal
    # This is a simplification - proper initialization would compute
    # E_internal from V, P, Q and machine reactance
    for (buf_idx, bus_idx) in enumerate(buf.gen_bus_indices)
        buf.E_internal[buf_idx] = V_bus[bus_idx]
    end
end

# ============================================================================
# Model-specific get_internal_voltage implementations
# ============================================================================

"""
    get_internal_voltage(ws, model::StructArray{<:GENCLS}, idx::Int) -> Complex

Get internal voltage for GENCLS generator.

E_internal = E_q * cis(delta)

where E_q is the internal voltage magnitude and delta is the rotor angle.
"""
function get_internal_voltage(ws, model::StructArray{<:GENCLS}, idx::Int)
    # GENCLS stores E_q (internal voltage magnitude) and delta (rotor angle)
    # Access from workspace var_vals
    E_q = ws.var_vals.E_q[idx]
    delta = ws.var_vals.delta[idx]
    return E_q * cis(delta)
end

"""
    get_internal_voltage(ws, model::StructArray{<:SubtransientMachine}, idx::Int) -> Complex

Get internal voltage for a GENROU, GENSAL or GENTPJ1 generator using V + Z*I formula.

The "internal voltage behind subtransient impedance" is computed from:
    E = V_terminal + Z_gen * I_stator

where:
- V_terminal = v * cis(theta) from algebraic variables
- I_stator = (Id - j*Iq) * cis(delta - π/2) transformed to network frame
- Z_gen = ra + j*Xd''

This approach uses values from the previous timestep, but is consistent with
the machine model equations.
"""
function get_internal_voltage(ws, model::StructArray{<:SubtransientMachine}, idx::Int)
    # Get algebraic variables (from previous timestep)
    v = ws.var_vals.v[idx]
    theta = ws.var_vals.theta[idx]
    Id = ws.var_vals.Id[idx]
    Iq = ws.var_vals.Iq[idx]
    delta = ws.var_vals.delta[idx]

    # Get impedance parameters
    ra = model.ra[idx]
    Xd2 = model.Xd2[idx]
    Z_gen = complex(ra, Xd2)

    # Terminal voltage in network frame
    V_terminal = v * cis(theta)

    # Transform stator current from d-q to network frame
    # From TDS init: I_dq = I_gen * cis(-delta + π/2), so I_gen = I_dq * cis(delta - π/2)
    # where I_dq = Id + j*Iq
    I_dq = complex(Id, Iq)
    I_stator = I_dq * cis(delta - π/2)

    # E_internal = V + Z*I
    E_internal = V_terminal + Z_gen * I_stator

    return E_internal
end

# ============================================================================
# Terminal Voltage Gathering (from TDS algebraic variables)
# ============================================================================

"""
    gather_terminal_voltages!(V_terminal, E_buf, sys_ws, dae)

Gather generator terminal voltages from TDS algebraic variables (v, theta).

This is called BEFORE the network solve to get V_terminal from the current
TDS state. The V_terminal values are used in the NC load solver and to compute
E_internal using the Y-bus network equations.

# Algorithm
For each generator model, gather (v, theta) from algebraic variables and compute:
    V_terminal = v * cis(theta)

# Arguments
- `V_terminal`: Output vector for terminal voltages (modified in place)
- `E_buf`: InternalVoltageBuffer with generator model map
- `sys_ws`: YBusSystemBuffers with workspaces
- `dae`: DAE with x, y vectors
"""
function gather_terminal_voltages!(V_terminal::AbstractVector{Complex{T}},
                                   E_buf::InternalVoltageBuffer{T},
                                   sys_ws, dae) where T
    x = hasfield(typeof(dae), :x) ? dae.x : T[]
    y = dae.y

    # Track which models have been gathered (call gather_variables! once per model)
    # Reuse pre-allocated Set to avoid allocation
    empty!(E_buf.gathered_models)

    for (E_idx, (model_name, model_idx)) in enumerate(E_buf.gen_model_map)
        # Gather variables once per model type
        if !(model_name in E_buf.gathered_models)
            ws = sys_ws.workspaces_by_name[model_name]
            meta = sys_ws.metadata_by_name[model_name]
            gather_variables!(ws, x, y)
            push!(E_buf.gathered_models, model_name)
        end

        # Get terminal voltage for this generator
        ws = sys_ws.workspaces_by_name[model_name]
        v = ws.var_vals.v[model_idx]
        theta = ws.var_vals.theta[model_idx]
        V_terminal[E_idx] = v * cis(theta)
    end
end

export InternalVoltageBuffer, setup_internal_voltage_buffer!
export gather_internal_voltages!, gather_internal_voltages_from_pf!
export gather_terminal_voltages!
export get_machine_impedance, populate_generator_impedances!
