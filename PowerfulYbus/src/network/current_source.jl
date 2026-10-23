# ============================================================================
# Current Source Gathering
# ============================================================================
#
# Gathers current injections I_src from current source models (e.g., REGCA1).
# Current sources inject I = (Ip + j*Iq) * exp(j*theta) at their terminal bus.
#
# This is parallel to InternalVoltageBuffer for generators.
#
# ============================================================================

using PowerfulTDS: AbstractREGCA1, REGFM_A1

"""
    CurrentSourceBuffer{T}

Buffer for storing gathered current components from current source models.

Current sources (e.g., REGCA1) inject I = (Ip - j*Iq) * exp(j*theta) where theta
is the terminal voltage angle (`source_current`). We store Ip and Iq separately so
the network Newton iteration uses the current iterate's voltage angle, not a stale
angle.

# Fields
- `Ip::Vector{T}`: Active current components (in-phase with voltage)
- `Iq::Vector{T}`: Reactive current components (positive = Q generation)
- `src_bus::Vector{Int}`: Internal bus index of each source (0 for an
  out-of-service source, which injects nothing)
- `src_model_map::Vector{Tuple{Symbol,UnitRange{Int}}}`: Model name to buffer index range
"""
struct CurrentSourceBuffer{T}
    Ip::Vector{T}
    Iq::Vector{T}
    src_bus::Vector{Int}
    src_model_map::Vector{Tuple{Symbol,UnitRange{Int}}}
end

"""
    CurrentSourceBuffer(::Type{T}, n_src::Int) where T

Create buffer for `n_src` current sources.
"""
function CurrentSourceBuffer(::Type{T}, n_src::Int) where T
    return CurrentSourceBuffer{T}(
        zeros(T, n_src),      # Ip
        zeros(T, n_src),      # Iq
        zeros(Int, n_src),
        Tuple{Symbol,UnitRange{Int}}[]
    )
end

"""
    setup_current_source_buffer!(buf::CurrentSourceBuffer, sys)

Setup buffer with current source model mapping.

Discovers current source models (`nc_load_type` is `CurrentSourceNCLoad`) and maps
each in-service source to its internal bus index.
"""
function setup_current_source_buffer!(buf::CurrentSourceBuffer{T}, sys) where T
    # Count current sources and build model map
    src_count = 0
    empty!(buf.src_model_map)

    for model_name in keys(sys.models)
        model = sys.models[model_name]
        model_type = eltype(model)
        network_partition_role(model_type) isa NonConformingLoad || continue
        nc_load_type(model_type) isa CurrentSourceNCLoad || continue

        n_model = length(model)
        if n_model > 0
            idx_range = (src_count + 1):(src_count + n_model)
            push!(buf.src_model_map, (model_name, idx_range))
            src_count += n_model
        end
    end

    # Resize buffers
    resize!(buf.Ip, src_count)
    resize!(buf.Iq, src_count)
    resize!(buf.src_bus, src_count)
    fill!(buf.Ip, zero(T))
    fill!(buf.Iq, zero(T))
    fill!(buf.src_bus, 0)

    # Populate the buses of the in-service sources
    for (model_name, idx_range) in buf.src_model_map
        src_bus = device_bus_indices(sys, model_name)
        for idx in nonlinear_instances(sys, model_name)  # instance idx is buffer entry idx_range[idx]
            buf.src_bus[idx_range[idx]] = src_bus[idx]
        end
    end

    return buf
end

"""
    gather_source_currents!(buf::CurrentSourceBuffer, sys_ws, dae)

Gather current components (Ip, Iq) from all current source models.

Uses the workspace approach:
1. Call gather_variables! to populate ws.var_vals from DAE vectors
2. Call get_source_current_components trait to get (Ip, Iq) from gathered values and
   the model's external parameters (system base)

The voltage angle theta is NOT used here - it comes from the network Newton
iterate. This ensures current sources use the latest
voltage angle, not a stale value from the previous timestep.

# Arguments
- `buf`: CurrentSourceBuffer to fill with Ip, Iq components
- `sys_ws`: YBusSystemBuffers with workspaces and models for each source
- `dae`: DAE with x, y vectors
"""
function gather_source_currents!(buf::CurrentSourceBuffer{T},
                                 sys_ws, dae) where T
    # Extract vectors (supports both full DAE and algebraic-only)
    x = hasfield(typeof(dae), :x) ? dae.x : T[]
    y = dae.y

    for (model_name, idx_range) in buf.src_model_map
        # Get workspace, metadata, and model from sys_ws
        ws = sys_ws.workspaces_by_name[model_name]
        meta = sys_ws.metadata_by_name[model_name]
        model = sys_ws.models_by_name[model_name]
        external = sys_ws.externals_by_name[model_name]

        # Gather variables into workspace using proper infrastructure
        gather_variables!(ws, x, y)

        # Use trait to get current components (Ip, Iq) from gathered values
        for (buf_idx, model_idx) in zip(idx_range, eachindex(model))
            Ip, Iq = get_source_current_components(ws, model, external, model_idx)
            buf.Ip[buf_idx] = Ip
            buf.Iq[buf_idx] = Iq
        end
    end
end

# ============================================================================
# Model-specific get_source_current_components implementations
# ============================================================================

"""
    get_source_current_components(ws, model::StructArray{<:Union{AbstractREGCA1, REGFM_A1}}, external, idx::Int) -> (Ip, Iq)

The current components (Ip, Iq) that a converter delivers at its terminal (`source_current`),
per unit on the system base: its active and reactive power injected into the bus, the algebraic
variables `Pe` and `Qe` on its own base (REGCA1 after its low-voltage active current and
high-voltage reactive current managements, REGFM_A1 behind its coupling reactance and current
limit), times `external.Mratio` (its base over the system base), over its terminal voltage
magnitude `v`, as `gather_variables!` gathered them, so that `V I* = Mratio (Pe + j Qe)` at that
voltage, its injection into the bus in the phasor TDS.
"""
function get_source_current_components(ws, model::StructArray{<:Union{AbstractREGCA1, REGFM_A1}}, external, idx::Int)
    scale = external.Mratio[idx] / ws.var_vals.v[idx]
    return (ws.var_vals.Pe[idx] * scale, ws.var_vals.Qe[idx] * scale)
end

export CurrentSourceBuffer, setup_current_source_buffer!
export gather_source_currents!, get_source_current_components
