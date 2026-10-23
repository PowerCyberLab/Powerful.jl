# ============================================================================
# Bus Voltage Inputs
# ============================================================================
#
# Controller variables that are bus voltages (foreign references to Bus.v and
# Bus.theta, e.g. SEXS.v) take their values from the bus voltages of the network
# solution, the one source of bus voltage in the Y-bus formulation. They are set
# after the controller gathers its variables from the DAE vectors, replacing the
# gathered values.
#
# ============================================================================

using PowerfulCore: ModelMetadata, get_base_model_type, get_group_ref

"""
    BusVoltageInput{A}

One controller variable that is a bus voltage.

# Fields
- `dest::A`: The variable's values in the controller workspace (`ws.var_vals`)
- `buses::Vector{Int}`: Internal bus index of each instance
- `angle::Bool`: `true` for the voltage angle (`Bus.theta`), `false` for the
  magnitude (`Bus.v`)
"""
struct BusVoltageInput{A<:AbstractVector}
    dest::A
    buses::Vector{Int}
    angle::Bool
end

"""
    bus_voltage_inputs(ws, metadata::ModelMetadata, sys, model_name) -> Vector{BusVoltageInput}

The variables of `model_name` that are foreign references to a network element with
no state in the Y-bus formulation (`FixedNetworkElement`, i.e. Bus), directly or
through a group (e.g. `:bus_group`), as resolved per instance by the system's
foreign references. Errors if a reference mixes network elements and other models,
or is to a quantity other than the bus voltage magnitude `v` or angle `theta`.
"""
function bus_voltage_inputs(ws, metadata::ModelMetadata, sys, model_name::Symbol)
    inputs = BusVoltageInput{eltype(values(ws.var_vals))}[]
    for var in metadata.vars
        var.source_model === nothing && continue
        link = get_group_ref(sys.order, model_name, var.ref_key)
        network = [network_partition_role(get_base_model_type(sys.models[m])) isa FixedNetworkElement
                   for m in link.target_type]
        any(network) || continue
        all(network) ||
            error("$model_name.$(var.name) references both Y-bus network elements and other models")
        var.source_name in (:v, :theta) ||
            error("$model_name.$(var.name) references $(var.source_model).$(var.source_name); a " *
                  "Y-bus network element provides only the bus voltage v and angle theta")
        push!(inputs, BusVoltageInput(getfield(ws.var_vals, var.name), link.target_idx,
                                      var.source_name === :theta))
    end
    return inputs
end

"""
    set_bus_voltages!(inputs::Vector{<:BusVoltageInput}, V_bus)

Set each bus voltage input from the complex bus voltages `V_bus`: `abs` for a
magnitude, `angle` for an angle.
"""
function set_bus_voltages!(inputs::Vector{<:BusVoltageInput}, V_bus::AbstractVector{<:Complex})
    for input in inputs
        dest, buses = input.dest, input.buses
        if input.angle
            @inbounds for i in eachindex(buses)
                dest[i] = angle(V_bus[buses[i]])
            end
        else
            @inbounds for i in eachindex(buses)
                dest[i] = abs(V_bus[buses[i]])
            end
        end
    end
    nothing
end

export BusVoltageInput, bus_voltage_inputs, set_bus_voltages!
