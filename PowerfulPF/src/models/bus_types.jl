"""
Bus classification for power flow.

Bus types are derived from in-service devices (Slack, PV), the same way the
power-flow residuals define them. `Bus.ide` from the data file is NOT used: a bus
marked PV whose generators are all out of service behaves as a PQ bus.
"""

"""
    _active_device_buses(sys, model) -> Set{Int}

Internal Bus indices hosting at least one in-service instance of `model`.
"""
function _active_device_buses(sys::SystemModel, model::Symbol)
    buses = Set{Int}()
    haskey(sys.models, model) || return buses
    devices = sys.models[model]
    isempty(devices) && return buses

    bus_idx = get_foreign_indices(sys, model, :i)
    status = getproperty(devices, get_status_field(devices))
    for k in eachindex(devices)
        status[k] && push!(buses, bus_idx[k])
    end
    return buses
end

"""
    get_slack_bus(sys::SystemModel) -> Int

Internal index of the reference bus: the bus hosting the in-service Slack
generator(s). Errors unless exactly one such bus exists.
"""
function get_slack_bus(sys::SystemModel)
    buses = _active_device_buses(sys, :Slack)
    length(buses) == 1 ||
        error("Expected exactly one bus with an in-service Slack generator, found $(length(buses))")
    return only(buses)
end

"""
    get_pv_buses(sys::SystemModel) -> Vector{Int}

Sorted internal indices of voltage-controlled buses: buses with an in-service PV
generator, excluding the slack bus.
"""
function get_pv_buses(sys::SystemModel)
    pv = _active_device_buses(sys, :PV)
    delete!(pv, get_slack_bus(sys))
    return sort!(collect(pv))
end

"""
    get_pq_buses(sys::SystemModel) -> Vector{Int}

Sorted internal indices of load buses: every bus that is neither PV nor slack.
"""
function get_pq_buses(sys::SystemModel)
    controlled = _active_device_buses(sys, :PV)
    push!(controlled, get_slack_bus(sys))
    return [i for i in 1:length(sys.models.Bus) if !(i in controlled)]
end

export get_pq_buses, get_pv_buses, get_slack_bus
