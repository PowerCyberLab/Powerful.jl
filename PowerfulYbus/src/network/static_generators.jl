# ============================================================================
# Static Generator Constraints
# ============================================================================
#
# Power flow generators left static in the TDS system (PV and Slack instances
# without a dynamic model) hold their bus as in the phasor TDS formulation. In
# the network solve they replace the KCL equations of their bus: the current
# they inject is whatever balances the bus, subject to
#   PV:    Re(V I*) = Σ pg   and   |V| = vs
#   Slack: V = vs∠va
#
# ============================================================================

using PowerfulPF: PV, Slack

"""
    StaticGeneratorBuses

Constraints of the in-service static generators (`StaticGenerator` role) at their
buses.

# Fields
- `pv_bus::Vector{Int}`: Internal indices of the PV-controlled buses
- `pv_p::Vector{Float64}`: Total active power `pg` of the PV generators at each
- `pv_vs::Vector{Float64}`: Their voltage magnitude setpoint `vs`
- `slack_bus::Vector{Int}`: Internal indices of the Slack buses
- `slack_V::Vector{ComplexF64}`: Their voltage `vs∠va`

A bus with a Slack generator holds its voltage; PV generators at the same bus then
set nothing further (their active power is absorbed like the Slack's).
"""
struct StaticGeneratorBuses
    pv_bus::Vector{Int}
    pv_p::Vector{Float64}
    pv_vs::Vector{Float64}
    slack_bus::Vector{Int}
    slack_V::Vector{ComplexF64}
end

Base.isempty(s::StaticGeneratorBuses) = isempty(s.pv_bus) && isempty(s.slack_bus)

"""
    static_generator_buses(sys) -> StaticGeneratorBuses

Collect the in-service static generators of `sys` by bus. Errors if generators at
one bus have different voltage setpoints `vs`.
"""
function static_generator_buses(sys)
    pv = Dict{Int, Tuple{Float64, Float64}}()       # bus => (Σ pg, vs)
    slack = Dict{Int, ComplexF64}()                  # bus => vs∠va
    for model_name in keys(sys.models)
        model = sys.models[model_name]
        network_partition_role(eltype(model)) isa StaticGenerator || continue
        isempty(model) && continue
        bus = device_bus_indices(sys, model_name)
        ext = getfield(sys.external, model_name)
        active = in_service(model)
        for idx in eachindex(model)
            active[idx] || continue
            _add_static_generator!(pv, slack, model, ext, idx, bus[idx])
        end
    end
    for (bus, V) in slack
        haskey(pv, bus) && _check_setpoint(bus, pv[bus][2], abs(V))
        delete!(pv, bus)
    end
    pv_bus = sort!(collect(keys(pv)))
    slack_bus = sort!(collect(keys(slack)))
    return StaticGeneratorBuses(pv_bus, [pv[b][1] for b in pv_bus], [pv[b][2] for b in pv_bus],
                                slack_bus, [slack[b] for b in slack_bus])
end

function _add_static_generator!(pv, slack, model::StructArray{<:PV}, ext, idx::Int, bus::Int)
    p, vs = get(pv, bus, (0.0, model.vs[idx]))
    _check_setpoint(bus, vs, model.vs[idx])
    pv[bus] = (p + model.pg[idx], vs)
end

# The angle reference `va` is the Slack's external parameter (its bus angle), as in
# the power flow Slack equations
function _add_static_generator!(pv, slack, model::StructArray{<:Slack}, ext, idx::Int, bus::Int)
    V = model.vs[idx] * cis(ext.va[idx])
    haskey(slack, bus) && !(slack[bus] ≈ V) &&
        error("Slack generators at internal bus $bus hold different voltages " *
              "$(slack[bus]) and $V")
    slack[bus] = V
end

_check_setpoint(bus, vs1, vs2) = vs1 ≈ vs2 ||
    error("Static generators at internal bus $bus have different voltage " *
          "setpoints vs = $vs1 and $vs2")

export StaticGeneratorBuses, static_generator_buses
