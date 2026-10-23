# ============================================================================
# YBus Event Processor
# ============================================================================
#
# Processes events and pre-computes all required Y-matrix topologies.
# Expands the events to primitives (BusFault by PowerfulTDS) and generates topology keys.
#
# ============================================================================

export process_ybus_events, check_ybus_events

using PowerfulCore: SystemModel, Disturbance
using PowerfulCore: SimConfig
using PowerfulTDS: FieldModification, SetpointChange, BusFault
using PowerfulTDS: expand_bus_faults!

"""
    process_ybus_events(events, sys, base_full_ybus, topology_cache, config) -> YBusEventSchedule

Process all events and pre-compute required Y-matrix topologies.

# Algorithm
1. Expand BusFault events to primitive FieldModifications
2. Sort all events by time
3. For each unique time point:
   - Apply cumulative modifications to full Y matrix
   - Generate TopologyKey for this state
   - Store its Y matrix in the cache

# Arguments
- `events::Vector{<:Disturbance}`: User-provided events
- `sys::SystemModel`: System model (may be modified to add Fault models)
- `base_full_ybus::FullYMatrix{T}`: Pre-fault full Y matrix
- `topology_cache::TopologyCache{T}`: Cache to store computed topologies
- `config::SimConfig`: The problem's config (fault impedance warning, see `expand_bus_faults!`)

# Returns
`YBusEventSchedule{T}` with times and topology keys for simulation.
"""
function process_ybus_events(
    events::Vector{<:Disturbance},
    sys::SystemModel,
    base_full_ybus::FullYMatrix{T},
    topology_cache::TopologyCache{T},
    config::SimConfig
) where T

    isempty(events) && return nothing

    # Step 1: Expand all events to FieldModifications
    # BusFaults become FaultApply/FaultClear (adding Fault models, as in TDS)
    check_ybus_events(sys, events)
    expanded_events = Union{FieldModification, SetpointChange}[e for e in expand_bus_faults!(sys, events, config)]
    check_ybus_events(sys, expanded_events)

    isempty(expanded_events) && return nothing

    # Step 2: Sort by time
    sort!(expanded_events, by=e -> e.time)

    # Step 3: Track cumulative state and generate topologies
    # We need to track the cumulative effect of all events up to each time point

    # Active state tracking
    active_faults = Dict{Tuple, Bool}()  # fault_key → is_active
    tripped_lines = Set{Tuple}()         # branch_key (from, to, ckt)

    times = T[]
    topology_keys = TopologyKey[]

    # Group events by time (multiple events can happen at same time)
    time_groups = Dict{T, Vector{Union{FieldModification, SetpointChange}}}()
    for event in expanded_events
        t = T(event.time)
        if !haskey(time_groups, t)
            time_groups[t] = Union{FieldModification, SetpointChange}[]
        end
        push!(time_groups[t], event)
    end

    # Warn about conflicting events at the same time
    _warn_conflicting_events(time_groups)

    # Process each time point in order
    for t in sort(collect(keys(time_groups)))
        events_at_t = time_groups[t]

        # Update cumulative state
        for event in events_at_t
            _update_cumulative_state!(active_faults, tripped_lines, event)
        end

        # Generate topology key for this state
        key = _make_topology_key(active_faults, tripped_lines, sys, T)

        # Pre-compute topology if not already in cache
        _store_topology!(topology_cache, key, base_full_ybus, sys)
        @debug "Topology $(key.name) at t=$t"

        push!(times, t)
        push!(topology_keys, key)
    end

    # The events that change a model's data or setpoints and no topology are applied at their times
    # (a topology event, in TOPOLOGY_MODELS, is applied through the topology key above instead)
    changes = [e for e in expanded_events if e isa SetpointChange || e.model ∉ TOPOLOGY_MODELS]
    prefault_key = TopologyKey(:prefault)
    return YBusEventSchedule{T}(times, topology_keys, prefault_key, changes)
end

# The models whose status a FieldModification may change: the topology events
const TOPOLOGY_MODELS = (:Branch, :Fault)

"""
    check_ybus_events(sys, events)

Throw an `ArgumentError` naming each event of `events` the Y-bus simulation cannot apply: an
event of a type other than `BusFault`, `FieldModification` and `SetpointChange`, and a
`FieldModification` of the status of a model other than a Branch or a Fault (it would change
which devices are in service, which the Y-bus buffers do not follow).
"""
function check_ybus_events(sys, events)
    refused = String[]
    for event in events
        if !(event isa Union{BusFault, FieldModification, SetpointChange})
            push!(refused, "$(nameof(typeof(event))) (the events it applies are BusFault, FieldModification and SetpointChange)")
        elseif event isa FieldModification && !(event.model in TOPOLOGY_MODELS) &&
               event.field === PowerfulCore.get_status_field(sys.models[event.model])
            push!(refused, "FieldModification of $(event.model).$(event.field) at t = $(event.time) (it changes the " *
                           "devices in service, which the Y-bus integrations do not follow; only the status of a " *
                           "Branch or a Fault, a topology event, changes)")
        end
    end
    isempty(refused) || throw(ArgumentError("The Y-bus simulation (PowerfulYbus) does not apply " * join(unique(refused), "; ")))
    return nothing
end

"""
Update cumulative state based on an event.
"""
function _update_cumulative_state!(
    active_faults::Dict{Tuple, Bool},
    tripped_lines::Set{Tuple},
    event::Union{FieldModification, SetpointChange}
)
    event isa SetpointChange && return nothing  # a setpoint changes no topology
    if event.model == :Fault
        # FaultApply: value=1, FaultClear: value=0
        active_faults[event.key] = (event.value == 1.0)
    elseif event.model == :Branch
        # BranchTrip: value=0, BranchClose: value=1
        if event.value == 0.0
            push!(tripped_lines, event.key)
        else
            delete!(tripped_lines, event.key)
        end
    end
end

"""
Warn about conflicting events at the same time (e.g., FaultApply + FaultClear for same fault).
"""
function _warn_conflicting_events(time_groups::Dict{T, <:Vector}) where T
    for (t, events) in time_groups
        length(events) <= 1 && continue

        # Group by (model, key), and the setpoint of a setpoint change, to find conflicts
        seen = Dict{Tuple{Symbol, Tuple, Union{Nothing, Symbol}}, Vector{Float64}}()
        for e in events
            k = (e.model, e.key, e isa SetpointChange ? e.setpoint : nothing)
            if !haskey(seen, k)
                seen[k] = Float64[]
            end
            push!(seen[k], e.value)
        end

        # Warn if same (model, key) has different values at same time
        for ((model, key, _), values) in seen
            if length(unique(values)) > 1
                @warn "Conflicting events at t=$t for $model $key: values=$values. Last event wins."
            end
        end
    end
end

"""
Generate a TopologyKey from cumulative state.
"""
function _make_topology_key(
    active_faults::Dict{Tuple, Bool},
    tripped_lines::Set{Tuple},
    sys::SystemModel,
    ::Type{T}
) where T
    # Collect active faults with their admittances
    faults_list = Tuple{Int, Complex{T}}[]
    for (key, is_active) in active_faults
        is_active && push!(faults_list, _fault_admittance(sys, get_index_by_key(sys, :Fault, key), T))
    end

    # Collect tripped line indices
    lines_out = Int[]
    for key in tripped_lines
        push!(lines_out, get_index_by_key(sys, :Branch, key))
    end

    # Generate descriptive name
    n_faults = count(values(active_faults))
    n_trips = length(tripped_lines)
    if n_faults == 0 && n_trips == 0
        name = :prefault
    elseif n_faults > 0 && n_trips == 0
        name = Symbol("fault_on_$(n_faults)")
    elseif n_faults == 0 && n_trips > 0
        name = Symbol("post_fault_$(n_trips)_lines_out")
    else
        name = Symbol("fault_$(n_faults)_lines_$(n_trips)")
    end

    return TopologyKey(name; lines_out=lines_out, faults=faults_list)
end



"""
    apply_changes!(prob, applied, t, runtime, models_by_name; residual = nothing) -> Bool

Apply the events of the schedule of `prob` that change a model's data or setpoints and no
topology (`YBusEventSchedule` `changes`: `FieldModification`s and `SetpointChange`s) that
are due at `t` and not yet applied (`applied`, the number applied so far, updated): the
fields to the problem's system, the setpoints to the buffers of every number type that
evaluates its equations (`apply_setpoint!`); and trace again what the automata of the
models they change read (`retrace_automata!`, the code of the phasor simulation's
`AutomatonSettle`) with `runtime` (an `AutomatonRuntime` or `nothing`) over the model arrays
`models_by_name` of its buffers. `residual` is the integration's `YBusResidual`, if it has
one yet (its buffers hold setpoints of their own). Returns whether any was due. Both
integrations call it where a topology event can fire, before they switch the topology. (An
event on a status field is refused when the problem is built, `check_ybus_events`.)
"""
function apply_changes!(prob, applied::Base.RefValue{Int}, t, runtime, models_by_name; residual = nothing)
    schedule = prob.event_schedule
    schedule === nothing && return false
    changes = schedule.changes
    first_due = applied[] + 1
    last_due = applied[]
    while last_due < length(changes) && changes[last_due + 1].time <= t
        last_due += 1
    end
    last_due < first_due && return false
    due = changes[first_due:last_due]
    sys = prob.base_prob.sys
    for event in due
        if event isa SetpointChange
            apply_setpoint!(prob, residual, event)
            continue
        end
        PowerfulTDS.apply!(sys, event)
    end
    applied[] = last_due
    invalidate!(prob.algebraic)
    runtime === nothing || retrace_automata!(runtime, models_by_name, _event_models(runtime, due, sys),
                                             _event_models(runtime, due, sys; setpoints = true))
    return true
end

"""
    apply_setpoint!(prob, residual, change::SetpointChange)

Apply `change` to the setpoint it names in every buffer that evaluates the equations of
`prob`: its own (`prob.sys_ws`, the value `:add` and `:multiply` start from), those of the
algebraic solve's Dual evaluations, and those of `residual` (a `YBusResidual` or `nothing`).
Buffers built later copy the setpoints of `prob.sys_ws` (`new_buffers`).
"""
function apply_setpoint!(prob, residual, change::SetpointChange)
    instance = get_index_by_key(prob.base_prob.sys, change.model, getfield(change, :key))
    current = setpoint_array(prob.sys_ws, change.model, change.setpoint)[instance]
    value = PowerfulTDS.modified_value(current, change)
    buffers = Any[prob.sys_ws]
    append!(buffers, (w.ws.sys_ws for w in values(prob.algebraic.workspaces)))
    if residual !== nothing
        push!(buffers, residual.float.sys_ws)
        append!(buffers, (w.sys_ws for w in values(residual.dual)))
    end
    for ws in buffers
        setpoint_array(ws, change.model, change.setpoint)[instance] = value
    end
    return nothing
end
