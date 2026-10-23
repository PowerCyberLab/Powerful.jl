# ============================================================================
# YBus Event Schedule and Time Segments
# ============================================================================
#
# Data structures for tracking topology changes during simulation.
# Implements time segment handling with dt adjustment.
#
# ============================================================================

export YBusEventSchedule, TimeSegment
export get_topology_at, has_topology_change, build_time_segments

# Note: TopologyKey is defined in ybus/topology_cache.jl and available in module scope

"""
    TimeSegment{T}

A time interval with uniform step size for integration.

Each segment between events has a (possibly adjusted) dt that divides the interval
evenly, ensuring the integration lands exactly on event times.

# Fields
- `t_start::T`: Start time of segment
- `t_end::T`: End time of segment
- `dt::T`: Time step for this segment (adjusted to divide evenly)
- `n_steps::Int`: Number of integration steps in this segment
"""
struct TimeSegment{T}
    t_start::T
    t_end::T
    dt::T
    n_steps::Int

    function TimeSegment{T}(t_start, t_end, dt, n_steps) where T
        t_end > t_start || error("t_end must be greater than t_start")
        dt > 0 || error("dt must be positive")
        n_steps >= 1 || error("n_steps must be at least 1")
        new{T}(T(t_start), T(t_end), T(dt), n_steps)
    end
end

TimeSegment(t_start::T, t_end::T, dt::T, n_steps::Int) where T = TimeSegment{T}(t_start, t_end, dt, n_steps)

"""
    YBusEventSchedule{T}

Tracks which topology is active at each time point.
Pre-computed before simulation starts.

# Fields
- `times::Vector{T}`: Sorted event times (topology changes at these times)
- `topology_keys::Vector{TopologyKey}`: Topology active AFTER each time
- `prefault_key::TopologyKey`: Initial topology (before any events)
- `changes::Vector{Union{FieldModification, SetpointChange}}`: The events that change a
  model's data or setpoints and no topology (a parameter of a model with automata, say),
  sorted by time; their times are
  among `times`, with the topology unchanged (`apply_changes!`)

# Usage
```julia
# Get topology at a specific time
key = get_topology_at(schedule, t)

# Check if topology changes in an interval
if has_topology_change(schedule, t, t_next)
    new_key = get_topology_at(schedule, t_next)
    switch_topology!(cache, new_key)
end
```
"""
struct YBusEventSchedule{T}
    times::Vector{T}
    topology_keys::Vector{TopologyKey}
    prefault_key::TopologyKey
    changes::Vector{Union{FieldModification, SetpointChange}}

    function YBusEventSchedule{T}(times, topology_keys, prefault_key, changes = FieldModification[]) where T
        length(times) == length(topology_keys) ||
            error("times and topology_keys must have same length")
        issorted(times) || error("times must be sorted")
        issorted(changes; by = e -> e.time) || error("changes must be sorted by time")
        new{T}(Vector{T}(times), topology_keys, prefault_key,
               Vector{Union{FieldModification, SetpointChange}}(changes))
    end
end

YBusEventSchedule(times::Vector{T}, topology_keys, prefault_key, changes = FieldModification[]) where T =
    YBusEventSchedule{T}(times, topology_keys, prefault_key, changes)

"""
    get_topology_at(schedule::YBusEventSchedule, t::Real) -> TopologyKey

Get the topology active at time `t`.

Returns prefault_key if t is before all events, otherwise returns
the topology that became active at the most recent event before or at t.
"""
function get_topology_at(schedule::YBusEventSchedule, t::Real)
    idx = searchsortedlast(schedule.times, t)
    return idx == 0 ? schedule.prefault_key : schedule.topology_keys[idx]
end

"""
    has_topology_change(schedule::YBusEventSchedule, t::Real, t_next::Real) -> Bool

Check if topology changes in the interval (t, t_next].

Returns true if there is at least one event time in (t, t_next].
"""
function has_topology_change(schedule::YBusEventSchedule, t::Real, t_next::Real)
    i1 = searchsortedlast(schedule.times, t)
    i2 = searchsortedlast(schedule.times, t_next)
    return i2 > i1
end

"""
    build_time_segments(tspan, event_times, nominal_dt; min_dt_ratio=0.5, warn_threshold=0.1)

Build time segments with adjusted dt plus protections.

Divides the simulation time span into segments based on event times.
Each segment has a uniform dt adjusted to divide the interval evenly,
ensuring we land exactly on event times.

# Arguments
- `tspan::Tuple{T,T}`: Simulation time span (t0, tf)
- `event_times::Vector{T}`: Sorted event times
- `nominal_dt::T`: User-specified nominal time step
- `min_dt_ratio::T`: Minimum allowed dt as fraction of nominal (default 0.5)
- `warn_threshold::T`: Warn if dt differs from nominal by more than this (default 0.1 = 10%)

# Returns
`Vector{TimeSegment{T}}` with one segment per interval between events.

# Example
```julia
tspan = (0.0, 5.0)
event_times = [1.0, 1.1]
segments = build_time_segments(tspan, event_times, 0.01)
# Returns 3 segments: [0,1], [1,1.1], [1.1,5]
# Each with dt adjusted to divide interval evenly
```

"""
function build_time_segments(
    tspan::Tuple{T,T},
    event_times::Vector{T},
    nominal_dt::T;
    min_dt_ratio::T = T(0.5),
    warn_threshold::T = T(0.1)
) where T
    segments = TimeSegment{T}[]
    min_dt = min_dt_ratio * nominal_dt

    # Build sorted list of all time points
    all_times = sort(unique([tspan[1]; event_times; tspan[2]]))

    # Filter out times outside tspan
    all_times = filter(t -> tspan[1] <= t <= tspan[2], all_times)

    i = 1
    while i < length(all_times)
        t_start = all_times[i]
        t_end = all_times[i + 1]
        interval = t_end - t_start

        # Protection: merge tiny segments with next
        if interval < min_dt && i + 1 < length(all_times)
            # Merge this interval into next by extending next segment backward
            all_times[i + 1] = t_start
            i += 1
            @debug "Merged tiny segment [$t_start, $t_end] (interval=$interval < min_dt=$min_dt)"
            continue
        end

        # Calculate steps and adjusted dt
        n_steps = max(1, round(Int, interval / nominal_dt))
        actual_dt = interval / n_steps

        # Warning if dt differs significantly from nominal
        dt_diff = abs(actual_dt - nominal_dt) / nominal_dt
        if dt_diff > warn_threshold
            @warn "Segment [$t_start, $t_end]: dt=$(round(actual_dt, sigdigits=4)) differs from nominal=$(nominal_dt) by $(round(dt_diff*100, digits=1))%"
        end

        push!(segments, TimeSegment{T}(t_start, t_end, actual_dt, n_steps))
        i += 1
    end

    return segments
end

# Convenience method without event times
function build_time_segments(tspan::Tuple{T,T}, nominal_dt::T; kwargs...) where T
    return build_time_segments(tspan, T[], nominal_dt; kwargs...)
end
