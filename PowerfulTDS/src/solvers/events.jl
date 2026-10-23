# ============================================================================
# Time-Domain Simulation Events/Disturbances
# ============================================================================
#
# This module provides a unified event system for TDS based on FieldModification.
# All event types (trips, faults, setpoint changes) reduce to the core primitive.
#
# Design Philosophy:
# - Single `apply!` function handles ALL models (no model-specific methods)
# - Typed wrappers (BranchTrip, LoadTrip) provide ergonomic API
# - Uses model's `get_key_name()` trait for composite key lookup
# - Integrates with OrdinaryDiffEq.jl via callbacks
#
# Dependencies:
# - Foundation: get_key_name, SystemModel
# - DAEInfra: IndexOrder (sys.order.key_indices)
# - DiffEqCallbacks: PresetTimeCallback, CallbackSet
# ============================================================================

using InlineStrings: String3
using PowerfulCore: ParameterSlice, gather_parameters, scatter_parameters!, refresh_external!

export Disturbance, FieldModification, SetpointChange
export BranchTrip, BranchClose
export apply!, to_callback, build_event_callbacks

# ============================================================================
# Core Primitive
# ============================================================================

# Note: `Disturbance` abstract type and `to_callback` function stub are defined
# in Foundation (core_types.jl) to allow model files to define specialized
# event types with their own to_callback methods.

"""
    FieldModification <: Disturbance

Universal primitive for modifying any model field at a specific time.
All event types reduce to this.

# Fields
- `time::Float64`: When to apply the modification (seconds)
- `model::Symbol`: Target model type (`:Branch`, `:PQ`, `:PV`, `:GENCLS`, etc.)
- `key::Tuple`: Composite key identifying the instance (from `get_key_name()`)
- `field::Symbol`: Field to modify (`:status`, `:p`, `:q`, `:x`, `:Vref`, etc.)
- `operation::Symbol`: How to modify - `:set`, `:add`, or `:multiply`
- `value::Float64`: Value for the operation

# Operations
- `:set` - Replace field with value: `field = value`
- `:add` - Add to field: `field = field + value`
- `:multiply` - Scale field: `field = field * value`

# Example
```julia
# Trip branch from bus 1 to bus 5, circuit "1" at t=1.0
FieldModification(1.0, :Branch, (1, 5, "1"), :status, :set, 0.0)

# Increase load P by 10 MW at bus 3, load id "1" at t=5.0
FieldModification(5.0, :PQ, (3, "1"), :p, :add, 10.0)

# Double branch impedance at t=2.0
FieldModification(2.0, :Branch, (4, 5, "1"), :x, :multiply, 2.0)
```
"""
struct FieldModification <: Disturbance
    time::Float64
    model::Symbol
    key::Tuple
    field::Symbol
    operation::Symbol
    value::Float64

    function FieldModification(time, model, key, field, operation, value)
        operation in (:set, :add, :multiply) ||
            error("Invalid operation :$operation. Must be :set, :add, or :multiply")
        new(Float64(time), model, key, field, operation, Float64(value))
    end
end

"""
    SetpointChange <: Disturbance

Change of a setpoint of a model at a specific time: the controller's reference value, as
set at initialization (`vref0` of an exciter, `pref0` of a governor), as a staged test
steps it (e.g. an AVR voltage reference step test).

# Fields
- `time::Float64`: When to apply the change (seconds)
- `model::Symbol`: Model of the controller (`:SEXS`, `:TGOV1`, ...)
- `key::Tuple`: Key of the instance (`(bus, id)` for a dynamic model)
- `setpoint::Symbol`: Setpoint to change (`:vref0`, `:pref0`, ...)
- `operation::Symbol`: `:set`, `:add` or `:multiply`, as for `FieldModification`
- `value::Float64`: Value for the operation

# Example
```julia
# Step the voltage reference of the exciter at bus 2 up by 2% at t = 1 s, back at t = 6 s
events = [SetpointChange(1.0, :SEXS, (2, "1"), :vref0, :add, 0.02),
          SetpointChange(6.0, :SEXS, (2, "1"), :vref0, :add, -0.02)]
```
"""
struct SetpointChange <: Disturbance
    time::Float64
    model::Symbol
    key::Tuple
    setpoint::Symbol
    operation::Symbol
    value::Float64

    function SetpointChange(time, model, key, setpoint, operation, value)
        operation in (:set, :add, :multiply) ||
            error("Invalid operation :$operation. Must be :set, :add, or :multiply")
        new(Float64(time), model, key, setpoint, operation, Float64(value))
    end
end

# ============================================================================
# Core Apply Function - Handles ALL Models
# ============================================================================

"""
    apply!(sys::SystemModel, mod::FieldModification) -> Nothing

Apply a field modification to the system.

This single function handles ALL models - no model-specific methods needed.
Uses existing IndexOrder infrastructure for O(1) key lookup. Errors if the model
or instance does not exist (an event that targets nothing is a setup bug).

# Algorithm
1. Look up instance index via `sys.order.key_indices[model][key]`
2. Get current field value
3. Apply operation (set/add/multiply)
4. Store new value

# Example
```julia
# Trip a branch
mod = FieldModification(1.0, :Branch, (1, 5, "1"), :status, :set, 0.0)
apply!(sys, mod)  # Branch 1-5-"1" is now offline
```
"""
function apply!(sys::SystemModel, mod::FieldModification)
    # Step 1: Look up instance index
    idx = get_index_by_key(sys, mod.model, mod.key)

    # Step 2: Get field array
    model_array = sys.models[mod.model]
    field_array = getproperty(model_array, mod.field)

    # Step 3: Apply operation and store the new value
    field_array[idx] = modified_value(field_array[idx], mod)

    # The external parameters that are copies of the field (read by the kernels of other
    # models) are copied again
    refresh_external!(sys, mod.model, mod.field)

    return nothing
end

"""
    modified_value(current, mod::Union{FieldModification, SetpointChange})

Value of the modified field or setpoint: `mod.value` (`:set`), `current + mod.value`
(`:add`) or `current * mod.value` (`:multiply`).
"""
function modified_value(current, mod::Union{FieldModification, SetpointChange})
    mod.operation === :set && return mod.value
    mod.operation === :add && return current + mod.value
    return current * mod.value  # :multiply (checked by the constructor)
end

"""
    check_modifiable(sys::SystemModel, mods::AbstractVector{FieldModification})

Throw an `ArgumentError` unless `mods` take effect in a running simulation. Each
modification, in time order (equal times in the order given) and on top of the ones
before it, must meet:

- the model, the instance and the field exist;
- the field is not one that the model's `ModelRating` declares, its own (a generator's
  `mbase`) or the linked device's (a machine's `mbase` is the replaced generator's): the
  per-unit coefficients and the derived ratios `Sn/Sb` (`RatingRatio`) are computed from the
  rating when the system is built and would be stale;
- the mass coefficients of the model stay as they are (e.g. a generator's inertia `H`
  sets the mass `2H` of its rotor speed), because the mass matrix is built once at setup;
- the modified instance satisfies its model's `check_parameter_domains` (the domains, the
  constraints that read parameters only, the divisors that must not be zero): generated with
  the model where it has any, declared empty otherwise. Every registered model declares one.

The modifications are applied to copies of the model arrays, in the order they will be
applied at run time, so an `:add` or `:multiply` that follows a `:set` is checked at its
result. See [`check_events`](@ref), which a problem calls when it is built.
"""
function check_modifiable(sys::SystemModel, mods::AbstractVector{FieldModification})
    modified = Dict{Symbol, Any}()
    for mod in sort(mods; by = m -> m.time, alg = Base.Sort.MergeSort)
        haskey(sys.models, mod.model) || throw(ArgumentError(
            "FieldModification of $(mod.model): no such model in the system"))
        model_array = sys.models[mod.model]
        declared = filter(r -> r.field === mod.field, model_metadata(get_base_model_type(model_array)).ratings)
        isempty(declared) || throw(ArgumentError(
            "FieldModification of $(mod.model).$(mod.field): $(mod.field) is the " *
            "$(nameof(typeof(only(declared).rating_type))) rating of $(mod.model)" *
            (only(declared).source isa LinkedDevice ? ", taken from the device it is linked to" : "") *
            ", from which the per-unit coefficients and the derived ratios Sn/Sb are computed when the system is " *
            "built; an event cannot change it. Set it before building the TDSProblem instead."))
        hasproperty(model_array, mod.field) || throw(ArgumentError(
            "FieldModification of $(mod.model).$(mod.field): $(mod.model) has no field " *
            "$(mod.field). Fields: $(join(propertynames(model_array), ", "))"))
        idx = get_index_by_key(sys, mod.model, mod.key)

        # The model with the modifications so far, and this one
        current = get!(() -> copy(model_array), modified, mod.model)
        field_array = getproperty(current, mod.field)
        field_array[idx] = modified_value(field_array[idx], mod)

        # Every mass coefficient of the model is what it is in the system
        for var in model_metadata(get_base_model_type(model_array)).vars
            mass_coeff = get_mass_coeff(var)
            isnothing(mass_coeff) && continue
            eval_mass_coeff(mass_coeff, current) == eval_mass_coeff(mass_coeff, model_array) && continue
            throw(ArgumentError(
                "FieldModification of $(mod.model).$(mod.field) changes the mass coefficient of " *
                "$(mod.model).$(var.name), but the mass matrix is fixed during a simulation. " *
                "Set $(mod.field) before building the TDSProblem instead."))
        end

        instance = current[idx]
        check_parameter_domains(instance,
            "FieldModification at t = $(mod.time) of $(mod.model)$(mod.key): $(mod.field) $(mod.operation) " *
            "$(mod.value), giving $(getproperty(instance, mod.field))")
    end
    return nothing
end

"""
    check_events(sys::SystemModel, events)

[`check_modifiable`](@ref) for the `FieldModification`s among `events` (which may hold
vectors of events, as the wrappers return). The one check of the events of a problem,
made when it is built, so an event that cannot take effect, or that sets a parameter
of a model outside its domain, is refused before any integration.
"""
function check_events(sys::SystemModel, events)
    mods = FieldModification[e for event in events for e in (event isa Vector ? event : [event])
                             if e isa FieldModification]
    isempty(mods) || check_modifiable(sys, mods)
    return nothing
end

"""
    check_constraint_modifiable(adapter, prob)

Refuse, before any integration, a `FieldModification` among the events of `prob` of a parameter that a reduced constraint
of the adapter reads (`StateConstraintReduction` `reads`): the model parameters that the equations it replaces or combines
have in them, found by tracing their kernels, whether the constraint holds throughout (the current balance of a bus joined
to the rest by inductors reads the bus's shunt conductance `g`, and a bridge's `g_neutral`, which would make it read the
voltages) or in some modes (REGCA1's `Volim`, which the equation of `dq` reads in the mode `clamp`, is the constant of the
constraint `Volim - vm = 0` of a measured state). The reduction holds that constraint on the states of the DAE (a
measurement does not jump when a setpoint changes, and a parameter that makes a balance read a voltage ends the
reduction), so an event cannot change what it reads. The error names the parameter, the constraint and the mode. (The
adapter finds the constraints and what they read when it is built, at the start of the solve; a problem alone has not
built them.)
"""
function check_constraint_modifiable(adapter::TDSAdapter, prob::TDSProblem)
    red = adapter.reduction
    red === nothing && return nothing
    mods = FieldModification[e for event in prob.events for e in (event isa Vector ? event : [event])
                             if e isa FieldModification]
    for event in mods, j in eachindex(red.rows)
        isempty(red.reads[j]) && continue
        instance = get_index_by_key(prob.sys, event.model, event.key)
        (event.model, event.field, instance) in red.reads[j] || continue
        gated = !isempty(red.gates[j].instance)
        held = gated ? "while $(red.gates[j].label)" : "throughout"
        throw(ArgumentError(
            "FieldModification at t = $(event.time) of $(event.model)$(event.key): $(event.field) is read by a constraint on " *
            "the states that the reduction holds $held: $(red.descriptions[j]) is a function of the states alone, and " *
            "$(event.field) enters it. A measured state does not jump when a setpoint changes, and a parameter that changes what a " *
            "constraint reads ends the reduction of the equation, so an event cannot change it; set $(event.field) before building " *
            "the TDSProblem instead." * (gated ? " (The mode is entered at the located crossing of the quantity it limits.)" : "")))
    end
    return nothing
end

"""
    apply!(sys::SystemModel, status::AbstractDict, mod::FieldModification)

Apply `mod` to `sys` and, when it changes the model's status field, to the solver's
status dict `status` too (the active mask of residual evaluation), so a device
switched by an event takes effect in the running simulation.
"""
function apply!(sys::SystemModel, status::AbstractDict, mod::FieldModification)
    apply!(sys, mod)
    if mod.field === get_status_field(sys.models[mod.model])
        idx = get_index_by_key(sys, mod.model, mod.key)
        status[mod.model][idx] = !iszero(getproperty(sys.models[mod.model], mod.field)[idx])
    end
    return nothing
end

# ============================================================================
# OrdinaryDiffEq.jl Integration
# ============================================================================

"""
    _field_changed!(jump, mod::FieldModification, instance::Int)

Extension point: a model field changed by an event, for an augmentation carried through
the integration beside the DAE state (`jump`), e.g. the derivatives of the field with
respect to the parameters it moves with. The default, for `jump === nothing` (the plain
simulation), does nothing; an augmentation that tracks the field adds a method for its own
`jump` type.
"""
_field_changed!(::Nothing, mod::FieldModification, instance::Int) = nothing

"""
    _network_change!(apply!, integrator, network)

Apply an event (`apply!()`) and, with `network = (adapter, changed, massless)`, set
`changed[]` to whether it changed the equations without mass, which the reinitialization
solves: the algebraic residuals and the rows of `x` whose mass is zero (`massless`), at
the integrator's state before and after it (read by `EventReinit`), as a settle of the
automata checks them (`settle_instant!`). With `network = nothing`, apply only.
"""
function _network_change!(apply!, integrator, network)
    network === nothing && return apply!()
    adapter, changed, massless = network
    nd = adapter.n_diff
    n = nd + adapter.n_alg
    before, after = zeros(n), zeros(n)
    _evaluate_packed!(before, adapter, view(integrator.u, 1:n), integrator.t, nothing)
    apply!()
    _evaluate_packed!(after, adapter, view(integrator.u, 1:n), integrator.t, nothing)
    changed[] = view(before, (nd + 1):n) != view(after, (nd + 1):n) || any(r -> before[r] != after[r], massless)
    return nothing
end

# The rows of `x` whose mass is zero (the diagonal of the adapter's mass matrix)
_massless_states(adapter) = [r for r in 1:adapter.n_diff if iszero(adapter.mass_matrix[r, r])]

"""
    to_callback(mod::FieldModification, sys::SystemModel, status::AbstractDict; initializealg=nothing, jump=nothing,
                network=nothing, automata=nothing) -> Callback

Convert a FieldModification to an OrdinaryDiffEq.jl callback.

Creates a PresetTimeCallback that triggers exactly at the specified time.

# Arguments
- `mod`: The field modification to apply
- `sys`: System model (for apply!)
- `status`: The solver's status dict, updated when `mod` switches a device
- `jump`: Receives the modification (see `_field_changed!`); `nothing` for the plain
  simulation
- `network`: `(adapter, changed, massless)` to record whether the modification changes
  the equations without mass (`_network_change!`, for `EventReinit`), or `nothing`
- `automata`: the simulation's `AutomatonRuntime` (or `nothing`, or one with no governed
  state, `AutomatonRuntime.governed`): `(mod.model, mod.field, instance)` is recorded in
  its `touched` before the automata settle runs, so a governed state whose right side
  reads this field can be repaired (`_repair_writes`, `governed_states.jl`)

# Returns
PresetTimeCallback that applies the modification at `mod.time`
"""
function to_callback(mod::FieldModification, sys::SystemModel, status::AbstractDict; initializealg=nothing,
                     jump=nothing, network=nothing, automata=nothing)
    # Type barrier: wrap SystemModel in Ref{Any} to prevent closure type explosion
    # SystemModel embeds all model types, causing massive closure type signatures.
    # Cost: one dynamic dispatch per event (~nanoseconds) vs compilation savings.
    sys_ref = Ref{Any}(sys)

    # affect!(integrator) - apply the modification
    # Closure type is now simple (Ref{Any}, FieldModification) instead of embedding SystemModel
    instance = get_index_by_key(sys, mod.model, mod.key)
    affect! = let m = mod, s_ref = sys_ref, st = status, jump = jump, instance = instance, network = network,
                  rt = automata
        (integrator) -> begin
            s = s_ref[]::SystemModel
            _network_change!(() -> apply!(s, st, m), integrator, network)
            _field_changed!(jump, m, instance)
            rt !== nothing && rt.governed && push!(rt.touched, (m.model, m.field, instance))
            @debug "Event at t=$(m.time): $(m.model)[$(m.key)].$(m.field) $(m.operation) $(m.value)"
            # Reduce dt after discontinuity for numerical stability
            DiffEq.set_proposed_dt!(integrator, min(integrator.dt, 0.001))
        end
    end

    # Use provided initializealg, or default to NoInit if none specified
    init_alg = isnothing(initializealg) ? DiffEq.NoInit() : initializealg

    return DiffEq.PresetTimeCallback(
        [mod.time], erased(affect!, Val(1));
        save_positions=(true, true),
        initializealg=erase_init(init_alg)
    )
end

"""
    _setpoint_changed!(jump, change::SetpointChange, instance::Int)

Extension point: a setpoint changed by an event, for an augmentation carried through the
integration beside the DAE state (`jump`), e.g. the derivatives of the setpoint with
respect to the parameters it moves with. The default, for `jump === nothing` (the plain
simulation), does nothing; an augmentation that tracks the setpoint adds a method for its
own `jump` type.
"""
_setpoint_changed!(::Nothing, change::SetpointChange, instance::Int) = nothing

"""
    to_callback(change::SetpointChange, sys::SystemModel, adapter; initializealg=nothing, jump=nothing,
                network=nothing, automata=nothing) -> Callback

PresetTimeCallback that applies `change` at `change.time` to the setpoint in every
workspace of `adapter` (the Float64 one the simulation reads and the Dual ones its
Jacobians read), and passes it to `jump` (see `_setpoint_changed!`). Errors at setup if
the model, the instance or the setpoint does not exist.

`automata`: as `to_callback(::FieldModification, ...)`'s own (`(change.model, change.setpoint, instance)` is
recorded in its `touched` before the automata settle runs): a setpoint is a column of the model like any other
(`_value_param_reads` traces it the same way), so a governed state whose right side reads one is repaired the
same way a `FieldModification` of an ordinary field is, not only a `FieldModification`. Unlike that one,
simultaneous `SetpointChange`s (or a `SetpointChange` simultaneous with a `FieldModification`) are not grouped
into one "last may repair": each pushes to `touched` and lets the settle after it repair, which is only
correct where no relation reads fields from more than one of them (see the header of `governed_states.jl`).
"""
function to_callback(change::SetpointChange, sys::SystemModel, adapter; initializealg = nothing,
                     jump = nothing, network = nothing, automata = nothing)
    cache = adapter.adcache
    instance = get_index_by_key(sys, change.model, change.key)
    slice = ParameterSlice(cache, change.model, change.setpoint, [instance], 1; source = :setpoint)

    affect! = let cache = cache, slice = [slice], change = change, instance = instance, jump = jump, network = network,
                  rt = automata
        (integrator) -> begin
            current = only(gather_parameters(cache, slice))
            _network_change!(() -> scatter_parameters!(cache, slice, [modified_value(current, change)]), integrator,
                             network)
            _setpoint_changed!(jump, change, instance)
            rt !== nothing && rt.governed && push!(rt.touched, (change.model, change.setpoint, instance))
            @debug "Event at t=$(change.time): $(change.model)[$(change.key)].$(change.setpoint) " *
                   "$(change.operation) $(change.value)"
            # Reduce dt after discontinuity for numerical stability
            DiffEq.set_proposed_dt!(integrator, min(integrator.dt, 0.001))
        end
    end

    return DiffEq.PresetTimeCallback(
        [change.time], erased(affect!, Val(1));
        save_positions = (true, true),
        initializealg = isnothing(initializealg) ? DiffEq.NoInit() : erase_init(initializealg)
    )
end

"""
    build_event_callbacks(events::Vector, sys::SystemModel, status::AbstractDict; address_manager=nothing, adapter=nothing, reevaluate_fn=nothing, initializealg=nothing, jump=nothing, network_changed=nothing) -> CallbackSet

Convert all events to an OrdinaryDiffEq.jl CallbackSet.

# Arguments
- `events`: Vector of Disturbance objects (FieldModification, events that reset states, etc.)
- `sys`: System model
- `status`: The solver's status dict (e.g. `adapter.status`); events that switch
  devices update it
- `address_manager`: AddressManager for events that need to modify state variables (optional)
- `adapter`: TDSAdapter for events that need to re-evaluate residuals after state modification (optional)
- `reevaluate_fn`: Closure `(integrator) -> ...` that re-evaluates residuals with current state (optional).
  Called after state modification to update adapter buffers before DAE reinitialization.
- `initializealg`: DAE initialization algorithm to use after events (default: NoInit)
  - Use `PowerFlowInit()` for robust power system reinitialization
  - Use `BrownFullBasicInit()` for generic DAE initialization
  - Use `NoInit()` to skip reinitialization (may cause solver failure)
- `jump`: What else a field or setpoint change changes, for an integration that carries
  more than the DAE state (see `_field_changed!`, `_setpoint_changed!`); `nothing` for
  the plain simulation
- `network_changed`: A `Ref{Bool}` the field and setpoint changes set to whether they
  change the equations without mass (the algebraic residuals and the rows of `x` whose
  mass is zero, with `adapter`), read by the reinitialization `EventReinit`; `nothing`
  not to record it
- `automata`: the simulation's `AutomatonRuntime` (or `nothing`); threaded to a
  `FieldModification`'s callback so it can record the field it touches (see
  `to_callback(::FieldModification, ...)`)

# Returns
CallbackSet containing all event callbacks, or nothing if no events

# Note
Some event types (those that reset states, e.g. a fault clearing) require `address_manager` to reset
state variables. These events define their own `to_callback` method that takes
`(event, sys, am, adapter, reevaluate_fn)` and is dispatched automatically.
"""
function build_event_callbacks(events::Vector, sys::SystemModel, status::AbstractDict;
                               address_manager=nothing, adapter=nothing, reevaluate_fn=nothing, initializealg=nothing,
                               jump=nothing, network_changed=nothing, automata=nothing)
    isempty(events) && return nothing
    network = network_changed === nothing ? nothing : (adapter, network_changed, _massless_states(adapter))

    # Flatten if a wrapper returned a Vector (e.g., LoadChange returns multiple FieldModifications)
    flattened = [e for event in events for e in (event isa Vector ? event : [event])]
    # The position of the last FieldModification at each exact time: only it may repair a governed state
    # (`AutomatonSettle` `repair`), from the fields all of them at that time have written by then; the others use
    # a copy of `initializealg` with repair turned off (`_with_repair`), everything else (the DAE reinitialization,
    # the ordinary transitions) unchanged
    last_at_time = Dict{Float64, Int}()
    for (p, e) in enumerate(flattened)
        e isa FieldModification && (last_at_time[e.time] = p)
    end

    callbacks = []

    for (p, e) in enumerate(flattened)
        if e isa FieldModification
            # Standard field modification - no state reset needed (checked when the problem was built)
            repair_here = last_at_time[e.time] == p
            alg = (!repair_here && initializealg isa AutomatonSettle) ? _with_repair(initializealg, false) :
                                                                        initializealg
            push!(callbacks, to_callback(e, sys, status; initializealg=alg, jump, network, automata))
        elseif e isa SetpointChange
            isnothing(adapter) && error("SetpointChange requires adapter but none provided")
            push!(callbacks, to_callback(e, sys, adapter; initializealg, jump, network, automata))
        elseif e isa Disturbance
            # Specialized Disturbance subtype (e.g., one that resets states)
            # These define to_callback(event, sys, am, adapter, reevaluate_fn) for state manipulation
            if isnothing(address_manager)
                error("Event type $(typeof(e)) requires address_manager but none provided")
            end
            if isnothing(adapter)
                error("Event type $(typeof(e)) requires adapter but none provided")
            end
            push!(callbacks, to_callback(e, sys, address_manager, adapter, reevaluate_fn))
        else
            error("Unknown event type: $(typeof(e)). Events must be FieldModification or Disturbance subtypes.")
        end
    end

    return isempty(callbacks) ? nothing : DiffEq.CallbackSet(callbacks...)
end

# ============================================================================
# Ergonomic Wrappers - Typed Constructors
# ============================================================================

# These use get_key_name() to document the expected key format.
# Example: Branch has get_key_name() = [:i, :j, :ckt], so BranchTrip takes (from, to, circuit)

"""
    BranchTrip(time, from_bus, to_bus, circuit) -> FieldModification

Trip (disconnect) a branch at the specified time.

# Arguments
- `time`: Event time in seconds
- `from_bus`: "From" bus number (Branch.i)
- `to_bus`: "To" bus number (Branch.j)
- `circuit`: Circuit identifier (Branch.ckt), typically "1" or " 1"

# Example
```julia
events = [
    BranchTrip(1.0, 1, 5, "1"),      # Trip line 1-5 circuit "1" at t=1.0s
    BranchTrip(1.5, 4, 5, " 1"),     # Trip line 4-5 at t=1.5s
]
prob = TDSProblem(setup_tds(pf_sol, dyn_models); events=events, tspan=(0.0, 10.0))
```

# Note
Branch key is `(i, j, ckt)` per `get_key_name(Branch) = [:i, :j, :ckt]`.
The circuit ID must match exactly (including spaces).
"""
function BranchTrip(time::Real, from_bus::Integer, to_bus::Integer, circuit::AbstractString)
    return FieldModification(time, :Branch, (Int32(from_bus), Int32(to_bus), String3(circuit)), :status, :set, 0.0)
end

# Convenience: allow passing circuit as Int (converts to "1", "2", etc.)
BranchTrip(time::Real, from_bus::Integer, to_bus::Integer, circuit::Integer) =
    BranchTrip(time, from_bus, to_bus, string(circuit))

"""
    BranchClose(time, from_bus, to_bus, circuit) -> FieldModification

Close (reconnect) a previously tripped branch.

# Example
```julia
events = [
    BranchTrip(1.0, 1, 5, "1"),    # Trip at t=1.0
    BranchClose(2.0, 1, 5, "1"),   # Reconnect at t=2.0
]
```
"""
function BranchClose(time::Real, from_bus::Integer, to_bus::Integer, circuit::AbstractString)
    return FieldModification(time, :Branch, (Int32(from_bus), Int32(to_bus), String3(circuit)), :status, :set, 1.0)
end

BranchClose(time::Real, from_bus::Integer, to_bus::Integer, circuit::Integer) =
    BranchClose(time, from_bus, to_bus, string(circuit))

# ============================================================================
# Fault Event Wrappers
# ============================================================================

"""
    FaultApply(time, bus, id="F1") -> FieldModification

Apply (turn on) a fault at the specified time.

Note: Sets status to 1 (will be converted to true for Bool field).

# Example
```julia
events = [FaultApply(1.0, 10)]  # Apply fault at bus 10 at t=1.0s
```
"""
function FaultApply(time::Real, bus::Integer, id::AbstractString="F1")
    return FieldModification(time, :Fault, (Int32(bus), String3(id)), :status, :set, 1)
end

"""
    FaultClear(time, bus, id="F1") -> FieldModification

Clear (turn off) a fault at the specified time.

Note: Sets status to 0 (will be converted to false for Bool field).

# Example
```julia
events = [
    FaultApply(1.0, 10),   # Apply fault at t=1.0s
    FaultClear(1.1, 10),   # Clear fault at t=1.1s
]
```
"""
function FaultClear(time::Real, bus::Integer, id::AbstractString="F1")
    return FieldModification(time, :Fault, (Int32(bus), String3(id)), :status, :set, 0)
end


# ============================================================================
# High-Level Fault Event
# ============================================================================

"""
    BusFault(time_on, time_off, bus, id="F1"; xf=0.003, rf=0.0, fault_type=0, rg=0.0, xg=0.0)

Fault at `bus` from `time_on` to `time_off`. `TDSProblem` adds the fault model of the
bus model holding `bus` (see `bus_fault`) and expands the event into its apply/clear
primitives: at a phasor bus, a `Fault` (balanced three-phase faults only,
`fault_type = 0`); a bus modeled in three phases takes every `fault_type`.

# Arguments
- `time_on::Real`: Fault application time (seconds)
- `time_off::Real`: Fault clearing time (seconds)
- `bus::Integer`: Bus number where fault occurs
- `id::AbstractString`: Fault identifier (default: "F1")
- `xf::Real`: Fault reactance in per-unit (default: 0.003)
- `rf::Real`: Fault resistance in per-unit (default: 0.0)
- `fault_type::Integer`: Fault type (default: 0, three-phase-to-ground); 1-3
  single-phase-to-ground on phase A/B/C, 4-6 line-to-line AB/BC/CA, 7-9 double
  line-to-ground AB/BC/CA, 10 three-phase without ground
- `rg::Real`, `xg::Real`: Fault-to-ground resistance and reactance (per-unit) of
  double line-to-ground faults (default: 0.0)

A second `BusFault` with the same bus and id reuses the fault model; its impedance
and type must match.

# Example
```julia
events = [
    BusFault(0.5, 0.6, 10, "F1"; xf=0.003, rf=0.0),
    BranchTrip(1.0, 1, 5, "1"),
]
tds_prob = TDSProblem(init; tspan=(0.0, 2.0), events=events)

# Single-phase-to-ground fault on phase A (at a bus modeled in three phases)
events = [BusFault(0.5, 0.6, 10; xf=0.1, fault_type=1)]
```
"""
struct BusFault <: Disturbance
    time_on::Float64
    time_off::Float64
    bus::Int32
    id::String3
    xf::Float64
    rf::Float64
    fault_type::Int8
    rg::Float64
    xg::Float64
end

function BusFault(time_on::Real, time_off::Real, bus::Integer, id::AbstractString="F1";
                  xf::Real=0.003, rf::Real=0.0, fault_type::Integer=0, rg::Real=0.0, xg::Real=0.0)
    time_on >= 0 || error("time_on must be non-negative")
    time_off > time_on || error("time_off must be greater than time_on")
    xf >= 0 || error("xf must be non-negative")
    rf >= 0 || error("rf must be non-negative")
    rg >= 0 && xg >= 0 || error("rg and xg must be non-negative")
    0 <= fault_type <= 10 || error("fault_type must be in 0:10, got $fault_type")
    return BusFault(Float64(time_on), Float64(time_off), Int32(bus), String3(id),
                    Float64(xf), Float64(rf), Int8(fault_type), Float64(rg), Float64(xg))
end

export BusFault

"""
    accepts_bus_fault(::Type{M}) -> Bool

Whether a `BusFault` attaches to the buses of model `M`, keyed by bus number. Such a
bus model also defines `bus_fault(::Val{name}, event)` for its name in the system.
"""
accepts_bus_fault(::Type{<:AbstractModel}) = false
accepts_bus_fault(::Type{<:Bus}) = true

"""
    bus_fault(::Val{vertex}, event::BusFault) -> (model, fault, apply, clear)

Fault model name, instance and apply/clear events for `event` at a bus of model
`vertex` (a model that `accepts_bus_fault`). Here the phasor `Fault` at a `:Bus`.
"""
function bus_fault(::Val{:Bus}, event::BusFault)
    event.fault_type == 0 && iszero(event.rg) && iszero(event.xg) ||
        error("BusFault at phasor bus $(event.bus): a phasor simulation represents balanced " *
              "three-phase faults only (fault_type = 0, no rg/xg), got fault_type = $(event.fault_type), " *
              "rg = $(event.rg), xg = $(event.xg). Unbalanced faults need the bus modeled in three phases.")
    fault = create_fault(event.bus; id = String(event.id), xf = event.xf, rf = event.rf, status = false)
    return (:Fault, fault, FaultApply(event.time_on, event.bus, event.id),
            FaultClear(event.time_off, event.bus, event.id))
end

# Name of the bus model holding bus number `bus`, among the models that accept bus faults; refuses a bus that none
# of them holds, naming the models that accept bus faults and how a bus model declares that it does
function _fault_bus_model(sys::SystemModel, bus)
    accepting = Symbol[]
    for (model_name, models) in pairs(sys.models)
        isempty(models) && continue
        accepts_bus_fault(get_base_model_type(models)) || continue
        push!(accepting, model_name)
        find_index(sys, model_name, (bus,)) === nothing || return model_name
    end
    throw(ArgumentError(
        "BusFault at bus $bus: no model of the system that accepts bus faults has bus $bus (those that do: " *
        (isempty(accepting) ? "none" : join(accepting, ", ")) * "). A bus model accepts them by declaring " *
        "PowerfulTDS.accepts_bus_fault(::Type{<:M}) = true and defining PowerfulTDS.bus_fault(::Val{name}, event) " *
        "for its name in the system"))
end

"""Fault reactance (pu) below which a fault warns unless loads are constant impedance."""
const LOW_FAULT_REACTANCE = 0.01

function _warn_low_fault_impedance(event::BusFault, config::SimConfig)
    event.xf < LOW_FAULT_REACTANCE || return nothing
    tds = get_solver_config(config, :TDS, TDSConfig)
    tds.convert_pq_to_zip && tds.p_zip == (1.0, 0.0, 0.0) && tds.q_zip == (1.0, 0.0, 0.0) && return nothing
    @warn """BusFault at bus $(event.bus) has xf = $(event.xf) < $(LOW_FAULT_REACTANCE) pu and loads are not \
    100% constant impedance, which may fail to converge during the fault. To convert them, pass
      SimConfig(solver = (TDS = TDSConfig(convert_pq_to_zip=true, p_zip=(1.0,0.0,0.0), q_zip=(1.0,0.0,0.0)),))
    as `config` to the setup of the simulation (e.g. `setup_tds`).""" maxlog=1
    return nothing
end

# Same fault parameters (every field but the status)
_same_fault(a, b) = all(f -> f === :status || getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

"""
    expand_bus_faults!(sys::SystemModel, events::Vector, config::SimConfig) -> Vector

Replace each `BusFault` in `events` by its apply/clear primitives and add its fault
model to `sys` (see `bus_fault`), once per (bus, id). Other events pass through.
Errors if a fault with the same (bus, id) already exists with different parameters.
Called by the `TDSProblem` constructors and the Y-bus event processing.
"""
function expand_bus_faults!(sys::SystemModel, events::Vector, config::SimConfig)
    expanded = Any[]
    new_faults = Dict{Symbol, Vector{Any}}()

    for event in events
        if !(event isa BusFault)
            push!(expanded, event)
            continue
        end
        model, fault, apply_event, clear_event = bus_fault(Val(_fault_bus_model(sys, event.bus)), event)
        key = (event.bus, event.id)

        pending = get!(new_faults, model, Any[])
        k = findfirst(f -> (f.i, f.id) == key, pending)
        idx = haskey(sys.models, model) ? find_index(sys, model, key) : nothing
        existing = !isnothing(k) ? pending[k] : !isnothing(idx) ? sys.models[model][idx] : nothing
        if isnothing(existing)
            push!(pending, fault)
        else
            _same_fault(existing, fault) ||
                error("BusFault at bus $(event.bus) id \"$(event.id)\" differs from the existing " *
                      "$model with that key; use a different id for a different fault.")
        end

        _warn_low_fault_impedance(event, config)
        push!(expanded, apply_event, clear_event)
    end

    # Add all new faults at once (one re-index)
    added = NamedTuple(model => StructArray([f for f in faults]) for (model, faults) in new_faults if !isempty(faults))
    isempty(added) || add_instances!(sys, added)

    return expanded
end

export expand_bus_faults!
