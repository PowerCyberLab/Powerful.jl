# TDS Event System

Disturbances in time-domain simulation are specified in `PowerfulTDS/src/solvers/events.jl`.
Branch trips and reclosures, load changes, and parameter and setpoint changes are all
expressed as a `FieldModification`, which a single `apply!` function carries out for every
model. Adding a model therefore requires no new event handler; the lookup uses the
`IndexOrder` of the system. Events that must also reset solver states (for example, clearing
a fault that leaves state behind) are separate `Disturbance` subtypes with their own
`to_callback` method.

## FieldModification

```julia
struct FieldModification <: Disturbance
    time::Float64      # When (seconds)
    model::Symbol      # What model (:Branch, :PQ, :GENCLS, ...)
    key::Tuple         # Which instance (composite key from get_key_name())
    field::Symbol      # Which field (:status, :pl, :x, ...)
    operation::Symbol  # How (:set, :add, :multiply)
    value::Float64     # Value for the operation
end
```

The constructor validates the operation and converts `time` and `value` to `Float64`. The
three operations are

| Operation | Effect |
|-----------|--------|
| `:set` | `field = value` |
| `:add` | `field = field + value` |
| `:multiply` | `field = field * value` |

Model fields hold per-unit values (see [Per-Unit & Multi-Island](miscellaneous.md)), so
values for power, impedance and similar fields are given in per-unit on the system base.

## Keys

Each model defines its composite key through the `get_key_name()` trait:

```julia
get_key_name(::Type{Bus}) = :i                    # (bus_number,)
get_key_name(::Type{T}) where T<:PQ = [:i, :id]   # (bus_number, load_id)
get_key_name(::Type{Branch}) = [:i, :j, :ckt]     # (from_bus, to_bus, circuit)
```

The key of an event has the same structure as the key of the model:

```julia
# Branch 1-5 circuit "1"
FieldModification(1.0, :Branch, (1, 5, "1"), :status, :set, 0.0)

# PQ load at bus 3, id "1": increase constant-power P by 0.1 pu (10 MW on a 100 MVA base)
FieldModification(5.0, :PQ, (3, "1"), :pl, :add, 0.1)
```

Keys are converted to the key field types (`Int32`, `String3`) by `normalize_key` when they
are looked up.

## The apply! Function

```julia
function apply!(sys::SystemModel, mod::FieldModification)
    # Step 1: Look up the instance (errors if it does not exist)
    idx = get_index_by_key(sys, mod.model, mod.key)

    # Step 2: Get the field array
    field_array = getproperty(sys.models[mod.model], mod.field)

    # Step 3: Apply the operation
    current = field_array[idx]
    new_val = if mod.operation == :set
        mod.value
    elseif mod.operation == :add
        current + mod.value
    elseif mod.operation == :multiply
        current * mod.value
    end

    field_array[idx] = new_val
    return nothing
end
```

A second method, `apply!(sys, status, mod)`, also updates the status dict of the solver when
the modified field is the status field of the model. The status dict holds the in-service
masks used by residual evaluation, so a device switched by an event takes effect in the
running simulation.

The function uses `sys.order.key_indices` for the lookup and `getproperty` for the field
access, so it contains no model-specific code. The residual workspaces hold references to the
same model arrays, and the change is visible at the next residual evaluation.

## Checks When the Problem Is Built

`TDSProblem` and `YBusTDSProblem` check every `FieldModification` of their events when they
are built (`check_events`). The events are checked in time order, each on top of the ones
before it. A check fails if

- the model, the instance or the field does not exist;
- the modification changes a mass coefficient (the mass matrix is fixed during a simulation);
- the modified instance violates `check_parameter_domains` of its model. This function covers
  the domains of the parameters, the constraints that read parameters only, and the divisors
  that must be nonzero. A model without the function is not checked.

The error names the model, the instance, the field, the operation, the value and the failed
check. Where the automata of a model govern a differential state, the admission of the
governed states is checked again when the event applies, for the models whose parameters it
changes.

## Typed Wrappers

Typed wrappers are constructors that return a `FieldModification`:

```julia
BranchTrip(t, from, to, ckt)   # FieldModification(t, :Branch, (from, to, ckt), :status, :set, 0.0)
BranchClose(t, from, to, ckt)  # FieldModification(t, :Branch, (from, to, ckt), :status, :set, 1.0)
FaultApply(t, bus, id="F1")    # FieldModification(t, :Fault, (bus, id), :status, :set, 1)
FaultClear(t, bus, id="F1")    # FieldModification(t, :Fault, (bus, id), :status, :set, 0)
```

A new wrapper is a one-line function, for example:

```julia
LoadTrip(t, bus, id) = FieldModification(t, :PQ, (bus, id), :status, :set, 0.0)
```

A wrapper may also return a vector of `FieldModification`s; `build_event_callbacks`
flattens vectors.

## DifferentialEquations.jl Integration

The path from a user-level event to the integrator is:

```
User API                    Core Primitive              Infrastructure
───────────────────────────────────────────────────────────────────────
BranchTrip(1.0, 1, 5, "1")
         │
         ▼
    FieldModification(...)  ──►  apply!(sys, status, mod)  ──►  get_index_by_key
         │                              │                             │
         │                              ▼                             ▼
         │                     status dict updated           getproperty(model, field)
         │                     (status fields only)                   │
         ▼                                                            ▼
    to_callback(...)                                       field_array[idx] = new
         │
         ▼
    PresetTimeCallback (initializealg = PowerFlowInit by default)
         │
         ▼
    DiffEq.solve(prob; callback = ...)
```

Each `FieldModification` becomes a `PresetTimeCallback`:

```julia
function to_callback(mod::FieldModification, sys::SystemModel, status::AbstractDict;
                     initializealg = nothing)
    sys_ref = Ref{Any}(sys)   # keeps the closure type small
    affect! = let m = mod, s_ref = sys_ref, st = status
        (integrator) -> begin
            apply!(s_ref[]::SystemModel, st, m)
            # Reduce the step after the discontinuity
            DiffEq.set_proposed_dt!(integrator, min(integrator.dt, 0.001))
        end
    end
    init_alg = isnothing(initializealg) ? DiffEq.NoInit() : initializealg
    return DiffEq.PresetTimeCallback([mod.time], affect!;
                                     save_positions = (true, true),
                                     initializealg = init_alg)
end
```

The `initializealg` of the callback reinitializes the algebraic variables after the event (see
"DAE Reinitialization: PowerFlowInit" below). Every callback that changes the trajectory
discontinuously saves the state before and after the event. This applies to these events and
to the located events of the automata. The solution then holds both states at the event time,
and the interpolation of the following step starts from the state after the event. An
interpolant that reads the first point of the step, as those of Rosenbrock and SDIRK methods
do, would otherwise start from the state before the event.

The TDS `solve()` builds the callbacks:

```julia
# Inside solve(prob::TDSProblem, ...)
if !isempty(prob.events)
    init_alg = _bind_adapter(initializealg, adapter)
    event_callbacks = build_event_callbacks(prob.events, prob.sys, adapter.status;
        address_manager = prob.address_manager, adapter = adapter,
        reevaluate_fn = reevaluate_fn, initializealg = init_alg)
end
```

Here `adapter.status` is the status dict of the solver. Other `Disturbance` subtypes (such as
`Fault3PhClearEvent`) are dispatched to
`to_callback(event, sys, address_manager, adapter, reevaluate_fn; initializealg)`, which can
modify state variables and re-evaluate residuals.

## User API

```julia
using Powerful, PowerfulCases

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
sol_pf = solve(ACPowerFlowProblem(sys), Newton())
dyr_data = parse_dynamics(case.dyr, sys)

init = setup_tds(sol_pf, dyr_data)
prob = TDSProblem(init;
    tspan = (0.0, 10.0),
    events = [
        BranchTrip(1.0, 1, 5, "1"),    # Trip at t = 1.0 s
        BranchClose(2.0, 1, 5, "1"),   # Reconnect at t = 2.0 s
    ],
)
sol = solve(prob)
```

An event without a wrapper is written directly as a `FieldModification`:

```julia
events = [
    # Set the damping of the generator at bus 1, id "1", to 2.0 at t = 3.0 s
    FieldModification(3.0, :GENCLS, (1, "1"), :D, :set, 2.0),
]
```

## Adding Support for New Models

A new model works with `FieldModification` without any change to the event code if it
satisfies three conditions:

1. It defines the `get_key_name()` trait.
2. Its fields are accessible through `getproperty`.
3. It is in the system, so that it has an entry in `sys.order.key_indices`.

A model that events should be able to switch also declares its `status_field` in the
metadata.

A modified field takes effect only if the residual kernels read it at every evaluation, so
a model must compute values derived from its parameters inside the kernel (as the fault does
with its admittance) and must not store them. The mass matrix is built once at setup, and
`build_event_callbacks` rejects, through `check_modifiable`, a `FieldModification` that
would change a mass coefficient (such as the `H` of a generator). Unknown models and fields
are rejected in the same way.

## Fault Model and Events

A bus fault is represented by a fault model added to the system (defined in
`PowerfulTDS/src/models/fault.jl`). The fault is a shunt admittance to ground at the faulted
bus that draws power only while its `status` is `true`.

```julia
@kwdef mutable struct Fault{Tv} <: AbstractModel
    i::Int32           # Bus number where the fault occurs
    id::String3        # Fault identifier (e.g. "F1")
    rf::Tv             # Fault resistance (pu)
    xf::Tv             # Fault reactance (pu)
    status::Bool       # false = cleared, true = applied
end
```

The fault admittance is `Yf = 1/Zf = Gf + jBf` with `Gf = Rf/(Rf² + Xf²)` and
`Bf = -Xf/(Rf² + Xf²)`. While the fault is applied it consumes `P = Gf·V²` and
`Q = -Bf·V²`. For an inductive fault (`rf ≈ 0`) the susceptance is `Bf ≈ -1/xf`, a large
negative value, and the absorbed reactive power depresses the bus voltage. The kernels
compute `gf` and `bf` from `rf` and `xf` at every evaluation (`fault_admittance`), so an
event that changes the fault impedance takes effect.

### BusFault Expansion

`BusFault(time_on, time_off, bus, id="F1"; xf=0.003, rf=0.0, fault_type=0, rg=0.0, xg=0.0)`
is the user-facing fault event. The `TDSProblem` constructors and the Y-bus event processing
call `expand_bus_faults!(sys, events, config)`, which processes each `BusFault` in four steps:

1. The fault model is chosen through the `accepts_bus_fault`/`bus_fault` extension point
   (`solvers/events.jl`), by the model type of the bus. The implementation for `Bus` is the
   phasor `Fault`, which rejects unbalanced `fault_type`s.
2. The model is added with `add_instances!` once per (bus, id). This updates the indices, the
   external parameters and `foreign_refs`. A repeated key must have the same parameters.
3. The event is replaced by the internal primitives that the chosen fault model declares
   (`FieldModification`s on the status of the phasor `Fault`, or the own events of another
   model).
4. A warning is issued if `xf < LOW_FAULT_REACTANCE` (0.01 pu) and the loads are not 100%
   constant impedance.

## DAE Reinitialization: PowerFlowInit

After topology changes (faults, trips), the algebraic variables must be reinitialized to
satisfy `g(x, y) = 0`. `PowerFlowInit` is defined in `PowerfulTDS/src/solvers/reinit.jl`.

### Current Mismatch Rows

Every term of a bus power balance carries the bus voltage magnitude
(`P_i = V_i Σ_j V_j (G_ij cos θ_ij + B_ij sin θ_ij) + ...`, and load and shunt powers vanish
with `V_i`), so `V_i = 0` satisfies both rows of bus `i` for any current flowing into it.
During a severe fault the voltage of the faulted bus is close to zero. A Newton iteration
started there after clearing converges to `V_i = 0`, the fault stays in place, and the
integrator stalls on the singular Jacobian (`∂P_i/∂θ = 0` at `V_i = 0`).
`reinitialize_algebraic!` therefore solves the bus rows as current mismatch (`P_i / V_i`,
`Q_i / V_i`, see `BusPowerRows`). These rows have the same roots with `V_i ≠ 0` and none at
`V_i = 0`, and the `abstol` of the reinitialization bounds the current mismatch of these
rows.

### Algorithm

```julia
struct PowerFlowInit <: DAEInitializationAlgorithm
    abstol::Union{Nothing, Float64}               # nothing: follow the integration
    maxiters::Int
    initial_guess::Symbol                         # :flat_start, :pre_fault, or :current
    adapter::Union{Nothing, TDSAdapter}           # Bound at solve() time
    pre_fault_v::Union{Nothing, Vector{Float64}}  # For :pre_fault mode
    verbose::Bool
end
```

The constructor is `PowerFlowInit(; abstol, maxiters, initial_guess)`. Without an `abstol`,
the reinitialization bound to a simulation stops at `REINIT_TOL_FRACTION` (1/100) of the
absolute tolerance of the integration. This is the fraction at which the Newton iterations of
the integrator stop: 1e-8 at the default `abstol = 1e-6` and 1e-12 at `1e-10`. The state
after an event is then consistent to the tolerance of the integration, as the check of
DiffEq after a callback (`CheckInit`) requires. An explicit `abstol` is kept.

An event reinitializes the network only where it changes the algebraic equations, that is,
where their residuals at the state differ before and after the event applies (`EventReinit`).
The automata follow the same rule for their transitions. A setpoint step of a controller
leaves the equations unchanged and the state as consistent as the integration keeps it. A
Newton solve there would move the state by less than the tolerance of the integration,
differently for nearby parameters, and would roughen the dependence of the trajectory on
them.

A reinitialization takes four steps:

1. The initial guess selected by `initial_guess` is applied (for a network event), or the
   iteration starts from the state of the integrator (for a transition of the automata, see
   below).
2. The Newton iteration solves `g(x, y) = 0` for the algebraic `y` with the differential `x`
   fixed, using the bus power balance rows as current mismatch. The iteration runs at the time
   `t` of the integrator, because a model whose residuals depend directly on time (a recording
   interpolated at `t`) needs it. Every evaluation states its time, and the ADCache has no
   default.
3. The result is checked against the DAE as the integration evaluates it at `t`. The rows
   without mass must be below the tolerance of the solve; otherwise the reinitialization
   fails like a solve that did not converge.
4. The state of the integrator is updated.

### Initial Guess Options

| Option | Voltage | Angle | Use |
|--------|---------|-------|-----|
| `:flat_start` | V = 1.0 | Current θ | Default. Converged after the short faults tested, and failed after a 2 s fault on ieee14 where `:current` converged; the outcome depends on the case. |
| `:pre_fault` | Pre-fault V | Current θ | Temporary faults after which the system returns to the pre-fault state. Unsuitable when the angles have moved far from their pre-fault values (after a long fault or a pole slip). |
| `:current` | Current V | Current θ | Small perturbations, load changes and setpoint modifications. Converged in the long-fault cases tested, and may still converge to the wrong solution from low voltages. |

The default comes from the `TDSConfig` of the problem: `:flat_start` when `flat_start = true`
(the default), and `:pre_fault` otherwise. The guess is selected with the `initializealg`
keyword of `solve`:

```julia
# Default (from TDSConfig; flat start unless flat_start = false)
sol = solve(prob)

# Pre-fault voltages (fault-on/fault-off scenarios)
sol = solve(prob; initializealg = PowerFlowInit(initial_guess = :pre_fault))

# Current voltages
sol = solve(prob; initializealg = PowerFlowInit(initial_guess = :current))
```

The guesses are applied as follows, and the pre-fault voltages are captured in `solve()`
before the simulation starts:

```julia
function _apply_initial_guess!(u, adapter, alg)
    if alg.initial_guess == :flat_start
        _apply_flat_start!(u, adapter)                      # V = 1.0, θ kept
    elseif alg.initial_guess == :pre_fault
        _apply_pre_fault_v!(u, adapter, alg.pre_fault_v)    # Pre-fault V, θ kept
    elseif alg.initial_guess == :current
        # Use the current values
    end
end
```

### Guess per Instant

The configured guess belongs to the network events (faults, trips, `FieldModification`).
These change the network itself, and the reinitialization that follows them (`EventReinit`,
run through the `initializealg` of DiffEq after the callback of the event) applies the guess.

The reinitializations at transitions of the automata start from the state of the integrator
(the guess `:current`, selected with `at_current_state`). This includes the settle right after
an event. Before a transition the state is consistent, and a transition changes the equations
at a crossing of its surface, so the state of the integrator is the nearest guess. A
transition reinitialized from a flat start inside a fault (bus voltages reset to 1 while the
fault holds them near 0.3) did not converge, and the current state converged at once.
Transitions reach the reinitialization through one function in
`PowerfulTDS/src/solvers/automata.jl`, which is the point where the two paths split.

The rounds at the start of a simulation also begin from the current state. This applies where
the equations of the initial modes are violated, at the point handed over by the
initialization, by more than the `abstol` of the `TDSConfig` of the problem (`StartResolve`),
and between the rounds of the start as at any other instant.

### Failure of the Reinitialization

`reinitialize_network!` returns whether the Newton iterations converged and whether the result
satisfies the DAE as the integration evaluates it at the instant. If either fails (a numerical
failure), a warning reports the time, the iterations, the residual and the worst rows (model,
variable, instance). The caller then ends the integration. All reinitializations share this
mechanism: `PowerFlowInit`, `EventReinit`, the settle of the automata (`AutomatonSettle`), the
`SensitivityInit` of the sensitivities, and the `HostReinit` of the Y-bus DAE. The solution has
`retcode == :ConvergenceFailure` and holds the trajectory up to the failed instant. The
wrappers stop at the first failure.

Errors in the models or the code throw, as does an inconsistent initial state (the check of
DiffEq at `t0`). An inconsistent start indicates a model or initialization bug.

Some reinitializations use `BrownFullBasicInit` of DifferentialEquations.jl. This applies when
the states to re-solve are zero-mass states that `PowerFlowInit` would otherwise hold fixed.

## Files

- `PowerfulTDS/src/solvers/events.jl`: `FieldModification`, `apply!`, `to_callback`, wrappers, `BusFault`
- `PowerfulTDS/src/solvers/reinit.jl`: `PowerFlowInit`, flat start, pre-fault capture
- `PowerfulTDS/src/solvers/solve.jl`: integration with `solve()`
- `PowerfulTDS/src/models/fault.jl`: phasor `Fault` model
- `Powerful/test/tds/`: tests
