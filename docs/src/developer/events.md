# TDS Event System

Disturbance specification for time-domain simulation
(`PowerfulTDS/src/solvers/events.jl`).

## Design Philosophy

**One primitive handles most events:**
- Branch trips and reclosures, load changes, parameter and setpoint changes
- All reduce to `FieldModification`
- A single `apply!` function for all models

**No model-specific methods:**
- Adding new models doesn't require new event handlers
- Uses the existing `IndexOrder` infrastructure

Events that must also reset solver states (e.g. clearing a fault that leaves state
behind) are separate `Disturbance` subtypes with their own `to_callback` method.

## Core Primitive: FieldModification

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

The constructor validates the operation and converts `time` and `value` to `Float64`.

**Operations:**
- `:set` → `field = value`
- `:add` → `field = field + value`
- `:multiply` → `field = field * value`

Model fields hold per-unit values (see [Per-Unit & Multi-Island](miscellaneous.md)), so
values for power, impedance and similar fields are in per-unit on the system base.

## How Keys Work

**Each model defines its composite key via the `get_key_name()` trait:**

```julia
get_key_name(::Type{Bus}) = :i                    # (bus_number,)
get_key_name(::Type{T}) where T<:PQ = [:i, :id]   # (bus_number, load_id)
get_key_name(::Type{Branch}) = [:i, :j, :ckt]     # (from_bus, to_bus, circuit)
```

**The event key must match the model's key structure:**

```julia
# Branch 1-5 circuit "1"
FieldModification(1.0, :Branch, (1, 5, "1"), :status, :set, 0.0)

# PQ load at bus 3, id "1": increase constant-power P by 0.1 pu (10 MW on a 100 MVA base)
FieldModification(5.0, :PQ, (3, "1"), :pl, :add, 0.1)
```

Keys are converted to the key field types (`Int32`, `String3`) by `normalize_key`
when they are looked up.

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

A second method, `apply!(sys, status, mod)`, also updates the solver's status dict when
the modified field is the model's status field. The status dict holds the in-service
masks used by residual evaluation, so a device switched by an event takes effect in
the running simulation.

**Why this works for all models:**
- Uses `sys.order.key_indices` (no new infrastructure)
- Uses `getproperty` for field access (generic)
- The residual workspaces hold references to the same model arrays, so the change is
  seen at the next residual evaluation
- No model-specific code paths

## Checks when the problem is built

`TDSProblem` (and `YBusTDSProblem`) check every `FieldModification` of its events when it
is built (`check_events`), in time order and each on top of the ones before it:

- the model, the instance and the field exist;
- the mass coefficients do not change (the mass matrix is fixed during a simulation);
- the modified instance satisfies `check_parameter_domains` of its model, where the model
  has it (generated models do): the domains, the constraints that read parameters only and
  the divisors that must not be zero. The error names the model, the instance, the field,
  the operation and the value, and the check. A hand-written model without the function is
  not checked.

Where the automata of a model govern a differential state, the admission of the governed
states is checked again when the event applies, for the models whose parameters it changes.

## Typed Wrappers

**Constructors that return a `FieldModification`:**

```julia
BranchTrip(t, from, to, ckt)   # FieldModification(t, :Branch, (from, to, ckt), :status, :set, 0.0)
BranchClose(t, from, to, ckt)  # FieldModification(t, :Branch, (from, to, ckt), :status, :set, 1.0)
FaultApply(t, bus, id="F1")    # FieldModification(t, :Fault, (bus, id), :status, :set, 1)
FaultClear(t, bus, id="F1")    # FieldModification(t, :Fault, (bus, id), :status, :set, 0)
```

**Adding a new wrapper** is a one-line function, for example:

```julia
LoadTrip(t, bus, id) = FieldModification(t, :PQ, (bus, id), :status, :set, 0.0)
```

A wrapper may also return a vector of `FieldModification`s; `build_event_callbacks`
flattens vectors.

## DifferentialEquations.jl Integration

**Each `FieldModification` becomes a `PresetTimeCallback`:**

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

The callback's `initializealg` reinitializes the algebraic variables after the event
(see "DAE Reinitialization: PowerFlowInit" below). Every callback that changes the
trajectory discontinuously (these events and the located events of the automata) saves the state before and after the event: the solution holds both at the
event time, and the interpolation of the step after it starts from the state after it
(an interpolant that reads the step's first point, as those of Rosenbrock and SDIRK
methods do, would otherwise start from the state before the event).

**TDS `solve()` builds the callbacks:**

```julia
# Inside solve(prob::TDSProblem, ...)
if !isempty(prob.events)
    init_alg = _bind_adapter(initializealg, adapter)
    event_callbacks = build_event_callbacks(prob.events, prob.sys, adapter.status;
        address_manager = prob.address_manager, adapter = adapter,
        reevaluate_fn = reevaluate_fn, initializealg = init_alg)
end
```

`adapter.status` is the solver's status dict. Other `Disturbance` subtypes (such as
`Fault3PhClearEvent`) are dispatched to
`to_callback(event, sys, address_manager, adapter, reevaluate_fn; initializealg)`, which
can modify state variables and re-evaluate residuals.

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

**Direct `FieldModification` for anything not wrapped:**

```julia
events = [
    # Set the damping of the generator at bus 1, id "1", to 2.0 at t = 3.0 s
    FieldModification(3.0, :GENCLS, (1, "1"), :D, :set, 2.0),
]
```

## Adding Support for New Models

**No event code changes are needed.** A new model works with `FieldModification` if:

1. It defines the `get_key_name()` trait
2. Its fields are accessible via `getproperty`
3. It is in the system, so it has an entry in `sys.order.key_indices`

To be switchable by events, the model also declares its `status_field` in metadata.

A modified field takes effect only if the residual kernels read it at every
evaluation, so a model must not store values derived from its parameters (compute
them in the kernel, as the fault does with its admittance). The mass matrix is built
once at setup: `build_event_callbacks` rejects, through `check_modifiable`, a
`FieldModification` that would change a mass coefficient (such as a generator's `H`),
as well as unknown models and fields.

## Architecture Summary

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

**Key point:** The event system reuses the existing `IndexOrder` infrastructure. No new
lookup mechanisms are needed.

## Fault Model and Events

A bus fault is represented by a fault **model** added to the system (defined in
`PowerfulTDS/src/models/fault.jl`), not by modifying other devices. The fault is a
shunt admittance to ground at the faulted bus that draws power only while its
`status` is `true`.

### Fault Model Structure

```julia
@kwdef mutable struct Fault{Tv} <: AbstractModel
    i::Int32           # Bus number where the fault occurs
    id::String3        # Fault identifier (e.g. "F1")
    rf::Tv             # Fault resistance (pu)
    xf::Tv             # Fault reactance (pu)
    status::Bool       # false = cleared, true = applied
end
```

**Physics:** Fault admittance `Yf = 1/Zf = Gf + jBf` with `Gf = Rf/(Rf² + Xf²)` and
`Bf = -Xf/(Rf² + Xf²)`. The fault consumes `P = Gf·V²` and `Q = -Bf·V²` while applied.
- For an inductive fault (rf ≈ 0), `Bf ≈ -1/xf` (large negative susceptance)
- The absorbed reactive power depresses the bus voltage

The kernels compute `gf` and `bf` from `rf` and `xf` at every evaluation
(`fault_admittance`), so an event that changes the fault impedance takes effect.

### BusFault Expansion

`BusFault(time_on, time_off, bus, id="F1"; xf=0.003, rf=0.0, fault_type=0, rg=0.0, xg=0.0)`
is the user-facing fault event. The `TDSProblem` constructors and the Y-bus event
processing call `expand_bus_faults!(sys, events, config)`, which for each `BusFault`:
1. Picks the fault model through the `accepts_bus_fault`/`bus_fault` extension point
   (`solvers/events.jl`), by the bus's own model type; `Bus`'s own implementation is the
   phasor `Fault`, which rejects unbalanced `fault_type`s
2. Adds the model with `add_instances!` once per (bus, id), which updates the indices,
   external params and `foreign_refs`; a repeated key must have the same parameters
3. Replaces the event by the internal primitives the chosen fault model declares
   (`FieldModification`s on the phasor `Fault`'s status, or another model's own events)
4. Warns if `xf < LOW_FAULT_REACTANCE` (0.01 pu) and loads are not 100% constant impedance

## DAE Reinitialization: PowerFlowInit

After topology changes (faults, trips), the algebraic variables must be reinitialized to
satisfy `g(x, y) = 0`. `PowerFlowInit` is defined in `PowerfulTDS/src/solvers/reinit.jl`.

### The Problem

Every term of a bus power balance carries the bus voltage magnitude
(`P_i = V_i Σ_j V_j (G_ij cos θ_ij + B_ij sin θ_ij) + ...`, and load and shunt powers
vanish with `V_i`), so `V_i = 0` satisfies both rows of bus `i` whatever current flows
into it. During a severe fault the faulted bus voltage is close to zero; Newton started
there after clearing converges to `V_i = 0`, the fault stays in place, and the
integrator stalls on the singular Jacobian (`∂P_i/∂θ = 0` at `V_i = 0`).
`reinitialize_algebraic!` therefore solves the bus rows as current mismatch
(`P_i / V_i`, `Q_i / V_i`, see `BusPowerRows`), which has the same roots with
`V_i ≠ 0` and none at `V_i = 0`; its `abstol` bounds the current mismatch of those rows.

### PowerFlowInit Algorithm

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

Construct it with `PowerFlowInit(; abstol, maxiters, initial_guess)`. Without an
`abstol`, the reinitialization bound to a simulation stops at `REINIT_TOL_FRACTION`
(1/100) of the integration's absolute tolerance, the fraction at which the integrator's
own Newton iterations stop: 1e-8 at the default `abstol = 1e-6`, 1e-12 at `1e-10`. The
state after an event is then consistent to the integration's tolerance, as DiffEq's check
after a callback (`CheckInit`) requires. An explicit `abstol` is kept.

An event reinitializes the network only where it changes the algebraic equations (their
residuals at the state differ before and after it applies; `EventReinit`), as the
automata do for their transitions. A setpoint step of a controller leaves them as they
were, and the state as consistent as the integration keeps it: a Newton solve there
would move it by less than the integration's tolerance, differently for nearby
parameters, and roughen the trajectory's dependence on them.

**Steps:**
1. **Apply the initial guess** selected by `initial_guess` (a network event), or start from
   the integrator's state (a transition of the automata, see below)
2. **Newton iteration**: Solve `g(x, y) = 0` for the algebraic `y` with the differential `x`
   fixed, bus power balance rows as current mismatch, at the integrator's time `t` (a model
   whose residuals depend on time directly, like a recording interpolated at `t`, needs it;
   every evaluation states its time, and the ADCache has no default)
3. **Check the result** against the DAE as the integration evaluates it at `t`
   (`_reinit_mismatch`): the rows without mass below the solve's tolerance, else the
   reinitialization fails like a solve that did not converge
4. Update the integrator state

### Initial Guess Options

| Option | Voltage | Angle | Use Case |
|--------|---------|-------|----------|
| `:flat_start` | V = 1.0 | Current θ | The default; converged after the short faults tested |
| `:pre_fault` | Pre-fault V | Current θ | Short faults, angles near their pre-fault values |
| `:current` | Current V | Current θ | Small disturbances; converged after the long faults tested |

The default comes from the problem's `TDSConfig`: `:flat_start` when `flat_start = true`
(the default), `:pre_fault` otherwise.

**Which instant uses which guess.** The configured guess belongs to the **network events**
(faults, trips, `FieldModification`): they change the network itself, and the
reinitialization that follows them (`EventReinit`, through DiffEq's `initializealg` after
the event's callback) applies it. The reinitializations at **transitions of the automata**
(`settle_instant!` → `_reinitialize!`, including those in the settle right after an event)
start from the integrator's state (`at_current_state`, the guess `:current`): before the
transitions the state is consistent, and a transition changes the equations at a crossing
of its surface, so the state is the nearest guess. A transition reinitialized from a flat
start inside a fault (bus voltages reset to 1 while the fault holds them near 0.3) did not
converge, where the current state converged at once. The two paths split in
`_reinitialize!` (automata.jl), the only caller of a reinitialization for transitions. The rounds of the start
use it too, where the equations of the initial modes are violated at the point the initialization handed over by more than the
problem's `TDSConfig` `abstol` (`StartResolve`), and between the rounds of the start as at any instant.

### Implementation

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

Pre-fault voltages are captured by `_bind_adapter()` in `solve()`, before the
simulation starts.

### Usage

```julia
# Default (from TDSConfig; flat start unless flat_start = false)
sol = solve(prob)

# Pre-fault voltages (good for fault-on/fault-off scenarios)
sol = solve(prob; initializealg = PowerFlowInit(initial_guess = :pre_fault))

# Current voltages
sol = solve(prob; initializealg = PowerFlowInit(initial_guess = :current))
```

**When to use each option:**
- `:flat_start` - The default; converged after the short faults tested, and did not after a
  2 s fault on ieee14 where `:current` did: which guess converges depends on the case
- `:pre_fault` - Temporary faults where the system should return to the pre-fault state;
  unsuitable when the angles have moved far from their pre-fault values (after a long fault
  or a pole slip)
- `:current` - Small perturbations, load changes, setpoint modifications; converged in the
  long-fault cases tested (it may still converge to the wrong solution from low voltages)

### Failure of the reinitialization

`reinitialize_network!` returns whether the Newton iterations converged and their result
satisfies the DAE as the integration evaluates it at the instant (`_reinit_mismatch`). Where
either fails (a numerical failure) it warns with the time, iterations, residual and worst rows
(model, variable, instance), and the caller ends the integration through
`_terminate_unless!`, the one mechanism of every reinitialization: `PowerFlowInit`,
`EventReinit`, the settle of the automata (`AutomatonSettle`), the sensitivities' `SensitivityInit`, and the Y-bus DAE's `HostReinit`.
The solution has `retcode == :ConvergenceFailure` and the trajectory up to the failed
instant. The wrappers stop at the first failure (`_running`). Errors in the models or the
code throw, as does an inconsistent initial state (DiffEq's own check at `t0`: an inconsistent
start is a model or initialization bug, not a physical outcome).

Some reinitializations use DifferentialEquations.jl's `BrownFullBasicInit` instead,
when the states to re-solve are zero-mass states that `PowerFlowInit` would otherwise
hold fixed.

## Files

- `PowerfulTDS/src/solvers/events.jl` - `FieldModification`, `apply!`, `to_callback`, wrappers, `BusFault`
- `PowerfulTDS/src/solvers/reinit.jl` - `PowerFlowInit`, flat start, pre-fault capture
- `PowerfulTDS/src/solvers/solve.jl` - Integration with `solve()`
- `PowerfulTDS/src/models/fault.jl` - Phasor `Fault` model
- `Powerful/test/tds/` - Tests
