# Time-Domain Simulation Architecture

This page describes the TDS infrastructure for developers. TDS extends the power flow architecture.

!!! note
    TDS internals are evolving. This page documents the current design; the interfaces described here may change.

## Relationship to Power Flow

TDS reuses the power flow infrastructure:

| Shared | TDS-Specific |
|--------|--------------|
| Autodiff cache (`create_adcache`) | `TDSProblem` |
| `SystemBuffers` | `TDSAdapter` |
| `ResidualBuffers` | `TDSInitialization`, `TDSInitData` |
| DAE vectors | Setpoints pipeline |

## Structure Hierarchy

```
TDSProblem                          ← Problem specification (user-owned)
├── sys: SystemModel                   (merged static + dynamic models)
├── address_manager: AddressManager
├── init: TDSInitData               ← Computed initial values
│   ├── dyn: initial states            (delta, omega, ...)
│   ├── setpoints: holding constants   (P_m0, E_q0, ...)
│   ├── pf: power flow snapshot        (v, theta)
│   └── network: non-power-flow network state (empty for phasor TDS)
├── events: disturbances
├── backend: ADBackend
└── config: SimConfig

        ↓ solve() creates solver infrastructure

TDSAdapter                          ← Created inside solve()
├── adcache                         ← Same cache type as power flow (create_adcache)
│   └── core.float_ws: SystemBuffers
│       └── workspaces_by_name      ← ResidualBuffers per model
│           └── setpoints           ← Populated from prob.init.setpoints
├── sys_ws: the Float64 SystemBuffers (direct reference for the hot path)
├── x_buf, y_buf: unpacking buffers
├── f_buf, g_buf: residual buffers
├── mass_matrix: for the DAE formulation
├── reinit_cache: lazily built Newton infrastructure for reinitialization
└── status: Dict{Symbol, Vector{Bool}} of in-service flags
```

Each level has one job. The autodiff cache, `SystemBuffers` and `ResidualBuffers` are shared
with power flow, and `TDSAdapter` is the only TDS-specific layer.

## Problem and Solver Infrastructure

A `TDSProblem` is a specification that user code creates and may reuse. The `TDSAdapter` and
the autodiff cache are working memory that `solve()` creates for each solve.

| | `TDSProblem` | `TDSAdapter` + autodiff cache |
|---|---|---|
| Created in | User code | `solve()` |
| Contains | Specification and initialization data | Working memory |
| Lifetime | Reusable | Per solve |
| Inspect | Yes | Internal |

## The Setpoints Flow

```
setup_tds(sol_pf, dyn_models) → TDSInitialization
    ↓ TDSProblem(init; ...) stores
prob.init.setpoints              ← Source of truth (in problem)
    ↓ solve() calls create_adcache(prob)
    ↓ workspaces now exist
populate_setpoints!(adcache, prob.init.setpoints)   (inside TDSAdapter)
    ↓ copies to
ResidualBuffers.setpoints        ← Working copy (in workspace)
```

The workspaces are created inside `solve()` and do not exist during problem construction. The
problem therefore stores the initial values, and each `solve()` creates fresh workspaces and
copies the setpoints from `prob.init`. The problem remains the canonical source for repeated
solves.

!!! warning "Algebraic Equation Form"
    Controller algebraic equations must use the explicit form `g = f - y` for Y-bus
    solver compatibility. See the "Explicit Algebraic Equation Form" warning in
    [Implementing Models: Residual Declarations](models.md#Residual-Declarations).

## Standard and Staged Construction

```julia
using Powerful, PowerfulCases
using PowerfulTDS: prepare_tds_system

case = load("ieee14")
sys_pf = SystemModel(parse_network(case.raw))
sol_pf = solve(ACPowerFlowProblem(sys_pf), Newton())
dyn_models = parse_dynamics(case.dyr, sys_pf)

# Standard API (most users)
init = setup_tds(sol_pf, dyn_models)
prob = TDSProblem(init; tspan = (0.0, 10.0))

# Staged API (debugging/customization)
sys_tds, links = prepare_tds_system(sys_pf, sol_pf, dyn_models)  # Inspect merged system
init = TDSInitialization(sys_tds, links, sol_pf)                  # Inspect initial states
prob = TDSProblem(init; tspan = (0.0, 10.0))
```

`TDSInitialization` has the fields `sys`, `dyn`, `setpoints`, `pf` and `config`.
`TDSProblem(init; ...)` takes `init.sys` as the system of the problem, packs `dyn`,
`setpoints` and `pf` into a `TDSInitData` (`prob.init`) and carries `init.config`. It also
expands `BusFault` events, which may add fault models to the system (see
[TDS Event System](events.md)).

## Key Functions

| Function | When Called | Purpose |
|----------|-------------|---------|
| `prepare_tds_system()` | `setup_tds` | Merge static and dynamic models, build `ModelLinks` |
| `TDSInitialization()` | `setup_tds` | Compute initial states |
| `create_adcache(prob)` | Inside `solve()` | Create the autodiff cache and its `SystemBuffers` |
| `TDSAdapter()` | Inside `solve()` | Wrap the cache for DifferentialEquations.jl |
| `populate_setpoints!()` | Inside `TDSAdapter()` | Copy setpoints into the workspaces |
| `init_tds_buffers!()` | Inside `TDSAdapter()` | Fill `x_buf`, `y_buf` from the initial conditions |
| `build_event_callbacks()` | Inside `solve()` | Turn events into callbacks |
| `build_automaton_callback()` | Inside `solve()` | Turn the automata of the models into a located callback |

## Callbacks in `solve`

`solve(prob::TDSProblem, alg; ...)` (`PowerfulTDS/src/solvers/solve.jl`) merges the following
callbacks into one `CallbackSet`:

1. The event callbacks from `build_event_callbacks`: one `PresetTimeCallback` per event,
   which reinitializes the algebraic variables with `initializealg` (default `PowerFlowInit`).
2. The automaton callback from `build_automaton_callback`: located events that switch the
   modes of the automata of the models (anti-windup among them) where their conditions change.
3. The recording callback, present when a `recorder` (`InjectionRecorder`) is passed: a
   `SavingCallback` that commits the recorded contributions at each saved step.
4. The `callback` keyword of the user, if given.

The recording callback is a `SavingCallback` because a `DiscreteCallback` interferes with
DAE reinitialization after events when the solver uses automatic differentiation (for
example `FBDF(autodiff=true)`), which leads to an `InitialFailure`.

The default algorithm is `QNDF(linsolve=KLUFactorization())`. The defaults of `abstol` and
`reltol` come from the `TDSConfig` of the problem, and explicit keywords take precedence.

### Rates Below the Integration Tolerance

A transition that fires on the zero of a value holds some other values at zero by continuity
(`Transition.continuity`) and decides them by their rate in the discrete state it enters. The
rate is computed at a state that the integration knows to its tolerance,
`tol_j = abstol_j + reltol_j |u_j|` for entry `j`. The rate `ġ` is therefore known to the
first-order bound

```math
B = \sum_j \left| \frac{\partial \dot g}{\partial u_j} \right| \mathrm{tol}_j ,
```

and its sign decides the value only where `|ġ| > B`. Otherwise the value stays at its zero
through the instant. It is kept there against the release of an instant that changes the
equations, and it is measured from where the instant leaves it.

The bound is large for a fast-oscillating state. A state oscillating at 10² to 10³ per second
has a rate orders of magnitude larger than a displacement of tolerance size suggests, so a
state error of 10⁻⁶ moves its rate by 10⁻⁴ of the scale of the oscillation. A clamp release
whose measured state is continuous across it (its rate there is zero) could then chatter up
to the budget of rounds at the tolerances where the computed sign came out positive.

The gradient of the rate with respect to the state is exact. It is computed from the Jacobian
`J` of the DAE at the state and at two points along the trajectory, from the gradient of the
value (all by Duals, `width` columns at a time), and from one adjoint solve with the
factorization of `F_ZZ`. With `p = ∇g` and `λ = F_ZZ⁻ᵀ p_Z`,

`∇ġ = d/dε ∇g(u + εu̇) + J_Wᵀ M_W⁻¹ (p_W - J_ZWᵀ λ) - d/dε J_Z(u + εu̇)ᵀ λ`.

The expression was derived with SymPy and tested against central differences in every
coordinate. It costs 0.6 ms per decision on a representative case, with 8 to 10 decisions per
run (0.2 % of a run of 2 to 4 s). A fixed-step integration has no tolerance and decides by the
raw rate.

### Known Gaps of the Decisions at an Instant

- A fixed-step integration decides the rate of a held value by its raw sign. The Y-bus
  trapezoidal integration (`solve_trapezoidal`, `FixedStepAutomata`) has no tolerance from
  which to bound the rate, so a value held at zero by continuity is released or kept by the
  sign of the computed rate, however small. Where that rate is zero in exact arithmetic (a
  release with the measured state continuous across it, as for an anti-windup clamp), the sign
  is rounding, and the release can chatter as it did in the adaptive integration before the
  bound.
- Y-bus hosts have no hold within the solve tolerance. After a reinitialization, the phasor
  TDS decides a value whose displacement from its anchor is within what the tolerance of the
  solve allows as if it were at the surface (bound `τ Σ_k |λ_k| s_k`). This requires the row
  scales of the convergence test of the reinitialization, which the Y-bus simulation does not
  provide. A value displaced by the solve alone is therefore decided by its sign there: the PI
  of REECA1 held at its limit can be released by a displacement of the size of the solve and
  held again an instant later.

### The QNDF Fork

The default solver is the QNDF method of OrdinaryDiffEqBDF, carried as a fork of that one
sublibrary. The root `Project.toml` pins the fork by commit in `[sources]`
(`cuihantao/OrdinaryDiffEq.jl`, branch `qndf-restart-v2.4.12`, `lib/OrdinaryDiffEqBDF`), and
every other OrdinaryDiffEq sublibrary comes from the registry. The branch is the v2.4.12
release plus three commits, each changing QNDF only:

1. The interpolation of a step that a callback truncates is rebased onto the length of the
   step (also proposed upstream as SciML/OrdinaryDiffEq.jl#4698).
2. A restart after an event is seeded from the vector field after the event, with
   `cache.order = 1`.
3. A Jacobian and iteration matrix are taken fresh at the first step after a discontinuity.

The dense output of an event-truncated step, the sensitivities read from it, the held states
of anti-windup limiters and the work after events depend on these commits. Tests gated on
`QNDF_REBASES_TRUNCATED_STEPS` assume them.

One known limitation is never reached from Powerful. The QNDF interpolation of a zero-length
interval whose stored derivatives are empty fails, and such an interval exists only where an
integration saves t0 twice. The start of the automaton callback (`build_automaton_callback`,
`initialize`) can move the state, by solving the algebraic point from the initial modes or by
re-solving the sensitivities after memory starts. It then replaces the saved t0 with the
settled state, so that DiffEq saves t0 only once. Code should not depend on a second save at t0.

To upgrade, cherry-pick the three commits onto the release of OrdinaryDiffEqBDF to move to,
run the tests of the BDF sublibrary, push the branch, and change `rev` (and the branch name)
in `[sources]`. Then run the Powerful test suite.

Packages installed from the registry or by URL (`Pkg.add`) do not read the root `[sources]`
and run the released OrdinaryDiffEqBDF. All functionality is available, with three
differences: a dense-output value inside a step that ends at a located event is off by about
1e-5 (the sampled and saved points are exact), restarts after events take more steps and
function evaluations, and a state that an anti-windup limiter holds drifts by the tolerance of
the integration, where the fork keeps it fixed. The fork is installed with
`Pkg.add(url = "https://github.com/cuihantao/OrdinaryDiffEq.jl", rev = "qndf-restart-v2.4.12", subdir = "lib/OrdinaryDiffEqBDF")`
in the environment.

## Augmenting a Simulation

A package can carry more state through a TDS integration than the plain DAE. A
parameter-sensitivity integration is an example: it integrates extra state alongside the DAE
and needs a map of what happens at each instant, what to read at a guard, and how to seed the
extra state at the start. The runtime offers these as a small set of extension points. Each
extension point dispatches on one parameter, conventionally named `jump`, which is `nothing`
for the plain simulation or a value of the type of the augmentation (a struct that holds the
state the augmentation needs). The augmentation adds methods on its own `jump` type and
changes nothing elsewhere. A `jump` of `nothing` costs only the dispatch already paid to
reach a no-op method.

### Extra State

The automata read the rates of the values of the DAE through an abstract `RateDAE` (defined
in `automaton_rates.jl`). Every subtype implements three functions:

- `_rate_host`: the host whose `dae_residuals!` evaluates the rates;
- `_dual_workspace`: the Dual-typed workspace in which the values of the automata are
  evaluated;
- `_plain_values!`: restores the plain values in that workspace before a rate evaluation.

`StateRates` is the `RateDAE` of the plain simulation, built from the DAE of the host of the
simulation. The function `_rate_dae(jump, rt, host, integrator)` returns the `RateDAE` to
use, and defaults to `StateRates` for `jump === nothing`. An augmentation that integrates
extra state (own rows, mass and Jacobian) defines its own `RateDAE` subtype over the combined
state, and a method of `_rate_dae` for its `jump` type that returns an instance of it.

### The Map at Each Instant

The following extension points are called at an instant (the crossing of a located
transition, an event, or the start). Each has a no-op default for `jump === nothing`.

| Extension point | Called | What it carries |
|---|---|---|
| `_instant_begins!(jump, integrator, located)` | Before the automata settle at the instant | The values located at their zero there |
| `_instant_ends!(jump, integrator, changed::Bool)` | After the settle | Whether any discrete state changed |
| `_memory_started!(jump, integrator, rt, t)` | Once every memory is at its initial value, before the modes start | None |
| `_transitions_firing!(jump, integrator, rt, steps, t)` | With the steps of each settle round, before their sets and mode changes apply | The steps about to fire |
| `_memory_seeds(jump)` | At the start, to decide whether it re-solves what `jump` carries | The own seed state of the augmentation, or `nothing` to always re-solve |
| `_memory_seeds_changed(jump, seeds)` | Compared against a fresh `_memory_seeds(jump)` after the start | Whether the seeds moved |
| `_states_written!(jump, integrator, rt, writes, t)` | When a mode governing a state writes it directly: an entry (a transition just entered that mode), or a repair (an event moved the right side of a state whose mode stayed active, see `governed_states.jl`) | The states and values written |
| `_field_changed!(jump, mod::FieldModification, instance)` | When an event changes a model field | The field changed |
| `_setpoint_changed!(jump, change::SetpointChange, instance)` | When an event changes a setpoint | The setpoint changed |

### Reinitialization

After an event or a transition that changes the algebraic equations, the DAE is reinitialized
through four extension points. The runtime passes `jump` to them explicitly wherever its own
code calls them. DiffEq calls one of them by itself, through the `initializealg` of a
callback that runs after an `affect!`, and this call has no parameter to carry `jump`. At
that place `jump` travels in a field of the initialization algorithm (`EventReinit`,
`AutomatonSettle`).

| Extension point | Called | Default for `jump === nothing` |
|---|---|---|
| `_bind_adapter(alg, adapter, jump; abstol, verbose)` | Once, binding the initialization algorithm of the simulation to its adapter | Binds `alg` and ignores `jump`, since binding never depends on it |
| `_event_reinit(alg, jump, changed)` | Building the reinitialization that an event triggers | Builds `EventReinit(alg, changed, jump)`, which carries `jump` to the call made by DiffEq |
| `at_current_state(alg, jump)` | Before the reinitialization of a transition, to start the Newton iterations from the state of the integrator | Delegates to the one-argument `at_current_state(alg)`, since `jump` never changes the guess that `alg` starts from |
| `_initialize_dae!(integrator, prob, alg, jump, ::Val{iip})` | To perform the reinitialization | Delegates to the four-argument form (the plain reinitialization of the simulation) |

An augmentation that must reinitialize its own extra state alongside the DAE (for example, to
make its sensitivities consistent with a reinitialized network) adds a method of
`_initialize_dae!` for the algorithm types it intercepts and its own `jump` type. The other
three extension points need an override only if the `jump` of an augmentation should change
which algorithm they build, and no augmentation has needed this so far.

### Governed States: the Repair After an Event

A mode that governs a state (`Mode.states`, see `governed_states.jl`) writes the state onto
its right side whenever a transition enters that mode. An event can also move the right side
of a state whose mode stays active during the whole instant, and the entry write never
touches such a state. The events concerned are a fixed-time `FieldModification` or
`SetpointChange`, both of which write a column of a model. At such an event the order is:

1. The data of the event is replaced.
2. In the first round of the settle, every guard is read on the retained state. Ordinary
   transitions fire in this round, including a release that already sees the new data.
3. In the same round, every state whose mode stayed active and whose right side was touched by
   the event is put back on its relation, through the same write as the entry.

The admissibility check of `check_governed` shows that a model field is the only channel of
this kind through which an event can move a right side. Tracing which fields each right side
reads (`g.param_reads`) and testing that set against the fields an event wrote is therefore
exact and needs no tolerance.

Only the last of several `FieldModification`s applied at one time may repair
(`AutomatonSettle.repair`), so a relation that reads fields from more than one of them is read
from all of them together. This grouping does not extend to a `SetpointChange`, and a
`Disturbance` that resets a differential state directly, bypassing a field write, never
reaches the repair. Both limits are documented in the header of `governed_states.jl`, and no
model in this codebase exercises either of them today.

The extension point `_states_written!` (see the table above) receives the same arguments for
both kinds of write, so it does not need to know which one wrote a state. A fixed-time event
locates no crossing in either case (`τ_p = 0`).

### Reading the Runtime at a Guard

An augmentation that evaluates something at the crossing of a guard (a correction to a kink
in what it integrates, for example) reads the structure of the automaton runtime directly.
The augmentation already depends on the package that defines the runtime, so the hooks above
are the only points where the runtime calls out, and no separate read-only extension point
exists for this purpose.

## See Also

- [Computation Architecture](workspaces.md): `SystemBuffers` and `ResidualBuffers` details
- [System Architecture Overview](overview.md): the power flow architecture that TDS builds on
- [TDS Event System](events.md): events and reinitialization
