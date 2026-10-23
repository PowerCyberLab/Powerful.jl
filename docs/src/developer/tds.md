# Time-Domain Simulation Architecture

TDS infrastructure for developers. TDS extends the power flow architecture.

!!! note
    TDS internals are evolving. This page documents the current design, not stable APIs.

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
├── init: TDSInitData               ← Computed values, not infrastructure
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

**Why the nesting?** Each layer has one job. The autodiff cache, `SystemBuffers` and
`ResidualBuffers` are shared with power flow; only `TDSAdapter` is TDS-specific.

## Problem vs. Solver Infrastructure

| | `TDSProblem` | `TDSAdapter` + autodiff cache |
|---|---|---|
| Created in | User code | `solve()` |
| Contains | Specification + initialization data | Working memory |
| Lifetime | Reusable | Per solve |
| Inspect? | Yes | Internal |

**Design rationale**: Problems are specifications. Solver infrastructure is an
implementation detail created for each `solve()`.

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

**Why the indirection?** Workspaces don't exist during problem construction; they are
created inside `solve()`. The problem stores what the initial values should be; the
solver creates buffers and fills them.

**Multiple solves**: Each `solve()` creates fresh workspaces and copies from
`prob.init`. The problem is the canonical source.

!!! warning "Algebraic Equation Form"
    Controller algebraic equations must use the explicit form `g = f - y` for Y-bus
    solver compatibility. See the "Explicit Algebraic Equation Form" warning in
    [Implementing Models: Residual Declarations](models.md#Residual-Declarations).

## Standard vs. Staged Construction

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

**Key point**: `TDSInitialization` has the fields `sys`, `dyn`, `setpoints`, `pf` and
`config`. `TDSProblem(init; ...)` takes `init.sys` as the problem's system, packs `dyn`,
`setpoints` and `pf` into a `TDSInitData` (`prob.init`) and carries `init.config`.
It also expands `BusFault` events, which may add fault models to the system (see
[TDS Event System](events.md)).

## Key Functions

| Function | When Called | Purpose |
|----------|-------------|---------|
| `prepare_tds_system()` | `setup_tds` | Merge static + dynamic models, build `ModelLinks` |
| `TDSInitialization()` | `setup_tds` | Compute initial states |
| `create_adcache(prob)` | Inside `solve()` | Create the autodiff cache and its `SystemBuffers` |
| `TDSAdapter()` | Inside `solve()` | Wrap the cache for DifferentialEquations.jl |
| `populate_setpoints!()` | Inside `TDSAdapter()` | Copy setpoints into the workspaces |
| `init_tds_buffers!()` | Inside `TDSAdapter()` | Fill `x_buf`, `y_buf` from the initial conditions |
| `build_event_callbacks()` | Inside `solve()` | Turn events into callbacks |
| `build_automaton_callback()` | Inside `solve()` | Turn the models' automata into a located callback |

## Callbacks in `solve`

`solve(prob::TDSProblem, alg; ...)` (`PowerfulTDS/src/solvers/solve.jl`) merges these
callbacks into one `CallbackSet`:

1. **Event callbacks** from `build_event_callbacks` - one `PresetTimeCallback` per event,
   reinitializing the algebraic variables with `initializealg` (default `PowerFlowInit`)
2. **Automaton callback** from `build_automaton_callback` - located events switching the
   modes of the models' automata (anti-windup among them) where their conditions change
3. **Recording callback** when a `recorder` (`InjectionRecorder`) is passed - a
   `SavingCallback` that commits the recorded contributions at each saved step
4. **The user's `callback` keyword**, if given

The recording callback is a `SavingCallback` rather than a `DiscreteCallback` because a
`DiscreteCallback` interferes with DAE reinitialization after events when the solver
uses automatic differentiation (e.g. `FBDF(autodiff=true)`), leading to an
`InitialFailure`.

The default algorithm is `QNDF(linsolve=KLUFactorization())`. `abstol` and `reltol`
default to the problem's `TDSConfig` values; explicit keywords take precedence.

### The sign of a rate the integration does not resolve

A transition that fires on the zero of a value holds some other values at zero by continuity (`Transition.continuity`) and
decides them by their rate in the discrete state it enters. That rate is computed at a state the integration knows to its tolerance,
`tol_j = abstol_j + reltol_j |u_j|` for entry `j`, so `ġ` is known to the first-order bound

```math
B = \sum_j \left| \frac{\partial \dot g}{\partial u_j} \right| \mathrm{tol}_j ,
```

and its sign decides only where `|ġ| > B` (`_resolved_rate`, `_rate_error_bound`); otherwise the value stays at its zero through the
instant (`zero`, kept against the release of an instant that changes the equations, and measured from where the instant leaves it).
A fast-oscillating state makes the bound large against its rate: a state oscillating at 10² to 10³ per second has a rate orders of
magnitude larger than a tolerance-scale displacement would suggest, so a state error of 10⁻⁶ moves its rate by 10⁻⁴ of the
oscillation's scale, and a clamp release whose measured state is continuous across it (its rate there zero) can chatter to the
budget of rounds at the tolerances where the computed sign came out positive.

The gradient of the rate with respect to the state is exact, from the Jacobian `J` of the DAE at the state and at two points along
the trajectory, the gradient of the value (all by Duals, `width` columns at a time) and one adjoint solve with the factorization of
`F_ZZ`: with `p = ∇g`, `λ = F_ZZ⁻ᵀ p_Z`, `∇ġ = d/dε ∇g(u + εu̇) + J_Wᵀ M_W⁻¹ (p_W - J_ZWᵀ λ) - d/dε J_Z(u + εu̇)ᵀ λ` (`_rate_gradient`;
derived with SymPy and tested against central differences in every coordinate). It costs 0.6 ms per decision on a representative
case, 8 to 10 decisions per run (0.2 % of a run of 2 to 4 s).

A fixed-step integration has no tolerance, and decides by the raw rate.

### Known gaps of the decisions at an instant

- **A fixed-step integration decides a held value's rate by its raw sign.** The Y-bus trapezoidal integration
  (`solve_trapezoidal`, `FixedStepAutomata`) has no tolerance from which to bound the rate, so a value held at zero by continuity is
  released or kept by the sign of the rate it computes, however small. Where that rate is zero in exact arithmetic (a release with
  the measured state continuous across it, as an anti-windup clamp's), its sign is rounding, and the release can chatter as it did
  in the adaptive integration before the bound.
- **Y-bus hosts have no hold within the solve tolerance.** After a reinitialization, the phasor TDS decides a value whose displacement
  from its anchor is within what the solve's tolerance allows as at the surface (`_hold_within_solve!`, bound `τ Σ_k |λ_k| s_k`).
  It needs the row scales of the reinitialization's convergence test (`_solve_row_scales`), which the Y-bus simulation does not
  give, so there a value displaced by the solve alone is decided by its sign: the PI of REECA1 held at its limit can be released by
  a displacement of the solve's size and held again an instant later.

### The QNDF fork

The default solver is QNDF of OrdinaryDiffEqBDF, carried as a fork of that one sublibrary:
the root `Project.toml` pins it by commit in `[sources]`
(`cuihantao/OrdinaryDiffEq.jl`, branch `qndf-restart-v2.4.12`, `lib/OrdinaryDiffEqBDF`), and every
other OrdinaryDiffEq sublibrary comes from the registry. The branch is the v2.4.12 release plus
three commits, each changing QNDF only:

1. the interpolation of a step that a callback truncates is rebased onto the step's length
   (also proposed upstream as SciML/OrdinaryDiffEq.jl#4698);
2. a restart after an event is seeded from the vector field after it, with `cache.order = 1`;
3. a Jacobian and iteration matrix are taken fresh at the first step after a discontinuity.

Powerful does not rely on upstream accepting them, so they are carried here. The dense
output of an event-truncated step, the sensitivities read from it, the held states of
anti-windup limiters and the work after events depend on them (tests gated on
`QNDF_REBASES_TRUNCATED_STEPS` assume them).

**Known, not reached from Powerful.** QNDF's interpolation of a zero-length interval whose stored
derivatives are empty fails; the interval exists only where an integration saves t0 twice. The
automaton callback's start (`build_automaton_callback`, `initialize`) can move the state (it solves the
algebraic point from the initial modes, or re-solves the sensitivities after memory starts), and then
replaces the saved t0 with the settled state instead of letting DiffEq save it again. Do not build on a
second save at t0.

**Upgrading.** Cherry-pick the three commits onto the release of OrdinaryDiffEqBDF to move to,
run the BDF sublibrary's tests, push the branch, and change `rev` (and the branch name) in
`[sources]`; then run the Powerful suite.

**Packages installed from the registry or by URL** (`Pkg.add`) do not read the root
`[sources]`: they run the released OrdinaryDiffEqBDF. Everything works, but the points above
differ: a dense-output value inside a step that ends at a located event is off by about 1e-5
(the sampled and saved points are not), restarts after events take more steps and function
evaluations, and a state an anti-windup limiter holds drifts by the integration's tolerance
instead of staying put. To get the fork, `Pkg.add(url = "https://github.com/cuihantao/OrdinaryDiffEq.jl",
rev = "qndf-restart-v2.4.12", subdir = "lib/OrdinaryDiffEqBDF")` in the environment.

## Augmenting a Simulation

A package can carry more state through a TDS integration than the plain DAE (a
parameter-sensitivity integration, for example): extra state integrated alongside it, a
map of what happens at each instant, what to read at a guard, and how to seed it at the
start. The runtime offers these as a small set of extension points, each dispatching on
one parameter, conventionally named `jump`: `nothing` for the plain simulation, or a value
of the augmentation's own type (a struct holding whatever
state the augmentation needs). The augmentation adds methods on its own `jump` type;
nothing elsewhere changes, and a `jump` of `nothing` costs nothing more than the one
dispatch already paid to reach a no-op method.

### Extra state

The automata read the rates of the DAE's values through an abstract `RateDAE` (defined in
`automaton_rates.jl`): every subtype implements `_rate_host` (the host whose
`dae_residuals!` evaluates it), `_dual_workspace` (the Dual-typed workspace the automata's
values are evaluated in) and `_plain_values!` (restores the plain values there before a
rate evaluation). `StateRates` is the plain simulation's `RateDAE`, built from the DAE of
the simulation's own host; `_rate_dae(jump, rt, host, integrator)` hands back the `RateDAE`
to use, defaulting to `StateRates` for `jump === nothing`. An augmentation that integrates
extra state (own rows, mass, Jacobian) defines its own `RateDAE` subtype over the combined
state and a method of `_rate_dae` for its own `jump` type that returns an instance of it.

### The map at each instant

Called at an instant (a transition's crossing located, an event, or the start), each with
a no-op default for `jump === nothing`:

| Extension point | Called | What it carries |
|---|---|---|
| `_instant_begins!(jump, integrator, located)` | Before the automata settle at the instant | The values located at their zero there |
| `_instant_ends!(jump, integrator, changed::Bool)` | After the settle | Whether any discrete state changed |
| `_memory_started!(jump, integrator, rt, t)` | Once every memory is at its initial value, before the modes start | — |
| `_transitions_firing!(jump, integrator, rt, steps, t)` | With the steps of each settle round, before their sets and mode changes apply | The steps about to fire |
| `_memory_seeds(jump)` | At the start, to decide whether it re-solves what `jump` carries | The augmentation's own seed state, or `nothing` to always re-solve |
| `_memory_seeds_changed(jump, seeds)` | Compared against a fresh `_memory_seeds(jump)` after the start | Whether the seeds moved |
| `_states_written!(jump, integrator, rt, writes, t)` | When a mode governing a state writes it directly: an entry (a transition just entered that mode), or a repair (an event moved the right side of a state whose mode stayed active, `_repair_writes`, `governed_states.jl`) | The states and values written |
| `_field_changed!(jump, mod::FieldModification, instance)` | When an event changes a model field | The field changed |
| `_setpoint_changed!(jump, change::SetpointChange, instance)` | When an event changes a setpoint | The setpoint changed |

### Reinitialization

After an event or a transition that changes the algebraic equations, the DAE is
reinitialized through four extension points; `jump` reaches them explicitly everywhere the
runtime's own code calls them, and through the initialization algorithm's own `jump` field
(`EventReinit`, `AutomatonSettle`) at the one place DiffEq calls in on its own, with no
parameter of its own to carry one (a callback's `initializealg`, run by DiffEq's own
callback machinery after an `affect!`).

| Extension point | Called | Default for `jump === nothing` |
|---|---|---|
| `_bind_adapter(alg, adapter, jump; abstol, verbose)` | Once, binding the simulation's initialization algorithm to its adapter | Binds `alg`, ignoring `jump`: binding never depends on it |
| `_event_reinit(alg, jump, changed)` | Building the reinitialization an event triggers | Builds `EventReinit(alg, changed, jump)`, carrying `jump` for the DiffEq-called boundary below |
| `at_current_state(alg, jump)` | Before a transition's reinitialization, to start the Newton iterations from the integrator's state | Delegates to `at_current_state(alg)`, the 1-argument form: `jump` never changes which guess `alg` starts from |
| `_initialize_dae!(integrator, prob, alg, jump, ::Val{iip})` | To actually reinitialize | Delegates to the 4-argument form (the plain simulation's reinitialization, unchanged) |

An augmentation that must reinitialize its own extra state alongside the DAE's (e.g. make
its sensitivities consistent with a reinitialized network) adds a method of
`_initialize_dae!` for the algorithm types it needs to intercept and its own `jump` type;
`_bind_adapter`, `_event_reinit` and `at_current_state` need no override unless an
augmentation's `jump` should change which algorithm these build, which none has needed so
far.

### Governed states: the repair after an event

A mode that governs a state (`Mode.states`, `governed_states.jl`) writes it onto its right side whenever a
transition enters that mode (`_enter_governed!`). An event (a fixed-time `FieldModification` or `SetpointChange`,
both a write to a model's own column) can also move the right side of a state whose mode stays active the whole
instant, which the entry alone never touches; the order at such an event is: the event's data is replaced, then
every guard is read on the retained state in the first round of the settle (ordinary transitions fire there,
including a release that already sees the new data), and in that same round every state whose mode stayed
active and whose right side the event touched is put back on its relation, through the entry's own write
(`_repair_writes`, `check_governed`'s admissibility is what proves a model field is the only channel of this
kind an event can move a right side through, so tracing which fields each right side reads, `g.param_reads`, and
testing that set against the fields an event actually wrote is exact, no tolerance). Only the last of several
`FieldModification`s applied at one time may repair (`AutomatonSettle.repair`, `_simulation_callbacks`), so a
relation that reads fields from more than one of them is read from all of them together; this grouping does not
extend to a `SetpointChange`, and a `Disturbance` that resets a differential state directly (bypassing a field
write) never reaches this repair at all -- both documented as open scope in `governed_states.jl`'s header, not
demonstrated against by any model in this codebase today. `_states_written!` (the table above) does not need to
know which of the two wrote a state: both call it with the same shape, and a fixed-time event locates no
crossing either way (`τ_p = 0`).

### Reading the runtime at a guard

An augmentation that evaluates something at a guard's crossing (a correction to a kink in
what it integrates, say) reads the automaton runtime's own structure directly: it already
depends on the package that defines it, so there is no separate read-only extension point
for this — only the hooks above, where the runtime calls out.

## See Also

- [Computation Architecture](workspaces.md) - SystemBuffers/ResidualBuffers details
- [System Architecture Overview](overview.md) - Power flow architecture (TDS builds on this)
- [TDS Event System](events.md) - Events and reinitialization
