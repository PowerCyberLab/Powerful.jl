"""
Time-Domain Simulation (TDS) system

This file contains the main TDS types and infrastructure for time-domain simulation.
Dynamic models replace static power flow models (PV/Slack) during TDS initialization.

# Design Philosophy
- Reuse `SystemModel` infrastructure (no new system struct)
- Add dynamic models to existing model collection
- Filter replaced PV/Slack instances
- Store links in problem metadata

# Module Structure
- tds.jl (this file): Core types (TDSConfig, TDSInitialization, TDSProblem)
- linking.jl: Model linking and system construction
- initialization.jl: State initialization system
- adapter.jl: Evaluation machinery (TDSAdapter)
- solve.jl: DifferentialEquations.jl integration
- validation.jl: Validation and solution analysis
"""

# Use abstract interface to decouple from PowerfulPF (enables parallel precompilation)
using PowerfulCore: AbstractSteadyStateSolution, get_system, get_address_manager, get_state_vector

# ============================================================================
# TDS Settings (Runtime Configuration)
# ============================================================================

"""Default maximum Newton iterations for algebraic reinitialization."""
const DEFAULT_REINIT_MAXITERS = 100

"""
    TDSConfig

Runtime configuration for Time-Domain Simulation.

These settings control solver behavior and can be loaded from TOML configuration files.

# Fields
- `convert_pq_to_zip::Bool`: Convert PQ loads to ZIP model during TDS setup (default: false)
- `p_zip::NTuple{3,Float64}`: Active power ZIP composition (Z%, I%, P%), must sum to 1.0 (default: (0.0, 0.0, 1.0) = constant power)
- `q_zip::NTuple{3,Float64}`: Reactive power ZIP composition (Z%, I%, P%), must sum to 1.0 (default: (0.0, 0.0, 1.0) = constant power)
- `flat_start::Bool`: Reinitialize algebraic variables after events from a flat start
  (V = 1, θ = 0) rather than the pre-fault voltages (default: true)
- `abstol::Float64`: Absolute tolerance for the DAE solver (default: 1e-6)
- `reltol::Float64`: Relative tolerance for the DAE solver (default: 1e-3)
- `maxiters::Int`: Maximum Newton iterations of the post-event reinitialization
  (default: $(DEFAULT_REINIT_MAXITERS))
- `adjust_limits::Bool`: When a limited variable (e.g. a governor valve position)
  initializes outside its limits, move the violated limit to the initial value and warn,
  so the simulation starts in steady state with the device at its limit; `false` throws
  an error listing the violations instead (default: true). See
  [`enforce_initial_limits!`](@ref).
- `locate_extrema::Bool`: Catch a limit or control-mode change that is reached and released
  again within a single solver time step, which the integrator's own sampling would
  otherwise step over entirely (default: true). This covers cases such as a governor valve
  that briefly hits its travel limit between two accepted steps, or a protection logic
  condition that turns on and back off within one step. Set to `false` to skip this check
  and let the solver step over such brief events unseen; doing so is cheaper (it avoids extra
  function evaluations to search inside each step) but risks missing short-lived limit or
  mode changes, which can leave the simulation in a state inconsistent with the physical
  device.

`solve(::TDSProblem)` uses these as defaults: explicit `abstol`/`reltol` keywords
and an explicit `initializealg` take precedence.

# ZIP Model
The ZIP model represents load as combination of:
- Z (constant impedance): P ∝ V², Q ∝ V²
- I (constant current): P ∝ V, Q ∝ V
- P (constant power): P = const, Q = const

# Example TOML
```toml
[solver.TDS]
convert_pq_to_zip = true
p_zip = [0.3, 0.3, 0.4]  # 30% Z, 30% I, 40% P for active power
q_zip = [0.2, 0.3, 0.5]  # 20% Z, 30% I, 50% P for reactive power
abstol = 1e-8
```

# Example Julia
```julia
settings = TDSConfig(
    convert_pq_to_zip = true,
    p_zip = (0.3, 0.3, 0.4),
    q_zip = (0.2, 0.3, 0.5)
)
```
"""
Base.@kwdef struct TDSConfig
    convert_pq_to_zip::Bool = false
    p_zip::NTuple{3,Float64} = (0.0, 0.0, 1.0)
    q_zip::NTuple{3,Float64} = (0.0, 0.0, 1.0)
    flat_start::Bool = true
    abstol::Float64 = 1e-6
    reltol::Float64 = 1e-3
    maxiters::Int = DEFAULT_REINIT_MAXITERS
    adjust_limits::Bool = true
    locate_extrema::Bool = true

    # Every TDSConfig is validated at construction (keyword or TOML)
    function TDSConfig(convert_pq_to_zip, p_zip, q_zip, flat_start, abstol, reltol, maxiters, adjust_limits,
                       locate_extrema)
        settings = new(convert_pq_to_zip, p_zip, q_zip, flat_start, abstol, reltol, maxiters, adjust_limits,
                       locate_extrema)
        validate(settings)
        return settings
    end
end

"""Tolerance on the sum of a (Z, I, P) load composition."""
const ZIP_SUM_TOL = 1e-6

"""
    validate(settings::TDSConfig)

Validate TDS settings. Throws ArgumentError if invalid. Called by the `TDSConfig`
constructor, so every `TDSConfig` in existence has passed it.

Checks, for both `p_zip` and `q_zip`:
- Each (Z, I, P) fraction lies in `[0, 1]` (NaN is rejected)
- The fractions sum to 1.0 (±$(ZIP_SUM_TOL)), so ZIP conversion preserves the load's
  consumption at the operating voltage
"""
function validate(settings::TDSConfig)
    _validate_zip_fractions(settings.p_zip, "p_zip")
    _validate_zip_fractions(settings.q_zip, "q_zip")
    return nothing
end

function _validate_zip_fractions(fractions::NTuple{3,Float64}, name::String)
    # Written as `all(in-range)` so NaN fails the check
    all(x -> 0.0 <= x <= 1.0, fractions) ||
        throw(ArgumentError("$name components (Z, I, P) must each lie in [0, 1], got $fractions"))
    s = sum(fractions)
    abs(s - 1.0) <= ZIP_SUM_TOL ||
        throw(ArgumentError("$name must sum to 1.0, got $s ($fractions)"))
    return nothing
end

export TDSConfig, validate

# TDS config is registered in PowerfulTDS.__init__() to survive precompilation
# See PowerfulTDS.jl for the register_solver_config!(:TDS, TDSConfig) call

# ============================================================================
# ZIP Load Conversion
# ============================================================================

"""
    apply_zip_conversion!(sys::SystemModel, sol_pf::AbstractSteadyStateSolution, settings::TDSConfig)

Convert the loads of `sys` (the models that are `is_zip_load`, such as `PQ`) to
the ZIP model while preserving their consumption at the operating voltage of the power
flow solution `sol_pf`. Refuses, before converting any load, a system holding a model that
does not declare whether it is a ZIP load (`check_zip_load_declared`).

This implements PSS/E's "Reconstruct loads using present voltage" operation: each
converted load consumes the **same power** at the solved voltage as before. Loads are
matched to the power flow by their (bus, id) key.

Only the loads of `sys` change. `setup_tds` converts the loads of the
systems it builds (copies of the power flow models), so the power flow system of
`sol_pf` keeps its loads and a later power flow on it is unaffected.

# ZIP Model
After conversion, load power consumption becomes voltage-dependent:
- P = pl·f(v) + ip*V + yp*V²  (constant power + current + impedance)
- Q = ql·f(v) + iq*V + yq*V²

with `f = constant_power_factor(v, vmin)` (1 above the constant-power breakpoint).

# Algorithm
For each load at operating voltage V_op, with total consumption P_total, Q_total:
```
yp_new = z% × P_total / V_op²   # Constant impedance (Z)
ip_new = i% × P_total / V_op    # Constant current (I)
pl_new = p% × P_total / f       # Constant power (P)
```

# Arguments
- `sys::SystemModel`: System whose loads are converted (in place)
- `sol_pf::AbstractSteadyStateSolution`: Solved power flow (bus voltages)
- `settings::TDSConfig`: Settings containing `p_zip` and `q_zip` tuples

# Example
```julia
sol = solve(ACPowerFlowProblem(sys), Newton())
settings = TDSConfig(
    convert_pq_to_zip = true,
    p_zip = (0.3, 0.3, 0.4),  # 30% Z, 30% I, 40% P
    q_zip = (0.2, 0.3, 0.5)   # 20% Z, 30% I, 50% P
)
apply_zip_conversion!(sys_tds, sol, settings)
```

# See Also
- [`TDSConfig`](@ref): Configuration struct with ZIP percentages
- [`validate`](@ref): Validates that ZIP percentages sum to 1.0
"""
function apply_zip_conversion!(sys::SystemModel, sol_pf::AbstractSteadyStateSolution, settings::TDSConfig)
    # Operating voltage of each power flow load, looked up by (bus, id) key
    sys_pf = get_system(sol_pf)
    v_addrs = get_addresses(get_address_manager(sol_pf), :PQ, :v)
    v_pq = get_state_vector(sol_pf)[v_addrs]
    vmin = sys_pf.settings.constant_power_vmin

    check_zip_load_declared(sys)
    for loads in values(sys.models)
        isempty(loads) && continue
        is_zip_load(get_base_model_type(loads)) || continue
        for k in eachindex(loads)
            pf_idx = get_index_by_key(sys_pf, :PQ, (loads.i[k], loads.id[k]))
            _redistribute_zip!(loads, k, v_pq[pf_idx], settings.p_zip, settings.q_zip, vmin)
        end
    end

    nothing
end

"""
    check_zip_load_declared(sys::SystemModel)

Throw an `ArgumentError` naming each model of `sys` that does not declare whether it is a ZIP
load (`PowerfulCore.is_zip_load`): the ZIP conversion cannot tell whether to convert it.
"""
function check_zip_load_declared(sys::SystemModel)
    refused = [String(name) for (name, model) in pairs(sys.models) if !declares_zip_load(get_base_model_type(model))]
    isempty(refused) || throw(ArgumentError(
        "The ZIP load conversion (convert_pq_to_zip) cannot tell whether these models are ZIP loads: " *
        join(refused, ", ") * " (no is_zip_load declared: a method of PowerfulCore.is_zip_load returning " *
        "true for a ZIP load and false otherwise)"))
    return nothing
end

"""
    _redistribute_zip!(loads, k, v, p_zip, q_zip, vmin)

Split load `k`'s active and reactive consumption at operating voltage `v` into
constant-impedance, constant-current and constant-power parts by the (Z, I, P)
fractions `p_zip` / `q_zip`, preserving consumption at `v`. Consumption follows
`calc_zip_load` with the constant-power breakpoint `vmin`; below it the constant-power
part draws `pl·f(v)`, so its share is divided by `f = constant_power_factor(v, vmin)`.
The coefficients have their phasor meaning for every ZIP load model (`is_zip_load`).
"""
@inline function _redistribute_zip!(loads, k::Int, v, p_zip, q_zip, vmin)
    z_p, i_p, p_p = p_zip  # Active power: (Z%, I%, P%)
    z_q, i_q, p_q = q_zip  # Reactive power: (Z%, I%, P%)
    v2 = v * v

    # Total power consumption at the operating voltage
    total = calc_zip_load(v, loads.pl[k], loads.ql[k], loads.ip[k], loads.iq[k],
                          loads.yp[k], loads.yq[k], vmin)
    f = constant_power_factor(v, vmin)

    # Redistribute to maintain same consumption at V_op
    loads.yp[k] = z_p * total.p / v2   # Constant impedance (Z)
    loads.ip[k] = i_p * total.p / v    # Constant current (I)
    loads.pl[k] = p_p * total.p / f              # Constant power (P)

    loads.yq[k] = z_q * total.q / v2
    loads.iq[k] = i_q * total.q / v
    loads.ql[k] = p_q * total.q / f
    nothing
end

export apply_zip_conversion!


# ============================================================================
# TDSInitialization - Initial conditions for TDS
# ============================================================================

"""
    TDSInitialization

The dynamic states and controller setpoints computed from a converged power flow
solution: the operating point a time-domain simulation starts from at `t = 0`. Build it
with [`setup_tds`](@ref), then pass it to [`TDSProblem`](@ref).

# Fields
- `dyn::NamedTuple` - Per-model initial states, e.g. `(GENROU=(delta=[], omega=[]), ...)`
  (rotor angle `delta` in rad, speed `omega` in pu, one entry per instance)
- `setpoints::NamedTuple` - Per-model controller references, e.g.
  `(GENROU=(P_m0=[], E_q0=[]), ...)` (mechanical power and field voltage setpoints
  computed from the power flow solution, so each machine starts in steady state)
- `pf::NamedTuple` - The power flow snapshot `(v=[], theta=[])` it was computed from

# Inspection
```julia
init = setup_tds(sol_pf, dyn_data)
init.dyn.GENROU.delta      # Initial rotor angles (rad)
init.dyn.GENROU.omega      # Initial speeds (should be ~1.0 pu)
init.setpoints.TGOV1.P_m0  # Mechanical power setpoints
```
"""
struct TDSInitialization
    sys::SystemModel
    dyn::NamedTuple
    setpoints::NamedTuple
    pf::NamedTuple
    config::SimConfig
end

"""
    TDSInitialization(sys::SystemModel, links::ModelLinks, pf_sol; config=SimConfig())

Compute initial conditions from power flow solution.

# Arguments
- `sys::SystemModel` - TDS system with dynamic models (from `prepare_tds_system`)
- `links::ModelLinks` - Model links mapping dynamic to static models
- `pf_sol` - Converged power flow solution
- `config::SimConfig` - Simulation config carried on to the `TDSProblem`

# Example
```julia
sys_tds, links = prepare_tds_system(sys_pf, sol_pf, dyr_data)
init = TDSInitialization(sys_tds, links, sol_pf)
```
"""
function TDSInitialization(sys::SystemModel, links::ModelLinks, pf_sol; config::SimConfig = SimConfig())
    result = enforce_initial_limits!(sys, get_solver_config(config, :TDS, TDSConfig).adjust_limits) do
        initialize_dynamic_models(sys, links, pf_sol)
    end
    return TDSInitialization(sys, result.dyn, result.setpoints, result.pf, config)
end

function Base.show(io::IO, init::TDSInitialization)
    println(io, "TDSInitialization:")
    for (model, states) in pairs(init.dyn)
        if !isempty(states)
            n = length(first(values(states)))
            fields = join(keys(states), ", ")
            println(io, "  $model ($n): $fields")
        end
    end
    if !isempty(init.setpoints)
        sp_models = join(keys(init.setpoints), ", ")
        print(io, "  Setpoints: $sp_models")
    end
end

export TDSInitialization

# ============================================================================
# Problem Formulation - TDS Problem Definition
# ============================================================================

"""
    AbstractSteadyState

The steady state a `TDSProblem` starts in: [`Equilibrium`](@ref) (every state at rest)
or [`SinusoidalSteadyState`](@ref) (instantaneous network states oscillating at the
system frequency). It sets the initial rates of change that `validate_tds_initialization`
checks and the default integrator of `solve`.
"""
abstract type AbstractSteadyState end

"""
    Equilibrium()

Steady state with every state at rest at the initial time (`M·ẋ = f = 0`, `g = 0`), as
a phasor simulation starts from a power flow.
"""
struct Equilibrium <: AbstractSteadyState end

"""
    SinusoidalSteadyState(rates)

Sinusoidal steady state: the instantaneous (three-phase) network states oscillate at the
system frequency. `rates` holds their rates of change at `t = 0`, the time of the phasors the
initialization builds the states from (with the states in `TDSInitData.network`), per model and
state in the layout of `TDSInitData.network`; states not listed are at rest. A problem that starts
at another time takes the steady state there (`init_tds_buffers!`, `init_tds_rates`: each listed
state `x(t₁) = x₀ cos ωt₁ + (ẋ₀/ω) sin ωt₁` and its rate, `ω` the system's `omega_base`).
"""
struct SinusoidalSteadyState <: AbstractSteadyState
    rates::NamedTuple
end

"""
    initial_rates(steady_state) -> NamedTuple

Rates of change at the initial time of the states that are not at rest, per model and
state (empty at an equilibrium).
"""
initial_rates(::Equilibrium) = NamedTuple()
initial_rates(s::SinusoidalSteadyState) = s.rates

export AbstractSteadyState, Equilibrium, SinusoidalSteadyState

"""
    TDSInitData(; dyn, setpoints, pf, network, steady_state)

States and setpoints a `TDSProblem` starts from.

# Fields
- `dyn::NamedTuple`: Per-model dynamic states, e.g. `(GENROU = (delta = [...], ...),)`
- `setpoints::NamedTuple`: Per-model setpoints, e.g. `(TGOV1 = (pref0 = [...],),)`
- `pf::NamedTuple`: Power flow snapshot `(v, theta)` of the phasor buses
- `network::NamedTuple`: Per-model values of network model variables beyond the phasor
  bus voltages of `pf`, e.g. `(Model = (var = [...], ...),)`: states and own algebraic
  variables of network models that a power flow does not solve; empty for a phasor
  simulation
- `steady_state::AbstractSteadyState`: The steady state these values form:
  [`Equilibrium`](@ref) (default) or, for oscillating network states,
  [`SinusoidalSteadyState`](@ref) with their initial rates of change
"""
Base.@kwdef struct TDSInitData
    dyn::NamedTuple = NamedTuple()
    setpoints::NamedTuple = NamedTuple()
    pf::NamedTuple = NamedTuple()
    network::NamedTuple = NamedTuple()
    steady_state::AbstractSteadyState = Equilibrium()
end

"""
    StructuralTraces

What the automata of a problem's models read, traced once from their kernels and shared by
every solve of the problem and its working copies (`TDSProblem` `traces`): the variables the
values of each model read, as the tracer path of `trace_value_reads` finds them, which holds
whatever the parameters, setpoints, state and time. Filled lazily by the first solve that
builds the automata of a model (`value_reads`), so that a problem with automata that cannot be
realized stays constructible. What depends on parameter or setpoint values (a kernel that
branches on them, `ValueReads` `branches`; the residual traces of `AlgebraicCheck`) is never
kept here.
"""
struct StructuralTraces
    reads::Dict{Symbol, Any}
    lock::ReentrantLock
end

StructuralTraces() = StructuralTraces(Dict{Symbol, Any}(), ReentrantLock())

"""
    TDSProblem{T, FType, B}

A time-domain simulation to run: a power system with its dynamic models (generators,
governors, exciters) already initialized at a steady-state operating point, a time span,
and a list of disturbance events. Pass it to `solve` to integrate the system's dynamics.

Build one from a [`TDSInitialization`](@ref) (the result of [`setup_tds`](@ref)):

```julia
prob = TDSProblem(init; tspan=(0.0, 5.0), events=[BusFault(0.5, 0.6, 10), BranchTrip(1.0, 1, 5, "1")])
sol = solve(prob)
```

# What You Can Read or Set
- `tspan::Tuple{T, T}`: simulation time span `(t_start, t_end)` in seconds, set as a
  keyword to `TDSProblem(init; tspan=...)`.
- `events::Vector`: the disturbances to apply during the run ([`BusFault`](@ref),
  [`BranchTrip`](@ref), [`BranchClose`](@ref), [`FieldModification`](@ref),
  [`SetpointChange`](@ref)), also set as a keyword.
- `sys::SystemModel`: the system being simulated, with dynamic models in place of the
  PV/Slack generators they replace; read this to inspect the model instances and their
  parameters.
- `config::SimConfig`: the solver settings ([`TDSConfig`](@ref)) that `solve` uses as
  defaults; explicit `solve` keywords override them.

# See Also
- [`setup_tds`](@ref): build the initialization this is constructed from
- [`BranchTrip`](@ref), [`BusFault`](@ref): event constructors
- [`validate_tds_initialization`](@ref): check initial residuals

# Example
```julia
using Powerful, PowerfulCases

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
sol_pf = solve(ACPowerFlowProblem(sys), Newton())

dyn_data = parse_dynamics(case.dyr, sys)
init = setup_tds(sol_pf, dyn_data)

events = [BusFault(0.5, 0.6, 10), BranchTrip(1.0, 1, 5, "1")]
prob = TDSProblem(init; tspan=(0.0, 5.0), events=events)

sol = solve(prob)
```

# Extended help
- `formulation::FType`: the DAE formulation, supported by every model of the system
  (checked by the constructor).
- `address_manager::AddressManager`: maps model variables and residuals to positions in
  the DAE's state and residual vectors.
- `backend::B`: the automatic-differentiation backend used to build the Jacobian
  (`ForwardDiffBackend` or `EnzymeBackend`).
- `traces::StructuralTraces`: which variables each model's operating-mode logic reads,
  discovered on the first solve of the problem and reused by subsequent solves and by
  working copies.
"""
struct TDSProblem{T, FType<:FormulationTrait, B<:ADBackend}
    sys::SystemModel
    formulation::FType
    address_manager::AddressManager
    tspan::Tuple{T, T}
    events::Vector{Any}  # Placeholder for fault events (Phase 9)

    # Initialization data (single source of truth): dyn states, setpoints,
    # pf (v, theta) snapshot and the values of other network models
    init::TDSInitData

    # Autodiff backend (ForwardDiffBackend or EnzymeBackend)
    backend::B

    # Simulation config: TDSConfig solver defaults for `solve`
    config::SimConfig

    # What the automata read, traced at the first solve and shared by all (see
    # `StructuralTraces`)
    traces::StructuralTraces
end

# Internal constructor - use TDSProblem(init::TDSInitialization; ...) for public API
function TDSProblem(
    sys::SystemModel;
    formulation::FormulationTrait = PowerInjection(),
    tspan::Tuple{Real, Real} = (0.0, 1.0),
    addressing::AddressingStrategy = NaturalAddressing(),
    init::TDSInitData = TDSInitData(),
    events::Vector = Any[],
    backend::ADBackend = ForwardDiffBackend(),
    config::SimConfig = SimConfig()
)
    check_formulation_support(sys, formulation)

    # Create AddressManager and setup addresses (like ACPowerFlowProblem)
    am = AddressManager(strategy=addressing)
    setup_with_formulation!(am, sys, formulation)
    apply_addressing_strategy!(am, sys, formulation)  # AMD/CustomPermutation reordering

    T = typeof(tspan[1])
    # Convert events to Vector{Any} to match struct type
    events_any = isempty(events) ? Any[] : Any[e for e in events]
    check_events(sys, events_any)  # (an event that cannot take effect is refused here, not mid-run)
    TDSProblem(sys, formulation, am, (T(tspan[1]), T(tspan[2])), events_any, init, backend, config,
               StructuralTraces())
end

"""
    create_adcache(prob::TDSProblem; backend=prob.backend)

Autodiff cache for `prob`.
"""
create_adcache(prob::TDSProblem; backend::ADBackend = prob.backend) =
    create_adcache(prob.address_manager, prob.sys, prob.formulation, backend)

"""
    working_copy(prob::TDSProblem) -> TDSProblem

`prob` on a copy of its system (see `copy(::SystemModel)`), for a solve to apply its
events to, or a problem derived from `prob` (Y-bus) to add its fault models to.
"""
working_copy(prob::TDSProblem) =
    TDSProblem(copy(prob.sys), prob.formulation, prob.address_manager, prob.tspan,
               prob.events, prob.init, prob.backend, prob.config, prob.traces)

"""
    TDSProblem(init::TDSInitialization; kwargs...)

Create a TDS problem from initialization (which contains the system).

# Arguments
- `init::TDSInitialization` - Complete initialization from `setup_tds`

# Keyword Arguments
- `tspan::Tuple` - Simulation time span (default: (0.0, 1.0))
- `formulation::FormulationTrait` - DAE formulation (default: PowerInjection()); every
  model of the system must support it, and the dynamic models (machines, exciters,
  governors, converters) are defined in the PowerInjection form only
- `addressing::AddressingStrategy` - Ordering of the algebraic addresses (default:
  NaturalAddressing(); AMDAddressing() or CustomPermutation(perm) reorder `y`/`g`,
  state addresses keep their natural order)
- `events::Vector` - Disturbance events (default: empty)
- `backend::ADBackend` - Autodiff backend (default: ForwardDiffBackend())

# Example
```julia
init = setup_tds(sol_pf, dyr_data)

# Debug: check initial conditions
@show init.dyn.GENCLS.omega  # Should be ~1.0

prob = TDSProblem(init; tspan=(0.0, 10.0))
sol = solve(prob)
```
"""
function TDSProblem(
    init::TDSInitialization;
    tspan::Tuple{Real, Real} = (0.0, 1.0),
    formulation::FormulationTrait = PowerInjection(),
    addressing::AddressingStrategy = NaturalAddressing(),
    events::Vector = Any[],
    backend::ADBackend = ForwardDiffBackend()
)
    check_formulation_support(init.sys, formulation)

    # The problem owns its system: a copy of init.sys, so the fault models added here
    # (BusFault events) do not reach init.sys or other problems built from it
    sys = copy(init.sys)
    expanded_events = expand_bus_faults!(sys, events, init.config)
    # Voltage playbacks start on the recording at tspan[1], their breaks as events
    expanded_events = expand_voltage_playback!(sys, expanded_events, tspan, init.pf)

    init_data = TDSInitData(dyn = init.dyn, setpoints = init.setpoints, pf = init.pf)

    return TDSProblem(
        sys;
        formulation = formulation,
        tspan = tspan,
        addressing = addressing,
        init = init_data,
        events = expanded_events,
        backend = backend,
        config = init.config
    )
end

"""
    setup_tds(sol_pf, dyn_models; config=SimConfig()) -> TDSInitialization

Prepare TDS system and compute initial conditions from power flow solution.

This is the recommended entry point for TDS simulation setup. It performs:
1. System preparation (merging dynamic models with power flow system)
2. Model linking (mapping dynamic → static models)
3. State initialization from power flow solution

# Arguments
- `sol_pf`: Power flow solution (contains system via `get_system(sol_pf)`)
- `dyn_models`: Dynamic models as NamedTuple (e.g., from `parse_dynamics`)

# Keyword Arguments
- `config::SimConfig`: Simulation configuration (default: `SimConfig()`)

# Returns
- `TDSInitialization`: Complete initialization containing system, states, setpoints, and pf snapshot

# Example
```julia
init = setup_tds(sol_pf, dyr_data)
prob = TDSProblem(init; tspan=(0.0, 10.0))
sol = solve(prob)

# Access system for inspection
init.sys.models  # View dynamic models

# With custom config (system-wide settings come from the SystemModel)
config = SimConfig(solver = (TDS = TDSConfig(abstol = 1e-8),))
init = setup_tds(sol_pf, dyr_data; config = config)
```
"""
function setup_tds(sol_pf, dyn_models::NamedTuple; config::SimConfig = SimConfig())
    sys_pf = get_system(sol_pf)
    sys_tds, links = prepare_tds_system(sys_pf, sol_pf, dyn_models; config = config)
    return TDSInitialization(sys_tds, links, sol_pf; config = config)
end


"""
    Base.show(io::IO, prob::TDSProblem)

Display a concise summary of the TDS problem.
"""
function Base.show(io::IO, prob::TDSProblem{T, F, B}) where {T, F, B}
    if get(io, :compact, false)
        print(io, "TDSProblem{$T, $(nameof(F))}(...)")
    else
        formulation_name = nameof(F)
        am = prob.address_manager

        # Count variables from next_idx
        n_state = get(am.next_idx, StateVar(), 1) - 1
        n_algebs = get(am.next_idx, AlgebVar(), 1) - 1

        # Count dynamic models (those with state variables)
        dyn_models = Symbol[]
        for (model_name, assignment) in am.assignments
            state_vars = get(assignment.by_type, StateVar(), Dict{Symbol, Vector{UInt}}())
            if !isempty(state_vars)
                push!(dyn_models, model_name)
            end
        end

        println(io, "TDSProblem{$T, $formulation_name}:")
        println(io, "  Time span: $(prob.tspan[1]) → $(prob.tspan[2]) s")
        println(io, "  State variables: $n_state")
        println(io, "  Algebraic variables: $n_algebs")
        println(io, "  Events: $(length(prob.events))")
        if !isempty(dyn_models)
            println(io, "  Dynamic models: $(join(sort(dyn_models), ", "))")
        end
    end
end

export TDSProblem, TDSInitData, setup_tds
