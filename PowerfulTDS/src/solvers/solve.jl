"""
TDS Solve Interface - Integration with OrdinaryDiffEq.jl

This module provides the solve() method for TDSProblem and helper functions
for creating ODE problems from TDS adapters.
"""

using DifferentiationInterface: AutoSparse
using ADTypes: AutoForwardDiff
using SparseMatrixColorings: GreedyColoringAlgorithm
using PowerfulCore: @vprint, get_float_ws, get_system, get_formulation, get_address_manager
using PowerfulCore: TDSSolution, max_address_in_space, StateVar, AlgebVar
using Printf: @printf

# ============================================================================
# Solve Interface - Integration with OrdinaryDiffEq.jl
# ============================================================================

"""
    create_diffeq_problem(adapter::TDSAdapter, tspan::Tuple; verbose) -> ODEProblem

Create OrdinaryDiffEq.jl ODEProblem with mass matrix formulation.

Uses ODEProblem with sparse mass matrix for BDF/Rosenbrock solvers.
This formulation: M * du/dt = f(u, t) where M is singular for DAE systems.

# Arguments
- `adapter::TDSAdapter` - Unified TDS adapter with cached buffers
- `tspan::Tuple` - Time span for simulation

# Returns
- `ODEProblem` configured with mass matrix, sparse Jacobian, and variable names
"""
function create_diffeq_problem(
    adapter::TDSAdapter,
    tspan::Tuple;
    verbose::Bool = false
)
    # Get initial conditions from adapter buffers
    n_diff = adapter.n_diff
    n_alg = adapter.n_alg
    n_total = n_diff + n_alg
    @vprint "\n  States: $n_diff differential, $n_alg algebraic ($n_total total)"

    u0 = pack_state(adapter)

    # Build variable names for plot labels
    am = get_address_manager(adapter.adcache)
    syms = Vector{Symbol}(undef, n_total)
    @inbounds for addr in 1:n_diff
        syms[addr] = variable_name(am, StateVar(), UInt(addr))
    end
    @inbounds for addr in 1:n_alg
        syms[n_diff + addr] = variable_name(am, AlgebVar(), UInt(addr))
    end

    # Temp buffers for AD computation
    du_temp = Vector{Float64}(undef, n_total)
    du_tgrad_temp = Vector{Float64}(undef, n_total)

    # ODE right-hand side using functor (stable type for precompilation)
    # ODERHSFunctor has stable type regardless of adapter contents
    ode_fn! = ODERHSFunctor(Ref{Any}(adapter))

    # Mutable time reference for Jacobian computation
    t_ref = Ref(tspan[1])

    # Mutable state reference for tgrad computation
    u_ref = Ref(copy(u0))

    # Jacobian closure functor (stable type for precompilation)
    jac_closure = JacobianClosureFunctor(Ref{Any}(adapter), t_ref)

    # Prepare sparse DAE Jacobian (includes differential states; done once)
    verbose && print("\n    Preparing Jacobian coloring... ")
    t_jac_prep = @elapsed jac_setup =
        prepare_sparse_jacobian(jac_closure, zero(du_temp), zero(u0), adapter.adcache; include_differential=true,
                                sparsity=tds_jacobian_sparsity(adapter))
    jac_prototype, jac_backend, jac_prep = jac_setup.J, jac_setup.backend, jac_setup.prep
    @vprint "done ($(round(t_jac_prep*1000, digits=1)) ms)"
    @vprint "  Jacobian: $(size(jac_prototype, 1)) x $(size(jac_prototype, 2)), nnz=$(nnz(jac_prototype))"

    # Sparse Jacobian functor (stable type for precompilation)
    sparse_jac! = SparseJacobianFunctor(
        Ref{Any}(jac_closure),
        Ref{Any}(jac_prep),
        Ref{Any}(jac_backend),
        t_ref,
        du_temp
    )

    # Tgrad closure functor (stable type for precompilation)
    tgrad_closure = TgradClosureFunctor(Ref{Any}(adapter), u_ref)

    # Prepare derivative computation (done once)
    tgrad_backend = AutoForwardDiff()
    verbose && print("    Preparing tgrad... ")
    t_tgrad_prep = @elapsed tgrad_prep = DifferentiationInterface.prepare_derivative(tgrad_closure, zero(du_temp), tgrad_backend, tspan[1])
    @vprint "done ($(round(t_tgrad_prep*1000, digits=1)) ms)"

    # Sparse tgrad functor (stable type for precompilation)
    sparse_tgrad! = SparseTgradFunctor(
        Ref{Any}(tgrad_closure),
        Ref{Any}(tgrad_prep),
        Ref{Any}(tgrad_backend),
        u_ref,
        du_tgrad_temp
    )

    # Wrap function with mass matrix using ODEFunction
    ode_func = DiffEq.ODEFunction(ode_fn!;
        mass_matrix = adapter.mass_matrix,
        sys = DiffEq.SymbolCache(syms),
        jac = sparse_jac!,
        jac_prototype = jac_prototype,
        tgrad = sparse_tgrad!,
    )
    return DiffEq.ODEProblem(ode_func, u0, tspan, nothing)
end

"""
    _default_tds_alg(prob::TDSProblem)

Default TDS integrator for the steady state the problem starts in
(`prob.init.steady_state`), with the KLU sparse direct solver.

At an [`Equilibrium`](@ref), QNDF, a quasi-constant stepsize NDF (orders up to 5):
excellent stability for stiff power system dynamics. In a [`SinusoidalSteadyState`](@ref),
FBDF limited to order 2: BDF methods above order 2 are not A-stable, and the lightly
damped high-frequency modes of a network oscillating at the system frequency
(inductances with bus capacitances) force them to tiny steps or make them fail; order
2 is A-stable, matching what such a fast-oscillating network needs.
"""
_default_tds_alg(prob::TDSProblem) = _default_tds_alg(prob.init.steady_state)
_default_tds_alg(::Equilibrium) = DiffEq.QNDF(linsolve=DiffEq.KLUFactorization())
_default_tds_alg(::SinusoidalSteadyState) =
    DiffEq.FBDF(max_order=Val(2), linsolve=DiffEq.KLUFactorization())

"""
    _merge_callbacks(cb1, cb2) -> Union{CallbackSet, Nothing}

Merge two callbacks (or callback sets) into a single CallbackSet.
Handles nothing values gracefully.
"""
function _merge_callbacks(cb1, cb2)
    if isnothing(cb1) && isnothing(cb2)
        return nothing
    elseif isnothing(cb1)
        return cb2
    elseif isnothing(cb2)
        return cb1
    else
        return DiffEq.CallbackSet(cb1, cb2)
    end
end

"""
    _default_initializealg(prob::TDSProblem) -> PowerFlowInit

Event reinitialization from the problem's `TDSConfig`: `maxiters` Newton
iterations, from a flat start if `flat_start` is set and from the pre-fault
voltages otherwise.
"""
function _default_initializealg(prob::TDSProblem)
    tds_config = get_solver_config(prob.config, :TDS, TDSConfig)
    return PowerFlowInit(maxiters = tds_config.maxiters,
                         initial_guess = tds_config.flat_start ? :flat_start : :pre_fault)
end

"""
    _simulation_callbacks(prob::TDSProblem, adapter, initializealg; abstol, jump=nothing, verbose=false)
        -> Union{CallbackSet, Nothing}

The callbacks of the simulation of `prob` with `adapter`, integrated with the absolute
tolerance `abstol`: its events, each followed by the reinitialization `initializealg`
bound to `adapter` and to `abstol` (`_bind_adapter`), and the automata of its models
(`adapter.automata`): their located callback (`build_automaton_callback`), and after
each event, once the DAE is reinitialized, a settle of the automata (`AutomatonSettle`). Where a transition changes
the algebraic equations, the bound `initializealg` reinitializes them. `jump` also
receives the field and setpoint changes of the events and the instants at which the
automata start and fire their transitions (see `_instant_begins!`).
"""
function _simulation_callbacks(prob::TDSProblem, adapter::TDSAdapter, initializealg;
                               abstol, jump = nothing, verbose::Bool = false)
    automata = adapter.automata
    # Bind adapter to initializealg (passthrough for non-PowerFlowInit algorithms)
    reinit_alg = isempty(prob.events) && automata === nothing ? nothing :
                 _bind_adapter(initializealg, adapter, jump; abstol, verbose)

    event_callbacks = nothing
    if !isempty(prob.events)
        # After an event, the network reinitialized where it changed the algebraic
        # equations (`EventReinit`), then the automata settle
        network_changed = Ref(true)
        init_alg = _event_reinit(reinit_alg, jump, network_changed)
        automata === nothing || (init_alg = AutomatonSettle(init_alg, automata, adapter, jump, prob.events, prob.sys))

        # Create re-evaluation closure for events that modify state variables.
        # This updates adapter's residual buffers with new state values, ensuring
        # consistency before DAE reinitialization or solver continuation.
        reevaluate_fn = let a = adapter
            (integrator) -> begin
                du = similar(integrator.u)
                evaluate!(du, a, integrator.u, integrator.t, ODEInterface())
            end
        end

        # Pass address_manager, adapter, and reevaluate_fn for events that need to
        # modify state variables and re-evaluate residuals (e.g., a fault clearing)
        event_callbacks = build_event_callbacks(prob.events, prob.sys, adapter.status;
            address_manager=prob.address_manager, adapter=adapter,
            reevaluate_fn=reevaluate_fn, initializealg=init_alg, jump, network_changed, automata)
    end
    automaton_callback = automata === nothing ? nothing :
        build_automaton_callback(automata, adapter, reinit_alg; jump,
                                 start_abstol = get_solver_config(prob.config, :TDS, TDSConfig).abstol,
                                 locate_extrema = get_solver_config(prob.config, :TDS, TDSConfig).locate_extrema)
    return _merge_callbacks(event_callbacks, automaton_callback)
end

"""
    solve(prob::TDSProblem, alg=_default_tds_alg(prob); initializealg, verbose=false, validate=true, recorder=nothing, kwargs...)

Solve time-domain simulation problem using data stored in `prob`.

Uses ODEProblem with mass matrix formulation (M * du/dt = f(u,t)) and sparse Jacobians.
Default solver is QNDF with KLUFactorization for efficient sparse linear solves, and
FBDF of order at most 2 for a `SinusoidalSteadyState` problem (see `_default_tds_alg`).

# Arguments
- `prob::TDSProblem` - TDS problem with initialization data
- `alg` - Solver algorithm (default: `_default_tds_alg(prob)`)
- `initializealg` - DAE initialization algorithm for events (default: `PowerFlowInit`
  with `maxiters` and `flat_start` from the problem's `TDSConfig`, converging to
  `REINIT_TOL_FRACTION` of `abstol`):
  - `PowerFlowInit()` - Flat start V=1, current theta (the default)
  - `PowerFlowInit(initial_guess=:pre_fault)` - Pre-fault V magnitudes, current theta
    (unsuitable when the angles have moved far from their pre-fault values: after a long
    fault or a pole slip it may not converge)
  - `PowerFlowInit(initial_guess=:current)` - Use current voltages (converged in the
    long-fault cases tested, where the other two did not; may converge to wrong solution
    from low voltages)
  Which guess converges depends on the case. A reinitialization that does not converge
  terminates the integration with `sol.retcode == :ConvergenceFailure`, with a warning, the
  trajectory kept up to that instant. An inconsistent initial state throws (DiffEq's check at `t0`)
  - `BrownFullBasicInit()` - DiffEq's generic initializer (may fail for severe faults)
  - `NoInit()` - Skip reinitialization (not recommended)
- `verbose` - Print progress information (default: false)
- `validate` - Validate initialization before solving (default: true). Checks that
  the DAE starts in steady state (see `validate_tds_initialization`).
- `recorder` - Optional `InjectionRecorder` for tracking per-model contributions to SharedRes
- `kwargs...` - Passed to OrdinaryDiffEq.solve(). `abstol` and `reltol` default to
  the problem's `TDSConfig` values.

# Example
```julia
# Basic solve (QNDF with KLUFactorization is default)
sol = solve(prob)

# With verbose output
sol = solve(prob; verbose=true)

# Use pre-fault voltages (good for fault-on/fault-off scenarios)
sol = solve(prob; initializealg=PowerFlowInit(initial_guess=:pre_fault))

# With injection recording (track model contributions to the bus power balance)
bus_keys = [(i,) for i in prob.sys.models.Bus.i]
monitor = InjectionMonitor(:Bus => [:p => bus_keys, :q => bus_keys])
recorder = InjectionRecorder(monitor, prob)  # Creates from TDSProblem
sol = solve(prob; recorder=recorder)
table = contributions_table(recorder, :Bus, :p, bus_keys[1])
```
"""
function solve(
    prob::TDSProblem,
    alg = _default_tds_alg(prob);
    initializealg = _default_initializealg(prob),
    verbose::Bool = false,
    validate::Bool = true,
    recorder = nothing,
    kwargs...
)
    if isempty(prob.init.dyn)
        error("TDSProblem must be initialized with dynamic states. Use setup_tds(sol_pf, dyn_models) first.")
    end

    # The solve owns its system: events modify a working copy, so the problem stays as
    # built and every solve of it starts from the same state
    prob = working_copy(prob)

    @vprint "="^60
    @vprint "Time Domain Simulation Setup"
    @vprint "  Time span: $(prob.tspan[1])s to $(prob.tspan[2])s"
    @vprint "="^60

    # Step 1: Create universal ADCache (shared infrastructure!)
    verbose && print("Creating ADCache... ")
    t_cache = @elapsed adcache = create_adcache(prob)
    @vprint "done ($(round(t_cache*1000, digits=1)) ms)"

    # Step 2: Create unified adapter (with optional recorder for type-stable hot path)
    verbose && print("Creating TDSAdapter... ")
    t_adapter = @elapsed adapter = TDSAdapter(adcache, prob; recorder=recorder)
    @vprint "done ($(round(t_adapter*1000, digits=1)) ms)"

    # Step 2b: Optional validation (reuses adapter - no duplication)
    if validate
        @vprint "Validating initialization..."
        result = validate_tds_initialization(adapter, prob.address_manager, prob.tspan[1])
        if !result.is_valid
            n_violations = length(result.non_zero_diff) + length(result.non_zero_alg)
            @printf("TDS initialization failed: %d violations (max_diff=%.2e, max_alg=%.2e)\n",
                    n_violations, result.max_diff_residual, result.max_alg_residual)
            println("  Run `print_validation_report(validate_tds_initialization(prob))` for details")
            @error "TDS initialization failed; simulation may freeze or diverge"
        end
        @vprint "  Initialization valid ✓"
    end

    # Step 3: Create ODE problem (always uses mass matrix formulation)
    verbose && print("Building ODEProblem (sparsity + coloring)... ")
    t_prob = @elapsed diffeq_prob = create_diffeq_problem(adapter, prob.tspan; verbose=verbose)
    @vprint "done ($(round(t_prob*1000, digits=1)) ms)"

    # Integrator keywords: the configured tolerances unless given, without the user
    # callback (merged with ours below)
    tds_config = get_solver_config(prob.config, :TDS, TDSConfig)
    filtered_kwargs = merge((abstol = tds_config.abstol, reltol = tds_config.reltol),
                            Base.structdiff(values(kwargs), NamedTuple{(:callback,)}))

    # Step 4: Build callbacks
    verbose && print("Building callbacks... ")
    t_callbacks = @elapsed begin
        # Event callbacks (faults, trips, etc.) and the callback of the models' automata
        simulation_callbacks = _simulation_callbacks(prob, adapter, initializealg; abstol = filtered_kwargs.abstol,
                                                     verbose)
    end
    n_events = length(prob.events)
    @vprint "done ($(round(t_callbacks*1000, digits=1)) ms, $n_events events)"

    # Recording callback (commits pending contributions at accepted steps)
    #
    # WARNING: Do NOT use DiscreteCallback here!
    # DiscreteCallback causes InitialFailure when combined with autodiff=true solvers
    # (e.g., FBDF(autodiff=true)) and event callbacks that trigger DAE reinitialization.
    # The DiscreteCallback interferes with the solver's internal state during the
    # initialization phase after events like faults.
    #
    # SavingCallback is designed for passive observation and doesn't affect solver state,
    # making it compatible with autodiff solvers.
    recording_callback = nothing
    if recorder !== nothing
        recording_callback = DiffEq.SavingCallback(
            (u, t, integrator) -> begin
                commit_contributions!(recorder, t)
                nothing  # Return nothing (we only want the side effect)
            end,
            DiffEq.SavedValues(Float64, Nothing)
        )
    end

    # Extract user callback from kwargs (if any) and merge with ours
    user_callback = get(kwargs, :callback, nothing)

    # Merge all callbacks (event + automaton + recording + user)
    callbacks = _merge_callbacks(simulation_callbacks, recording_callback)
    callbacks = _merge_callbacks(callbacks, user_callback)

    t_setup = t_cache + t_adapter + t_prob + t_callbacks
    @vprint ""
    @vprint "Setup complete ($(round(t_setup*1000, digits=1)) ms total)"
    @vprint "  Algorithm: $(nameof(typeof(alg)))"
    @vprint "="^60
    @vprint "Solving..."

    # Step 5: Solve using OrdinaryDiffEq.jl
    t_solve = @elapsed begin
        if isnothing(callbacks)
            sol = DiffEq.solve(diffeq_prob, alg; filtered_kwargs...)
        else
            sol = DiffEq.solve(diffeq_prob, alg; callback=callbacks, filtered_kwargs...)
        end
    end

    @vprint "="^60
    @vprint "TDS Complete"
    @vprint "  Solve time: $(round(t_solve, digits=3))s"
    @vprint "  Return code: $(sol.retcode)"
    @vprint "  Time points: $(length(sol.t))"
    @vprint "="^60

    # Wrap in TDSSolution for unified API (zero data copying)
    wrapped = _wrap_tds_solution(sol, prob)
    _check_solution(wrapped)
    return wrapped
end

"""
    check_solution(::Type{M}, models, sol::TDSSolution)

Diagnostics a model type adds after a solve: warnings about the solution (a device current beyond its declared limit, for
one). `models` are the instances of type `M` of the solved system. A check reads the solution and changes nothing; the default
is none.
"""
check_solution(::Type, models, sol) = nothing

# The checks of every model type of the solved system (after the solve: the integration is untouched)
function _check_solution(sol::TDSSolution)
    for models in values(sol.sys.models)
        isempty(models) && continue
        check_solution(get_base_model_type(models), models, sol)
    end
    return nothing
end

"""
    _wrap_tds_solution(ode_sol, prob::TDSProblem) -> TDSSolution

Wrap an ODESolution in TDSSolution for unified API access.

Zero data copying - stores reference to ODESolution directly.
Enables `get_series(sol, :Model, :var, key)` API.
"""
function _wrap_tds_solution(ode_sol, prob::TDSProblem)
    sys = prob.sys
    am = prob.address_manager
    n_diff = Int(max_address_in_space(am, StateVar()))
    n_alg = Int(max_address_in_space(am, AlgebVar()))

    # Convert retcode to Symbol for storage
    retcode = Symbol(ode_sol.retcode)

    return TDSSolution(
        ode_sol,      # Keep for interpolation + data access
        nothing,      # Use sol.t via accessor
        nothing,      # Use sol.u via accessor
        retcode,
        n_diff,
        n_alg,
        sys,
        am
    )
end
