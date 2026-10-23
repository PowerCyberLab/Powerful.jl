# ============================================================================
# Y-Bus TDS Solve Integration
# ============================================================================
#
# Integrates the Y-bus DAE residual (YBusResidual) with OrdinaryDiffEq.jl in
# mass-matrix form, with either network treatment. Topology events switch the
# network and make the algebraic variables consistent with the states
# (consistent_algebraics!), as is the initial state. The automata of the models run on
# this DAE with the runtime of the phasor simulation (PowerfulTDS): located transitions,
# and a settle after each event.
#
# ============================================================================

import CommonSolve: solve

using DifferentiationInterface: DifferentiationInterface
using PowerfulCore: StateVar, AlgebVar, variable_name, max_address_in_space
using PowerfulCore: TDSSolution, jacobian_sparsity_pattern, create_di_backend, ForwardDiffBackend
using PowerfulCore.DiffEq
using PowerfulCore.DiffEq: PresetTimeCallback, CallbackSet, u_modified!, NoInit
using SparseArrays: nnz

"""
    solve(prob::YBusTDSProblem[, alg]; network=SimultaneousNetwork(), maxiters=50,
          network_abstol=1e-8, verbose=false, kwargs...) -> TDSSolution

Solve the Y-bus TDS problem with OrdinaryDiffEq.jl as a mass-matrix DAE
([`YBusResidual`](@ref)).

# Arguments
- `prob::YBusTDSProblem`: Problem with Y-bus formulation (and events)
- `alg`: OrdinaryDiffEq algorithm; default by `network`:
  - `SimultaneousNetwork()`: `QNDF(linsolve=KLUFactorization())` with the sparse
    ForwardDiff Jacobian (pattern from the model structure, coloring)
  - `EliminatedNetwork()`: `QNDF()` with OrdinaryDiffEq's dense ForwardDiff Jacobian
- `network::NetworkTreatment`: [`SimultaneousNetwork`](@ref) (default) or
  [`EliminatedNetwork`](@ref)
- `maxiters`, `network_abstol`: network solve limits (`EliminatedNetwork` evaluations
  and [`consistent_algebraics!`](@ref))
- `verbose::Bool`: Print setup information
- `kwargs...`: Passed to OrdinaryDiffEq.solve (`initializealg` defaults to `NoInit()`:
  the initial state and the state after each event are made consistent by
  `consistent_algebraics!`)

Each model's own limiter or controller mode logic (e.g. the anti-windup limits of SEXS
and TGOV1) starts at its initial mode, and the solver catches the instant each limit or
mode is reached, as it does for a network event. Where a mode change affects an
algebraic variable or a controller state with zero time constant, those are re-solved
to stay consistent with it, and after any event or mode change the solver resolves every
other mode change it triggers before continuing.

# Returns
`TDSSolution` wrapping the ODESolution (`retcode` as a Symbol). Where the algebraic solve
of the initial state or after an event or a transition does not converge (a numerical failure,
reported by a warning with the time, the iterations, the residual and its worst rows), the
integration stops with `retcode = :ConvergenceFailure` (the trajectory so far, or the initial
state alone if it is the initial one).

# Example
```julia
ybus_prob = YBusTDSProblem(tds_prob; events=[BusFault(0.1, 0.2, 10, "F1")])
sol = solve(ybus_prob)                              # sparse, network in the DAE
sol = solve(ybus_prob; network=EliminatedNetwork())  # dense, network solved inside
```
"""
solve(prob::YBusTDSProblem; network::NetworkTreatment = SimultaneousNetwork(), kwargs...) =
    solve(prob, default_algorithm(network); network, kwargs...)

"""
    default_algorithm(network::NetworkTreatment)

OrdinaryDiffEq algorithm used by `solve(::YBusTDSProblem)` for `network`.
"""
default_algorithm(::SimultaneousNetwork) = DiffEq.QNDF(linsolve=DiffEq.KLUFactorization())
default_algorithm(::EliminatedNetwork) = DiffEq.QNDF()

function solve(
    prob::YBusTDSProblem,
    alg;
    network::NetworkTreatment = SimultaneousNetwork(),
    maxiters::Int = 50,
    network_abstol::Real = 1e-8,
    verbose::Bool = false,
    kwargs...
)
    tspan = prob.base_prob.tspan
    @vprint "="^60
    @vprint "Y-Bus TDS Simulation (DiffEq.jl, $(nameof(typeof(network))))"
    @vprint "  Time span: $(tspan[1])s to $(tspan[2])s"
    @vprint "="^60

    residual = YBusResidual(prob, network; maxiters, abstol=network_abstol)
    n_diff, n_total = residual.n_diff, residual.n_diff + residual.n_alg
    @vprint "  States: $(residual.n_diff) differential, $(residual.n_alg) algebraic"
    automata = _automaton_runtime(residual)
    reinit = HostReinit((u, t) -> consistent_algebraics!(u, residual, t))
    # (the tolerance the initialization is validated to, below which the start does not solve the point again)
    start_abstol = PowerfulCore.get_solver_config(prob.base_prob.config, :TDS, PowerfulTDS.TDSConfig).abstol

    # Initial state: pre-fault topology, the automata started (every evaluation reads the
    # discrete state they start with), algebraic variables consistent with the states
    _switch_to!(prob, _topology_at(prob, tspan[1]))
    u0 = zeros(Float64, n_total)
    init_tds_buffers!(view(u0, 1:n_diff), view(u0, (n_diff + 1):n_total), prob.base_prob)
    if automata !== nothing
        # (the algebraic variables are made consistent between the rounds of the start, from the modes entered)
        instant = FixedInstant(u0, Float64(tspan[1]), (u, t) -> consistent_algebraics!(u, residual, t),
                               () -> error("No rates are evaluated at the start of the automata"), false)
        start!(automata, view(u0, 1:n_diff), view(u0, (n_diff + 1):n_total), tspan[1];
               resolve = StartResolve(residual, instant, nothing, start_abstol))
        instant.failed && return TDSSolution(nothing, [tspan[1]], reshape(u0, :, 1),
                                             Symbol(DiffEq.ReturnCode.ConvergenceFailure), residual.n_diff, residual.n_alg,
                                             prob.base_prob.sys, prob.base_prob.address_manager)
    end
    consistent_algebraics!(u0, residual, tspan[1]) ||
        return TDSSolution(nothing, [tspan[1]], reshape(u0, :, 1), Symbol(DiffEq.ReturnCode.ConvergenceFailure),
                           residual.n_diff, residual.n_alg, prob.base_prob.sys, prob.base_prob.address_manager)

    ode_func = _ode_function(network, residual, u0, tspan[1], _variable_names(prob); verbose)
    diffeq_prob = DiffEq.ODEProblem(ode_func, u0, tspan, nothing)

    event_callback, event_tstops = create_event_callbacks(prob, residual, automata, reinit)
    isnothing(event_callback) || @vprint "  Events: $(length(event_tstops)) topology changes"
    automaton_callback = automata === nothing ? nothing : build_automaton_callback(automata, residual, reinit; start_abstol)

    # Merge with user-provided tstops and callback (if any)
    all_tstops = sort(unique([event_tstops; get(kwargs, :tstops, Float64[])]))
    callbacks = filter(!isnothing, Any[event_callback, automaton_callback, get(kwargs, :callback, nothing)])
    callback = isempty(callbacks) ? nothing : length(callbacks) == 1 ? only(callbacks) : CallbackSet(callbacks...)
    solve_kwargs = merge((initializealg = NoInit(),),
                         Base.structdiff(values(kwargs), NamedTuple{(:tstops, :callback)}))

    @vprint "  Algorithm: $(nameof(typeof(alg)))"
    @vprint "Solving..."
    sol = DiffEq.solve(diffeq_prob, alg; callback, tstops=all_tstops, solve_kwargs...)
    @vprint "Done! $(sol.stats.naccept) accepted steps"

    return TDSSolution(sol, nothing, nothing, Symbol(network_failure_retcode(residual, sol.retcode)), residual.n_diff,
                       residual.n_alg, prob.base_prob.sys, prob.base_prob.address_manager)
end

# The automata of the models of the residual's problem, run on its Float64 buffers with
# the runtime of PowerfulTDS (`nothing` when no model declares any)
function _automaton_runtime(r::YBusResidual)
    sys_ws = r.float.sys_ws
    return automaton_runtime(sys_ws, r.prob.base_prob.sys, sys_ws.status, r.n_diff, r.n_alg,
                             () -> sync_automaton_memory!(r), Vector(diag(r.prob.M)))
end

# Variable names (plot labels) in the DAE address layout
function _variable_names(prob::YBusTDSProblem)
    am = prob.base_prob.address_manager
    n_diff = Int(max_address_in_space(am, StateVar()))
    n_alg = Int(max_address_in_space(am, AlgebVar()))
    return [[variable_name(am, StateVar(), UInt(a)) for a in 1:n_diff];
            [variable_name(am, AlgebVar(), UInt(a)) for a in 1:n_alg]]
end

# Eliminated network: OrdinaryDiffEq computes the dense Jacobian with ForwardDiff
function _ode_function(::EliminatedNetwork, residual::YBusResidual, u0, t0, syms; verbose)
    return DiffEq.ODEFunction(residual; mass_matrix=residual.prob.M, sys=DiffEq.SymbolCache(syms))
end

# The pattern of the Jacobian of the Y-bus DAE with the network simultaneous, from the model
# structure (the rates of the automata of a fixed-step integration are evaluated along its
# trajectory, `FixedStepAutomata`)
function _jacobian_pattern(residual::YBusResidual)
    bp = residual.prob.base_prob
    am = bp.address_manager
    return jacobian_sparsity_pattern(am, bp.sys, PowerInjection(); include_differential=true,
                                     models=filter(m -> haskey(am.assignments, m), keys(bp.sys.models)))
end

# Simultaneous network: sparse ForwardDiff Jacobian with coloring. The pattern of the
# model structure (`jacobian_sparsity_pattern` of the models with DAE addresses, in the
# TDS system's address layout) contains the Y-bus residual's: machine and controller
# residuals are the same model residuals, and a Bus p/q row depends on the voltages of
# the buses connected to it by a branch and on the variables of the devices at the bus,
# which includes every entry of Y and every generator internal voltage and injection.
# Fault models added by events have no addresses: a fault changes only the diagonal
# of Y, already in the pattern.
function _ode_function(::SimultaneousNetwork, residual::YBusResidual, u0, t0, syms; verbose)
    pattern = _jacobian_pattern(residual)
    closure = YBusJacobianClosure(residual, Ref(Float64(t0)))
    backend = create_di_backend(ForwardDiffBackend(), pattern)
    du = similar(u0)
    prep = DifferentiationInterface.prepare_jacobian(closure, du, backend, u0)
    jac = YBusSparseJacobian(closure, prep, backend, du)
    jac_prototype = Float64.(pattern)
    @vprint "  Jacobian: $(size(jac_prototype, 1)) x $(size(jac_prototype, 2)), nnz=$(nnz(jac_prototype))"
    return DiffEq.ODEFunction(residual; mass_matrix=residual.prob.M, jac, jac_prototype,
                              sys=DiffEq.SymbolCache(syms))
end

# du = F(u) at the time held in `t`, for DifferentiationInterface
struct YBusJacobianClosure{R<:YBusResidual}
    residual::R
    t::Base.RefValue{Float64}
end
(c::YBusJacobianClosure)(du, u) = c.residual(du, u, nothing, c.t[])

# OrdinaryDiffEq Jacobian J = ∂F/∂u at (u, t): sparse ForwardDiff with coloring
struct YBusSparseJacobian{C<:YBusJacobianClosure, PR, B}
    closure::C
    prep::PR
    backend::B
    du::Vector{Float64}
end
function (j::YBusSparseJacobian)(J, u, p, t)
    j.closure.t[] = t
    DifferentiationInterface.jacobian!(j.closure, j.du, J, j.prep, j.backend, u)
    nothing
end

"""
    rate_dual_type(r::YBusResidual{SimultaneousNetwork}, width) -> Type

The Dual type of the Jacobian of the Y-bus DAE with the network simultaneous
(`YBusJacobianClosure`'s tag), in whose buffers of `r` the rates of the values of the automata
are evaluated (PowerfulTDS `StateRates`), as the phasor simulation does with its Jacobian's:
the integration builds them at its first Jacobian, before any step and so before any located
event, and the rates are evaluated at events, never within an evaluation of the Jacobian, so
no Dual of the one is a perturbation of the other's. With the network eliminated the
Jacobian is OrdinaryDiffEq's dense one, taken one column per evaluation (one partial), and the
rates keep Duals of their own tag and `width` partials: `F_ZZ` is dense there too.
"""
function rate_dual_type(r::YBusResidual{SimultaneousNetwork}, ::Int)
    for T in keys(r.dual)
        T <: (ForwardDiff.Dual{ForwardDiff.Tag{F, Float64}, Float64} where F <: YBusJacobianClosure) && return T
    end
    throw(ErrorException(
        "The rates of the automata are evaluated in the Dual buffers of the integration's Jacobian, which this " *
        "integration has not built (no Jacobian of `YBusJacobianClosure` was evaluated before the located event " *
        "that needs them)."))
end

# Topology active at time t (pre-fault before the first event)
_topology_at(prob::YBusTDSProblem, t) =
    isnothing(prob.event_schedule) ? TopologyKey(:prefault) : get_topology_at(prob.event_schedule, t)

_switch_to!(prob::YBusTDSProblem, key::TopologyKey) =
    switch_topology!(prob.network_cache.topology_cache, key)

"""
    create_event_callbacks(prob::YBusTDSProblem, residual::YBusResidual, automata, reinit)
        -> (callback, tstops)

Create DiffEq callbacks for topology switching at event times.

At each event time the topology of the schedule (`get_topology_at`, as in the
trapezoidal solver) becomes active and the algebraic variables are made consistent
with the states under it ([`consistent_algebraics!`](@ref)); the automata (`automata`,
an `AutomatonRuntime` or `nothing`) then settle at that state (`settle_instant!`,
reinitializing with `reinit` where a transition changes the algebraic equations).

# Returns
- `(callback, tstops)` tuple if events exist
- `(nothing, Float64[])` if no events
"""
function create_event_callbacks(prob::YBusTDSProblem{T}, residual::YBusResidual, automata, reinit) where T
    event_schedule = prob.event_schedule
    if isnothing(event_schedule) || isempty(event_schedule.times)
        return (nothing, T[])
    end

    applied = Ref(0)
    models = residual.float.sys_ws.models_by_name
    affect! = function (integrator)
        t = integrator.t
        apply_changes!(prob, applied, t, automata, models; residual)
        new_key = get_topology_at(event_schedule, t)
        @debug "DiffEq callback" t new_topology=new_key.name before=prob.network_cache.topology_cache.current[].name
        _switch_to!(prob, new_key)
        if !consistent_algebraics!(integrator.u, residual, t)
            DiffEq.terminate!(integrator, DiffEq.ReturnCode.ConvergenceFailure)
            return nothing
        end
        automata === nothing || settle_instant!(automata, residual, integrator, reinit; event = true)
        u_modified!(integrator, true)
    end

    # (the closure behind an `ErasedCall`: `PresetTimeCallback` builds its callback and function from it)
    callback = PresetTimeCallback(event_schedule.times, erased(affect!, Val(1)))
    return (callback, copy(event_schedule.times))
end
