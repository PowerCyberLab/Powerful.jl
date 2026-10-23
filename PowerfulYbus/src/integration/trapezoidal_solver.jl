# ============================================================================
# Trapezoidal Solver for Y-Bus TDS
# ============================================================================
#
# Simple fixed-step trapezoidal integration:
# - Predictor: x_{k+1} = x_k + h * dx_k
# - Corrector: x_{k+1} = x_k + h/2 * (dx_k + dx_{k+1})
#
# Features:
# - Time segments with dt adjustment to hit event times exactly
# - Topology switching at segment boundaries (O(1) pre-computed Y matrix swap)
# - Much faster than DiffEq.jl for simple simulations

export solve_trapezoidal

using LinearAlgebra
using SparseArrays
using PowerfulCore: TDSSolution


"""
    FixedStepHost

The Y-bus DAE with the network simultaneous (`YBusResidual`) as the host of the automata of a
fixed-step integration (`FixedStepAutomata`): its residuals for the check of the algebraic
equations and its Dual workspaces for the rates of the values their transitions hold at zero
by continuity. Built when first needed, which a run whose transitions neither change the
algebraic equations nor hold values at zero never does.
"""
mutable struct FixedStepHost{P}
    prob::P
    maxiters::Int
    abstol::Float64
    residual::Union{Nothing, YBusResidual}
end

FixedStepHost(prob, maxiters, abstol) = FixedStepHost(prob, maxiters, Float64(abstol), nothing)

function _residual(h::FixedStepHost)
    h.residual === nothing && (h.residual = YBusResidual(h.prob, SimultaneousNetwork(); maxiters = h.maxiters,
                                                         abstol = h.abstol))
    return h.residual
end

dae_residuals!(du, h::FixedStepHost, u, t) = dae_residuals!(du, _residual(h), u, t)
dual_workspace(h::FixedStepHost, ::Type{T}) where T = dual_workspace(_residual(h), T)
algebraic_equations_changed!(::FixedStepHost) = nothing

"""
    solve_trapezoidal(prob::YBusTDSProblem; dt=0.01, maxiters=50, abstol=1e-8)

Solve a Y-bus TDS problem with a fixed-step trapezoidal integrator: a simple,
predictable alternative to the adaptive-step OrdinaryDiffEq.jl solvers, often faster
for straightforward simulations.

Uses time segments with a slightly adjusted `dt` within each one, so that the
integration lands exactly on every event time rather than overshooting it and
interpolating back. The network topology switches automatically at each event time,
using the pre-computed, already-factored Y matrix for the new topology.

# Arguments
- `prob`: YBusTDSProblem to solve
- `dt`: Nominal time step (default 0.01s) - adjusted per segment to divide evenly
- `maxiters`: Max Newton steps of each network solve (default 50)
- `abstol`: Network solve tolerance on the max KCL current mismatch, in pu
  (default 1e-8)
- `verbose`: Print progress info (default false)

# Failure
If a solve of the network or of the algebraic rows does not converge within `maxiters` (a
numerical failure of a valid problem, e.g. a fault the network cannot carry), the simulation
stops with a warning (the time, the iterations, the residual and its worst rows) and returns the
trajectory up to the last completed step with `retcode = :ConvergenceFailure` (SciML's
`ReturnCode.ConvergenceFailure`). Errors in the models or the code (a block whose rows determine
no variable, a missing field) throw.

# Time Segment Behavior
When events are present, the time span is divided into segments at event times.
Each segment has a uniform dt adjusted to divide the interval evenly:
- If nominal dt=0.01 and segment is 0.095s, we use dt=0.095/9 ≈ 0.01056
- This ensures we land exactly on event times without tiny final steps

Segments smaller than 0.5×nominal_dt are merged with adjacent segments.
Warnings are issued if dt varies more than 10% from nominal.

# Limiters and Controller Modes
Each model's own limiter or controller mode logic (e.g. the anti-windup limits of SEXS
and TGOV1) starts at its initial mode. Where a limit or mode condition is reached inside
a step, that step is cut short at the instant it is reached (found to within
`1e-10 max(1, |t_end|)` seconds), the algebraic variables are re-solved there for
consistency with the new mode, and the step then continues to its original end time.
The same re-solve happens after every topology switch. The output keeps the time grid
of the segments (not the shortened sub-steps).
"""
function solve_trapezoidal(
    prob::YBusTDSProblem{T};
    dt::T = T(0.01),
    maxiters::Int = 50,
    abstol::T = T(1e-8),
    verbose::Bool = false
) where T

    # Extract problem info
    am = prob.base_prob.address_manager
    sys = prob.base_prob.sys
    tspan = prob.base_prob.tspan
    event_schedule = prob.event_schedule

    n_diff = Int(max_address_in_space(am, StateVar()))
    n_alg = Int(max_address_in_space(am, AlgebVar()))
    n_total = n_diff + n_alg

    # Build time segments with dt adjustment
    event_times = isnothing(event_schedule) ? T[] : event_schedule.times
    segments = build_time_segments(tspan, event_times, dt)

    # Calculate total steps for output allocation
    total_steps = sum(seg.n_steps for seg in segments) + 1

    verbose && println("Trapezoidal solver: $(length(segments)) segments, " *
                       "$(length(event_times)) events, $(total_steps) total steps")

    # Allocate solution storage
    t_out = Vector{T}(undef, total_steps)
    u_out = zeros(T, n_total, total_steps)

    # Initialize state vector from prob.init
    u = zeros(T, n_total)
    init_tds_buffers!(view(u, 1:n_diff), view(u, (n_diff + 1):n_total), prob.base_prob)

    # Store initial condition
    t_out[1] = tspan[1]
    u_out[:, 1] .= u

    # Allocate work arrays
    du = zeros(T, n_total)

    # Views for differential and algebraic parts
    x = view(u, 1:n_diff)
    f = view(du, 1:n_diff)
    f_pred = zeros(T, n_diff)

    # Mass matrix inverse for DAE: dx/dt = f * M_inv (pre-computed in prob)
    M_inv_diag = prob.M_inv_diag

    # Track current topology
    topology_cache = prob.network_cache.topology_cache
    current_topology_key = isnothing(event_schedule) ?
        TopologyKey(:prefault) : event_schedule.prefault_key

    # Ensure we start with prefault topology
    if haskey(topology_cache.topologies, current_topology_key)
        switch_topology!(topology_cache, current_topology_key)
    end

    # The automata of the problem's buffers, started again at the initial state; the last
    # network solution of a step or algebraic solve they run, to report a failure
    y = view(u, (n_diff+1):n_total)
    automata = prob.automata
    last_net = Ref{Any}(nothing)
    last_t = Ref(tspan[1])
    step! = (u, t, h) -> (net = _trapezoidal_step!(u, du, f_pred, prob, t, h, maxiters, abstol);
                          last_net[] = net; last_t[] = t + h; net.converged)
    reinit! = (u, t) -> (net = solve_algebraic!(u, du, prob, t, maxiters, abstol);
                         last_net[] = net; last_t[] = t; net.converged)
    time_tol = 1e-10 * max(1.0, abs(tspan[2]))
    # The automata run on the Y-bus DAE with the network simultaneous (its residuals for the
    # check of the algebraic equations, and the rates of the values their transitions hold at
    # zero, along its trajectory), whose Dual workspaces get their memory from the runtime
    work = nothing
    host = nothing
    if automata !== nothing
        host = FixedStepHost(prob, maxiters, abstol)
        prob.memory_sync[] = () -> host.residual === nothing || sync_automaton_memory!(host.residual)
        work = FixedStepAutomata(automata, host, () -> sparse(Float64.(_jacobian_pattern(_residual(host)))), reinit!)
        # (the algebraic variables are solved again between the rounds of the start, from the modes entered)
        instant = work.instant
        instant.u, instant.t, instant.failed = u, tspan[1], false
        start!(automata, x, y, tspan[1]; resolve = StartResolve(host, instant, nothing,
                                                                  PowerfulCore.get_solver_config(prob.base_prob.config, :TDS, PowerfulTDS.TDSConfig).abstol))
        instant.failed && return _network_not_converged(prob, last_net[], tspan[1], t_out, u_out, 1)
    end

    # The algebraic variables at the initial state, and its derivatives (~0 at equilibrium)
    net = solve_algebraic!(u, du, prob, tspan[1], maxiters, abstol)
    net.converged || return _network_not_converged(prob, net, tspan[1], t_out, u_out, 1)
    if verbose
        println("[Trapezoidal solver] Initial derivatives (should be ~0 at equilibrium):")
        max_idx = argmax(abs.(f))
        println("  f_max = ", maximum(abs, f), " at index ", max_idx)
        println("  f[1:10] (GENCLS delta/omega) = ", f[1:min(10,n_diff)])
        if n_diff > 10
            println("  f[11:20] (TGOV1?) = ", f[11:min(20,n_diff)])
        end
        if n_diff > 20
            println("  f[21:30] (SEXS?) = ", f[21:min(30,n_diff)])
        end
    end

    # Main loop over segments
    k_out = 1  # Output index (starts at 1, already stored initial)

    changes_applied = Ref(0)
    for seg in segments
        segment_dt = seg.dt

        # The events that change a model's data are applied (and the automata trace what they
        # read again) where the segment starts, before the topology switches
        changed = apply_changes!(prob, changes_applied, seg.t_start, automata, prob.sys_ws.models_by_name;
                                 residual = host === nothing ? nothing : host.residual)

        # Check for topology change at segment start
        switched = false
        if !isnothing(event_schedule)
            new_key = get_topology_at(event_schedule, seg.t_start)
            if new_key != current_topology_key
                switch_topology!(topology_cache, new_key)
                current_topology_key = new_key
                switched = true
                verbose && println("  t=$(seg.t_start): switched to topology $(new_key.name)")
            end
        end
        # The automata settle at the state an event makes consistent
        if (switched || changed) && automata !== nothing
            (reinit!(u, seg.t_start) && settle_fixed!(work, u, seg.t_start; event = true)) ||
                return _network_not_converged(prob, last_net[], last_t[], t_out, u_out, k_out)
        end

        # Integrate this segment
        for step in 1:seg.n_steps
            t = seg.t_start + (step - 1) * segment_dt
            t_next = t + segment_dt

            # One trapezoidal step; with automata, cut where a condition changes its truth
            if automata === nothing
                net = _trapezoidal_step!(u, du, f_pred, prob, t, segment_dt, maxiters, abstol)
                net.converged || return _network_not_converged(prob, net, t_next, t_out, u_out, k_out)
            else
                step_automata!(work, step!, u, t, segment_dt; time_tol) ||
                    return _network_not_converged(prob, last_net[], last_t[], t_out, u_out, k_out)
            end

            # Check for NaN - fail fast with useful error
            if any(isnan, u)
                nan_idx = findfirst(isnan, u)
                error("NaN detected in state vector at t=$t_next, index $nan_idx. " *
                      "Check model equations and initial conditions.")
            end

            # Store solution
            k_out += 1
            t_out[k_out] = t_next
            u_out[:, k_out] .= u
        end
    end

    return TDSSolution(
        nothing,    # No ODESolution for Ybus
        t_out,      # Time vector
        u_out,      # State matrix
        :Success,
        n_diff,
        n_alg,
        sys,        # SystemModel reference
        am          # AddressManager reference
    )
end

"""
    _trapezoidal_step!(u, du, f_pred, prob, t, h, maxiters, abstol) -> NetworkSolution

One step of length `h` from `t` of the trapezoidal integration of `u = [x; y]`. Every solve is
`solve_algebraic!`, which leaves the derivatives `f` at its result in `du`: the derivatives at
`t` (scaled by the inverse mass), the Euler predictor of `x` and the solve at `t + h` with the
derivatives there, the trapezoidal corrector of `x`, and the solve again. Returns the network
solution of the last solve, or of the first that did not converge (the step stops there).
"""
function _trapezoidal_step!(u::Vector{T}, du::Vector{T}, f_pred::Vector{T}, prob::YBusTDSProblem{T}, t, h,
                            maxiters::Int, abstol::T) where T
    n_diff = length(f_pred)
    x = view(u, 1:n_diff)
    f = view(du, 1:n_diff)
    M_inv_diag = prob.M_inv_diag
    t_next = t + h

    # Derivatives at the current state, scaled by M_inv. The algebraic solve leaves the
    # residuals of the DAE at its result in `du`: with `y` consistent (it is, from the step
    # before) it takes no step and `f` is that of the state.
    # DAE: M * dx/dt = f  →  dx/dt = f * M_inv
    net = solve_algebraic!(u, du, prob, t, maxiters, abstol)
    net.converged || return net
    @. f = f * M_inv_diag

    # Save f_k for trapezoidal corrector
    @. f_pred = f

    # Predictor: x_pred = x + dt * f_k (Euler step for differential states only)
    @. x = x + h * f

    # Solve algebraic equations to get consistent y for the predicted x; `f` is then the
    # derivatives at the predicted state
    net = solve_algebraic!(u, du, prob, t_next, maxiters, abstol)
    net.converged || return net
    @. f = f * M_inv_diag

    # Corrector: x = x + dt/2 * (f_{k+1} - f_k)
    @. x = x + h / 2 * (f - f_pred)

    # Re-solve algebraic equations for corrected x
    return solve_algebraic!(u, du, prob, t_next, maxiters, abstol)
end

"""
    _network_not_converged(prob, net_sol, t, t_out, u_out, k_out) -> TDSSolution

Warn about a network solve that did not converge at `t` and return the
trajectory up to output index `k_out` (the last completed step) with
`retcode = :ConvergenceFailure`.
"""
function _network_not_converged(prob, net_sol::NetworkSolution, t, t_out, u_out, k_out::Int)
    @warn network_failure_message(net_sol, t) * "; stopping the simulation"
    am = prob.base_prob.address_manager
    return TDSSolution(
        nothing,
        t_out[1:k_out],
        u_out[:, 1:k_out],
        Symbol(ReturnCode.ConvergenceFailure),
        Int(max_address_in_space(am, StateVar())),
        Int(max_address_in_space(am, AlgebVar())),
        prob.base_prob.sys,
        am
    )
end

"""
    solve_algebraic!(u, du, prob, t, maxiters, abstol; iterations_limit = maxiters) -> NetworkSolution

Solve the algebraic equations g(x, y) = 0 for y given x, the differential states frozen, by the
Newton iteration on the models' algebraic rows and the network equations together
(`AlgebraicNewton`, on the Newton engine of PowerfulCore): the bus voltages (`v` and `theta` of every Bus) are unknowns beside the
models' algebraic variables, the network's KCL equations are rows beside theirs, and the
sparse Jacobian of the system couples them. The iteration runs until the largest residual of
the algebraic rows (the models' and the network's) is below `abstol`, and leaves the
residuals of the DAE at the result in `du`.

The Jacobian is kept from one solve to the next until the equations change by something other
than the state (a topology, the data of a model, a mode of an automaton) or an iteration
with it contracts the residual too little. Where the iteration does not converge in
`iterations_limit` iterations (the limit of the network solve `maxiters` by default), or finds no
decrease of the residual along the step of a Jacobian taken at the iterate, the returned
solution is not `converged` and its `detail` names the time, the iterations, the residual and its
worst rows: the callers stop with a failing retcode, as a solve that stopped and went on would
give results that do not satisfy the equations. A block whose rows determine no variable (a
singular Jacobian) throws `ArgumentError`.

The returned `NetworkSolution` is that of the buffers at the result.
"""
function solve_algebraic!(
    u::Vector{T},
    du::Vector{T},
    prob::YBusTDSProblem{T},
    t::T,
    maxiters::Int,
    abstol::T;
    iterations_limit::Int = maxiters
) where T
    return algebraic_solve!(prob.algebraic, prob, u, du, t, abstol, iterations_limit)
end
