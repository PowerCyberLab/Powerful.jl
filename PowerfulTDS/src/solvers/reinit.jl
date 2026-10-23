"""
Algebraic Reinitialization for TDS

This module provides algebraic variable reinitialization after topology changes
(faults, line trips, etc.) during time-domain simulation.

Key components:
- `PowerFlowInit`: DAE initialization algorithm for OrdinaryDiffEq callbacks
- `reinitialize_algebraic!`: Newton solver for the rows without mass (g(x,y)=0, and the
  rows of the states whose mass is zero)

The reinitialization solves network equations with the states with mass fixed,
which is more robust than generic DAE initializers for power systems.
"""

using SparseArrays: SparseMatrixCSC, sparse, nonzeros, nzrange
using KLU: KLU
using LinearAlgebra: SingularException
using DifferentiationInterface
using DiffEqBase: DAEInitializationAlgorithm
using PowerfulCore: PrecomputedSparsityDetector, @vprint, AlgebVar, AlgebRes, StateRes, get_addresses,
    build_reverse_map!, lookup_address,
    get_float_ws, get_address_manager, get_backend, prepare_sparse_jacobian,
    ForwardDiffBackend, EnzymeBackend, get_tds_dual_buffers,
    NewtonPolicy, NewtonBuffers, BUS_VOLTAGE_STEP_LIMIT, KLULinearSolver, newton_solve!, reset_factorization!, jacobian_sparsity_pattern,
    factorize!, linear_solve!, StateVar,
    get_system, get_formulation

# Import _initialize_dae! from OrdinaryDiffEqCore to extend it
import OrdinaryDiffEqCore: _initialize_dae!

# ============================================================================
# Backend-Specific Algebraic Variable Access
# ============================================================================
#
# ForwardDiff and Enzyme have different memory layout requirements for autodiff:
#
# - **ForwardDiff**: Works with SubArrays (views). Uses Dual number propagation
#   where derivatives are tracked via type system. Views into the integrator's
#   u vector work correctly because Dual types preserve memory addresses.
#
# - **Enzyme**: Requires contiguous Vector (not views). Enzyme performs AD at
#   the LLVM level using shadow arrays for derivatives. SubArrays introduce
#   pointer indirection that breaks Enzyme's memory layout analysis.
#   See: https://enzyme.mit.edu/julia/stable/faq/#faq-views-subarrays
#
# Performance Impact:
# - ForwardDiff: Zero-copy (view into u)
# - Enzyme: ~100ns copy in + ~100ns copy out for typical n_alg=1000
# - This overhead is negligible compared to Jacobian evaluation (~1-100ms)
#
# The `needs_copyback` flag ensures Enzyme results are written back to u.
# ============================================================================

"""
    _get_y_working_buffer(u, adapter, ::ForwardDiffBackend) -> (y_buf, needs_copyback)

ForwardDiff: Use a view into u directly (zero-copy).
Returns (view, false) - no copy-back needed since view modifies u in place.

ForwardDiff's Dual number propagation works correctly with SubArrays because
derivatives are tracked via the type system, not memory layout.
"""
@inline function _get_y_working_buffer(u::AbstractVector, adapter::TDSAdapter, ::ForwardDiffBackend)
    n_diff = adapter.n_diff
    y_view = @view u[n_diff+1:end]
    return (y_view, false)
end

"""
    _get_y_working_buffer(u, adapter, ::EnzymeBackend) -> (y_buf, needs_copyback)

Enzyme: Copy to adapter.y_buf (Enzyme requires contiguous Vector, not SubArray).
Returns (y_buf, true) - must copy back to u after convergence.

Enzyme performs AD at the LLVM level using shadow arrays. SubArrays introduce
pointer indirection that breaks Enzyme's static analysis of memory locations.
The copy overhead (~100ns for n_alg=1000) is negligible vs Jacobian cost.
"""
@inline function _get_y_working_buffer(u::AbstractVector, adapter::TDSAdapter, ::EnzymeBackend)
    n_diff = adapter.n_diff
    n_alg = adapter.n_alg
    y_buf = adapter.y_buf
    @inbounds for i in 1:n_alg
        y_buf[i] = u[n_diff + i]
    end
    return (y_buf, true)
end

"""
    _copyback_y!(u, y_buf, adapter, needs_copyback)

Copy y_buf back to u if needed (only for Enzyme backend).

This ensures Enzyme's computed values are written back to the integrator's u vector.
ForwardDiff doesn't need this because it operates on a view that modifies u in place.
"""
@inline function _copyback_y!(u::AbstractVector, y_buf, adapter::TDSAdapter, needs_copyback::Bool)
    if needs_copyback
        n_diff = adapter.n_diff
        n_alg = adapter.n_alg
        @inbounds for i in 1:n_alg
            u[n_diff + i] = y_buf[i]
        end
    end
    nothing
end

# ============================================================================
# Constants
# ============================================================================

# DEFAULT_REINIT_MAXITERS is defined in tds.jl, next to the TDSConfig field it defaults

"""Default convergence tolerance of `reinitialize_algebraic!` called on its own."""
const DEFAULT_REINIT_TOL = 1e-8

"""
    REINIT_TOL_FRACTION

The reinitialization after an event stops at this fraction of the integration's absolute
tolerance (`PowerFlowInit` without an `abstol` of its own): the fraction κ = 1/100 at
which the integrator's own Newton iterations stop (OrdinaryDiffEq's `NLNewton`), so that
the state after an event is as consistent as the integrator leaves its steps. It bounds
the current mismatch of the bus rows (`BusPowerRows`), so the power residual that DiffEq
checks against `abstol` after the event (`CheckInit`, root mean square) stays below it
for bus voltages up to 100 pu. With the default `abstol = 1e-6`, 1e-8.
"""
const REINIT_TOL_FRACTION = 1 / 100

"""Default Jacobian update frequency for Modified Newton method.

Setting this to 1 gives standard Newton (Jacobian every iteration).
Higher values use Modified Newton (reuse Jacobian for multiple iterations).
Typical values: 3-5 for good balance between convergence and speed.
"""
const DEFAULT_JAC_UPDATE_FREQ = 5

# ============================================================================
# PowerFlowInit - DAE Initialization Algorithm (Public API)
# ============================================================================

"""
    PowerFlowInit <: DAEInitializationAlgorithm

Re-solve the network after an event (a fault, a breaker trip, a setpoint change) while
holding the dynamic states (rotor angle, speed, flux, controller integrators — any state
with a time constant) fixed at their value just before the event. This is the standard way
to bring a dynamic simulation back to a consistent operating point after a discrete change
in the network or in a limiter's mode, without perturbing the machine and controller states
that evolved smoothly through it.

Internally this solves the power-flow-like algebraic equations (and any state whose time
constant is zero, which the solver treats like an algebraic variable) with Newton
iteration, which is more robust for power system network equations than generic DAE
initializers.

# Initial Guess Options
Pass `initial_guess` to control the starting point of the Newton iteration:
- `:flat_start` (the default) — start from V = 1 pu, current bus angles. Converges in most
  short-fault cases.
- `:pre_fault` — start from the voltage magnitudes just before the event, current bus
  angles. Good for temporary, short faults; a poor guess once the angles have swung far
  from their pre-event values (a long fault or a pole slip), since the pre-event
  magnitudes no longer match the post-event angles.
- `:current` — start from whatever voltages the integrator currently holds. Tends to
  converge after long or severe faults where the other two guesses fail, but a Newton
  iteration started from depressed voltages can also converge to a spurious low-voltage
  solution.

No single guess is the most robust for every case; if a reinitialization fails to converge
with the default, retry with one of the others. (Example: ieee14 with GENROU, a 2 s fault
at bus 3, `xf = 0.003`: `:flat_start` and `:pre_fault` fail to converge at clearing, and
`:current` succeeds. Short faults at buses 4 and 10, 0.1 s, converge with every guess.)

This initial guess applies only to reinitialization after a network event. When a limiter
or controller changes mode on its own (for example a governor valve reaching its travel
limit), the solver instead resumes from the state the integrator already holds, since that
state was already consistent just before the mode change.

# Reinitialization Failure
A reinitialization that does not converge within `maxiters`, or whose result does not
satisfy the DAE at that instant, is treated as a numerical failure of the simulation (not a
code error): it warns with the time, iteration count, residual and worst-offending rows,
and the simulation stops with `sol.retcode == :ConvergenceFailure`, keeping the trajectory
up to that point. An inconsistent *initial* state at the start of a simulation (`t0`)
still throws, since that indicates a model or setup error rather than a physical outcome.

# Usage
```julia
sol = solve(prob, FBDF(); initializealg=PowerFlowInit())                       # flat start (default)
sol = solve(prob, FBDF(); initializealg=PowerFlowInit(initial_guess=:pre_fault))
sol = solve(prob, FBDF(); initializealg=PowerFlowInit(initial_guess=:current))
```

# Extended help
- `abstol`: convergence tolerance on the current mismatch; `nothing` (the default) ties it
  to the simulation's own absolute tolerance once bound to a solve, 1/100th of it.
- `maxiters`: maximum Newton iterations (default: $(DEFAULT_REINIT_MAXITERS)).
- `adapter`, `pre_fault_v`: bound automatically by `solve()`; not set by users.
"""
struct PowerFlowInit <: DAEInitializationAlgorithm
    abstol::Union{Nothing, Float64}
    maxiters::Int
    initial_guess::Symbol  # :flat_start, :pre_fault, or :current
    adapter::Union{Nothing, TDSAdapter}
    pre_fault_v::Union{Nothing, Vector{Float64}}  # Stored V magnitudes for :pre_fault
    verbose::Bool  # Print progress during reinitialization
end

"""
    PowerFlowInit(; abstol, maxiters, initial_guess) -> PowerFlowInit

Configure how the network is re-solved after an event. See [`PowerFlowInit`](@ref) for
the initial-guess choices and what happens when reinitialization fails.

# Arguments
- `abstol`: convergence tolerance on the current mismatch (default: `nothing`, which ties
  it to 1/100th of the simulation's own absolute tolerance — 1e-8 at the default
  `abstol = 1e-6`).
- `maxiters`: maximum Newton iterations (default: $(DEFAULT_REINIT_MAXITERS)).
- `initial_guess`: `:flat_start` (default), `:pre_fault`, or `:current`.

# Example
```julia
sol = solve(prob; initializealg=PowerFlowInit())                        # flat start (default)
sol = solve(prob; initializealg=PowerFlowInit(initial_guess=:pre_fault))
sol = solve(prob; initializealg=PowerFlowInit(initial_guess=:current))
sol = solve(prob; initializealg=PowerFlowInit(abstol=1e-10, maxiters=30))
```
"""
function PowerFlowInit(;
    abstol::Union{Nothing, Real} = nothing,
    maxiters::Int = DEFAULT_REINIT_MAXITERS,
    initial_guess::Symbol = :flat_start
)
    if initial_guess ∉ (:flat_start, :pre_fault, :current)
        error("PowerFlowInit: initial_guess must be :flat_start, :pre_fault, or :current")
    end
    PowerFlowInit(abstol === nothing ? nothing : Float64(abstol), maxiters, initial_guess, nothing, nothing, false)
end

export PowerFlowInit

# ============================================================================
# Bus Rows as Current Mismatch
# ============================================================================

"""
    BusPowerRows

Algebraic rows of the phasor bus power balance (`Bus` `p` and `q` residuals of the
PowerInjection formulation) and the address of each bus voltage magnitude `v`.

Every term of a bus power balance carries the bus voltage magnitude
(`P_i = V_i Σ_j V_j (G_ij cos θ_ij + B_ij sin θ_ij) + ...`, and load and shunt powers
vanish with `V_i`), so `V_i = 0` satisfies both rows of bus `i` whatever current flows
into it, and there `∂P_i/∂θ = ∂Q_i/∂θ = 0`. Newton started near that point (a faulted
bus after clearing, whose fault-on voltage is near zero) converges to it: the "cleared"
fault persists and the integrator then stalls on the singular Jacobian.
[`reinitialize_algebraic!`](@ref) therefore solves these rows divided by `V_i`, the
current mismatch, which has the same roots with `V_i ≠ 0` and none at `V_i = 0`.
"""
struct BusPowerRows
    p::Vector{Int}
    q::Vector{Int}
    v::Vector{Int}
end

"""
    bus_power_rows(am::AddressManager) -> Union{BusPowerRows, Nothing}

The bus power balance rows of `am`, or `nothing` in a formulation without them
(CurrentInjection balances currents, which do not vanish with the voltage).
"""
function bus_power_rows(am)
    haskey(am.addresses, (AlgebRes(), :Bus, :p)) || return nothing
    return BusPowerRows(Int.(get_addresses(am, AlgebRes(), :Bus, :p)),
                        Int.(get_addresses(am, AlgebRes(), :Bus, :q)),
                        Int.(get_addresses(am, AlgebVar(), :Bus, :v)))
end

"""
    to_current_mismatch!(g, y, rows, bus_status)

Divide the power balance rows of each in-service bus in `g` by the bus voltage
magnitude in `y` (see [`BusPowerRows`](@ref)). Rows of out-of-service buses hold
their variables and are left as they are.
"""
@inline function to_current_mismatch!(g, y, rows::BusPowerRows, bus_status)
    @inbounds for k in eachindex(rows.v)
        bus_status[k] || continue
        inv_v = inv(y[rows.v[k]])
        g[rows.p[k]] *= inv_v
        g[rows.q[k]] *= inv_v
    end
    return g
end
to_current_mismatch!(g, y, ::Nothing, bus_status) = g

# ============================================================================
# Core Newton Solver
# ============================================================================

"""
    reinitialize_algebraic!(u::AbstractVector, adapter::TDSAdapter, t; kwargs...)

Reinitialize the variables without mass of `u`, the state at time `t`, after a topology
change (fault, trip, etc) or a transition of an automaton: the
algebraic variables `y`, and the differential states whose mass is zero (a zero inertia,
a bus without capacitance), which the DAE `M u̇ = F(u)` constrains as algebraic.

Solves the rows without mass of the DAE (`g(x, y) = 0`, and the rows of the states
without mass, as the integration evaluates them: a row the reduction of state-only
constraints replaces by its derivative in that form) for those variables, with the
states with mass fixed, using Modified Newton:

    z_{k+1} = z_k - J^{-1} * r(z_k)

where the Jacobian J = ∂r/∂z is only recomputed every `jac_update_freq` iterations.
Without states without mass, `z = y` and `r = g`.
The bus power balance rows are solved as current mismatch (divided by the bus
voltage magnitude, see `BusPowerRows`), so `abstol` bounds the current
mismatch of those rows.

# Arguments
- `u`: Combined state vector [x; y] - modified in place
- `adapter`: TDSAdapter with residual evaluation capability
- `t`: The time of the state `u` (no default: a model's residuals can depend on it
  directly, so a solve at another time solves other equations)
- `maxiters`: Maximum Newton iterations (default: $(DEFAULT_REINIT_MAXITERS))
- `abstol`: Convergence tolerance on ||r|| (default: $(DEFAULT_REINIT_TOL))
- `jac_update_freq`: Jacobian update frequency (default: $(DEFAULT_JAC_UPDATE_FREQ))
  - `1` = standard Newton (Jacobian every iteration)
  - `>1` = Modified Newton (reuse Jacobian for multiple iterations)
- `verbose`: Print per-iteration progress (default: false)

# Returns
`true` if converged, `false` otherwise

# Performance
Uses Modified Newton by default: the Jacobian is only recomputed every `jac_update_freq`
iterations, significantly reducing computation time for large systems. The KLU
factorization is also reused when the Jacobian hasn't changed.

For a 70k-bus system with ~117ms per Jacobian computation, this can reduce
reinitialization time from ~2.5s (14 iterations × 160ms) to ~1s (3 Jacobian
updates + 11 cheap iterations).

**Cached infrastructure**: The closure, Jacobian prep, matrix, and backend are
cached in TDSAdapter after the first call, saving ~50-200ms per subsequent event.
"""
function reinitialize_algebraic!(u::AbstractVector, adapter::TDSAdapter, t::Real; kwargs...)
    failure = _reinitialize_algebraic(u, adapter, t; kwargs...)
    failure === nothing && return true
    @warn "Algebraic reinit: " * _reinit_failure_message(failure) * _unreduced_constraint_note(adapter)
    return false
end

# The rate of the start for a value on a surface (`_start_algebraic_rate`): the variables without mass `z` (the
# algebraic variables and the states whose mass is zero) at the rate the linear system of the reinitialization gives,
# `F_z ż = -(F_x ẋ + F_t)` on the rows without mass, the other states at `rates`. `F_x ẋ + F_t` is one central
# difference of the reinitialization's residual over `±h` along `(ẋ, 1)` at the start's `z` (the residual it leaves
# there cancels); `F_z` is its Jacobian at the start, factorized afresh (with pivoting) by its solver
function _start_algebraic_rate(adapter::TDSAdapter, x::AbstractVector, y::AbstractVector, t::Real,
                               rates::AbstractVector, h::Real)
    n_diff = adapter.n_diff
    u = vcat(x, y)
    adapter.reinit_cache === nothing && (adapter.reinit_cache = _reinit_cache(adapter, u, t))
    cache = adapter.reinit_cache
    udot = vcat(rates, zeros(length(y)))
    cache.n_unknowns == 0 && return udot
    z, _ = _reinit_unknowns(u, adapter, cache, get_backend(adapter.adcache))
    moving = trues(n_diff)
    moving[cache.massless] .= false
    along = function (out, s)
        @inbounds for k in 1:n_diff
            adapter.x_buf[k] = moving[k] ? x[k] + s * rates[k] : x[k]
        end
        cache.t[] = t + s
        cache.closure(out, z)
        return out
    end
    driven = (along(zeros(cache.n_unknowns), h) .- along(zeros(cache.n_unknowns), -h)) ./ (2h)
    along(cache.newton_bufs.g, 0.0)  # (the closure at the start, for its Jacobian)
    J = cache.newton_bufs.J
    DifferentiationInterface.jacobian!(cache.closure, cache.newton_bufs.g_trial, J, cache.jac_prep, cache.jac_backend, z)
    solver = cache.linear_solver
    reset_factorization!(solver)
    try
        factorize!(solver, J)
    catch e
        e isa SingularException || rethrow()
        return "the Jacobian of the rows without mass with respect to their variables (the reinitialization's) " *
               "is singular at the start: KLU found it so" * _start_singular_rows(J, adapter, cache)
    end
    κ = KLU.condest(solver.fact)
    κ * eps() >= 1 && return "the Jacobian of the rows without mass with respect to their variables (the " *
        "reinitialization's) is singular at the start: KLU estimates its condition number at $κ" *
        _start_singular_rows(J, adapter, cache)
    zdot = similar(driven)
    linear_solve!(solver, zdot, driven)
    n_z = length(cache.massless)
    for (k, a) in enumerate(cache.massless)
        udot[a] = -zdot[k]
    end
    for k in 1:adapter.n_alg
        udot[n_diff + k] = -zdot[n_z + k]
    end
    return udot
end

# The variables of the columns of `J` (the reinitialization's Jacobian) without a nonzero entry, named, for the
# message of a singular one
function _start_singular_rows(J, adapter::TDSAdapter, cache::ReinitCache)
    am = get_address_manager(adapter.adcache)
    isnothing(am.reverse_map) && build_reverse_map!(am)
    n_z = length(cache.massless)
    empty = [c for c in 1:size(J, 2) if all(iszero, view(nonzeros(J), nzrange(J, c)))]
    isempty(empty) && return ""
    name(c) = c <= n_z ? lookup_address(am, StateVar(), UInt(cache.massless[c])) :
                         lookup_address(am, AlgebVar(), UInt(c - n_z))
    return "; no row reads " * join((string(name(c)) for c in empty[1:min(5, end)]), ", ") *
           " (model, variable, instance) there"
end

# `reinitialize_algebraic!` reporting what failed: `nothing` where the solve converged, else
# the failure (`_reinit_failure`) for the caller to report and stop on
function _reinitialize_algebraic(
    u::AbstractVector{T},
    adapter::TDSAdapter,
    t::Real;
    maxiters::Int = DEFAULT_REINIT_MAXITERS,
    abstol::Real = DEFAULT_REINIT_TOL,
    jac_update_freq::Int = DEFAULT_JAC_UPDATE_FREQ,
    verbose::Bool = false
) where T
    n_diff = adapter.n_diff

    # Copy x from u into adapter.x_buf (stable reference for cached closure)
    @inbounds for i in 1:n_diff
        adapter.x_buf[i] = u[i]
    end

    # Lazy init ReinitCache (created on first call, reused for all subsequent events)
    adapter.reinit_cache === nothing && (adapter.reinit_cache = _reinit_cache(adapter, u, t))
    cache = adapter.reinit_cache
    cache.n_unknowns == 0 && return nothing  # Nothing without mass
    cache.t[] = t   # the closure evaluates the residuals at the time of the state

    # The unknowns: a view of y in u (no states without mass, ForwardDiff), else a copy;
    # the iterations behind a function barrier, specialized on which
    z, needs_copyback = _reinit_unknowns(u, adapter, cache, get_backend(adapter.adcache))
    return _reinit_newton!(u, z, needs_copyback, adapter, cache; maxiters, abstol, jac_update_freq, verbose)
end

# The Newton iterations of `reinitialize_algebraic!` on the unknowns `z`, and the copy back
# of the result into `u`
function _reinit_newton!(u::AbstractVector, z::AbstractVector, needs_copyback::Bool, adapter::TDSAdapter,
                         cache::ReinitCache; maxiters::Int, abstol::Real, jac_update_freq::Int, verbose::Bool)
    # Modified Newton (shared engine): refresh the Jacobian on the first iteration
    # (topology may have changed since the last event), every jac_update_freq
    # iterations, or when |r| drops by less than 10%. Steps are limited and keep every bus
    # voltage magnitude positive by the rule of the Y-bus algebraic solve
    # (`BUS_VOLTAGE_STEP_LIMIT`, `positive`) and are taken without backtracking. The KLU
    # factorization is reused across events.
    residual! = (r, z) -> cache.closure(r, z)
    jacobian! = (Jm, z) -> DifferentiationInterface.jacobian!(
        cache.closure, cache.newton_bufs.g_trial, Jm, cache.jac_prep, cache.jac_backend, z)
    policy = NewtonPolicy(; jacobian_every=jac_update_freq, stall_ratio=0.9, max_step=BUS_VOLTAGE_STEP_LIMIT, positive=cache.voltages, max_tries=1)

    result = newton_solve!(residual!, jacobian!, z, cache.newton_bufs, cache.linear_solver, policy;
                           abstol, maxiters, catch_linear_failure=true, verbose)

    if result.converged
        @debug "Algebraic reinit converged" result.iterations result.residual_norm abstol
        _copyback_unknowns!(u, z, adapter, cache, needs_copyback)
        return nothing
    end

    # The residual at the point the solve stopped at, for the worst rows
    r = cache.newton_bufs.g
    cache.closure(r, z)
    result.status === :singular || _copyback_unknowns!(u, z, adapter, cache, needs_copyback)
    return _reinit_failure(result, r, abstol, adapter, cache)
end

# What failed in a reinitialization: the Newton `status` (`:maxiters`, `:diverged`, or
# `:singular` for a Jacobian that could not be factorized), the iterations, the largest residual
# (a current mismatch at the bus rows) against `abstol`, and the worst rows named by
# (model, variable, instance)
function _reinit_failure(result, r::AbstractVector, abstol, adapter::TDSAdapter, cache::ReinitCache)
    return (; status = result.status, iterations = result.iterations,
            residual = norm(r, Inf), abstol, rows = _worst_reinit_rows(r, adapter, cache))
end

# The five largest of the rows without mass `r` (the states without mass, then g), named by
# (model, variable, instance) through the reverse maps
function _worst_reinit_rows(r::AbstractVector, adapter::TDSAdapter, cache::ReinitCache)
    am = get_address_manager(adapter.adcache)
    isnothing(am.reverse_map) && build_reverse_map!(am)
    n_massless = length(cache.massless)
    row_name(k) = k <= n_massless ? lookup_address(am, StateRes(), UInt(cache.massless[k])) :
                                    lookup_address(am, AlgebRes(), UInt(k - n_massless))
    worst = sortperm(abs.(r); rev = true)[1:min(5, length(r))]
    return join(("$(row_name(k)): $(abs(r[k]))" for k in worst), "; ")
end

# The constraints on the states alone that the reduction found in the modes the automata are in and does not make
# index 1 (a singular reinitialization is likely theirs), as a sentence for the failure's message
function _unreduced_constraint_note(adapter::TDSAdapter)
    reasons = _unreduced_constraints(adapter)
    isempty(reasons) && return ""
    return ". The reduction found a constraint on the states alone in this discrete state that one differentiation " *
           "does not make index 1: " * join(reasons, "; ")
end

function _reinit_failure_message(f)
    f.status === :mismatch && return "the reinitialized state does not satisfy the DAE as the integration " *
        "evaluates it at this time (the reinitialization solved other equations): max residual " *
        "$(f.residual) against $(f.abstol), worst rows (model, variable, instance) $(f.rows)"
    why = f.status === :singular ? "the Jacobian is singular" :
          f.status === :diverged ? "diverged" : "did not converge"
    return "$why after $(f.iterations) iterations, max residual $(f.residual) against $(f.abstol), " *
           "worst rows (model, variable, instance) $(f.rows)"
end

"""
    _reinit_mismatch(u, adapter, t, abstol) -> Union{Nothing, NamedTuple}

The check of a reinitialization's result: the residual of the DAE as the integration
evaluates it (`_evaluate_packed!`, with the reduction of state-only constraints) at the
reinitialized state `u` and the integrator's time `t`, at the rows without mass (the
states without mass and `g`, the bus power balance rows as current mismatch, as the solve
measures them). `nothing` where the largest is below the solve's tolerance `abstol` (the
solve's own criterion), else the failure (status `:mismatch`), its worst rows named. It
catches a solve whose equations are not the integration's: before `t` was passed to the
reinitialization, it solved a time-dependent model's residuals at `t = 0` instead of the
integration's actual time, and reported success.
"""
function _reinit_mismatch(u::AbstractVector, adapter::TDSAdapter, t::Real, abstol::Real)
    cache = adapter.reinit_cache
    cache.n_unknowns == 0 && return nothing
    du, r = cache.check_buf, cache.check_rows
    _evaluate_packed!(du, adapter, u, t, nothing)
    n_diff, n_z = adapter.n_diff, length(cache.massless)
    @inbounds begin
        for (k, a) in enumerate(cache.massless)
            r[k] = du[a]
        end
        for k in 1:adapter.n_alg
            r[n_z + k] = du[n_diff + k]
        end
    end
    cache.rows === nothing || to_current_mismatch!(view(r, (n_z + 1):cache.n_unknowns), view(u, (n_diff + 1):length(u)),
                                                   cache.rows, adapter.status[:Bus])
    residual = norm(r, Inf)
    residual < abstol && return nothing
    return (; status = :mismatch, iterations = 0, residual, abstol, rows = _worst_reinit_rows(r, adapter, cache))
end

# The states of `x` without mass (indices into `x`), from the diagonal of the mass matrix
_states_without_mass(adapter::TDSAdapter) =
    (mass = diag(adapter.mass_matrix); [k for k in 1:adapter.n_diff if iszero(mass[k])])

# The reinitialization's closure, Jacobian and Newton buffers for the adapter's unknowns:
# y alone where no state lacks mass (the residual `g`, evaluated into the closure's output
# directly), else [x without mass; y] (the packed DAE residual at their rows). The Jacobian
# is prepared on the unknowns as the iterations hold them (a view of y in `u`, or a copy),
# behind a function barrier on the closure of the layout
function _reinit_cache(adapter::TDSAdapter, u::AbstractVector, t::Real)
    massless = _states_without_mass(adapter)
    n_w = length(massless) + adapter.n_alg
    # The time the closure evaluates at: the time of each reinitialization's state, set by it
    # (here that of the first, at which the preparation of the Jacobian evaluates the closure)
    tref = Ref(Float64(t))
    n_w == 0 && return ReinitCache(nothing, nothing, nothing, nothing, nothing, nothing, massless, 0, Float64[], Int[],
                                   tref, nothing, Float64[], Float64[])
    rows = bus_power_rows(get_address_manager(adapter.adcache))
    if isempty(massless)
        closure, pattern = _reinit_closure_y(adapter, rows, tref)
        z = first(_get_y_working_buffer(u, adapter, get_backend(adapter.adcache)))
        return _prepared_reinit_cache(closure, pattern, z, rows, adapter, massless, n_w, Float64[], tref)
    end
    closure, pattern = _reinit_closure_z(adapter, rows, massless, tref)
    unknowns = zeros(n_w)
    return _prepared_reinit_cache(closure, pattern, unknowns, rows, adapter, massless, n_w, unknowns, tref)
end

# The residual g on y at adapter.x_buf (stable across calls; x is copied into x_buf at the
# start of each reinitialize_algebraic! call) and at the time `tref[]` (set by each call to
# the time of its state), and its sparsity
function _reinit_closure_y(adapter::TDSAdapter, rows, tref::Base.RefValue{Float64})
    adcache, status, n_diff = adapter.adcache, adapter.status, adapter.n_diff
    closure = let f_buf=adapter.f_buf, x_buf=adapter.x_buf, adcache=adcache, n_diff=n_diff, status=status, rows=rows,
                  tref=tref
        (g_out, y_in) -> begin
            Tg = eltype(g_out)
            if Tg == Float64
                adcache(f_buf, g_out, x_buf, y_in, status, tref[])
            else
                # Cached Dual-typed f buffer for the ForwardDiff path
                f_dual = get_tds_dual_buffers(adcache, Tg, n_diff, length(g_out)).f_buf
                adcache(f_dual, g_out, x_buf, y_in, status, tref[])
            end
            rows === nothing || to_current_mismatch!(g_out, y_in, rows, status[:Bus])
            nothing
        end
    end
    return closure, jacobian_sparsity_pattern(get_address_manager(adcache), get_system(adcache), get_formulation(adcache))
end

# The packed DAE residual (`_evaluate_packed!`'s, the reduced rows in their reduced form) on
# z = [x without mass; y] at adapter.x_buf with its states without mass from z, at the time
# `tref[]`, read at those states' rows and at g, and its sparsity
function _reinit_closure_z(adapter::TDSAdapter, rows, massless::Vector{Int}, tref::Base.RefValue{Float64})
    adcache, status = adapter.adcache, adapter.status
    n_diff, n_alg, n_z = adapter.n_diff, adapter.n_alg, length(massless)
    n_w = n_z + n_alg
    red = adapter.reduction
    work = Dict{DataType, Tuple{Vector, Vector}}()
    closure = let x_buf=adapter.x_buf, adcache=adcache, n_diff=n_diff, n_alg=n_alg, n_z=n_z, n_w=n_w,
                  massless=massless, status=status, rows=rows, red=red, work=work, tref=tref
        (r_out, z_in) -> begin
            Tz = eltype(r_out)
            uw, duw = get!(() -> (Vector{Tz}(undef, n_diff + n_alg), Vector{Tz}(undef, n_diff + n_alg)),
                           work, Tz)::Tuple{Vector{Tz}, Vector{Tz}}
            @inbounds begin
                for k in 1:n_diff
                    uw[k] = x_buf[k]
                end
                for (k, a) in enumerate(massless)
                    uw[a] = z_in[k]
                end
                for k in 1:n_alg
                    uw[n_diff + k] = z_in[n_z + k]
                end
            end
            adcache(view(duw, 1:n_diff), view(duw, (n_diff + 1):(n_diff + n_alg)), view(uw, 1:n_diff),
                    view(uw, (n_diff + 1):(n_diff + n_alg)), status, tref[])
            red === nothing || _apply_reduction!(duw, red, uw, tref[])
            @inbounds begin
                for (k, a) in enumerate(massless)
                    r_out[k] = duw[a]
                end
                for k in 1:n_alg
                    r_out[n_z + k] = duw[n_diff + k]
                end
            end
            rows === nothing || to_current_mismatch!(view(r_out, (n_z + 1):n_w), view(z_in, (n_z + 1):n_w),
                                                    rows, status[:Bus])
            nothing
        end
    end
    full = red === nothing ?
        jacobian_sparsity_pattern(get_address_manager(adcache), get_system(adcache), get_formulation(adcache);
                                  include_differential=true) :
        red.pattern
    rows_without_mass = [massless; n_diff .+ (1:n_alg)]
    return closure, full[rows_without_mass, rows_without_mass]
end

# The cache of a closure on unknowns held as `z`: its Jacobian (a current-mismatch row also
# depends on its bus voltage magnitude) prepared on them, and the Newton buffers
function _prepared_reinit_cache(closure, pattern, z, rows, adapter::TDSAdapter, massless::Vector{Int}, n_w::Int,
                                unknowns::Vector{Float64}, tref::Base.RefValue{Float64})
    n_z = length(massless)
    if rows !== nothing
        pattern = pattern .| sparse(n_z .+ vcat(rows.p, rows.q), n_z .+ vcat(rows.v, rows.v), true, n_w, n_w)
    end
    r_buf = zeros(n_w)
    (; J, backend, prep) = prepare_sparse_jacobian(closure, r_buf, z, adapter.adcache; sparsity=pattern)

    # Newton buffers: z/r are the iterate; trial arrays and step are owned here
    bufs = NewtonBuffers(r_buf, similar(r_buf), Vector{Float64}(undef, n_w), similar(r_buf), J)
    voltages = rows === nothing ? Int[] : n_z .+ rows.v
    return ReinitCache(closure, prep, J, backend, KLULinearSolver(), bufs, massless, n_w, unknowns, voltages,
                       tref, rows, zeros(adapter.n_diff + adapter.n_alg), zeros(n_w))
end

# The unknowns of the Newton iterations at u, and whether they are copied back to u after
# them: y as a view of u without states without mass (ForwardDiff), else a copy
function _reinit_unknowns(u::AbstractVector, adapter::TDSAdapter, cache::ReinitCache, backend)
    isempty(cache.massless) && return _get_y_working_buffer(u, adapter, backend)
    z = cache.unknowns
    n_diff, n_z = adapter.n_diff, length(cache.massless)
    @inbounds begin
        for (k, a) in enumerate(cache.massless)
            z[k] = u[a]
        end
        for k in 1:adapter.n_alg
            z[n_z + k] = u[n_diff + k]
        end
    end
    return (z, true)
end

function _copyback_unknowns!(u::AbstractVector, z, adapter::TDSAdapter, cache::ReinitCache, needs_copyback::Bool)
    isempty(cache.massless) && return _copyback_y!(u, z, adapter, needs_copyback)
    n_diff, n_z = adapter.n_diff, length(cache.massless)
    @inbounds begin
        for (k, a) in enumerate(cache.massless)
            u[a] = z[k]
        end
        for k in 1:adapter.n_alg
            u[n_diff + k] = z[n_z + k]
        end
    end
    nothing
end

export reinitialize_algebraic!

# ============================================================================
# Initial Guess Strategies (Internal)
# ============================================================================

"""
    bus_voltage_addresses(sys::SystemModel, am::AddressManager) -> Vector{UInt}

Algebraic addresses of the voltage magnitude `v` of every bus model of `sys` (the models
implementing the `:bus_group` group), in model order; empty when the formulation has no
`v` (CurrentInjection).
"""
function bus_voltage_addresses(sys::SystemModel, am::AddressManager)
    addrs = UInt[]
    for (model_name, models) in pairs(sys.models)
        isempty(models) && continue
        model_type = get_base_model_type(models)
        implements_group(model_type) === :bus_group || continue
        key = (AlgebVar(), model_name, resolve_group_var_name(:bus_group, :v, model_type))
        haskey(am.addresses, key) && append!(addrs, am.addresses[key])
    end
    return addrs
end

"""
    _apply_flat_start!(u::AbstractVector, adapter::TDSAdapter)

Apply flat start initial conditions to voltage magnitudes (V=1.0).
Theta is NOT reset to preserve the delta-theta coupling in generator models.
"""
function _apply_flat_start!(u::AbstractVector, adapter::TDSAdapter)
    am = get_address_manager(adapter.adcache)
    n_diff = adapter.n_diff

    # PowerInjection formulation: Reset v of every bus model only (theta preserved)
    for addr in adapter.bus_v_addrs
        # Addresses are in algebraic space, offset by n_diff in u vector
        u[n_diff + addr] = 1.0
    end

    # WARNING: Do NOT reset theta here. Generator models use (delta - theta),
    # so resetting theta=0 while delta is fixed breaks the coupling.

    # CurrentInjection formulation: Reset v_r=1.0 and v_i=0.0
    if haskey(am.addresses, (AlgebVar(), :Bus, :v_r))
        v_r_addrs = get_addresses(am, :Bus, :v_r)
        for addr in v_r_addrs
            u[n_diff + addr] = 1.0
        end
    end

    if haskey(am.addresses, (AlgebVar(), :Bus, :v_i))
        v_i_addrs = get_addresses(am, :Bus, :v_i)
        for addr in v_i_addrs
            u[n_diff + addr] = 0.0
        end
    end

    @debug "Applied flat start: V=1.0 for all buses"
    nothing
end

"""
    _apply_pre_fault_v!(u::AbstractVector, adapter::TDSAdapter, pre_fault_v::Vector{Float64})

Apply pre-fault voltage magnitudes as initial guess.
Theta is NOT reset to preserve the delta-theta coupling in generator models.
"""
function _apply_pre_fault_v!(u::AbstractVector, adapter::TDSAdapter, pre_fault_v::Vector{Float64})
    am = get_address_manager(adapter.adcache)
    n_diff = adapter.n_diff

    # PowerInjection formulation: Set v of every bus model to pre-fault values (theta
    # preserved); pre_fault_v is in the order of adapter.bus_v_addrs
    for (i, addr) in enumerate(adapter.bus_v_addrs)
        u[n_diff + addr] = pre_fault_v[i]
    end

    # WARNING: Do NOT reset theta here. Generator models use (delta - theta),
    # so resetting theta=0 while delta is fixed breaks the coupling.

    # CurrentInjection formulation: Set v_r to pre-fault V, v_i to 0
    if haskey(am.addresses, (AlgebVar(), :Bus, :v_r))
        v_r_addrs = get_addresses(am, :Bus, :v_r)
        for (i, addr) in enumerate(v_r_addrs)
            u[n_diff + addr] = pre_fault_v[i]  # V * cos(0) = V
        end
    end

    if haskey(am.addresses, (AlgebVar(), :Bus, :v_i))
        v_i_addrs = get_addresses(am, :Bus, :v_i)
        for addr in v_i_addrs
            u[n_diff + addr] = 0.0  # V * sin(0) = 0
        end
    end

    @debug "Applied pre-fault voltage magnitudes"
    nothing
end

"""
    _apply_initial_guess!(u, adapter, alg)

Apply initial guess based on `alg.initial_guess` setting:
- `:flat_start` - V=1.0, theta kept
- `:pre_fault` - pre-fault V magnitudes, theta kept
- `:current` - no change (use current voltages)
"""
function _apply_initial_guess!(u::AbstractVector, adapter::TDSAdapter, alg)
    if alg.initial_guess == :flat_start
        _apply_flat_start!(u, adapter)
    elseif alg.initial_guess == :pre_fault
        if isnothing(alg.pre_fault_v)
            error("PowerFlowInit: pre_fault_v not captured. This should not happen.")
        end
        _apply_pre_fault_v!(u, adapter, alg.pre_fault_v)
    elseif alg.initial_guess == :current
        @debug "Using current voltages as initial guess"
        # Do nothing - use current values
    end
    nothing
end

"""
    _capture_voltage_magnitudes(adapter::TDSAdapter) -> Vector{Float64}

Capture current voltage magnitudes from the adapter's algebraic buffer.
Works for PowerInjection (`v` of every bus model, see `bus_voltage_addresses`) and
CurrentInjection (v_r, v_i of `Bus`) formulations.
"""
function _capture_voltage_magnitudes(adapter::TDSAdapter)
    am = get_address_manager(adapter.adcache)
    y_buf = adapter.y_buf

    # PowerInjection formulation: directly use v
    isempty(adapter.bus_v_addrs) || return [y_buf[addr] for addr in adapter.bus_v_addrs]

    # CurrentInjection formulation: compute V = sqrt(v_r^2 + v_i^2)
    if haskey(am.addresses, (AlgebVar(), :Bus, :v_r))
        v_r_addrs = get_addresses(am, :Bus, :v_r)
        v_i_addrs = get_addresses(am, :Bus, :v_i)
        return [sqrt(y_buf[v_r_addrs[i]]^2 + y_buf[v_i_addrs[i]]^2) for i in eachindex(v_r_addrs)]
    end

    error("_capture_voltage_magnitudes: unknown formulation (no bus model with :v, and no :v_r in Bus)")
end

# ============================================================================
# OrdinaryDiffEq Integration (Internal)
# ============================================================================

"""
    _bind_adapter(alg::PowerFlowInit, adapter::TDSAdapter, jump; abstol, verbose=false) -> PowerFlowInit

Bind a TDSAdapter to the PowerFlowInit configuration, for an integration with the
absolute tolerance `abstol` (a number, or one per state: the smallest counts). Called
by `solve()` when processing events. Without an `abstol` of its own, the
reinitialization stops at `REINIT_TOL_FRACTION` of the integration's.

If `initial_guess == :pre_fault`, captures the current V magnitudes from the adapter
(which represent the pre-fault steady-state solution).

`jump` (an augmentation carried through the integration beside the DAE state, `nothing`
for the plain simulation; see `_initialize_dae!`): unused here, accepted so every
caller passes it uniformly; binding `alg` to `adapter` does not depend on it.
"""
function _bind_adapter(alg::PowerFlowInit, adapter::TDSAdapter, jump; abstol, verbose::Bool = false)
    pre_fault_v = nothing

    if alg.initial_guess == :pre_fault
        # Capture pre-fault voltage magnitudes from adapter's initial state
        pre_fault_v = _capture_voltage_magnitudes(adapter)
    end

    tol = something(alg.abstol, REINIT_TOL_FRACTION * minimum(abstol))
    PowerFlowInit(tol, alg.maxiters, alg.initial_guess, adapter, pre_fault_v, verbose)
end

# Passthrough for other init algorithms (BrownFullBasicInit, NoInit, etc.)
_bind_adapter(alg, ::TDSAdapter, jump; abstol, verbose::Bool = false) = alg

"""
    reinitialize_network!(u, t, alg::PowerFlowInit) -> Bool

Re-solve the variables without mass of `u = [x; y]` at time `t` (the algebraic
variables and the states whose mass is zero, `reinitialize_algebraic!`) with the states
with mass held, as `alg` configures (initial guess, tolerance, iterations), after restoring
the state-only constraints an event may have broken. Returns `true` when it converged. Where
the Newton iterations do not (a numerical failure: the state after a long fault may have no
solution) it warns, with the time, the iterations, the residual and the worst rows, and
returns `false`, the callers terminating the integration (`_terminate_unless!`). A converged
result is checked against the DAE as the integration evaluates it at `t` (`_reinit_mismatch`),
and fails the same way where it does not satisfy it. Errors in the models or the code throw.
"""
function reinitialize_network!(u::AbstractVector, t, alg::PowerFlowInit)
    adapter = alg.adapter
    if isnothing(adapter) || isnothing(alg.abstol)
        error("PowerFlowInit: not bound to a simulation. This should not happen - please report a bug.")
    end

    # Extract verbose for @vprint macro
    verbose = alg.verbose

    # Restore state-only constraints an event may have broken (e.g. a status switch that
    # changes which inductor currents meet at a bus)
    enforce_state_constraints!(u, adapter, t)

    # Apply initial guess based on setting
    _apply_initial_guess!(u, adapter, alg)

    # Call Newton solver on algebraic variables with timing, at the time of the state
    t_reinit = @elapsed failure = _reinitialize_algebraic(
        u, adapter, t;
        maxiters=alg.maxiters, abstol=alg.abstol, verbose=verbose
    )
    # The result checked against the DAE the integration evaluates, at the same time
    failure === nothing && (failure = _reinit_mismatch(u, adapter, t, alg.abstol))

    t_ms = round(t_reinit * 1000, digits=1)
    if failure !== nothing
        @warn "PowerFlowInit: t=$(t)s " * _reinit_failure_message(failure) * _unreduced_constraint_note(adapter) *
              "; stopping the simulation"
        return false
    end
    @vprint "PowerFlowInit: t=$(t)s converged ($(t_ms) ms)"
    return true
end

"""
    EventReinit(network::PowerFlowInit, changed::Base.RefValue{Bool}, jump) <: DAEInitializationAlgorithm

The reinitialization of the network after an event of `events.jl`: `network` where the
event changed the equations without mass, nothing otherwise. The event sets `changed` (the
algebraic residuals, or the residuals of the states without mass, at the state differ
before and after it applies, see `_network_change!`), as the automata do for their transitions (`settle_instant!`): an
event that leaves the algebraic equations as they were (a setpoint of a controller)
leaves the state as consistent as the integration keeps it, and a Newton solve there
would move it by less than the integration's tolerance, differently for nearby
parameters. The flag is read once and returns to `true`, so a later reinitialization at
the same instant (a transition of an automaton that changes the algebraic equations)
re-solves the network.

`jump` (`nothing` for the plain simulation) is carried alongside `network` so that
`_initialize_dae!` has it where DiffEq calls in through the callback's `initializealg`,
with no parameter of its own to pass one through (see `_initialize_dae!`).
"""
struct EventReinit{J} <: DAEInitializationAlgorithm
    network::PowerFlowInit
    changed::Base.RefValue{Bool}
    jump::J
end

function reinitialize_network!(u::AbstractVector, t, alg::EventReinit)
    changed, alg.changed[] = alg.changed[], true
    changed || return true
    algebraic_equations_changed!(alg.network.adapter)
    return reinitialize_network!(u, t, alg.network)
end

"""
    at_current_state(alg) -> DAEInitializationAlgorithm
    at_current_state(alg, jump) -> DAEInitializationAlgorithm

The reinitialization of an instant of transitions of the automata (`_reinitialize!`, from
`settle_instant!`), from the simulation's `alg`: its Newton iterations start from the
integrator's state (`initial_guess = :current`). Before the transitions the state is
consistent, and a transition changes the equations without mass at a crossing of its
surface, so the state is the guess nearest the solution; the configured guess (`:pre_fault`
or `:flat_start`, which reset the bus voltage magnitudes) belongs to network events, which
change the network itself (faults, trips, `FieldModification`), and they keep it
(`EventReinit`, `_initialize_dae!`). A transition inside a fault reinitialized from the
default flat start (every bus voltage magnitude reset to 1 while the fault holds them near
0.3) did not converge (the converter at bus 2 entering its current limit at t ≈ 1.5 in the
calibration replay e1 (3, 3), round 9e; `:pre_fault` failed the same way), where the current
state converges at once.

- `PowerFlowInit`: the same tolerance, iterations and adapter, the guess `:current`.
- `EventReinit`: the transition variant of its network's (an instant of transitions always
  reinitializes: its algebraic check found the equations changed).
- The others start from the integrator's state already (the Y-bus DAE's `HostReinit`
  solves from `u`; DiffEq's algorithms start from the integrator's state; `NoInit` solves
  nothing) and are returned as they are.

The 2-argument form is the extension point: `jump` (an augmentation carried through the
integration beside the DAE state, `nothing` for the plain simulation) never changes which
guess `alg` starts from, so the default delegates to the 1-argument form for every `alg`;
an augmentation adds a method only where it needs a different `alg` transformation.
"""
at_current_state(alg::PowerFlowInit) =
    PowerFlowInit(alg.abstol, alg.maxiters, :current, alg.adapter, alg.pre_fault_v, alg.verbose)
at_current_state(alg::EventReinit) = at_current_state(alg.network)
at_current_state(alg::DAEInitializationAlgorithm) = alg
at_current_state(alg, jump) = at_current_state(alg)

"""
    algebraic_equations_changed!(host)

The equations without mass of the DAE of `host` changed at the instant being settled: an
event changed them (`EventReinit`, which reinitializes only then: the event's network
check, `_network_change!`, found the residuals without mass at the same state differ
before and after it, or the event does not check), or a transition of the automata did
(their algebraic check, `settle_instant!`). The numerical structure of the
reinitialization's Jacobian may have changed with them (a mode can move the nonzero of a
row off its own variable: REGCA1's high-voltage clamp holds `Vt = Volim` in the equation
of `dq`, which has no `dq` in it), so the pivots of an earlier factorization are not those
of this one. Recorded on the reinitialization's linear solver (`reset_factorization!` of
the `ReinitCache`'s `KLULinearSolver`), and consumed by its next `factorize!`, a full
factorization with pivoting rather than a numeric refactorization with the old pivots;
a matrix that is singular throws there all the same. Before the first reinitialization
there is nothing to invalidate: its first factorization is a full one.
"""
function algebraic_equations_changed!(adapter::TDSAdapter)
    cache = adapter.reinit_cache
    (cache === nothing || cache.linear_solver === nothing) && return nothing
    reset_factorization!(cache.linear_solver)
    return nothing
end

"""
    _initialize_dae!(integrator, prob, alg, jump, ::Val{iip})

The extension point of an integration's reinitialization: every caller we control passes
`jump` (an augmentation carried through the integration beside the DAE state, `nothing`
for the plain simulation) explicitly, but DiffEq calls a callback's `initializealg` with
no parameter of its own to carry one (`alg::EventReinit`, `alg::AutomatonSettle`: their
own `jump` field is where that boundary reads it, then forwards to this 5-argument form).
The default for `jump === nothing` is the integration's 4-argument reinitialization,
unchanged whatever `alg` is; an augmentation adds a method for the `alg` types it needs to
intercept (e.g. `Union{PowerFlowInit, EventReinit}`, reinitializing the DAE and then the
extra state it adds).
"""
_initialize_dae!(integrator, prob, alg, ::Nothing, ::Val{iip}) where {iip} =
    _initialize_dae!(integrator, prob, alg, Val(iip))

# DiffEq calls a callback's `initializealg` with no parameter for `jump`: the `EventReinit`
# it is holding carries its own (`EventReinit`), read here and forwarded to the extension
# point above
_initialize_dae!(integrator, prob, alg::EventReinit, ::Val{iip}) where {iip} =
    _initialize_dae!(integrator, prob, alg, alg.jump, Val(iip))

function _initialize_dae!(integrator, prob, alg::EventReinit, ::Nothing, ::Val{iip}) where {iip}
    _terminate_unless!(integrator, reinitialize_network!(integrator.u, integrator.t, alg)) &&
        _restart_history!(integrator)
    return nothing
end

# The reinitialization of the events of a simulation reinitializing with `alg` (see
# `EventReinit`, which carries `jump` for `_initialize_dae!`); another algorithm (e.g.
# `BrownFullBasicInit`) as it is. `jump` never changes which algorithm this builds, only
# whether it is there for `EventReinit` to carry.
_event_reinit(alg::PowerFlowInit, jump, changed::Base.RefValue{Bool}) = EventReinit(alg, changed, jump)
_event_reinit(alg, jump, changed::Base.RefValue{Bool}) = alg

"""
    _restart_history!(integrator)

Copy the reinitialized state into the previous states of a multistep method.
"""
function _restart_history!(integrator)
    if hasproperty(integrator, :uprev)
        copyto!(integrator.uprev, integrator.u)
    end
    if hasproperty(integrator, :uprev2)
        copyto!(integrator.uprev2, integrator.u)
    end
    return nothing
end

# Implement the _initialize_dae! interface for OrdinaryDiffEq
# This gets called by the callback system after an affect! modifies the state
function _initialize_dae!(
    integrator,
    prob,  # prob: unused but required by interface
    alg::PowerFlowInit,
    ::Val{iip}  # isinplace: handle both true/false
) where {iip}
    _terminate_unless!(integrator, reinitialize_network!(integrator.u, integrator.t, alg)) &&
        _restart_history!(integrator)
    return nothing
end
