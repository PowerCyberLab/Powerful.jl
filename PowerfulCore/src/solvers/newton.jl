"""
Newton-Raphson engine shared by power flow, contingency analysis and TDS
algebraic reinitialization.

Callers provide the residual and Jacobian as closures, preallocated buffers, a
linear solver (factorization cache) and a `NewtonPolicy`. Differences between
solvers (Jacobian reuse, trust region, line search) are policy settings, not
separate loops.
"""

using KLU: KLU, klu, klu!
using LinearAlgebra: ldiv!, norm

# ============================================================================
# Policy
# ============================================================================

"""
    NewtonPolicy(; jacobian_every=1, stall_ratio=Inf, stall_reused_only=false, reuse_jacobian=false,
                   max_step=Inf, positive=Int[], positive_fraction=0.001, damping=1.0, max_tries=8,
                   backtrack_factor=0.5, merit_norm=Inf, sufficient_decrease=0.0, strict_search=false)

Iteration policy for `newton_solve!`.

- `jacobian_every`: refresh the Jacobian every `k` iterations (1 = full Newton,
  >1 = "dishonest"/modified Newton reusing the factorization in between).
  The Jacobian is always refreshed on the first iteration, unless `reuse_jacobian`.
- `stall_ratio`: also refresh when `‖g‖∞ > stall_ratio · ‖g_prev‖∞` (an iteration that
  contracts too little); with `stall_reused_only`, only after an iteration that reused a Jacobian
  taken at an earlier iterate (a chord iteration), not after one with the Jacobian at its own iterate.
- `reuse_jacobian`: the first iteration of a call uses the factorization the solver holds
  from an earlier call (a chord Newton across calls) and refreshes only where it holds
  none. The caller discards it (`reset_factorization!`) where the equations changed.
- `max_step`: trust region; the Newton step is scaled so `‖Δy‖∞ ≤ max_step`.
- `positive`: indices of variables that stay positive (bus voltage magnitudes): a trial point
  is clipped so that each of them remains above `positive_fraction` of its value before the
  step (a variable that crosses zero reaches another solution of equations polar in it, with
  the opposite sign). Applies to variables that are positive at the start of the iteration.
- `damping`: constant factor applied to every step.
- `max_tries`: backtracking line search trials (step halved by `backtrack_factor`
  each time); a trial is accepted when `‖g‖` falls below `(1 - sufficient_decrease · α) ‖g_prev‖`
  (α the fraction of the step taken) in the norm `merit_norm` (`Inf`, or `2`, the merit
  function along which the Newton step is a descent direction). `max_tries=1` takes the
  (damped, limited) step as is.
- `strict_search`: where no trial is accepted, the iterate stays and the Jacobian is
  taken again at the next iteration (status `:stalled` if it was taken at this iterate);
  otherwise the last trial is accepted regardless.
"""
Base.@kwdef struct NewtonPolicy
    jacobian_every::Int = 1
    stall_ratio::Float64 = Inf
    stall_reused_only::Bool = false
    reuse_jacobian::Bool = false
    max_step::Float64 = Inf
    positive::Vector{Int} = Int[]
    positive_fraction::Float64 = 0.001
    damping::Float64 = 1.0
    max_tries::Int = 8
    backtrack_factor::Float64 = 0.5
    merit_norm::Float64 = Inf
    sufficient_decrease::Float64 = 0.0
    strict_search::Bool = false
end

"""
    BUS_VOLTAGE_STEP_LIMIT

The step limit (`NewtonPolicy` `max_step`, 0.5 pu or rad) of the Newton iterations on bus voltages
in polar form (the reinitialization of the phasor simulation and the algebraic solve of the Y-bus
integrations), which also keep the voltage magnitudes positive (`positive`, see
[`NewtonPolicy`](@ref)). At a bus whose voltage is near zero the column of the Jacobian for its angle
is proportional to the voltage, and an unlimited Newton step in the angle is whole turns long, to a
root that is the same modulo 2π; a voltage that crosses zero reaches the other solution of the polar
equations, `(-v, θ + π)`.
"""
const BUS_VOLTAGE_STEP_LIMIT = 0.5

# ============================================================================
# Linear solvers (factorization caches)
# ============================================================================

"""
    AbstractNewtonLinearSolver

Factorization cache for `newton_solve!`. Implementations provide
`factorize!(solver, J)` (full factorization on first use, numeric refactorization
afterwards), `linear_solve!(solver, dy, rhs)` and `reset_factorization!(solver)`.
"""
abstract type AbstractNewtonLinearSolver end

"""
    KLULinearSolver(ordering=:AMD)

Sparse KLU factorization cache with fill-reducing `ordering`: `:AMD` (default),
`:COLAMD`, or `:natural` (no permutation, for pre-ordered matrices).
"""
mutable struct KLULinearSolver <: AbstractNewtonLinearSolver
    ordering::Symbol
    fact::Union{Nothing, KLU.KLUFactorization{Float64, Int}}

    function KLULinearSolver(ordering::Symbol = :AMD)
        ordering in (:AMD, :COLAMD, :natural) ||
            error("KLU ordering must be :AMD, :COLAMD or :natural, got :$ordering")
        new(ordering, nothing)
    end
end

"""
    reset_factorization!(solver::AbstractNewtonLinearSolver) -> solver

Discard the cached factorization, so the next [`factorize!`](@ref) is a full one: a fresh
symbolic analysis and a numeric factorization with pivoting. Needed when the sparsity
pattern of the matrix changes, and where its numerical structure may have (its equations
changed, so that a nonzero of a row can move off the entry an earlier factorization
pivoted on): a refactorization reuses the pivots of the last full factorization, which
can then be zero or unstable, and a matrix that is singular throws there all the same.
"""
function reset_factorization! end

"""
    is_factorized(solver::AbstractNewtonLinearSolver) -> Bool

Whether `solver` holds a factorization that [`linear_solve!`](@ref) can use.
"""
function is_factorized end

"""
    factorize!(solver::AbstractNewtonLinearSolver, J) -> solver

Factorize the sparse matrix `J` into `solver`. The first call, and the first after
[`reset_factorization!`](@ref), performs the symbolic analysis and the numeric
factorization with pivoting; later calls reuse the symbolic analysis and the pivots and
only refactorize numerically, so `J` must keep the same sparsity pattern and numerical
structure (otherwise call [`reset_factorization!`](@ref) first). A singular `J` throws.
"""
function factorize! end

"""
    linear_solve!(solver::AbstractNewtonLinearSolver, dy, rhs) -> dy

Solve `J * dy = rhs` in place with the factorization from [`factorize!`](@ref).
"""
function linear_solve! end

reset_factorization!(s::KLULinearSolver) = (s.fact = nothing; s)
is_factorized(s::KLULinearSolver) = s.fact !== nothing

function factorize!(s::KLULinearSolver, J)
    if s.fact === nothing
        s.fact = s.ordering === :AMD ? klu(J) : _klu_with_ordering(J, s.ordering)
    else
        klu!(s.fact, J)  # Numeric refactorization, symbolic analysis and pivots reused
    end
    return s
end

# Ordering must be set before symbolic analysis, so build the factorization manually.
# KLU's Common.ordering codes: 0 = AMD, 1 = COLAMD, 2 = natural.
function _klu_with_ordering(J, ordering::Symbol)
    fact = KLU.KLUFactorization(size(J, 1), J.colptr .- 1, J.rowval .- 1, J.nzval)
    fact.common.ordering = ordering === :AMD ? 0 : ordering === :COLAMD ? 1 : 2
    KLU.klu_analyze!(fact)
    KLU.klu_factor!(fact)
    return fact
end

function linear_solve!(s::KLULinearSolver, dy, rhs)
    copyto!(dy, rhs)
    ldiv!(s.fact, dy)
    return dy
end

"""
    gpu_linear_solver() -> AbstractNewtonLinearSolver

GPU sparse linear solver for `newton_solve!`. Has no method here; a GPU package defines
it.
"""
function gpu_linear_solver end

# ============================================================================
# Engine
# ============================================================================

"""
    NewtonBuffers(g, g_trial, y_trial, dy, J)

Preallocated arrays for `newton_solve!`: residual at the iterate, residual and
point for line-search trials (also used as Jacobian scratch), Newton step, and
the Jacobian (with its sparsity pattern).
"""
struct NewtonBuffers{V<:AbstractVector, W<:AbstractVector, D<:AbstractVector, M}
    g::V
    g_trial::V
    y_trial::W
    dy::D
    J::M
end

_newton_result(converged, status, iterations, steps, residual_norm, refresh_due) =
    (; converged, status, iterations, steps, residual_norm, refresh_due)

# The norm of the line search's merit function: `Inf`, or 2 (a plain sum of squares, as the
# line search evaluates it at every trial; BLAS's overflow-safe norm costs more than the residual)
_merit_norm(g, p::Float64) = p == Inf ? norm(g, Inf) : p == 2 ? sqrt(sum(abs2, g)) : norm(g, p)

"""
    newton_solve!(residual!, jacobian!, y, bufs, solver, policy; abstol, maxiters,
                  catch_linear_failure=false, verbose=false, history=nothing,
                  fresh_jacobian=false, refresh_due=false)
        -> (converged, status, iterations, steps, residual_norm, refresh_due)

Solve `g(y) = 0` in place starting from `y`, with `residual!(g, y)` and
`jacobian!(J, y)`. Convergence is `‖g‖∞ < abstol` at the start of an iteration, or after the
last step. The residual last evaluated is that at `y` where the solve converged.

`status` is `:converged`, `:maxiters`, `:diverged` (non-finite residual or iterate),
`:stalled` (`policy.strict_search`: no decrease along the step of a Jacobian taken at the
iterate) or `:singular` (linear solve threw; only with `catch_linear_failure=true`,
otherwise the exception propagates). `iterations` counts the iterations begun (the one
that finds the residual below `abstol` included), `steps` the steps taken.

The solver's factorization is reused across calls where `policy.reuse_jacobian`. A step that
reused the Jacobian and contracted the residual by less than `policy.stall_ratio` leaves
`refresh_due` true: the Jacobian is taken again at the next iteration, of this call or, passed
as `refresh_due` to the next call, of that one. `fresh_jacobian`: the caller took the
factorization held at `y` (a Jacobian just taken is not taken again for a stall). Call
`reset_factorization!(solver)` for a fresh symbolic analysis. `history`, a vector, receives
the residual norm at the start and after each accepted step.
"""
function newton_solve!(residual!, jacobian!, y::AbstractVector, bufs::NewtonBuffers,
                       solver::AbstractNewtonLinearSolver, policy::NewtonPolicy;
                       abstol::Real, maxiters::Int,
                       catch_linear_failure::Bool = false, verbose::Bool = false,
                       history::Union{Nothing, Vector{Float64}} = nothing,
                       fresh_jacobian::Bool = false, refresh_due::Bool = false)
    (; g, g_trial, y_trial, dy, J) = bufs
    residual!(g, y)
    res_norm = norm(g, Inf)
    history === nothing || push!(history, res_norm)
    last_refresh = 0
    steps = 0
    fresh = fresh_jacobian   # whether the Jacobian of the last step was taken at its iterate
    due = refresh_due        # whether the Jacobian is to be taken again

    for iter in 1:maxiters
        @vprint "Iteration $iter: ‖g‖∞ = $res_norm"
        res_norm < abstol && return _newton_result(true, :converged, iter, steps, res_norm, due)
        isfinite(res_norm) || return _newton_result(false, :diverged, iter, steps, res_norm, due)

        refresh = due || !is_factorized(solver) ||
                  (iter == 1 ? !policy.reuse_jacobian : iter - last_refresh >= policy.jacobian_every)
        due = false

        # Newton step: J Δy = -g
        try
            if refresh
                jacobian!(J, y)
                factorize!(solver, J)
                last_refresh = iter
            end
            linear_solve!(solver, dy, g)
            dy .*= -1
        catch e
            catch_linear_failure || rethrow()
            @vprint "  Linear solve failed: $e"
            return _newton_result(false, :singular, iter, steps, res_norm, due)
        end
        fresh = refresh || fresh && iter == 1

        # Trust region
        step = norm(dy, Inf)
        step > policy.max_step && (dy .*= policy.max_step / step)

        # Backtracking line search; the last trial is accepted regardless unless strict
        alpha = 1.0
        accepted = false
        trial_norm = res_norm
        merit = policy.merit_norm == Inf ? res_norm : _merit_norm(g, policy.merit_norm)
        for trial in 1:policy.max_tries
            @. y_trial = y + (alpha * policy.damping) * dy
            @inbounds for a in policy.positive
                v = y[a]
                v > 0 && (y_trial[a] = max(y_trial[a], policy.positive_fraction * v))
            end
            residual!(g_trial, y_trial)
            trial_norm = norm(g_trial, Inf)
            trial_merit = policy.merit_norm == Inf ? trial_norm : _merit_norm(g_trial, policy.merit_norm)
            accepted = trial_merit < (1 - policy.sufficient_decrease * alpha) * merit
            (accepted || trial == policy.max_tries) && break
            alpha *= policy.backtrack_factor
        end
        steps += 1
        if !accepted && policy.strict_search
            @vprint "  no decrease along the step (Jacobian $(fresh ? "taken here" : "reused"))"
            fresh && return _newton_result(false, :stalled, iter, steps, res_norm, due)
            due = true
            continue
        end
        copyto!(y, y_trial)
        copyto!(g, g_trial)
        # (a Jacobian kept from an earlier iterate that contracts too little is taken again)
        due = !(fresh && policy.stall_reused_only) && trial_norm > policy.stall_ratio * res_norm
        fresh = false
        res_norm = trial_norm
        history === nothing || push!(history, res_norm)
        @vprint "  step: ‖Δy‖∞ = $step, alpha = $alpha$(refresh ? "" : " (Jacobian reused)")"

        any(!isfinite, y) && return _newton_result(false, :diverged, iter, steps, res_norm, due)
    end

    # The residual after the last step
    res_norm < abstol && return _newton_result(true, :converged, maxiters + 1, steps, res_norm, due)
    return _newton_result(false, :maxiters, maxiters, steps, res_norm, due)
end
