# ============================================================================
# State-Only Constraint Reductions
# ============================================================================
#
# A zero-mass row of the mass-matrix DAE M du/dt = F(u) is an algebraic equation.
# The integrators solve such rows for the zero-mass variables, which requires each
# row to depend on them. A row that constrains differential states only does not;
# the DAE has index 2. Its time derivative depends on the zero-mass variables
# through the state equations, so replacing the row by it restores index 1. The
# rows are found once, when the adapter is built (`install_state_reduction!`,
# `adapter.jl`), by whichever algorithm a solver config registers through the
# `reduction_algorithm` extension point; plain TDS registers none. The functions
# here are generic over the result (`adapter.reduction::Union{Nothing,
# StateConstraintReduction}`, `adapter.jl`), regardless of which algorithm, if any,
# produced it.
#
# The original constraint A·x = 0 is linear and holds initially, and linear
# multistep and Runge-Kutta methods preserve linear invariants, so it stays
# satisfied. Events that make the states jump break it; `enforce_state_constraints!`
# restores it.
# ============================================================================

export enforce_state_constraints!, constraint_values

"""
    tds_jacobian_sparsity(adapter) -> Union{Nothing, SparseMatrixCSC}

Sparsity of the DAE Jacobian as the adapter evaluates it: the reduced pattern when
rows are replaced by their derivative, `nothing` (the model pattern) otherwise.
"""
tds_jacobian_sparsity(adapter::TDSAdapter) =
    adapter.reduction === nothing ? nothing : adapter.reduction.pattern

"""
    constraint_values(adapter, u, t) -> Vector{Float64}

The value `c(t)·x + b` of each constraint of the adapter's reduction (`reduction.rows`) at the state `u` and time `t`, read
from the residual as the models define it (the weighted sum of the rows of its combination); `NaN` for one whose modes
the automata are not in (its gate is closed, its row the original equation). Empty without a reduction.
"""
function constraint_values(adapter::TDSAdapter, u::AbstractVector, t::Real)
    red = adapter.reduction
    red === nothing && return Float64[]
    du = zeros(length(u))
    _evaluate_unreduced!(du, adapter, u, t, nothing, adapter.status)
    return [gate_open(red.gates[j]) ? Float64(_constraint_value(red, du, j)) : NaN for j in eachindex(red.rows)]
end

"""
    enforce_state_constraints!(u, adapter, t; atol=1e-10) -> Bool

Restore the state-only constraints of the adapter's reduction at `u`, the state at time `t`, after an event
that made states jump, such as a fault current reset when a fault clears at a bus
whose every connection is an inductor, or a mode that holds a constraint is entered (the constraints of
the modes the automata are in, `gate_open`). Returns whether `u` changed.

The reduced rows keep `A(t)·x + b` at its value, so a jump that leaves `A(t)·x + b = c ≠ 0` would persist
(`A(t)` the constraint coefficients at `t`, `constraint_matrix`; `c` the value of the constraint
read from the residual, `_constraint_value`).
The correction is the one that conserves the flux of every inductor loop: the change
`Δx` minimizes `Δxᵀ M Δx` subject to `A(x + Δx) = 0`, i.e.
`Δx = -M⁻¹Aᵀ (A M⁻¹ Aᵀ)⁻¹ c` with `A = A(t)`, over the states of in-service instances (the states of
out-of-service instances, e.g. a cleared fault's reset currents, keep their values).
This is `project_states!` with the `MinimumNorm` policy and the inverse masses as weights, which moves every state
of a constraint; a mode that governs a state uses the `AlongStates` policy (governed_states.jl).
"""
function enforce_state_constraints!(u::AbstractVector, adapter::TDSAdapter, t::Real; atol::Real = 1e-10)
    red = adapter.reduction
    red === nothing && return false
    n_diff = adapter.n_diff
    # (the constraints that hold now: a combination's gate is open in the modes it exists in)
    active = [j for j in eachindex(red.rows) if gate_open(red.gates[j])]
    isempty(active) && return false

    du = zeros(length(u))
    _evaluate_unreduced!(du, adapter, u, t, nothing, adapter.status)
    c = [_constraint_value(red, du, j) for j in active]
    maximum(abs, c) <= atol && return false

    mass = diag(adapter.mass_matrix)
    inv_mass = [!iszero(mass[k]) && (red.owner_model[k] === :none ||
                adapter.status[red.owner_model[k]][red.owner_index[k]]) ? 1 / mass[k] : 0.0
                for k in 1:n_diff]
    return project_states!(u, constraint_matrix(red, t)[active, :], c, MinimumNorm(inv_mass))
end

"""
    enforce_state_constraints!(integrator, adapter) -> Bool

Restore the state-only constraints at the integrator's state after an event, and mark
the state as modified when it changed.
"""
function enforce_state_constraints!(integrator, adapter::TDSAdapter)
    changed = enforce_state_constraints!(integrator.u, adapter, integrator.t)
    changed && DiffEq.u_modified!(integrator, true)
    return changed
end
