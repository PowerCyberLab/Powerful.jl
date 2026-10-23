# ============================================================================
# Projection of Differential States onto Linear Constraints
# ============================================================================
#
# One routine, `project_states!`, for the two places where an event leaves states off a linear
# relation they must satisfy, and puts them back. The policies differ in which states move:
#
# - `MinimumNorm`: the state-only constraints of a reduced DAE (index_reduction.jl) after an
#   event that made states jump. With the constraint `A x = b`, its residual `c = A x - b` at the
#   current `x` and a weight `w_k ≥ 0` for every state (zero: it does not move), the correction
#   of least `Σ_k d_k² / w_k` that satisfies it is
#
#       d = -W Aᵀ (A W Aᵀ)⁻¹ c,        W = diag(w)
#
#   (derived with SymPy: `A (x + d) = b` for any `A` and `w`). With `w_k = 1 / M_k`, the inverse
#   mass, it conserves the flux of every inductor loop.
#
# - `AlongStates`: the states a mode governs (automaton_runtime.jl, `Mode.states`) where a
#   transition enters the mode. Only the designated states move, each of its constraint's rows
#   holding one, and the others' contribution is in the right side `b`: the weights of the
#   formula above are 1 for the designated states and 0 for the others, so `A W Aᵀ = A_S A_Sᵀ`
#   (`A_S` the columns of the designated states) and `d = -A_S⁻¹ c`, i.e. the designated states
#   become `x_S = A_S⁻¹ b`. They are solved from `b` rather than moved by `d`: for a relation
#   `x_g = r(x)` the state is then `r` itself, not `x_g - (x_g - r)`, which rounds.
# ============================================================================

using SparseArrays: AbstractSparseMatrix

"""
    ProjectionPolicy

Which states `project_states!` moves: `MinimumNorm` or `AlongStates`.
"""
abstract type ProjectionPolicy end

"""
    MinimumNorm(weights)

The policy of the least weighted correction `d = -W Aᵀ (A W Aᵀ)⁻¹ c` over every state with a nonzero weight
(`weights[k]` for state `k`).
"""
struct MinimumNorm{V <: AbstractVector{Float64}} <: ProjectionPolicy
    weights::V
end

"""
    AlongStates(states)

The policy that moves the designated `states` (indices into `x`) alone, one per row of the constraint, to
the solution of `A_S x_S = b`.
"""
struct AlongStates <: ProjectionPolicy
    states::Vector{Int}
end

"""
    project_states!(x, A, c, policy::MinimumNorm) -> Bool

Add to the differential states `x` (its first `size(A, 2)` entries) the correction `d = -W Aᵀ (A W Aᵀ)⁻¹ c`
that satisfies the constraint `A x = b` whose residual at `x` is `c = A x - b`, moving only the states with a
nonzero weight (see the header). The rows of `A` must be independent over them: a singular `A W Aᵀ` throws.
Returns whether any entry of `x` changed.

    project_states!(x, A, b, policy::AlongStates) -> Bool

Set the designated states of `policy` to the solution of `A x_S = b`: `A` is square, a row for each designated
state and its column; what the other states contribute to the constraint is in `b` (a relation `x_g = r(x)`:
`A = 1` and `b = r`). Returns whether any entry of `x` changed.
"""
function project_states!(x::AbstractVector, A::AbstractSparseMatrix, c::AbstractVector, policy::MinimumNorm)
    isempty(c) && return false
    AW = A * Diagonal(policy.weights)
    lambda = (AW * transpose(A)) \ c
    d = -(transpose(AW) * lambda)
    changed = false
    @inbounds for k in eachindex(d)
        iszero(d[k]) && continue
        x[k] += d[k]
        changed = true
    end
    return changed
end

function project_states!(x::AbstractVector, A::AbstractSparseMatrix, b::AbstractVector, policy::AlongStates)
    isempty(b) && return false
    size(A) == (length(b), length(policy.states)) || throw(DimensionMismatch(
        "A governed-state projection solves one row per designated state: A is $(size(A)) for " *
        "$(length(policy.states)) states and a right side of $(length(b))"))
    solved = A \ b
    changed = false
    @inbounds for (k, j) in enumerate(policy.states)
        x[j] == solved[k] && continue
        x[j] = solved[k]
        changed = true
    end
    return changed
end
