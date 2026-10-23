"""
Shared Utilities for Linearization and Small-Signal Analysis

This module provides common infrastructure for:
- Eigenvalue analysis (small-signal stability)
- Linearization/Bode analysis (frequency response)

Both analyses share:
- State name collection from AddressManager
- Zero-mass state reduction (DAE to ODE conversion)
- Jacobian sub-matrix extraction
"""

using LinearAlgebra: Diagonal, lu
using PowerfulCore: StateVar, variable_name, prepare_sparse_jacobian

# ============================================================================
# Result Types
# ============================================================================

"""
    JacobianSubMatrices{T, M<:AbstractMatrix{T}}

Sub-matrices extracted from full DAE Jacobian.

Parametric on matrix type M to preserve sparsity when extracted from sparse Jacobian.

# Fields
- `fx::M` - ∂f/∂x (n_diff × n_diff)
- `fy::M` - ∂f/∂y (n_diff × n_alg)
- `gx::M` - ∂g/∂x (n_alg × n_diff)
- `gy::M` - ∂g/∂y (n_alg × n_alg)
"""
struct JacobianSubMatrices{T, M<:AbstractMatrix{T}}
    fx::M
    fy::M
    gx::M
    gy::M
end

"""
    ZeroMassReductionResult{T}

Result of reducing the DAE `M ẋ = f(x, y, u)`, `0 = g(x, y, u)` to the ODE
`ẋ_r = A x_r + B u` of the states with nonzero mass (`x_r`).

The eliminated quantities `w = [y; z]` (the algebraic variables `y` in address
order, then the zero-mass states `z` in `zero_mass_indices` order) follow from
`x_r` and `u` as `w = -(W_inv_Wx x_r + W_inv_Wu u)`.

# Fields
- `A::Matrix{T}` - Reduced state matrix
- `B::Union{Nothing, Matrix{T}}` - Reduced input matrix (nothing if not computed)
- `retained_indices::Vector{Int}` - State addresses of `x_r`, in order
- `zero_mass_indices::Vector{Int}` - State addresses with zero mass coefficient (in `w`)
- `reduced_state_names::Vector{Symbol}` - State names of `x_r`
- `W_inv_Wx::Matrix{T}` - `W⁻¹ W_x`, the negated sensitivity of `w` to `x_r`
- `W_inv_Wu::Union{Nothing, Matrix{T}}` - `W⁻¹ W_u`, the negated sensitivity of `w`
  to `u` (nothing if not computed)
"""
struct ZeroMassReductionResult{T}
    A::Matrix{T}
    B::Union{Nothing, Matrix{T}}
    retained_indices::Vector{Int}
    zero_mass_indices::Vector{Int}
    reduced_state_names::Vector{Symbol}
    W_inv_Wx::Matrix{T}
    W_inv_Wu::Union{Nothing, Matrix{T}}
end

# ============================================================================
# State Name Collection
# ============================================================================

"""
    collect_state_names(am, n_diff::Int) -> Vector{Symbol}

Collect state variable names from AddressManager in address order.

# Arguments
- `am`: AddressManager containing variable metadata
- `n_diff`: Number of differential states

# Returns
Vector of state names in address order.
"""
function collect_state_names(am, n_diff::Int)
    names = Vector{Symbol}(undef, n_diff)
    for addr in 1:n_diff
        names[addr] = variable_name(am, StateVar(), UInt(addr))
    end
    return names
end

# ============================================================================
# Jacobian Sub-matrix Extraction
# ============================================================================

"""
    extract_jacobian_submatrices(J_full, n_diff) -> JacobianSubMatrices

Extract sub-matrices from full Jacobian of combined [x; y] system.

# Arguments
- `J_full`: Full Jacobian of size (n_diff + n_alg) × (n_diff + n_alg)
- `n_diff`: Number of differential states

# Returns
`JacobianSubMatrices` containing fx, fy, gx, gy blocks.
"""
function extract_jacobian_submatrices(J_full::M, n_diff::Int) where {T, M<:AbstractMatrix{T}}
    fx = J_full[1:n_diff, 1:n_diff]
    fy = J_full[1:n_diff, n_diff+1:end]
    gx = J_full[n_diff+1:end, 1:n_diff]
    gy = J_full[n_diff+1:end, n_diff+1:end]
    # Let Julia infer the matrix type from the slices
    return JacobianSubMatrices(fx, fy, gx, gy)
end

# ============================================================================
# Zero-Mass State Reduction
# ============================================================================

"""
    dae_jacobian(adapter::TDSAdapter, t) -> (J, u0)

Sparse Jacobian ∂[f;g]/∂[x;y] of the DAE at the adapter's current state
`u0 = pack_state(adapter)` and time `t`, the time of that state (the operating point for
linearization; for the state a `TDSAdapter` is built with, the start of its problem,
`prob.tspan[1]`). `t` has no default: a model's residuals can depend on it directly.
"""
function dae_jacobian(adapter::TDSAdapter, t::Real)
    u0 = pack_state(adapter)
    du = zeros(length(u0))
    closure = JacobianClosureFunctor(Ref{Any}(adapter), Ref(Float64(t)))
    (; J, backend, prep) = prepare_sparse_jacobian(closure, du, u0, adapter.adcache; include_differential=true,
                                                   sparsity=tds_jacobian_sparsity(adapter))
    DifferentiationInterface.jacobian!(closure, du, J, prep, backend, u0)
    return J, u0
end

"""
    reduce_zero_mass_states(fx, fy, gx, gy, M_diag, state_names; f_u=nothing, g_u=nothing)

Reduce the linearized DAE to an ODE in the states with nonzero mass.

States with zero mass coefficient (`M_diag[i] == 0`) are algebraic in effect: with
the algebraic variables `y` they form the eliminated quantities `w = [y; z]`,
solved from `[g; f_z] = 0`:

    W w = -(W_x x_r + W_u u),   W = [gy  gx_z; fy_z  fx_zz]

so `w = -(W⁻¹W_x x_r + W⁻¹W_u u)`, and

    A = M_r⁻¹ (fx_rr - f_rw W⁻¹W_x),   B = M_r⁻¹ (f_u,r - f_rw W⁻¹W_u)

where `f_rw = [fy_r  fx_rz]` are the rows of the retained states.

# Arguments
- `fx, fy, gx, gy`: Jacobian sub-matrices from ∂[f;g]/∂[x;y]
- `M_diag`: Diagonal of mass matrix for differential states
- `state_names`: Names for each differential state

# Keyword Arguments
- `f_u, g_u`: Input Jacobians (optional). If provided, `B` and `W_inv_Wu` are computed.

# Returns
`ZeroMassReductionResult`
"""
function reduce_zero_mass_states(
    fx, fy, gx, gy,
    M_diag::AbstractVector,
    state_names::AbstractVector{Symbol};
    f_u::Union{Nothing, AbstractMatrix} = nothing,
    g_u::Union{Nothing, AbstractMatrix} = nothing
)
    compute_B = !isnothing(f_u) && !isnothing(g_u)

    zero_mass_idx = findall(iszero, M_diag)
    retained_idx = findall(!iszero, M_diag)
    isempty(retained_idx) && error("All states have zero mass: the reduced system has no states")

    r, z = retained_idx, zero_mass_idx
    M_inv = Diagonal(inv.(view(M_diag, r)))
    fx_rr, f_rw, W, W_x = _reduction_blocks(fx, fy, gx, gy, r, z)

    # Only the right-hand side of the sparse solve needs to be dense
    W_fact = lu(W)
    W_inv_Wx = W_fact \ Matrix(W_x)
    A = M_inv * (fx_rr - f_rw * W_inv_Wx)

    W_inv_Wu, B = if compute_B
        f_u_r, W_u = isempty(z) ? (f_u, g_u) : (f_u[r, :], [g_u; f_u[z, :]])
        W_inv_Wu = W_fact \ Matrix(W_u)
        (W_inv_Wu, M_inv * (f_u_r - f_rw * W_inv_Wu))
    else
        (nothing, nothing)
    end

    return ZeroMassReductionResult(A, B, r, z, state_names[r], W_inv_Wx, W_inv_Wu)
end

"""
    _reduction_blocks(fx, fy, gx, gy, r, z) -> (fx_rr, f_rw, W, W_x)

Blocks of the eliminated system of `reduce_zero_mass_states` for the retained states
`r` and the zero-mass states `z`: `fx_rr`, `f_rw = [fy_r fx_rz]`,
`W = [gy gx_z; fy_z fx_zz]` and `W_x = [gx_r; fx_zr]`. Without zero-mass states they are
the Jacobian blocks themselves (no copies). The blocks are linear in the Jacobian, so
the same function gives the blocks of a Jacobian derivative.
"""
_reduction_blocks(fx, fy, gx, gy, r, z) = isempty(z) ? (fx, fy, gy, gx) :
    (fx[r, r], [fy[r, :] fx[r, z]], [gy gx[:, z]; fy[z, :] fx[z, z]], [gx[:, r]; fx[z, r]])

"""
    reduced_state_matrix_derivative(red, blocks, W_fact, M_diag, d_jac, dM_diag) -> Matrix

Derivative `dA` of the state matrix `A = M_r⁻¹ (fx_rr - f_rw W⁻¹ W_x)` of
`reduce_zero_mass_states` (`red`) along a direction in which the Jacobian blocks change
by `d_jac = (fx, fy, gx, gy)` and the mass diagonal by `dM_diag`:

    dX = W⁻¹ (dW_x - dW X),   X = W⁻¹ W_x
    dA = -M_r⁻¹ dM_r A + M_r⁻¹ (dfx_rr - df_rw X - f_rw dX)

`blocks` are the blocks of the Jacobian (`_reduction_blocks`) and `W_fact` a
factorization of their `W`. Errors if a state with zero mass gains mass along the
direction, which changes the structure of the reduction.
"""
function reduced_state_matrix_derivative(red::ZeroMassReductionResult, blocks, W_fact,
                                         M_diag::AbstractVector, d_jac, dM_diag::AbstractVector)
    r, z = red.retained_indices, red.zero_mass_indices
    all(iszero, view(dM_diag, z)) ||
        error("A state with zero mass gains mass: the reduction to the retained states changes")
    _, f_rw, _, _ = blocks
    dfx_rr, df_rw, dW, dW_x = _reduction_blocks(d_jac..., r, z)

    X = red.W_inv_Wx
    dX = W_fact \ (Matrix(dW_x) - dW * X)
    M_inv = Diagonal(inv.(view(M_diag, r)))
    return -(M_inv * Diagonal(view(dM_diag, r))) * red.A + M_inv * (dfx_rr - df_rw * X - f_rw * dX)
end
