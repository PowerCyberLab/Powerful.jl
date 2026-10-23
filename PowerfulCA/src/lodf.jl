"""
    lodf.jl

Line Outage Distribution Factor (LODF) computation.

LODF[m, l] = fraction of outaged line l's flow that redistributes to line m.

Mathematical formulation:
    H = PTDF * A'           # L x L matrix
    h = diag(H)             # Self-sensitivities
    LODF = H ./ (1 - h')    # Element-wise division
    LODF[diag] = -1         # Outaged line loses all flow

Post-outage flow estimation:
    P_new[m] = P_base[m] + LODF[m, l] * P_base[l]
"""

using LinearAlgebra: diag

"""
    LODF{T}

Line Outage Distribution Factor matrix.

# Fields
- `matrix::Matrix{T}`: Dense L x L LODF matrix
- `branch_indices::Vector{Tuple{Symbol, Int}}`: Branch model mapping

# Usage
```julia
lodf = compute_lodf(ptdf, A)

# Post-outage flow on branch m when branch l is outaged
P_new_m = P_base_m + lodf.matrix[m, l] * P_base_l
```
"""
struct LODF{T<:Real}
    matrix::Matrix{T}
    branch_indices::Vector{Tuple{Symbol, Int}}
end

"""
    compute_lodf(ptdf::PTDF, A::SparseMatrixCSC{Int8}) -> LODF

Compute LODF matrix from PTDF and incidence matrix.

# Arguments
- `ptdf::PTDF`: Pre-computed PTDF matrix
- `A::SparseMatrixCSC{Int8}`: Branch-bus incidence matrix

# Returns
`LODF` struct containing L x L dense matrix.

# Algorithm
```
H = PTDF * A'           # L x L: sensitivity of flow l to injection pattern of branch m
h = diag(H)             # Self-sensitivity
LODF[m, l] = H[m, l] / (1 - H[l, l])  for m ≠ l
LODF[l, l] = -1         # Outaged line loses all flow
```

# Example
```julia
ptdf = compute_ptdf(sys)
A, _ = build_incidence_matrix(sys)
lodf = compute_lodf(ptdf, A)

# Check: when line 3 outages, how much flows to line 5?
redistribution_fraction = lodf.matrix[5, 3]
```
"""
function compute_lodf(ptdf::PTDF{T}, A::SparseMatrixCSC{Int8}) where T
    # H = PTDF * A' (L x L)
    # H[m, l] = sum_k PTDF[m, k] * A[l, k]
    #         = PTDF[m, from(l)] - PTDF[m, to(l)]
    # This is the sensitivity of flow on m to the injection pattern of branch l
    H = ptdf.matrix * A'

    L = size(H, 1)
    lodf_matrix = zeros(T, L, L)

    for l in 1:L
        _lodf_column!(view(lodf_matrix, :, l), view(H, :, l), l)
    end

    return LODF(lodf_matrix, ptdf.branch_indices)
end

"""
    compute_lodf(sys::SystemModel) -> LODF

Convenience method: compute LODF directly from system.
"""
function compute_lodf(sys::SystemModel)
    ptdf = compute_ptdf(sys)
    A, _ = build_incidence_matrix(sys)
    return compute_lodf(ptdf, A)
end

"""
    compute_lodf_column(lptdf::LazyPTDF, outaged_branch::Int) -> Vector{Float64}

Compute a single LODF column (flow redistribution when branch is outaged).

More efficient than computing full LODF when only a few outages are of interest.

# Arguments
- `lptdf::LazyPTDF`: Virtual PTDF for on-demand computation
- `outaged_branch::Int`: Index of branch being outaged

# Returns
Vector of length L: LODF[:, outaged_branch]
"""
function compute_lodf_column(lptdf::LazyPTDF{T}, outaged_branch::Int) where T
    L = size(lptdf.A, 1)
    lodf_col = zeros(T, L)

    # Get endpoints of outaged branch (use pre-computed indices for O(1) access)
    from_idx = lptdf.from_indices[outaged_branch]
    to_idx = lptdf.to_indices[outaged_branch]

    # H[:, l] = PTDF[:, from(l)] - PTDF[:, to(l)]
    # Get PTDF columns for from and to buses
    ptdf_from = get_ptdf_column(lptdf, from_idx)
    ptdf_to = get_ptdf_column(lptdf, to_idx)

    # H column for outaged branch
    return _lodf_column!(lodf_col, ptdf_from - ptdf_to, outaged_branch)
end

"""
    _lodf_column!(lodf_col, H_col, l) -> lodf_col

LODF column for outaged branch `l` from `H_col = H[:, l]` (H = PTDF · Aᵀ):
`H[m, l] / (1 - H[l, l])` for m ≠ l and -1 for the outaged branch itself. A
radial branch (1 - H[l, l] ≈ 0, outage disconnects the system) gets NaN.
"""
function _lodf_column!(lodf_col::AbstractVector{T}, H_col::AbstractVector, l::Int) where T
    divisor = one(T) - H_col[l]
    if abs(divisor) < 1e-10
        lodf_col .= T(NaN)
    else
        lodf_col .= H_col ./ divisor
    end
    lodf_col[l] = -one(T)  # Outaged line loses all flow
    return lodf_col
end

# Post-outage flows: base flows plus the outaged branch's flow redistributed by its LODF column
_post_outage_flows(base_flows::AbstractVector, lodf_col::AbstractVector, outaged_branch::Int) =
    base_flows .+ lodf_col .* base_flows[outaged_branch]

"""
    estimate_post_outage_flows(base_flows::Vector{T}, lodf::LODF{T}, outaged_branch::Int) -> Vector{T}

Estimate branch flows after an outage using LODF.

# Arguments
- `base_flows::Vector{T}`: Pre-outage branch flows (MW)
- `lodf::LODF{T}`: Pre-computed LODF matrix
- `outaged_branch::Int`: Index of branch being outaged

# Returns
Estimated post-outage flows for all branches.

# Example
```julia
# Get base case flows from power flow solution
base_flows = [get branch flows from solution]

# Estimate flows after outaging branch 5
post_outage = estimate_post_outage_flows(base_flows, lodf, 5)

# Check for overloads
for (i, flow) in enumerate(post_outage)
    if abs(flow) > rating[i]
        println("Branch \$i overloaded: \$flow MW")
    end
end
```
"""
function estimate_post_outage_flows(
    base_flows::Vector{T},
    lodf::LODF{T},
    outaged_branch::Int
) where T
    return _post_outage_flows(base_flows, view(lodf.matrix, :, outaged_branch), outaged_branch)
end

"""
    estimate_post_outage_flows(base_flows::Vector{T}, lptdf::LazyPTDF{T}, outaged_branch::Int) -> Vector{T}

Estimate post-outage flows using virtual PTDF (computes LODF column on demand).
"""
function estimate_post_outage_flows(
    base_flows::Vector{T},
    lptdf::LazyPTDF{T},
    outaged_branch::Int
) where T
    return _post_outage_flows(base_flows, compute_lodf_column(lptdf, outaged_branch), outaged_branch)
end
