"""
    ptdf.jl

Power Transfer Distribution Factor (PTDF) computation.

PTDF[l, k] = change in flow on branch l per unit injection at bus k (with slack absorbing).

Mathematical formulation:
    PTDF = D * A * B^(-1)

where:
- D = diag(1/x_l) is the diagonal susceptance matrix (L x L)
- A = incidence matrix (L x N)
- B = bus susceptance matrix (N x N), with slack row/column removed

CRITICAL: Never compute B^(-1) explicitly. Use sparse factorization and solves.
"""

using PowerfulPF: get_slack_bus

# Scatter a slack-reduced vector into full bus coordinates (slack entry = 0)
function _expand_reduced!(full::AbstractVector, reduced::AbstractVector, non_slack::Vector{Int})
    fill!(full, zero(eltype(full)))
    for (i, orig_idx) in enumerate(non_slack)
        full[orig_idx] = reduced[i]
    end
    return full
end

# DC branch flows: flows[l] = b[l] * (theta[from(l)] - theta[to(l)])
function _dc_branch_flows!(flows::AbstractVector, b, from_indices, to_indices, theta)
    @inbounds for l in eachindex(flows)
        flows[l] = b[l] * (theta[from_indices[l]] - theta[to_indices[l]])
    end
    return flows
end

"""
    PTDF{T}

Power Transfer Distribution Factor matrix: `matrix[l, k]` is the change in flow on
branch `l` per unit (1 pu, i.e. 1 MW at the system base) of injection at bus `k`, with
the slack bus absorbing the balance.

# Fields
- `matrix::Matrix{T}`: dense L x N PTDF matrix
- `branch_indices::Vector{Tuple{Symbol, Int}}`: branch model mapping
- `slack_bus::Int`: slack bus internal index (reference)

Both `l` and `k` are *internal* array indices (1-based positions in the branch and bus
arrays), not the bus numbers or branch keys from the raw case data. A bus's internal
index can only be read back from a bus number with
`PowerfulCore.get_index_by_key(sys, :Bus, (bus_number,))`.

# Usage
```julia
ptdf = compute_ptdf(sys)
flow_change = ptdf.matrix[branch_idx, bus_idx] * injection
```
"""
struct PTDF{T<:Real}
    matrix::Matrix{T}
    branch_indices::Vector{Tuple{Symbol, Int}}
    slack_bus::Int
end

"""
    compute_ptdf(sys::SystemModel; slack_bus::Union{Int, Nothing}=nothing) -> PTDF

Compute the full PTDF matrix for a system.

# Arguments
- `sys::SystemModel`: System with topology analyzed
- `slack_bus::Union{Int, Nothing}`: Slack bus internal index (default: from topology)

# Returns
`PTDF` struct containing L x N dense matrix.

# Algorithm
1. Build B matrix from incidence: B = A' * D * A
2. Remove slack bus row/column to make B invertible
3. Factorize reduced B using KLU sparse solver
4. For each bus k, solve B * theta_k = e_k
5. Compute PTDF[:, k] = D * A * theta_k

# Example
```julia
sys = SystemModel(parse_network(load("ieee14").raw))
ptdf = compute_ptdf(sys)

# Flow change on branch 3 from 1 MW injection at bus 5
delta_flow = ptdf.matrix[3, 5] * 1.0
```
"""
function compute_ptdf(sys::SystemModel; slack_bus::Union{Int, Nothing}=nothing)
    # Use provided slack or get from system (uses get_slack_bus helper for consistency)
    slack_idx = if isnothing(slack_bus)
        get_slack_bus(sys)
    else
        slack_bus
    end

    # DC network: incidence, susceptances and factorized slack-reduced B
    (; A, branch_indices, b, from_indices, to_indices, non_slack, F) = dc_network(sys, slack_idx)
    L = length(b)
    N = size(A, 2)

    # Compute PTDF column by column
    # PTDF[:, k] = D * A * theta_k where B * theta_k = e_k
    ptdf_matrix = zeros(Float64, L, N)

    # Temporary buffers
    e_k = zeros(Float64, N - 1)
    theta_k_reduced = zeros(Float64, N - 1)
    theta_k_full = zeros(Float64, N)

    for k in 1:N
        if k == slack_idx
            # PTDF column for slack is zero (reference bus)
            continue
        end

        # Map k to reduced index
        k_reduced = k < slack_idx ? k : k - 1

        # Set up unit injection at bus k
        fill!(e_k, 0.0)
        e_k[k_reduced] = 1.0

        # Solve B_reduced * theta = e_k
        ldiv!(theta_k_reduced, F, e_k)

        # PTDF column: D * A * theta (slack angle = 0)
        _expand_reduced!(theta_k_full, theta_k_reduced, non_slack)
        _dc_branch_flows!(view(ptdf_matrix, :, k), b, from_indices, to_indices, theta_k_full)
    end

    return PTDF(ptdf_matrix, branch_indices, slack_idx)
end


"""
    LazyPTDF{T}

Lazy/virtual PTDF that computes rows or columns on demand.

For large systems, storing the full PTDF (L x N) is expensive.
LazyPTDF stores only the factorization and computes entries as needed.

# Fields
- `F`: KLU factorization of reduced B matrix
- `A::SparseMatrixCSC{Int8, Int}`: Incidence matrix
- `b::Vector{T}`: Branch susceptances (1/x)
- `slack_bus::Int`: Slack bus index
- `non_slack::Vector{Int}`: Non-slack bus indices
- `branch_indices::Vector{Tuple{Symbol, Int}}`: Branch model mapping
- `row_cache::Dict{Int, Vector{T}}`: Cached PTDF rows
- `col_cache::Dict{Int, Vector{T}}`: Cached PTDF columns

# Example
```julia
lptdf = LazyPTDF(sys)

# Get PTDF row for branch 5 (on-demand computation)
row = get_ptdf_row(lptdf, 5)

# Get PTDF column for bus 10
col = get_ptdf_column(lptdf, 10)
```
"""
mutable struct LazyPTDF{T<:Real}
    F::Any  # KLU factorization
    A::SparseMatrixCSC{Int8, Int}
    b::Vector{T}
    slack_bus::Int
    non_slack::Vector{Int}
    branch_indices::Vector{Tuple{Symbol, Int}}
    from_indices::Vector{Int}  # Pre-computed: from bus for each branch
    to_indices::Vector{Int}    # Pre-computed: to bus for each branch
    row_cache::Dict{Int, Vector{T}}
    col_cache::Dict{Int, Vector{T}}
    max_cache_size::Int
end

"""
    LazyPTDF(sys::SystemModel; slack_bus=nothing, max_cache_size=100) -> LazyPTDF

Create a virtual PTDF that computes entries on demand.
"""
function LazyPTDF(sys::SystemModel; slack_bus::Union{Int, Nothing}=nothing, max_cache_size::Int=100)
    # Use provided slack or get from system (uses get_slack_bus helper for consistency)
    slack_idx = if isnothing(slack_bus)
        get_slack_bus(sys)
    else
        slack_bus
    end

    (; A, branch_indices, b, from_indices, to_indices, non_slack, F) = dc_network(sys, slack_idx)

    return LazyPTDF(
        F, A, b, slack_idx, non_slack, branch_indices,
        from_indices, to_indices,
        Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(),
        max_cache_size
    )
end

"""
    LazyPTDF(F, A, b, branch_indices, from_indices, to_indices, slack_idx, non_slack; max_cache_size=100)

Internal constructor for creating LazyPTDF from pre-computed components INCLUDING factorization.
Used by screen_contingencies_dc to avoid redundant B matrix construction/factorization.

This is the fast path - F is already computed, no rebuild needed.
"""
function LazyPTDF(
    F,  # Pre-computed KLU factorization
    A::SparseMatrixCSC{Int8, Int},
    b::Vector{T},
    branch_indices::Vector{Tuple{Symbol, Int}},
    from_indices::Vector{Int},
    to_indices::Vector{Int},
    slack_idx::Int,
    non_slack::Vector{Int};
    max_cache_size::Int = 100
) where T
    return LazyPTDF(
        F, A, b, slack_idx, non_slack, branch_indices,
        from_indices, to_indices,
        Dict{Int, Vector{T}}(),
        Dict{Int, Vector{T}}(),
        max_cache_size
    )
end

"""
    get_ptdf_row(lptdf::LazyPTDF, branch_idx::Int) -> Vector{Float64}

Get PTDF row for a specific branch (sensitivity of branch flow to all bus injections).
"""
function get_ptdf_row(lptdf::LazyPTDF{T}, branch_idx::Int) where T
    # Check cache
    if haskey(lptdf.row_cache, branch_idx)
        return lptdf.row_cache[branch_idx]
    end

    N = size(lptdf.A, 2)
    row = zeros(T, N)

    # PTDF row l = b[l] * (e_from - e_to)' * B^(-1)
    # This requires solving B' * x = b[l] * (e_from - e_to)
    from_idx = lptdf.from_indices[branch_idx]
    to_idx = lptdf.to_indices[branch_idx]

    # For row computation, we need B^(-T) which equals B^(-1) since B is symmetric
    # Build RHS: b[l] * (e_from - e_to) for reduced system
    rhs = zeros(T, length(lptdf.non_slack))
    bl = lptdf.b[branch_idx]

    for (i, orig_idx) in enumerate(lptdf.non_slack)
        if orig_idx == from_idx
            rhs[i] = bl
        elseif orig_idx == to_idx
            rhs[i] = -bl
        end
    end

    # Solve B_reduced * x = rhs; slack column is zero
    _expand_reduced!(row, lptdf.F \ rhs, lptdf.non_slack)

    # Cache with eviction
    if length(lptdf.row_cache) >= lptdf.max_cache_size
        # Simple eviction: remove a random entry
        delete!(lptdf.row_cache, first(keys(lptdf.row_cache)))
    end
    lptdf.row_cache[branch_idx] = row

    return row
end

"""
    get_ptdf_column(lptdf::LazyPTDF, bus_idx::Int) -> Vector{Float64}

Get PTDF column for a specific bus (sensitivity of all branch flows to injection at bus).
"""
function get_ptdf_column(lptdf::LazyPTDF{T}, bus_idx::Int) where T
    # Check cache
    if haskey(lptdf.col_cache, bus_idx)
        return lptdf.col_cache[bus_idx]
    end

    L = size(lptdf.A, 1)
    col = zeros(T, L)

    # Slack bus: all sensitivities are zero
    if bus_idx == lptdf.slack_bus
        return col
    end

    # Solve B * theta = e_k
    bus_reduced_idx = findfirst(==(bus_idx), lptdf.non_slack)
    if isnothing(bus_reduced_idx)
        error("Bus $bus_idx not found in non-slack buses")
    end

    e_k = zeros(T, length(lptdf.non_slack))
    e_k[bus_reduced_idx] = one(T)

    theta_full = _expand_reduced!(zeros(T, size(lptdf.A, 2)), lptdf.F \ e_k, lptdf.non_slack)

    # Compute column: PTDF[l, k] = b[l] * (theta[from] - theta[to])
    _dc_branch_flows!(col, lptdf.b, lptdf.from_indices, lptdf.to_indices, theta_full)

    # Cache with eviction
    if length(lptdf.col_cache) >= lptdf.max_cache_size
        delete!(lptdf.col_cache, first(keys(lptdf.col_cache)))
    end
    lptdf.col_cache[bus_idx] = col

    return col
end

"""
    get_ptdf_element(lptdf::LazyPTDF, branch_idx::Int, bus_idx::Int) -> Float64

Get a single PTDF element (efficient for sparse access patterns).
"""
function get_ptdf_element(lptdf::LazyPTDF{T}, branch_idx::Int, bus_idx::Int) where T
    # Check if row or column is cached
    if haskey(lptdf.row_cache, branch_idx)
        return lptdf.row_cache[branch_idx][bus_idx]
    elseif haskey(lptdf.col_cache, bus_idx)
        return lptdf.col_cache[bus_idx][branch_idx]
    end

    # Compute single element without caching full row/column
    if bus_idx == lptdf.slack_bus
        return zero(T)
    end

    # Solve for theta at bus_idx only
    bus_reduced_idx = findfirst(==(bus_idx), lptdf.non_slack)
    e_k = zeros(T, length(lptdf.non_slack))
    e_k[bus_reduced_idx] = one(T)

    theta_reduced = lptdf.F \ e_k

    # Get theta at from/to buses (using pre-computed indices)
    from_idx = lptdf.from_indices[branch_idx]
    to_idx = lptdf.to_indices[branch_idx]

    theta_from = if from_idx == lptdf.slack_bus
        zero(T)
    else
        from_reduced = findfirst(==(from_idx), lptdf.non_slack)
        theta_reduced[from_reduced]
    end

    theta_to = if to_idx == lptdf.slack_bus
        zero(T)
    else
        to_reduced = findfirst(==(to_idx), lptdf.non_slack)
        theta_reduced[to_reduced]
    end

    return lptdf.b[branch_idx] * (theta_from - theta_to)
end
