"""
Automatic Jacobian sparsity pattern detection from model structure.

This module analyzes the system's address manager and model metadata to
determine which residuals depend on which variables, producing the sparsity
pattern of the Jacobian matrix ∂g/∂y.
"""

using SparseArrays
using SparseMatrixColorings
using ADTypes: AbstractSparsityDetector
import ADTypes

export jacobian_sparsity_pattern, instance_residual_rows, matrix_colors, PrecomputedSparsityDetector

"""
    jacobian_sparsity_pattern(am, sys, formulation; include_differential=false,
                              models=keys(sys.models)) -> SparseMatrixCSC{Bool}

Compute Jacobian sparsity pattern from model connectivity (conservative estimate).

# Arguments
- `am`: AddressManager with variable/residual positions
- `sys`: SystemModel
- `formulation`: PowerInjection or CurrentInjection
- `include_differential`: Include StateVar/StateRes for full DAE Jacobian
- `models`: Names of the models whose residuals enter the pattern (all by default)
"""
function jacobian_sparsity_pattern(
    am::AddressManager, sys::SystemModel, formulation::FormulationTrait;
    include_differential::Bool = false,
    models = keys(sys.models)
)
    # Get dimensions
    n_diff = include_differential ? Int(max_address_in_space(am, StateVar())) : 0
    n_alg = Int(max_address_in_space(am, AlgebVar()))
    m = n_diff + n_alg  # Total matrix size

    # Per model: the residual rows of each instance and the variable columns of each
    # instance; every residual of instance i depends on every variable of instance i
    blocks = Tuple{Matrix{Int}, Vector{Vector{Int}}}[]
    nnz_total = 0

    for model_name in models
        model_vec = sys.models[model_name]
        # Skip empty models
        length(model_vec) == 0 && continue

        form_metadata = get_metadata_for_formulation(
            model_metadata(get_base_model_type(model_vec)), formulation)

        # Collect variable addresses for this model with offset mapping
        # Each instance i has addresses for its variables
        var_addr_lists = Vector{Vector{Int}}()

        for var in form_metadata.vars
            # Determine offset based on address type
            if var.address_type isa StateVar
                include_differential || continue  # Skip for power flow
                offset = 0  # StateVar: columns 1:n_diff
            elseif var.address_type isa AlgebVar
                offset = n_diff  # AlgebVar: columns n_diff+1:n_total
            else
                continue
            end

            addrs = item_addresses(am, model_name, var)

            # Apply offset mapping
            mapped_addrs = [Int(a) + offset for a in addrs]
            push!(var_addr_lists, mapped_addrs)
        end

        rows = _instance_residual_rows(am, model_name, length(model_vec), form_metadata,
                                       n_diff, include_differential)
        for var_addrs in var_addr_lists
            nnz_total += size(rows, 1) * min(length(var_addrs), size(rows, 2))
        end
        push!(blocks, (rows, var_addr_lists))
    end

    # (row, col) pairs, allocated once at their final length
    I = Vector{Int}(undef, nnz_total)  # Residual indices (rows)
    J = Vector{Int}(undef, nnz_total)  # Variable indices (columns)
    k = 0
    for (rows, var_addr_lists) in blocks
        k = _fill_pattern_pairs!(I, J, k, rows, var_addr_lists)
    end

    # Build sparse boolean matrix (an instance not delivered through a residual has no row in it)
    # Note: sparse() automatically handles duplicate (i,j) pairs
    resize!(I, k)
    resize!(J, k)
    return sparse(I, J, trues(k), m, m)
end

# Write the (row, col) pairs of one model after position k of I and J: each residual row
# of instance i with each variable column of instance i. Returns the last position.
function _fill_pattern_pairs!(I::Vector{Int}, J::Vector{Int}, k::Int, rows::Matrix{Int},
                              var_addr_lists::Vector{Vector{Int}})
    @inbounds for i in axes(rows, 2), r in 1:size(rows, 1), var_addrs in var_addr_lists
        if i <= length(var_addrs) && rows[r, i] != NO_ROW
            k += 1
            I[k] = rows[r, i]
            J[k] = var_addrs[i]
        end
    end
    return k
end

"""
    instance_residual_rows(am, sys, formulation, model_name; include_differential=false)
        -> Matrix{Int}

Rows of the DAE residual `[f; g]` that each instance of `model_name` writes: column `i`
holds those of instance `i` (`0` where the instance is not delivered through the residual),
one row per residual of the model: its differential
residuals (rows `1:n_diff`, only with `include_differential`) and its algebraic
residuals (offset by `n_diff`), including shared residuals such as bus balances.

A quantity read only by the equations of one instance (its parameters and setpoints)
enters exactly these rows.
"""
function instance_residual_rows(
    am::AddressManager, sys::SystemModel, formulation::FormulationTrait, model_name::Symbol;
    include_differential::Bool = false
)
    model_vec = sys.models[model_name]
    n = length(model_vec)
    n == 0 && return Matrix{Int}(undef, 0, 0)

    n_diff = include_differential ? Int(max_address_in_space(am, StateVar())) : 0
    form_metadata = get_metadata_for_formulation(
        model_metadata(get_base_model_type(model_vec)), formulation)
    return _instance_residual_rows(am, model_name, n, form_metadata, n_diff, include_differential)
end

# Whether a residual has rows in the pattern: algebraic always, differential with
# `include_differential`
_in_pattern(residual, include_differential::Bool) =
    (include_differential && residual.address_type isa StateRes) || residual.address_type isa AlgebRes

# `instance_residual_rows` for the `n` instances of a model with formulation metadata
# `form_metadata`, the differential rows taking the first `n_diff` rows
function _instance_residual_rows(am::AddressManager, model_name::Symbol, n::Int, form_metadata,
                                 n_diff::Int, include_differential::Bool)
    n_rows = count(r -> _in_pattern(r, include_differential), form_metadata.residuals)
    rows = Matrix{Int}(undef, n_rows, n)
    k = 0
    for residual in form_metadata.residuals
        _in_pattern(residual, include_differential) || continue
        k += 1
        offset = residual.address_type isa StateRes ? 0 : n_diff
        res_addrs = item_addresses(am, model_name, residual)
        length(res_addrs) == n || error("Residual $(residual.name) of $model_name has " *
            "$(length(res_addrs)) addresses for $n instances")
        _fill_rows!(rows, k, res_addrs, offset)
    end
    return rows
end

# The row of an instance in a residual it is not delivered through (`UNROUTED`)
const NO_ROW = 0

# Row k of `rows`: the addresses shifted by `offset` (a function barrier: the residuals
# of a model have different types)
function _fill_rows!(rows::Matrix{Int}, k::Int, addrs::AbstractVector, offset::Int)
    @inbounds for i in eachindex(addrs)
        rows[k, i] = addrs[i] == UNROUTED ? NO_ROW : Int(addrs[i]) + offset
    end
    return nothing
end

"""
    matrix_colors(sparsity::SparseMatrixCSC) -> Vector{Int}

Compute graph coloring for sparse Jacobian matrix.

Uses SparseMatrixColorings.jl to compute a distance-2 coloring of the sparsity
pattern, which enables efficient compressed ForwardDiff Jacobian computation.

# Arguments
- `sparsity::SparseMatrixCSC`: Boolean or numeric sparsity pattern

# Returns
- `Vector{Int}`: Color assignment for each column (1-indexed)

# Example
```julia
sparsity = jacobian_sparsity_pattern(sys)
colors = matrix_colors(sparsity)  # Column colors for compression
```
"""
function matrix_colors(sparsity::SparseMatrixCSC)
    # Use greedy coloring algorithm for column coloring
    # This is the same algorithm SparseDiffTools used
    result = coloring(sparsity, ColoringProblem(), GreedyColoringAlgorithm())
    return column_colors(result)
end

# ==============================================================================
# PrecomputedSparsityDetector - Fast Sparsity Integration with DifferentiationInterface
# ==============================================================================

"""
    PrecomputedSparsityDetector(pattern::AbstractSparseMatrix)

Custom sparsity detector for DifferentiationInterface.jl that uses a pre-computed
sparsity pattern instead of automatic tracing.

This provides a ~38x speedup over `TracerSparsityDetector` by using the fast
`jacobian_sparsity_pattern()` function which analyzes model structure directly.

# Performance Comparison (ACTIVSg 70k bus system)
- `TracerSparsityDetector()`: ~4.5 seconds
- `PrecomputedSparsityDetector()`: ~0.12 seconds (38x faster!)

# Example
```julia
using DifferentiationInterface
using DifferentiationInterface: AutoSparse
using ADTypes: AutoForwardDiff
using SparseMatrixColorings: GreedyColoringAlgorithm

# Fast sparsity detection from model structure
sparsity = jacobian_sparsity_pattern(sys)

# Create custom detector
detector = PrecomputedSparsityDetector(sparsity)

# Use with DifferentiationInterface
backend = AutoSparse(
    AutoForwardDiff();
    sparsity_detector = detector,  # Use our fast detector!
    coloring_algorithm = GreedyColoringAlgorithm()
)

# Prepare and compute Jacobian as usual
prep = prepare_jacobian(residual!, zero(g), backend, zero(y))
jacobian!(residual!, g, J, prep, backend, y)
```

# See Also
- `jacobian_sparsity_pattern`: Computes sparsity from model structure (fast)
- `TracerSparsityDetector`: Automatic tracing (slow but fully automatic)
"""
struct PrecomputedSparsityDetector{M<:AbstractMatrix} <: AbstractSparsityDetector
    pattern::M  # SparseMatrixCSC{Bool}, or a GPU sparse matrix for a GPU system (see device_sparse)
end

"""
    ADTypes.jacobian_sparsity(f!, y, x, detector::PrecomputedSparsityDetector)

ADTypes interface implementation for in-place functions `f!(y, x)`.
Returns the pre-computed sparsity pattern without tracing.
"""
function ADTypes.jacobian_sparsity(_f!, _y, _x, detector::PrecomputedSparsityDetector)
    return detector.pattern
end

"""
    ADTypes.jacobian_sparsity(f, x, detector::PrecomputedSparsityDetector)

ADTypes interface implementation for out-of-place functions `y = f(x)`.
Returns the pre-computed sparsity pattern without tracing.
"""
function ADTypes.jacobian_sparsity(_f, _x, detector::PrecomputedSparsityDetector)
    return detector.pattern
end
