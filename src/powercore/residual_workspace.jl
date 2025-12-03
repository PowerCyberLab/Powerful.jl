"""
Residual calculation workspace with zero-allocation design.

This module provides pre-allocated workspaces for efficient residual evaluation
in power system DAE models. Uses NamedTuples for type-stable, dictionary-like
access with zero runtime overhead.
"""

export ResidualWorkspace
export gather_variables!, zero_residuals!, scatter_residuals!

# ============================================================================
# Helper Functions for Constructor
# ============================================================================

"""
    _build_variable_storage(::Type{T}, metadata, am, array_type) where T

Build variable address and value storage for a model.

Returns two vectors of pairs: (var_addr_pairs, var_val_pairs)
- Addresses are transferred to GPU if array_type <: CuVector
- Values are allocated with the specified array type
"""
function _build_variable_storage(
    ::Type{T},
    metadata::ModelMetadata,
    am::AddressManager,
    array_type::Type{<:AbstractVector}
) where T<:Number
    AT = array_type{T}
    model_name = metadata.name

    var_addr_pairs = Pair{Symbol, Any}[]
    var_val_pairs = Pair{Symbol, AT}[]

    for var in metadata.vars
        var_name = var.name

        # Get addresses from address manager
        addrs_cpu = if is_internal(var)
            am.assignments[model_name].by_type[var.address_type][var_name]
        else
            am.retrieved[model_name].by_type[var.address_type][var_name]
        end

        # Transfer to GPU if needed (one-time cost)
        addrs = if AT <: CUDA.CuVector
            CUDA.CuArray(addrs_cpu)
        else
            addrs_cpu
        end

        # Allocate value storage
        push!(var_addr_pairs, var_name => addrs)
        push!(var_val_pairs, var_name => AT(zeros(T, length(addrs))))
    end

    return var_addr_pairs, var_val_pairs
end

"""
    _build_residual_storage(::Type{T}, metadata, am, array_type) where T

Build residual address and value storage for a model.

Returns two vectors of pairs: (res_addr_pairs, res_val_pairs)
- Addresses are transferred to GPU if array_type <: CuVector
- Values are allocated with the specified array type
"""
function _build_residual_storage(
    ::Type{T},
    metadata::ModelMetadata,
    am::AddressManager,
    array_type::Type{<:AbstractVector}
) where T<:Number
    AT = array_type{T}
    model_name = metadata.name

    res_addr_pairs = Pair{Symbol, Any}[]
    res_val_pairs = Pair{Symbol, AT}[]

    for res in metadata.residuals
        res_name = res.name

        # Get addresses from address manager
        addrs_cpu = if is_internal(res)
            am.assignments[model_name].by_type[res.address_type][res_name]
        else
            am.retrieved[model_name].by_type[res.address_type][res_name]
        end

        # Transfer to GPU if needed (one-time cost)
        addrs = if AT <: CUDA.CuVector
            CUDA.CuArray(addrs_cpu)
        else
            addrs_cpu
        end

        # Allocate value storage
        push!(res_addr_pairs, res_name => addrs)
        push!(res_val_pairs, res_name => AT(zeros(T, length(addrs))))
    end

    return res_addr_pairs, res_val_pairs
end

# ============================================================================
# Core Workspace Structure
# ============================================================================

"""
    ResidualWorkspace{T, AT, VarAddrs, VarVals, ResAddrs, ResVals}

Pre-allocated workspace for a single model type's residual calculation.

Uses NamedTuples for type-stable access to variable and residual arrays.
All memory is allocated during construction - zero allocation during calculation.

# Type Parameters
- `T`: Float type (Float64, Float32, etc.)
- `AT`: Array type for data storage (Vector{T} for CPU, CuVector{T} for GPU)
- `VarAddrs`: NamedTuple type for variable addresses
- `VarVals`: NamedTuple type for variable values
- `ResAddrs`: NamedTuple type for residual addresses
- `ResVals`: NamedTuple type for residual values

# Fields
- `var_addrs`: NamedTuple mapping variable names to address vectors (always CPU)
- `var_vals`: NamedTuple mapping variable names to value vectors (can be CPU/GPU)
- `res_addrs`: NamedTuple mapping residual names to address vectors (always CPU)
- `res_vals`: NamedTuple mapping residual names to value vectors (can be CPU/GPU)
- `active_mask`: Array indicating which instances are in-service (matches AT)
- `n_instances`: Number of model instances

# Example
For a PQ load model on CPU:
```julia
ws = ResidualWorkspace{Float64, Vector{Float64}, ...}(...)
v_values = ws.var_vals.v  # Vector{Float64}
```

For GPU:
```julia
ws = ResidualWorkspace{Float64, CuVector{Float64}, ...}(...)
v_values = ws.var_vals.v  # CuVector{Float64}
```
"""
struct ResidualWorkspace{T<:Number,
                         AT<:AbstractVector{T},
                         MaskT<:AbstractVector{Bool},
                         VarAddrs<:NamedTuple,
                         VarVals<:NamedTuple,
                         ResAddrs<:NamedTuple,
                         ResVals<:NamedTuple}
    var_addrs::VarAddrs  # Always CPU (just UInt64 indices)
    var_vals::VarVals    # Can be CPU or GPU
    res_addrs::ResAddrs  # Always CPU
    res_vals::ResVals    # Can be CPU or GPU
    active_mask::MaskT   # Matches backend (Vector{Bool} or CuVector{Bool})
    n_instances::Int
end

# ============================================================================
# Constructor from ModelMetadata
# ============================================================================

"""
    ResidualWorkspace(::Type{T}, metadata::ModelMetadata,
                           am::AddressManager, model_vec::AbstractModel,
                           array_type::Type{<:AbstractVector}) where T

Construct a workspace for a model with pre-allocated arrays.

The array type is explicitly specified (Vector for CPU, CuVector for GPU).
Typically called from SystemResidualWorkspace which gets array_type from sys.array_type.

# Arguments
- `T`: Float type for calculations
- `metadata`: ModelMetadata describing variables and residuals
- `am`: AddressManager containing address assignments
- `model_vec`: Model vector instance (e.g., PQVec)
- `array_type`: Array type constructor (Vector, CuVector, etc.)

# Returns
- `ResidualWorkspace`: Fully allocated workspace ready for calculations

# Example
```julia
# Typically used internally by SystemResidualWorkspace
ws = ResidualWorkspace(Float64, metadata, am, model_vec, Vector)
ws_gpu = ResidualWorkspace(Float64, metadata, am, model_vec_gpu, CuVector)
```
"""
function ResidualWorkspace(
    ::Type{T},
    metadata::ModelMetadata,
    am::AddressManager,
    model_vec::StructArray{M},
    array_type::Type{<:AbstractVector}
) where {T<:Number, M<:AbstractModel}

    n_instances = length(model_vec)
    AT = array_type{T}

    # Build variable and residual storage using helper functions
    var_addr_pairs, var_val_pairs = _build_variable_storage(T, metadata, am, array_type)
    res_addr_pairs, res_val_pairs = _build_residual_storage(T, metadata, am, array_type)

    # Create active mask (all instances active by default)
    active_mask_base = hasproperty(model_vec, :status) ? model_vec.status : trues(n_instances)

    # Transfer active mask to GPU if needed
    active_mask = AT <: CUDA.CuVector ? CUDA.CuArray(active_mask_base) : active_mask_base

    # Construct NamedTuples
    var_addrs_nt = NamedTuple(var_addr_pairs)
    var_vals_nt = NamedTuple(var_val_pairs)
    res_addrs_nt = NamedTuple(res_addr_pairs)
    res_vals_nt = NamedTuple(res_val_pairs)

    # Return workspace with full type parameters
    return ResidualWorkspace{T, AT, typeof(active_mask),
                             typeof(var_addrs_nt), typeof(var_vals_nt),
                             typeof(res_addrs_nt), typeof(res_vals_nt)}(
        var_addrs_nt,
        var_vals_nt,
        res_addrs_nt,
        res_vals_nt,
        active_mask,
        n_instances
    )
end

# ============================================================================
# Generic Operations (Zero Allocation)
# ============================================================================

"""
    gather_variables!(ws::ResidualWorkspace{T}, dae::DAE{T}) where T

Gather variable values from global DAE vector into local workspace buffers.

Reads from `dae.y[addr]` and writes to `ws.var_vals[name][i]` for all variables
defined in the workspace. Zero allocation - operates on pre-allocated arrays.

# Arguments
- `ws`: Workspace with pre-allocated variable buffers
- `dae`: DAE system containing global variable vector

# Notes
- This function is `@inline` for compiler optimization
- Uses `@inbounds` and `@simd` for performance
- The compiler unrolls NamedTuple iteration at compile time
"""
@inline function gather_variables!(
    ws::ResidualWorkspace{T, Vector{T}, <:AbstractVector{Bool}},
    dae::DAE{<:Vector{T}}
) where T
    # CPU version: Direct indexed access
    for (var_name, addrs) in pairs(ws.var_addrs)
        vals = getfield(ws.var_vals, var_name)

        # Vectorized gather operation
        @inbounds @simd for i in eachindex(addrs)
            vals[i] = dae.y[addrs[i]]
        end
    end

    nothing
end

"""
GPU gather kernel - zero allocation variable gathering.

Each thread copies one element from the global DAE vector to the local workspace
using the address mapping.
"""
function gather_kernel!(dest, src, indices, n)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= n
        @inbounds dest[i] = src[indices[i]]
    end
    return nothing
end

@inline function gather_variables!(
    ws::ResidualWorkspace{T, <:CUDA.CuVector{T}, <:CUDA.CuVector{Bool}},
    dae::DAE{<:CUDA.CuVector{T}}
) where T
    # GPU version: Use custom kernel to avoid broadcast allocation
    threads = 256

    for (var_name, addrs) in pairs(ws.var_addrs)
        vals = getfield(ws.var_vals, var_name)
        n = length(addrs)
        blocks = cld(n, threads)

        # Launch gather kernel (async, no sync needed - kernels execute in order)
        CUDA.@cuda threads=threads blocks=blocks gather_kernel!(vals, dae.y, addrs, n)
    end

    # No synchronize - let kernels execute async on default stream
    nothing
end

"""
    zero_residuals!(ws::ResidualWorkspace{T}) where T

Zero out local residual buffers before calculation.

Sets all values in `ws.res_vals` to zero. Should be called before
`calc_residuals!` to ensure clean slate.

# Arguments
- `ws`: Workspace with pre-allocated residual buffers
"""
@inline function zero_residuals!(ws::ResidualWorkspace{T}) where T
    # Iterate over all residuals in the NamedTuple
    for (res_name, vals) in pairs(ws.res_vals)
        fill!(vals, zero(T))
    end

    nothing
end

"""
    scatter_residuals!(ws::ResidualWorkspace{T}, dae::DAE{T},
                      metadata::ModelMetadata) where T

Scatter (accumulate) residual values from local workspace into global DAE vector.

Reads from `ws.res_vals[name][i]` and writes to `dae.g[addr]` for all residuals.
Respects `active_mask` to skip inactive instances. Zero allocation.

# Residual Access Patterns
- `SharedRes`: Accumulates into shared address space (e.g., bus power balance)
  Uses `+=` to add contributions from multiple models
- `PrivateRes`: Direct write to private address space (e.g., branch flow limits)
  Uses `=` for direct assignment

# Arguments
- `ws`: Workspace with computed residual values
- `dae`: DAE system with global residual vector
- `metadata`: Model metadata specifying residual access patterns
"""
@inline function scatter_residuals!(
    ws::ResidualWorkspace{T, Vector{T}, <:AbstractVector{Bool}},
    dae::DAE{<:Vector{T}},
    metadata::ModelMetadata
) where T
    # CPU version: Direct indexed access
    for res in metadata.residuals
        res_name = res.name
        addrs = getfield(ws.res_addrs, res_name)
        vals = getfield(ws.res_vals, res_name)

        if res.access isa SharedRes
            # SharedRes: accumulate contributions (multiple models write to same address)
            @inbounds @simd for i in eachindex(addrs)
                if ws.active_mask[i]
                    dae.g[addrs[i]] += vals[i]
                end
            end
        else  # PrivateRes
            # PrivateRes: direct write (only this model writes to these addresses)
            @inbounds @simd for i in eachindex(addrs)
                if ws.active_mask[i]
                    dae.g[addrs[i]] = vals[i]
                end
            end
        end
    end

    nothing
end

"""
GPU scatter kernel for SharedRes (with atomic add)
"""
function scatter_shared_kernel!(dest, vals, addrs, active, n)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= n && active[i]
        @inbounds CUDA.@atomic dest[addrs[i]] += vals[i]
    end
    return nothing
end

"""
GPU scatter kernel for PrivateRes (direct write, no atomics needed)
"""
function scatter_private_kernel!(dest, vals, addrs, active, n)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= n && active[i]
        @inbounds dest[addrs[i]] = vals[i]
    end
    return nothing
end

@inline function scatter_residuals!(
    ws::ResidualWorkspace{T, <:CUDA.CuVector{T}, <:CUDA.CuVector{Bool}},
    dae::DAE{<:CUDA.CuVector{T}},
    metadata::ModelMetadata
) where T
    # GPU version: Direct scatter using CUDA kernels (zero allocation)
    for res in metadata.residuals
        res_name = res.name
        addrs = getfield(ws.res_addrs, res_name)  # Already CuArray!
        vals = getfield(ws.res_vals, res_name)
        n = length(addrs)

        # Launch configuration
        threads = 256
        blocks = cld(n, threads)

        if res.access isa SharedRes
            # SharedRes: Use atomic add for thread-safe accumulation
            CUDA.@cuda threads=threads blocks=blocks scatter_shared_kernel!(
                dae.g, vals, addrs, ws.active_mask, n
            )
        else  # PrivateRes
            # PrivateRes: direct write (no atomics needed)
            CUDA.@cuda threads=threads blocks=blocks scatter_private_kernel!(
                dae.g, vals, addrs, ws.active_mask, n
            )
        end
    end

    # No synchronize - let kernels execute async on default stream
    nothing
end
