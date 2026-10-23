"""
Generic gather/scatter kernels for CPU and GPU.

These kernels are used throughout the codebase for efficient data movement:
- residual_workspace.jl: Gather variables, scatter residuals
- parameter_system.jl: Gather/scatter optimization parameters

All kernels are zero-allocation and use efficient SIMD/GPU parallelism.
"""

# ============================================================================
# Gather Operations: dest[i] = src[indices[i]]
# ============================================================================

"""
    gather_cpu!(dest, src, indices)

CPU implementation of gather using SIMD vectorization.
dest[i] = src[indices[i]] for all i.
"""
@inline function gather_cpu!(dest::AbstractVector, src::AbstractVector, indices::AbstractVector)
    @inbounds @simd for i in eachindex(indices)
        dest[i] = src[indices[i]]
    end
    nothing
end

"""
    gather!(dest, src, indices)

Generic gather operation with automatic CPU/GPU dispatch.
dest[i] = src[indices[i]] for all i.

Zero allocation, uses SIMD on CPU.

Handles SubArray (from view) and mixed element types (e.g., Vector{Float64} and Vector{Number}).
"""
@inline gather!(dest::AbstractVector, src::AbstractVector, indices::AbstractVector) =
    gather_cpu!(dest, src, indices)

# ============================================================================
# Scatter Operations: dest[indices[i]] = src[i]
# ============================================================================

"""
    scatter_cpu!(dest, src, indices, active, mode::Symbol)

CPU implementation of scatter with optional masking.
- mode = :private: Direct write (dest[indices[i]] = src[i]) for PrivateRes
- mode = :add: Accumulate (dest[indices[i]] += src[i]) for SharedRes(:add)
- mode = :set: Direct write (dest[indices[i]] = src[i]) for SharedRes(:set)

Only writes if active[i] is true.
"""
@inline function scatter_cpu!(
    dest::AbstractVector,
    src::AbstractVector,
    indices::AbstractVector,
    active::AbstractVector{Bool},
    mode::Symbol
)
    eltype(dest) == eltype(src) ||
        error("scatter: element types differ (dest $(eltype(dest)), src $(eltype(src)))")

    if mode === :private || mode === :set
        # Direct write for PrivateRes and setter SharedRes(:set)
        # SIMD is safe here - no loop-carried dependencies
        @inbounds @simd for i in eachindex(indices)
            if active[i]
                dest[indices[i]] = src[i]
            end
        end
    else  # :add (default SharedRes)
        # Accumulate for adder SharedRes(:add)
        # NOTE: @simd removed because with duplicate indices, multiple iterations
        # may write to the same address, creating loop-carried dependencies
        # that @simd cannot handle safely
        @inbounds for i in eachindex(indices)
            if active[i]
                dest[indices[i]] += src[i]
            end
        end
    end
    nothing
end

"""
    PartialAddresses(addrs)

The addresses of a residual that only some instances are delivered through (see
[`injection_routes`](@ref)): an instance with the address [`UNROUTED`](@ref) is not scattered.
Residuals every instance is delivered through keep a plain address vector and the loops above.
"""
struct PartialAddresses <: AbstractVector{UInt}
    addrs::Vector{UInt}
end
Base.size(a::PartialAddresses) = size(a.addrs)
Base.getindex(a::PartialAddresses, i::Int) = a.addrs[i]

@inline function scatter_cpu!(
    dest::AbstractVector,
    src::AbstractVector,
    indices::PartialAddresses,
    active::AbstractVector{Bool},
    mode::Symbol
)
    eltype(dest) == eltype(src) ||
        error("scatter: element types differ (dest $(eltype(dest)), src $(eltype(src)))")
    addrs = indices.addrs
    adding = mode === :add
    @inbounds for i in eachindex(addrs)
        (active[i] && addrs[i] != UNROUTED) || continue
        if adding
            dest[addrs[i]] += src[i]
        else
            dest[addrs[i]] = src[i]
        end
    end
    nothing
end

"""
    scatter!(dest, src, indices, active, mode::Symbol)

Generic scatter operation with automatic CPU/GPU dispatch.
- mode = :private: Direct write
- mode = :add: Accumulate
- mode = :set: Direct write (setter mode)

Only writes if active[i] is true. Zero allocation.
"""
@inline scatter!(
    dest::AbstractVector,
    src::AbstractVector,
    indices::AbstractVector,
    active::AbstractVector{Bool},
    mode::Symbol
) = scatter_cpu!(dest, src, indices, active, mode)

# ============================================================================
# Simplified scatter for unmasked operations (for parameter system)
# ============================================================================

"""
    scatter_cpu_unmasked!(dest, src, indices)

CPU scatter without active mask checking.
"""
@inline function scatter_cpu_unmasked!(
    dest::AbstractVector,
    src::AbstractVector,
    indices::AbstractVector
)
    @inbounds @simd for i in eachindex(indices)
        dest[indices[i]] = src[i]
    end
    nothing
end

"""
    scatter_unmasked!(dest, src, indices)

Generic unmasked scatter with automatic CPU/GPU dispatch.
Used when all elements are active (e.g., parameter optimization).

Handles SubArray (from view) and mixed element types (e.g., Vector{Number} and Vector{Float64}).
"""
@inline scatter_unmasked!(dest::AbstractVector, src::AbstractVector, indices::AbstractVector) =
    scatter_cpu_unmasked!(dest, src, indices)
