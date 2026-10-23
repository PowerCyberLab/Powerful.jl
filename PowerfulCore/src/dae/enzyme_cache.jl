"""
Enzyme-based autodiff cache with Float64-only path.

This module provides the Enzyme backend implementation:
- EnzymeCache: Uses Float64 throughout - Enzyme handles derivatives via LLVM

Key difference from ForwardDiff:
- No Dual number types
- No workspace type proliferation
- Zero allocation immediately (no warmup needed)

Usage with DifferentiationInterface:
- Use `AutoEnzyme(; mode=Enzyme.Forward, function_annotation=Enzyme.Duplicated)`
- This tells DI to properly shadow the functor's mutable state
"""

export EnzymeCache

# ============================================================================
# EnzymeCache
# ============================================================================

"""
    EnzymeCache{F} <: AbstractADCache{F}

Enzyme-based autodiff cache with Float64-only path.

Always uses Float64 - Enzyme handles derivatives via LLVM transformation
and shadow arrays, not Dual numbers.

**Zero allocation**: Immediately (no warmup needed)

# Type Parameters
- `F`: FormulationTrait type (PowerInjection or CurrentInjection)

# Fields
- `core::ADCacheCore{F}`: Shared state (sys, formulation, address_manager, float_ws)
- `empty_f::Vector{Float64}`: Pre-allocated empty vector for power flow (avoids allocation in differentiated code)
- `empty_x::Vector{Float64}`: Pre-allocated empty vector for power flow (avoids allocation in differentiated code)

# Example
```julia
# Create via factory function (recommended)
adcache = create_adcache(am, sys, formulation, EnzymeBackend())

# Use as callable
adcache(g, y)        # Power flow
adcache(f, g, x, y)  # TDS
```

# Usage with DifferentiationInterface
```julia
using DifferentiationInterface: AutoSparse
using ADTypes: AutoEnzyme
using Enzyme

# IMPORTANT: Use function_annotation=Duplicated to properly shadow the functor
backend = AutoSparse(
    AutoEnzyme(; mode=Enzyme.Forward, function_annotation=Enzyme.Duplicated);
    sparsity_detector = PrecomputedSparsityDetector(sparsity),
    coloring_algorithm = GreedyColoringAlgorithm()
)

prep = DifferentiationInterface.prepare_jacobian(adcache, g, backend, y)
DifferentiationInterface.jacobian!(adcache, g, J, prep, backend, y)
```

# Implementation Notes
Enzyme works at the LLVM level:
- Forward mode: Uses (x, dx) pairs where both are Float64
- The workspace stays Float64 - no SystemBuffers{Dual} needed
- DifferentiationInterface handles shadow array management
- **Important**: Empty vectors for f/x are pre-allocated because Enzyme cannot
  differentiate through array allocations inside the target function
"""
struct EnzymeCache{F<:FormulationTrait} <: AbstractADCache{F}
    core::ADCacheCore{F}
    empty_f::Vector{Float64}
    empty_x::Vector{Float64}
end

# ============================================================================
# Factory Function
# ============================================================================

"""
    create_adcache(am, sys, formulation, ::EnzymeBackend)

Create EnzymeCache for Enzyme backend.
"""
function create_adcache(am::AddressManager, sys, formulation::F,
                        ::EnzymeBackend) where {F<:FormulationTrait}
    core = ADCacheCore(am, sys, formulation)
    # Pre-allocate empty vectors for power flow (Enzyme can't allocate in differentiated code)
    empty_f = Vector{Float64}(undef, 0)
    empty_x = Vector{Float64}(undef, 0)
    EnzymeCache{F}(core, empty_f, empty_x)
end

"""
    get_backend(::EnzymeCache)

Returns EnzymeBackend() for EnzymeCache.
"""
get_backend(::EnzymeCache) = EnzymeBackend()

# Enzyme differentiates the Float64 workspace, which holds the memory of the automata:
# there is no other workspace to copy it into (see the ForwardDiffCache method)
sync_automaton_memory!(::EnzymeCache) = nothing

# ============================================================================
# Callable Interfaces
# ============================================================================

"""
    (cache::EnzymeCache)(g::AbstractVector{Float64}, y::AbstractVector{Float64}, status, t)

Evaluate algebraic residuals (power flow interface) at time `t`, the time of the state `y`.

Enzyme always uses Float64 - no Dual dispatch needed.
"""
function (cache::EnzymeCache)(g::AbstractVector{Float64}, y::AbstractVector{Float64},
                              status::Dict{Symbol, <:AbstractVector{Bool}},
                              t::Float64)
    # Use pre-allocated empty vectors (zero allocation)
    cache(cache.empty_f, g, cache.empty_x, y, status, t)
    return nothing
end

"""
    (cache::EnzymeCache)(f, g, x, y, status, t; recorder=nothing)

Evaluate DAE residuals (both differential and algebraic) at time `t`, the time of the state
`(x, y)`.

Enzyme always uses Float64 - derivatives computed via LLVM shadow arrays.
"""
function (cache::EnzymeCache)(
    f::AbstractVector{Float64},
    g::AbstractVector{Float64},
    x::AbstractVector{Float64},
    y::AbstractVector{Float64},
    status::Dict{Symbol, <:AbstractVector{Bool}},
    t::Float64;
    recorder = nothing
)
    # Always Float64 path - Enzyme handles derivatives externally
    evaluate_system_residuals!(cache.core.float_ws, f, g, x, y, status, t; recorder=recorder)
    return nothing
end

# ============================================================================
# TDS Dual Buffers - Not Applicable for Enzyme
# ============================================================================

"""
    get_tds_dual_buffers(cache::EnzymeCache, ::Type{Float64}, n_diff::Int, n_alg::Int)

EnzymeCache does not use Dual number buffers. For Float64, TDS buffers are
managed by the TDSAdapter, not the ADCache.

This method exists for API compatibility but should not be called in normal usage.
If called, it indicates a code path error where Enzyme is being used like ForwardDiff.
"""
function get_tds_dual_buffers(::EnzymeCache, ::Type{Float64}, n_diff::Int, n_alg::Int)
    error("""
        get_tds_dual_buffers should not be called on EnzymeCache.

        EnzymeCache uses Float64 throughout - TDS buffers are managed by TDSAdapter.
        This error indicates a code path that incorrectly treats Enzyme like ForwardDiff.

        Check that the calling code dispatches on backend type correctly.
        """)
end

"""
    get_tds_dual_buffers(cache::EnzymeCache, ::Type{T}, n_diff::Int, n_alg::Int) where T

EnzymeCache does not support Dual number types.

Enzyme performs differentiation at the LLVM level using shadow arrays,
not via Dual number propagation. If this method is called with a Dual type,
it indicates a backend configuration error.
"""
function get_tds_dual_buffers(::EnzymeCache, ::Type{T}, n_diff::Int, n_alg::Int) where T
    error("""
        EnzymeCache does not support Dual number types.

        Got type: $T

        Enzyme performs differentiation at the LLVM level using shadow arrays,
        not via Dual number propagation like ForwardDiff.

        If you need Dual number support, use ForwardDiffBackend() instead.
        """)
end
