"""
Autodiff backend infrastructure for dual ForwardDiff/Enzyme support.

This module provides the shared infrastructure for both backends:
- Backend trait types (ForwardDiffBackend, EnzymeBackend)
- ADCacheCore (shared state across all backends)
- AbstractADCache (abstract type for polymorphism)
- create_di_backend (DifferentiationInterface backend factory)
"""

export ADBackend, ForwardDiffBackend, EnzymeBackend
export AbstractADCache, ADCacheCore
export create_adcache, create_di_backend, get_backend
export get_float_ws, get_formulation

using DifferentiationInterface: AutoSparse
using ADTypes: AutoForwardDiff, AutoEnzyme
using Enzyme: Enzyme
using SparseMatrixColorings: GreedyColoringAlgorithm

# ============================================================================
# Backend Trait Types
# ============================================================================

"""
    ADBackend

Abstract type for autodiff backend selection.

Subtypes:
- `ForwardDiffBackend`: Use ForwardDiff.jl with Dual number caching
- `EnzymeBackend`: Use Enzyme.jl with LLVM-level AD

# Example
```julia
prob = ACPowerFlowProblem(sys; backend=EnzymeBackend())
```
"""
abstract type ADBackend end

"""
    ForwardDiffBackend <: ADBackend

ForwardDiff.jl backend with Dual number workspace caching.

Zero allocation after warmup (first call per unique Dual type).
"""
struct ForwardDiffBackend <: ADBackend end

"""
    EnzymeBackend <: ADBackend

Enzyme.jl backend with LLVM-level automatic differentiation.

Zero allocation immediately (no warmup needed) - uses Float64 throughout.
"""
struct EnzymeBackend <: ADBackend end

# ============================================================================
# ADCacheCore - Shared State
# ============================================================================

"""
    ADCacheCore{F}

Shared state for all autodiff cache backends.

Contains the common fields needed by both ForwardDiffCache and EnzymeCache:
- System model reference
- Formulation trait
- Address manager
- Float64 workspace (used by all backends)

# Type Parameters
- `F`: FormulationTrait type (PowerInjection or CurrentInjection)

# Type Erasure Note
The `sys` field is typed as `Any` rather than `SystemModel` to prevent excessive
compile-time type specialization. Since `SystemModel` contains model-specific types
in its `models::NamedTuple`, specializing on sys would create separate compiled versions
of ADCacheCore for every unique model combination. The runtime cost of one dynamic
dispatch is negligible (~nanoseconds) compared to avoiding combinatorial compilation times.

All accessor functions use type barriers (e.g., `::SystemModel`) to restore type stability
at the call site.
"""
struct ADCacheCore{F<:FormulationTrait}
    sys::Any  # SystemModel reference (type-erased to break type chain - see note above)
    formulation::F
    address_manager::AddressManager
    float_ws::SystemBuffers{Float64, F}

    function ADCacheCore(sys, formulation::F, am::AddressManager, float_ws::SystemBuffers{Float64, F}) where F
        new{F}(sys, formulation, am, float_ws)
    end
end

# Core with a freshly built Float64 workspace (shared start of every create_adcache)
ADCacheCore(am::AddressManager, sys, formulation::FormulationTrait) =
    ADCacheCore(sys, formulation, am, SystemBuffers(Float64, am, sys, formulation))

# ============================================================================
# AbstractADCache - Interface
# ============================================================================

"""
    AbstractADCache{F}

Abstract type for autodiff caches.

All cache types must implement:
- Callable interface: `cache(g, y, status, t)` for power flow
- Callable interface: `cache(f, g, x, y, status, t)` for TDS

`t` is the time of the state and has no default: a model's residuals can depend on it
directly (not only through the state), so every caller states it.

Concrete subtypes:
- `ForwardDiffCache`: ForwardDiff with Dual workspace caching
- `EnzymeCache`: Enzyme with Float64-only path
"""
abstract type AbstractADCache{F<:FormulationTrait} end

"""
    (cache::AbstractADCache)(dae::DAE{Vector{Float64}}, status::Dict{Symbol, <:AbstractVector{Bool}}, t)

Evaluate residuals for Float64 DAE (zero-copy path) at time `t`, the time of the state in `dae`.
"""
function (cache::AbstractADCache)(dae::DAE{Vector{Float64}}, status::Dict{Symbol, <:AbstractVector{Bool}}, t::Real)
    evaluate_system_residuals!(cache.core.float_ws, dae, status, t)
    return nothing
end

# ============================================================================
# create_di_backend - DifferentiationInterface Backend Factory
# ============================================================================

"""
    create_di_backend(::ForwardDiffBackend, sparsity)

Create DifferentiationInterface backend for ForwardDiff with sparse Jacobian support.
"""
function create_di_backend(::ForwardDiffBackend, sparsity)
    AutoSparse(
        AutoForwardDiff();
        sparsity_detector = PrecomputedSparsityDetector(sparsity),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
end

"""
    create_di_backend(::EnzymeBackend, sparsity)

Create DifferentiationInterface backend for Enzyme with sparse Jacobian support.

Configuration:
- `set_runtime_activity(Forward)`: Enable runtime activity analysis for complex code paths
  where Enzyme cannot statically determine which memory locations are active.
- `function_annotation=Duplicated`: Tell Enzyme the closure captures mutable workspace
  state that affects derivatives. Without this, Enzyme defaults to Const.

See: https://enzyme.mit.edu/julia/stable/faq/#faq-runtime-activity
"""
function create_di_backend(::EnzymeBackend, sparsity)
    AutoSparse(
        AutoEnzyme(
            mode = Enzyme.set_runtime_activity(Enzyme.Forward),
            function_annotation = Enzyme.Duplicated
        );
        sparsity_detector = PrecomputedSparsityDetector(sparsity),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
end

# ============================================================================
# create_adcache - Factory Function (implemented by each backend)
# ============================================================================

"""
    create_adcache(am::AddressManager, sys, formulation, backend::ADBackend)

Factory function to create the appropriate ADCache for the given backend.

# Arguments
- `am`: AddressManager for variable/residual positions
- `sys`: SystemModel reference
- `formulation`: FormulationTrait (PowerInjection or CurrentInjection)
- `backend`: ADBackend (ForwardDiffBackend or EnzymeBackend)

# Returns
Backend-specific ADCache (ForwardDiffCache or EnzymeCache)
"""
function create_adcache end

# ============================================================================
# Accessor Functions
# ============================================================================
#
# All accessors use @inline for guaranteed inlining across module boundaries.
# Type barriers (::Type annotations) restore type stability for type-erased fields.
# ============================================================================

"""
    set_hold!(cache::AbstractADCache, y::AbstractVector)
    set_hold!(cache::AbstractADCache, x::AbstractVector, y::AbstractVector)

Hold the own algebraic variables of out-of-service instances at their values in `y`,
and their own zero-mass states at their values in `x`, for every workspace of the
cache (they share the held values).
"""
System.set_hold!(cache::AbstractADCache, y::AbstractVector) = System.set_hold!(cache.core.float_ws, y)
System.set_hold!(cache::AbstractADCache, x::AbstractVector, y::AbstractVector) =
    System.set_hold!(cache.core.float_ws, x, y)

"""
    get_float_ws(cache::AbstractADCache) -> SystemBuffers{Float64}

Get the Float64 SystemBuffers workspace from any ADCache type.
"""
@inline get_float_ws(cache::AbstractADCache) = cache.core.float_ws

"""
    get_system(cache::AbstractADCache) -> SystemModel

Get the SystemModel reference from any ADCache type.

Uses type barrier to restore type stability (sys is stored as Any).
"""
@inline get_system(cache::AbstractADCache)::SystemModel = cache.core.sys::SystemModel

"""
    get_formulation(cache::AbstractADCache) -> FormulationTrait

Get the formulation trait from any ADCache type.
"""
@inline get_formulation(cache::AbstractADCache) = cache.core.formulation

"""
    get_address_manager(cache::AbstractADCache) -> AddressManager

Get the AddressManager from any ADCache type.
"""
@inline get_address_manager(cache::AbstractADCache) = cache.core.address_manager

"""
    get_backend(cache::AbstractADCache) -> ADBackend

Get the backend type from any ADCache.

Implemented by each concrete cache type.
"""
function get_backend end
