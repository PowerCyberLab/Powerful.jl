"""
Cached workspace for efficient autodiff evaluations.

This module provides workspace caching for automatic differentiation with
support for both ForwardDiff and Enzyme backends.

## Architecture

The autodiff cache system has three layers:

1. **ad_backend.jl**: Shared infrastructure
   - ADBackend trait types (ForwardDiffBackend, EnzymeBackend)
   - ADCacheCore (shared state)
   - AbstractADCache (interface)
   - create_di_backend (DI backend factory)

2. **forwarddiff_cache.jl**: ForwardDiff backend
   - ForwardDiffCache with Dual workspace caching
   - Zero allocation after warmup

3. **enzyme_cache.jl**: Enzyme backend
   - EnzymeCache with Float64-only path
   - Zero allocation immediately

## Usage

```julia
# Power flow with ForwardDiff (default)
prob = ACPowerFlowProblem(sys)
adcache = create_adcache(am, sys, formulation, ForwardDiffBackend())

# Power flow with Enzyme
adcache = create_adcache(am, sys, formulation, EnzymeBackend())

# Both backends implement the same callable interface
adcache(g, y, status, t)        # Power flow
adcache(f, g, x, y, status, t)  # TDS (t, the time of the state, has no default)
```

## Backend Comparison

| Backend | Mechanism | Zero Alloc After |
|---------|-----------|------------------|
| ForwardDiff | Dict-cached SystemBuffers{Dual} | Warmup |
| Enzyme | Single SystemBuffers{Float64} | Immediately |
"""

# Include shared infrastructure first
using DifferentiationInterface: DifferentiationInterface

include("ad_backend.jl")

# Include parallel backend implementations
include("forwarddiff_cache.jl")
include("enzyme_cache.jl")

# Parameters written into the workspaces by element type
include("parameters.jl")

"""
    prepare_sparse_jacobian(f!, y, x, adcache, contexts...; include_differential=false, strict=Val(true),
                            sparsity=nothing) -> (J, backend, prep)

Set up a sparse Jacobian of `f!(y, x, contexts...)` over the DAE variables of `adcache`
(`contexts`: DifferentiationInterface contexts, e.g. `Constant(t)` for the time of the state,
passed to `jacobian!` in the same types):
the sparsity pattern from the model structure (`jacobian_sparsity_pattern`, or the
given `sparsity`), a `Float64` Jacobian with that pattern, the DifferentiationInterface
backend for the cache's AD backend, and the prepared DI extras. For a GPU system the
pattern and Jacobian are device sparse matrices (`device_sparse`), so the Jacobian
is decompressed on the GPU.

`y` and `x` are passed to `prepare_jacobian` as given; their types must match the
arrays used later in `jacobian!` (e.g. a view stays a view). Pass `sparsity` to
share one pattern across several preparations (e.g. per-thread workspaces).
"""
function prepare_sparse_jacobian(f!, y, x, adcache::AbstractADCache, contexts...;
                                 include_differential::Bool = false,
                                 strict = Val(true),
                                 sparsity = nothing)
    pattern = isnothing(sparsity) ?
        jacobian_sparsity_pattern(get_address_manager(adcache), get_system(adcache),
                                  get_formulation(adcache); include_differential) :
        sparsity
    array_type = get_system(adcache).array_type
    backend = create_di_backend(get_backend(adcache), device_sparse(array_type, pattern))
    prep = DifferentiationInterface.prepare_jacobian(f!, y, backend, x, contexts...; strict)
    return (J = device_sparse(array_type, Float64.(pattern)), backend = backend, prep = prep)
end

# Re-export all symbols
export prepare_sparse_jacobian
