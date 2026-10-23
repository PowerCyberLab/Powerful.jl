# Model Registry - Central registration system for model types

using .Foundation: AbstractModel

export @register_model, @init_model_registry
export MODEL_REGISTRY
export get_registered_models, get_model_type
export all_models
export register_pending_models!


"""
Central registry for model types and their metadata
"""
struct ModelRegistry{T<:Tuple}
    models::T
end

"""
    MODEL_REGISTRY

Global `Ref` holding the registry of all model types registered with
[`@register_model`](@ref) across the loaded sub-packages. Sub-packages add their
models at load time with [`register_pending_models!`](@ref); query it with
[`all_models`](@ref) and [`get_registered_models`](@ref).
"""
const MODEL_REGISTRY = Ref{ModelRegistry}(ModelRegistry(()))

"""
Type-stable helper to add a model type to registry tuple
"""
function _add_model(registry::ModelRegistry{T}, ::Type{M}) where {T, M<:AbstractModel}
    # Create new tuple type with added model
    new_models = (registry.models..., M)
    ModelRegistry{typeof(new_models)}(new_models)
end

"""
    register_pending_models!(models)

Register a collection of model types to the shared MODEL_REGISTRY.

Called from sub-package `__init__()` functions to register models that were
collected during precompilation. This is the standard pattern for cross-package
model registration in Julia.

# Example
```julia
# In sub-package module file:
const _PENDING_MODELS = Type[]

# ... includes with @register_model calls ...

function __init__()
    PowerfulCore.register_pending_models!(_PENDING_MODELS)
end
```
"""
function register_pending_models!(models)
    for M in models
        MODEL_REGISTRY[] = _add_model(MODEL_REGISTRY[], M)
    end
end

"""
    @init_model_registry

Initialize model registration infrastructure for a sub-package.

Call this macro once at the top of a sub-package module (after `using PowerfulCore`)
to set up the `_PENDING_MODELS` vector and `_EXPORTED_MODEL_NAMES` for tracking.
Then use `@register_model` in individual model files, and call
`register_pending_models!(_PENDING_MODELS)` in `__init__()`.

After all models are registered (after all includes), copy `_EXPORTED_MODEL_NAMES`
to `EXPORTED_MODELS` for use by the umbrella package.

# Example
```julia
module PowerfulPF
using PowerfulCore
@init_model_registry  # Creates _PENDING_MODELS and _EXPORTED_MODEL_NAMES

include("models/bus.jl")  # Contains @register_model Bus
include("models/pq.jl")   # Contains @register_model PQ

# After all includes, create EXPORTED_MODELS for umbrella package
const EXPORTED_MODELS = copy(_EXPORTED_MODEL_NAMES)

function __init__()
    register_pending_models!(_PENDING_MODELS)
end
end
```
"""
macro init_model_registry()
    esc(quote
        const _PENDING_MODELS = Type[]
        const _EXPORTED_MODEL_NAMES = Symbol[]
    end)
end

"""
Macro to register a model type.

When called from a sub-package that defines `_PENDING_MODELS::Vector{Type}`,
the model type is added to that vector during precompilation. The sub-package's
`__init__()` should then call `register_pending_models!(_PENDING_MODELS)` to
transfer them to the shared registry at runtime.

Also tracks the model name in `_EXPORTED_MODEL_NAMES` for automatic export
by the umbrella package (Powerful.jl).

When called from PowerfulCore itself, directly updates MODEL_REGISTRY.
"""
macro register_model(model_type::Symbol)
    mod = __module__
    quote
        if isdefined($mod, :_PENDING_MODELS)
            # Sub-package: collect for later registration in __init__()
            push!($mod._PENDING_MODELS, $(esc(model_type)))
            # Also track name for export by umbrella package
            if isdefined($mod, :_EXPORTED_MODEL_NAMES)
                push!($mod._EXPORTED_MODEL_NAMES, $(QuoteNode(model_type)))
            end
        else
            # Direct registration (PowerfulCore or interactive use)
            MODEL_REGISTRY[] = _add_model(
                MODEL_REGISTRY[], $(esc(model_type))
            )
        end
    end
end

# Type-stable registry queries
"""
    get_registered_models(registry::ModelRegistry) -> Tuple{Vararg{Type}}

Return tuple of all model types in the registry.

For the global registry, prefer `all_models()` which handles lazy loading.
"""
get_registered_models(registry::ModelRegistry{T}) where T = registry.models

"""
    get_model_type(name::Symbol) -> Type{<:AbstractModel}

Get model type from global registry by name.

# Example
```julia
model_type = get_model_type(:Bus)    # → Bus
model_type = get_model_type(:GENCLS) # → GENCLS
```
"""
function get_model_type(name::Symbol)
    for M in MODEL_REGISTRY[].models
        if nameof(M) === name
            return M
        end
    end
    error("Model $name not found in registry")
end

"""
    all_models() -> Tuple{Vararg{Type}}

Return all registered model types. Preferred over direct registry access.
"""
all_models() = MODEL_REGISTRY[].models
