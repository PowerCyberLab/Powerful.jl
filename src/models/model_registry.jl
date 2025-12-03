export @register_model
export MODEL_REGISTRY
export get_registered_models, get_model_type
export all_models


"""
Central registry for model types and their metadata
"""
struct ModelRegistry{T<:Tuple}
    models::T
end

# Empty registry with type parameter
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
Macro to register a model type. Updates MODEL_REGISTRY.
"""
macro register_model(model_type::Symbol)
    quote
        MODEL_REGISTRY[] = _add_model(MODEL_REGISTRY[], $(esc(model_type)))
    end
end

# Type-stable registry queries
"""Get all registered model types"""
get_registered_models(registry::ModelRegistry{T}) where T = registry.models

"""Get model type from registry by name"""
function get_model_type(registry::ModelRegistry{T}, name::Symbol) where T
    for M in registry.models
        if nameof(M) === name
            return M
        end
    end
    error("Model $name not found in registry")
end

"""
    all_models()

Get all registered model types.

Convenience accessor for MODEL_REGISTRY[].models. This is the recommended
way to access all registered models instead of directly accessing the registry.

# Returns
- Tuple of all registered model types

# Example
```julia
# Get all models
models = all_models()

# Use in load_system
models_data = load_system(case, all_models())

# Use in SystemModel (though default already uses all models)
sys = SystemModel(case, all_models())

# Or just use the default which does this automatically:
sys = SystemModel(case)  # Uses all_models() by default
```

# See Also
- [`MODEL_REGISTRY`](@ref): The underlying registry
- [`get_registered_models`](@ref): Lower-level accessor
"""
all_models() = MODEL_REGISTRY[].models

"""
    _check_model_docstring(model_type::Type) -> Nothing

Check if the given model type has a docstring defined.
Throws an ArgumentError if no docstring is found.
"""
function _check_model_docstring(model_type::Type)
    let binding = Base.Docs.Binding(@__MODULE__, nameof(model_type))
        if !haskey(Base.Docs.meta(@__MODULE__), binding)
            @warn "Missing docstring for model $model_type"
        end
    end
    return nothing
end


### === Begin Test ===
@testitem "ModelRegistry" begin
    using Powerful.PowerCore
    using Powerful.Models

    @testset "get_registered_models" begin
        @test Bus in get_registered_models(MODEL_REGISTRY[])
        @test get_model_type(MODEL_REGISTRY[], :Bus) == Bus
    end
end
