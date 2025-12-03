"""
Interface for model metadata.

All models must implement model_metadata(::Type{ModelType}) to describe
their variables, residuals, and memory layout strategy.
"""

"""
    model_metadata(::Type{M}; layout::LayoutStrategy=ContiguousVariables()) where {M<:AbstractModel}

Get metadata for a model type. All models must implement this method.

The metadata describes the model's variables (both internal and external),
residual equations, output variables, and memory layout strategy for the
address manager.

# Arguments
- Model type (e.g., `PQ`, `Bus`, `Branch`)
- `layout`: Optional layout strategy for address assignment (default: ContiguousVariables())

# Returns
- `ModelMetadata` object describing the model's structure

# Example
```julia
# Get metadata for PQ load model
metadata = model_metadata(PQ)

# With custom layout
metadata = model_metadata(Bus; layout=ContiguousInstances())

# Access metadata properties
metadata.vars      # Vector of ModelVar
metadata.residuals # Vector of ModelResidual
metadata.layout    # Layout strategy
```

# Implementation
Each model must define this method in its source file:
```julia
function model_metadata(::Type{MyModel}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :MyModel,
        vars = [...],
        output_vars = [...],
        residuals = [...],
        layout = layout
    )
end
```
"""
function model_metadata(model::Type; layout::LayoutStrategy=ContiguousVariables())
    _ = layout  # Unused in stub, but part of interface
    error("No `model_metadata` function defined for model $model. " *
          "Each model must implement model_metadata(::Type{$model}) to describe its structure.")
end

"""
    model_metadata(model_vec::AbstractModel; kwargs...)

Convenience method to get metadata from a model vector instance.

# Arguments
- `model_vec`: Model vector instance (e.g., PQVec, BusVec)
- `kwargs...`: Additional keyword arguments passed to model_metadata(Type)

# Returns
- `ModelMetadata` for the base model type

# Example
```julia
pq_vec = parse_model(PQ, case)

# Old way
metadata = model_metadata(PQ)

# New way (more convenient)
metadata = model_metadata(pq_vec)

# With custom layout
metadata = model_metadata(pq_vec; layout=ContiguousInstances())
```
"""
function model_metadata(model_vec::AbstractModel; kwargs...)
    base_type = get_base_model_type(model_vec)
    return model_metadata(base_type; kwargs...)
end

"""
    calc_residuals!(ws::ResidualWorkspace, model_vec::AbstractModel)

Calculate residuals for a model.

This is a generic interface that must be implemented by each model type.
Model-specific implementations compute the residual equations and store
the results in the workspace's residual buffers.

# Arguments
- `ws`: Pre-allocated workspace with gathered variables
- `model_vec`: Model vector instance (e.g., PQVec, BusVec)

# Implementation
Each model must implement this method:
```julia
function calc_residuals!(ws::ResidualWorkspace, model_vec::MyModelVec)
    # Access variables from ws.var_vals
    # Compute residuals
    # Store in ws.res_vals
end
```

# See Also
- [`ResidualWorkspace`](@ref): Workspace structure
- [`evaluate_system_residuals!`](@ref): System-level evaluation
"""
# Fallback for StructArray-wrapped models
function calc_residuals!(ws::ResidualWorkspace{T, AT}, model_vec::StructArray{M}) where {T, AT, M<:AbstractModel}
    # Fallback: do nothing (appropriate for models like Bus that don't compute residuals)
    # Models that need to compute residuals should implement their own method
    nothing
end

export model_metadata, calc_residuals!
