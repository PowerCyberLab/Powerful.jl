"""
Metadata query predicates for filtering model metadata.

Use these predicates with Julia's `filter()`, `any()`, `count()`, etc.

# Quick Reference
```julia
filter(is_state_var, metadata.vars)
filter(is_algeb_var, metadata.vars)
filter(is_internal, metadata.vars)
filter(is_shared_res(:add), metadata.residuals)
filter(v -> is_internal(v) && is_state_var(v), metadata.vars)
```
"""

# ============================================================================
# Address Type Predicates
# ============================================================================

"""Check if item.address_type isa StateVar"""
is_state_var(item) = item.address_type isa StateVar

"""Check if item.address_type isa AlgebVar"""
is_algeb_var(item) = item.address_type isa AlgebVar

"""Check if item.address_type isa StateRes"""
is_state_res(item) = item.address_type isa StateRes

"""Check if item.address_type isa AlgebRes"""
is_algeb_res(item) = item.address_type isa AlgebRes

# `is_internal` / `is_external` are defined for ModelVar and ModelResidual in model_var.jl

# ============================================================================
# Access Pattern Predicates
# ============================================================================

"""Check if item.access isa PrivateRes"""
is_private_res(item) = item.access isa PrivateRes

"""
    is_shared_res(mode::Symbol=:any) -> Function

Predicate builder: returns predicate checking SharedRes with optional mode.

# Examples
```julia
filter(is_shared_res(), metadata.residuals)        # Any SharedRes
filter(is_shared_res(:add), metadata.residuals)    # Only SharedRes(:add)
filter(is_shared_res(:set), metadata.residuals)    # Only SharedRes(:set)
```
"""
function is_shared_res(mode::Symbol = :any)
    if mode == :any
        item -> item.access isa SharedRes
    else
        item -> item.access isa SharedRes && item.access.mode === mode
    end
end

# ============================================================================
# Convenience Queries
# ============================================================================

"""
    variable_limits(metadata::ModelMetadata) -> Vector{VariableLimit}

The bounds of the model's limited variables, as declared: its `VariableLimit`s
(`metadata.variable_limits`, the limits its automata enforce).
"""
variable_limits(metadata::ModelMetadata) = metadata.variable_limits

# ============================================================================
# Formulation-Specific Queries (not predicates)
# ============================================================================

"""
    get_metadata_for_formulation(metadata::ModelMetadata, formulation::FormulationTrait) -> ModelMetadata

Copy of `metadata` whose `vars` and `residuals` are the formulation-specific
overrides (from `vars_by_formulation` / `residuals_by_formulation`) when the
model declares them, and the base lists otherwise. All other fields are kept.

```julia
current_meta = get_metadata_for_formulation(metadata, CurrentInjection())
```
"""
function get_metadata_for_formulation(metadata::ModelMetadata, formulation::FormulationTrait)
    F = typeof(formulation)
    vars = get(metadata.vars_by_formulation, F, metadata.vars)
    residuals = get(metadata.residuals_by_formulation, F, metadata.residuals)
    # The positional constructor, the fields in their order (this runs for every model of
    # every problem built)
    return ModelMetadata(ntuple(Val(fieldcount(ModelMetadata))) do i
        name = fieldname(ModelMetadata, i)
        name === :vars ? vars : name === :residuals ? residuals : getfield(metadata, i)
    end...)
end

# Export all predicates
export is_state_var, is_algeb_var, is_state_res, is_algeb_res
export is_external, is_internal
export is_private_res, is_shared_res
export variable_limits
export get_metadata_for_formulation
