"""
Utility functions for working with model types and instances.

TODO: After the StructArray migration, these functions can be merged to another module.
"""

"""
    get_base_model_type(model_vec::AbstractModel) -> Type

Get the base model type from a model vector instance.


# Arguments
- `model_vec`: Model vector instance (e.g., PQVec, BusVec)

# Returns
- Base model type (e.g., PQ, Bus)

# Example
```julia
pq_vec = parse_model(PQ, case)  # Type: PQVec{Float64}
base_type = get_base_model_type(pq_vec)  # Returns: PQ

# Can then get metadata
metadata = model_metadata(base_type)
```

# Implementation Notes
For StructArray-wrapped models, extracts the base type from the element type.
"""
function get_base_model_type(model_vec::StructArray{T}) where T <: AbstractModel
    # For StructArray, the eltype IS the base type
    # Return the unparameterized base type (e.g., Bus instead of Bus{Float64})
    return T.name.wrapper
end
