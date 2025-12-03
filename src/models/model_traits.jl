import Base: length

# model_metadata interface is defined in PowerCore.model_metadata
# Individual models implement model_metadata(::Type{ModelType}) in their source files
# NOTE: We do NOT export model_metadata here to avoid conflict with PowerCore export

"""
    length(model::AbstractModel)

Returns the number of devices represented by this model by checking the length of 
any of its data fields. All fields in a model must have equal lengths.
"""
function length(model::AbstractModel)
    fields = fieldnames(typeof(model))
    isempty(fields) && return 0
    return length(getfield(model, first(fields)))
end
