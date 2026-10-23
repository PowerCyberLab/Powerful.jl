"""
    with_parameters(dyn_models, parameters, values) -> NamedTuple

A copy of the dynamic models `dyn_models` with the parameters `parameters`,
`(model, field, instance)` tuples, set to `values`: the point of an optimizer's
iteration, to simulate or differentiate again.
"""
function with_parameters(dyn_models::NamedTuple, parameters::AbstractVector, values::AbstractVector{<:Real})
    length(parameters) == length(values) ||
        throw(ArgumentError("$(length(parameters)) parameters but $(length(values)) values"))
    models = map(copy, dyn_models)
    for ((model, field, instance), value) in zip(parameters, values)
        getproperty(models[model], field)[instance] = value
    end
    return models
end

export with_parameters
