using PowerFlowData
using Powerful.PowerCore: parse_model, PSSE, supports_format
using Powerful.PowerCore: AbstractModel

"""
$(SIGNATURES)

Load raw data for specified model types
"""
function load_system(raw_data::PowerFlowData.Network, model_types::Vector{T}) where T
    return NamedTuple(nameof(T) => parse_model(T, raw_data) for T in model_types)
end

function load_system(raw_data::PowerFlowData.Network, model_types::Tuple{Vararg{Type}}) 
    return NamedTuple(nameof(T) => parse_model(T, raw_data) for T in model_types)
end

export load_system, parse_model