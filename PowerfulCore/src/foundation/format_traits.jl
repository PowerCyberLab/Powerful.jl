export FormatSupport
export supports_format, parse_model, is_format_supported, validate_raw
export format_records, declares_format_records, record_counts

import PowerFlowData
export load_models

"""
Trait system for data format conversions
"""

"""
    FormatSupport{F}()

Trait value returned by [`supports_format`](@ref) when a model can be parsed from data
format `F` (a [`DataFormat`](@ref)); [`parse_model`](@ref) dispatches on it.
"""
struct FormatSupport{F} end

"""
    NoFormatSupport()

Trait value returned by [`supports_format`](@ref) when a model does not support a data
format. [`parse_model`](@ref) throws an `ArgumentError` for it.
"""
struct NoFormatSupport end

"""
    supports_format(::Type{M}, format::Type{F}) where {M,F}

Returns the format support trait for a given model type and format.
Default implementation returns NoFormatSupport.
"""
supports_format(::Type, ::Type) = NoFormatSupport()

"""
    parse_model(::Type{M}, data, trait) -> StructArray{M}

Convert raw parsed data to a model StructArray based on format trait dispatch.

# Arguments
- `M`: Target model type (e.g., `Bus`, `Branch`)
- `data`: Raw data from file parser
- `trait`: Format trait from `supports_format(M, F)`

Implement for specific model types and formats.
"""
function parse_model end

# Default fallback with informative error
function parse_model(::Type{M}, data, ::NoFormatSupport) where M
    throw(ArgumentError(
        "No conversion implemented for model type $M. " *
        "Please implement `supports_format` and `parse_model` for this type."
    ))
end

"""
    parse_model(::Type{M}, raw_data, ::Type{F}=PSSE) -> StructArray{M}

Convert raw data to model StructArray using automatic format dispatch.

Looks up the format trait via `supports_format(M, F)` and dispatches.
"""
function parse_model(::Type{M}, raw_data, ::Type{F}=PSSE) where {M,F}
    trait = supports_format(M, F)
    return parse_model(M, raw_data, trait)
end

"""
    is_format_supported(::Type{M}, ::Type{F}) -> Bool

Check if model type `M` supports conversion from format `F` (e.g., PSSE).
"""
is_format_supported(::Type{M}, ::Type{F}) where {M,F} = !(supports_format(M, F) isa NoFormatSupport)

"""
    format_records(::Type{M}, ::Type{F}) -> Tuple{Vararg{Symbol}}

The record categories of data format `F` that model `M` reads: for PSS/E RAW, the record
tables of a `PowerFlowData.Network` (e.g. `(:buses,)` for Bus). Every model that supports a
format (`supports_format`) declares them; the RAW reader warns of the categories with records
that no model it loads reads, and refuses a supported model that declares none.
"""
function format_records end

"""
    declares_format_records(::Type{M}, ::Type{F}) -> Bool

Whether model `M` declares the records of format `F` it reads (a method of `format_records`).
"""
declares_format_records(::Type{M}, ::Type{F}) where {M,F} = hasmethod(format_records, Tuple{Type{M}, Type{F}})

"""
    record_counts(raw_data::PowerFlowData.Network) -> Vector{Pair{Symbol, Int}}

The record categories of a parsed PSS/E RAW case (its `PowerFlowData.Records` tables), each
with its number of records.
"""
function record_counts(raw_data::PowerFlowData.Network)
    counts = Pair{Symbol, Int}[]
    for name in fieldnames(typeof(raw_data))
        records = getfield(raw_data, name)
        records isa PowerFlowData.Records && push!(counts, name => length(records))
    end
    return counts
end

"""
Helper to validate raw data before conversion
"""
function validate_raw(::Type{M}, data, format::Type{F}) where {M,F}
    if !is_format_supported(M, format)
        throw(ArgumentError("Model type $M does not support conversion from format $F"))
    end
    # TODO: Additional validation logic...
end

"""
    load_models(raw_data, model_types)

Load models from parsed network data for specified model types.

This is a stub function that must be implemented by downstream packages
(e.g., PowerfulPF for PSS/E RAW format).

# Arguments
- `raw_data`: Parsed network data (e.g., `PowerFlowData.Network`)
- `model_types`: Tuple or vector of model types to load

# Returns
- `NamedTuple` of model StructArrays, keyed by model name
"""
function load_models end

# Stub for SystemModel convenience constructor - implemented in system_builder.jl
function _build_system_from_namedtuple end
