export FormatSupport
export supports_format, parse_model, is_format_supported, validate_raw

"""
Trait system for data format conversions
"""

# Trait types
struct FormatSupport{F} end
struct NoFormatSupport end

"""
    supports_format(::Type{M}, format::Type{F}) where {M,F}

Returns the format support trait for a given model type and format.
Default implementation returns NoFormatSupport.
"""
supports_format(::Type, ::Type) = NoFormatSupport()

"""
    parse_model(::Type{M}, data, trait)

Convert raw data to model instance based on trait dispatch.
Must be implemented for specific model types and formats.
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
Load components for any model type that supports the format
"""
function parse_model(::Type{M}, raw_data, ::Type{F}=PSSE) where {M,F}
    trait = supports_format(M, F)
    return parse_model(M, raw_data, trait)
end

"""
Check if a model type supports conversion from a format
"""
is_format_supported(::Type{M}, ::Type{F}) where {M,F} = !(supports_format(M, F) isa NoFormatSupport)

"""
Helper to validate raw data before conversion
"""
function validate_raw(::Type{M}, data, format::Type{F}) where {M,F}
    if !is_format_supported(M, format)
        throw(ArgumentError("Model type $M does not support conversion from format $F"))
    end
    # TODO: Additional validation logic...
end
