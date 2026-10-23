export Owner, OwnerVec
"""
    $TYPEDEF

PSS/E allows users to identify ownership of facilities, equipment, and loads
through owner assignments. Major network elements can have up to four different
owners, enabling interpretation and reporting of results based on ownership.

# Fields

$TYPEDFIELDS
"""
@kwdef struct Owner <: AbstractModel
    """
    Owner number (1 through the maximum number of owners at the current size level)
    """
    i::Int16

    """
    Alphanumeric identifier assigned to owner i. The name may contain up to
    twelve characters and must be enclosed in single quotes. Can contain any
    combination of blanks, uppercase letters, numbers and special characters.
    Set to twelve blanks by default.
    """
    owname::String15
end

@register_model Owner

# Type alias for StructArray-based vectorized type
"""
    OwnerVec

Alias for `StructArray{Owner}`, the struct-of-arrays container holding all
[`Owner`](@ref) instances of a system.
"""
const OwnerVec = StructArray{Owner}

get_key_name(::Type{T}) where T<:Owner = [:i]

supports_format(::Type{Owner}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Owner}, ::Type{PSSE}) = (:owners,)

function parse_model(::Type{Owner}, raw, ::FormatSupport{PSSE})
    owners = raw.owners

    # Create StructArray directly from named tuple of vectors
    return StructArray{Owner}((
        i = owners.i,
        owname = owners.owname,
    ))
end

function model_metadata(::Type{Owner})
    ModelMetadata(
        name = :Owner,
        vars = ModelVar[],
        observable_vars = Symbol[],
        residuals = ModelResidual[],
    )
end

# No-op for PowerInjection formulation (Owner has no residuals)
@inline function calc_residuals!(
    _ws::ResidualBuffers{T},
    _model::StructArray{Owner},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T)
) where {T}
    nothing
end

# No-op for CurrentInjection formulation (Owner has no residuals)
@inline function calc_residuals!(
    _ws::ResidualBuffers{T},
    _model::StructArray{Owner},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T)
) where {T}
    nothing
end

# No-op for DCPowerFlow formulation (Owner has no residuals)
@inline function calc_residuals!(
    _ws::ResidualBuffers{T},
    _model::StructArray{Owner},
    _external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T)
) where {T}
    nothing
end

