export AreaInterchange, AreaInterchangeVec

"""
    $TYPEDEF

Area interchange data records.

# Fields

$TYPEDFIELDS
"""
@kwdef struct AreaInterchange{Tv} <: AbstractModel
    """
    Area number (1 through the maximum number of areas at the current size level)
    """
    i::Int16

    """
    Area slack bus number for interchange control, given by its bus number or
    extended bus name enclosed in single quotes. The bus must be a generator
    (type two) bus in the specified area. Any area containing a system swing bus
    (type three) must have either that swing bus or a bus number of zero
    specified for its area slack bus number. `isw` = 0 by default.
    """
    isw::Int32

    """
    Desired net interchange leaving the area (export) in MW.
    """
    pdes::Tv

    """
    Interchange tolerance bandwidth; entered in MW.
    `ptol` = 10.0 by default.
    """
    ptol::Tv

    """
    Alphanumeric identifier assigned to area I.
    The name may contain up to twelve characters.
    `arname` is set to twelve blanks by default.
    """
    arname::String15
end

@register_model AreaInterchange

# Type aliases for StructArray-based vectorized types
"""
    AreaInterchangeVec{T}

Alias for `StructArray{AreaInterchange{T}}`, the struct-of-arrays container holding all
[`AreaInterchange`](@ref) instances of a system.
"""
const AreaInterchangeVec{T} = StructArray{AreaInterchange{T}}

get_key_name(::Type{T}) where T<:AreaInterchange = [:i]

supports_format(::Type{AreaInterchange}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{AreaInterchange}, ::Type{PSSE}) = (:area_interchanges,)

function parse_model(::Type{AreaInterchange}, raw, ::FormatSupport{PSSE})
    area_interchanges = raw.area_interchanges

    # Create StructArray directly from named tuple of vectors
    return StructArray{AreaInterchange{Float64}}((
        i = area_interchanges.i,
        isw = area_interchanges.isw,
        pdes = area_interchanges.pdes,
        ptol = area_interchanges.ptol,
        arname = area_interchanges.arname,
    ))
end

function model_metadata(::Type{AreaInterchange})
    ModelMetadata(
        name = :AreaInterchange,
        vars = [],
        observable_vars = [],
        residuals = [],

    )
end

# No-op for PowerInjection formulation (AreaInterchange has no residuals)
@inline function calc_residuals!(
    _ws::ResidualBuffers{T},
    _model::StructArray{<:AreaInterchange},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T)
) where {T}
    nothing
end

# No-op for CurrentInjection formulation (AreaInterchange has no residuals)
@inline function calc_residuals!(
    _ws::ResidualBuffers{T},
    _model::StructArray{<:AreaInterchange},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T)
) where {T}
    nothing
end

# No-op for DCPowerFlow formulation (AreaInterchange has no residuals)
@inline function calc_residuals!(
    _ws::ResidualBuffers{T},
    _model::StructArray{<:AreaInterchange},
    _external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T)
) where {T}
    nothing
end

