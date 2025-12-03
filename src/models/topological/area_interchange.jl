export AreaInterchange, AreaInterchangeVec
export AreaInterchangeNum, AreaInterchangeNumVec

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
const AreaInterchangeVec{T} = StructArray{AreaInterchange{T}}

get_key_name(::Type{T}) where T<:AreaInterchange = [:i]

supports_format(::Type{AreaInterchange}, ::Type{PSSE}) = FormatSupport{PSSE}()

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

function model_metadata(::Type{AreaInterchange}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :AreaInterchange,
        vars = [],
        output_vars = [],
        residuals = [],
        layout = layout
    )
end

## === Begin Tests === ###

@testitem "AreaInterchange" begin
    using PowerFlowData
    using Powerful
    using Powerful.Models
    using Powerful.PowerCore
    using InlineStrings

    case = PowerFlowData.parse_network(
        joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    )

    area_interchange = parse_model(AreaInterchange, case)
    @test area_interchange isa AreaInterchangeVec
    @test all(area_interchange.i .== [1, 2])
    @test all(area_interchange.arname .== String15["AREA1", "AREA2"])
end
