export Zone, ZoneVec

"""
    $TYPEDEF

All buses and loads can be assigned to reside in a zone of the network. Zones
enable users to develop reports and check results on a zonal basis, allowing
highly specific reporting and interpretation of analytical results. Zone
assignments are made by specifying zone numbers in the data records for buses
and loads.

# Fields

$TYPEDFIELDS
"""
@kwdef struct Zone <: AbstractModel
    "Zone number (1 through the maximum number of zones at the current size level)"
    i::Int16

    """
    Alphanumeric identifier assigned to zone i. The name may contain up to
    twelve characters and must be enclosed in single quotes. Can contain any
    combination of blanks, uppercase letters, numbers and special characters.
    Set to twelve blanks by default.
    """
    zoname::String15
end

@register_model Zone

# Type alias for StructArray-based vectorized type
const ZoneVec = StructArray{Zone}

get_key_name(::Type{T}) where T<:Zone = [:i]

supports_format(::Type{Zone}, ::Type{PSSE}) = FormatSupport{PSSE}()

function parse_model(::Type{Zone}, raw, ::FormatSupport{PSSE})
    zones = raw.zones

    # Create StructArray directly from named tuple of vectors
    return StructArray{Zone}((
        i = zones.i,
        zoname = zones.zoname,
    ))
end

function model_metadata(::Type{Zone}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :Zone,
        vars = [],
        output_vars = [],
        residuals = [],
        layout = layout
    )
end

### === Begin Tests === ###

@testitem "Zone" begin
    using PowerFlowData
    using Powerful
    using Powerful.Models
    using Powerful.PowerCore
    using InlineStrings

    case = PowerFlowData.parse_network(
        joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    )

    zone = parse_model(Zone, case)
    @test zone isa ZoneVec
    @test all(zone.i .== [1, 2])
    @test all(zone.zoname .== String15["ZONE1", "ZONE2"])
end
