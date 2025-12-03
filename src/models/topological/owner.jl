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
const OwnerVec = StructArray{Owner}

get_key_name(::Type{T}) where T<:Owner = [:i]

supports_format(::Type{Owner}, ::Type{PSSE}) = FormatSupport{PSSE}()

function parse_model(::Type{Owner}, raw, ::FormatSupport{PSSE})
    owners = raw.owners

    # Create StructArray directly from named tuple of vectors
    return StructArray{Owner}((
        i = owners.i,
        owname = owners.owname,
    ))
end

function model_metadata(::Type{Owner}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :Owner,
        vars = [],
        output_vars = [],
        residuals = [],
        layout = layout
    )
end

### === Begin Tests === ###

@testitem "Owner" begin
    using PowerFlowData
    using Powerful
    using Powerful.Models
    using Powerful.PowerCore
    using InlineStrings

    case = PowerFlowData.parse_network(
        joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    )

    owner = parse_model(Owner, case)
    @test owner isa OwnerVec
    @test all(owner.i .== [1, 2])
    @test all(owner.owname .== String15["OWNER1", "OWNER2"])
end