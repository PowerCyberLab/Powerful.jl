export Bus, BusVec, BusNum, BusNumVec

"""
    $TYPEDEF

Network bus data records.

# Fields

$TYPEDFIELDS
"""
Base.@kwdef mutable struct Bus{Tv} <: AbstractBus{Tv}
    "Bus number (1 to 999997)."
    i::Int32

    """
    Alphanumeric identifier assigned to bus 'I'.

    The name may be up to twelve characters and must be enclosed in single
    quotes. NAME may contain any combination of blanks, uppercase letters,
    numbers and special characters, but the first character must not be a minus
    sign.
    """
    name::String15
    
    "Bus base voltage in kV."
    basekv::Tv

    """
    Bus type code:
    - 1 - load bus or other bus without any generator boundary condition. 
    - 2 - generator or plant bus either regulating voltage or with a fixed reactive power (Mvar). A generator that reaches its reactive power limit will no longer control voltage but rather hold reactive power at its limit. 
    - 3 - swing bus or slack bus. It has no power or reactive limits and regulates voltage at a fixed reference angle. 
    - 4 - disconnected or isolated bus.
    """
    ide::Int8

    "Area number. 1 through the maximum number of areas at the current size level."
    area::Int16

    "Zone number. 1 through the maximum number of zones at the current size level."
    zone::Int16

    "Owner number. 1 through the maximum number of owners at the current size level."
    owner::Int16

    "Bus voltage magnitude in pu."
    vm::Tv

    "Bus voltage angle in degrees."
    va::Tv

    "Normal voltage magnitude high limit in pu."
    nvhi::Tv

    "Normal voltage magnitude low limit in pu."
    nvlo::Tv

    "Emergency voltage magnitude high limit in pu."
    evhi::Tv

    "Emergency voltage magnitude low limit in pu."
    evlo::Tv
end

@register_model Bus

# Type aliases for StructArray-based vectorized types
const BusVec{T} = StructArray{Bus{T}}
const BusNum{T} = Bus{T}  # No longer needed, but kept for compatibility
const BusNumVec{T} = StructArray{Bus{T}}  # Same as BusVec

get_key_name(::Type{Bus}) = :i

# Define Bus model variables
const BUS1PH_VARS = [
    ModelVar(
        :theta, 
        AlgebVar(),
        description="Bus voltage angle",
        units="rad",
        bounds=(-Inf, Inf)
    ),

    ModelVar(
        :v, 
        AlgebVar(),
        description="Bus voltage magnitude",
        units="pu",
        bounds=(0.0, 2.0)
    )
]

const BUS1PH_OUTPUT_VARS = [:theta, :v]

# Define residual variables similar to model variables
const BUS1PH_RESIDUALS = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes();
        description="Active power balance equation",
    ),
    ModelResidual(
        :q, 
        AlgebRes(),
        SharedRes();
        description="Reactive power balance equation",
    )
]

function model_metadata(::Type{Bus}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :Bus,
        vars = BUS1PH_VARS,
        output_vars = BUS1PH_OUTPUT_VARS,
        residuals = BUS1PH_RESIDUALS,
        layout = layout
    )
end


### === Format Support === ###

supports_format(::Type{Bus}, ::Type{PSSE}) = FormatSupport{PSSE}()

function parse_model(::Type{Bus}, raw, ::FormatSupport{PSSE})
    buses = raw.buses

    # Create StructArray directly from named tuple of vectors
    return StructArray{Bus{Float64}}((
        i = buses.i,
        name = buses.name,
        basekv = buses.basekv,
        ide = buses.ide,
        area = buses.area,
        zone = buses.zone,
        owner = buses.owner,
        vm = buses.vm,
        va = buses.va,
        nvhi = buses.nvhi,
        nvlo = buses.nvlo,
        evhi = buses.evhi,
        evlo = buses.evlo,
    ))
end


@testitem "Bus" begin
    using PowerFlowData
    using Powerful
    using Powerful.Models
    using Powerful.PowerCore
    using StructArrays

    case = PowerFlowData.parse_network(
        joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    )

    bus = parse_model(Bus, case)
    # @test bus isa BusVec

    # struct_array = to_struct_array(bus)
    # @test struct_array isa StructArray
end
