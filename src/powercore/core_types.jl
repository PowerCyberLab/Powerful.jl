abstract type AbstractModel end

### ============================= ###
abstract type VarOwnership end
struct OwnVar <: VarOwnership end
struct ForeignVar{Ti<:Integer,T<:AbstractVector{Ti}} <: VarOwnership
    model_name::Symbol
    var_name::Symbol
    indexer::T
end

# === Variable Traits === #
abstract type VarProperty end

struct Description <: VarProperty
    value::String
end

struct Units <: VarProperty
    value::String
end

struct Bounds <: VarProperty
    min::Float64
    max::Float64
end

# Type-safe property collection
const PropertyDict = Dict{Type{<:VarProperty},VarProperty}

abstract type VarTrait end
struct Internal <: VarTrait end
struct External <: VarTrait end

abstract type AddressableType end
abstract type AddressableVar <: AddressableType end
abstract type AddressableRes <: AddressableType end

struct AlgebVar <: AddressableVar end
struct StateVar <: AddressableVar end
struct ObservedVar <: AddressableVar end
struct AlgebRes <: AddressableRes end
struct StateRes <: AddressableRes end

"""
    ModelVar{T<:VarTrait, AT, SM, IN, Props} represents a variable in a power system model

Type parameters:
- T: Variable trait (Internal/External)
- AT: Addressable type (e.g., Algeb, State)
- SM: Source model type or Nothing
- IN: Indexer type or Nothing
- Props: Type of properties dictionary
"""
struct ModelVar{
    T<:VarTrait,
    AT<:AddressableVar,
    SM<:Union{Symbol,Nothing},
    IN<:Union{Symbol,Nothing},
    Props<:PropertyDict
}
    name::Symbol
    address_type::AT
    source_model::SM
    source_name::Union{Symbol,Nothing}
    indexer::IN
    properties::Props
end

abstract type ResAccess end
struct SharedRes <: ResAccess end
struct PrivateRes <: ResAccess end

struct ModelResidual{
    T<:VarTrait,
    AT<:AddressableRes,
    RA<:ResAccess,
    SM<:Union{Symbol,Nothing},
    SR<:Union{Symbol,Nothing},
    IN<:Union{Symbol,Nothing}
}
    name::Symbol
    address_type::AT
    access::RA
    source_model::SM
    source_name::SR
    indexer::IN
    description::String
end

export AddressableType, AddressableVar, AddressableRes
export AlgebVar, StateVar, ObservedVar
export AlgebRes, StateRes
export VarTrait, Internal, External
export VarOwnership, OwnVar, ForeignVar
export VarProperty, Description, Units, Bounds, PropertyDict
export ModelVar

export ResTrait, PartialRes, FullRes, SharedRes
export ModelResidual
### ============================= ###

### ========== Address Manager ========== ###

"""
Abstract type for model-specific addresses
Generated based on model requirements
"""
abstract type ModelAddresses end


"""
Tracks address assignment status for a specific model type
"""
mutable struct ModelAddressAssignment
    model_name::Symbol
    # maps addressable type => variable name => range
    by_type::Dict{AddressableType,Dict{Symbol,Vector{UInt}}}
    is_complete::Bool  # true if all required variables are allocated
end

"""
Track retrieved addresses for a specific model instance
"""
struct ModelRetrievedAddresses
    model_name::Symbol
    by_type::Dict{AddressableType,Dict{Symbol,Vector{UInt}}}
end

"""
Enhanced address manager with explicit address assignment tracking
"""
Base.@kwdef struct AddressManager
    addresses::Dict{Tuple{AddressableType,Symbol,Symbol},Vector{UInt}} = Dict()
    next_idx::Dict{AddressableType,Int} = Dict()

    # Address assignment tracking
    assignments::Dict{Symbol,ModelAddressAssignment} = Dict()
    retrieved::Dict{Symbol,ModelRetrievedAddresses} = Dict()
end

"""
Base types for layout strategy and metadata
"""
abstract type LayoutStrategy end
struct ContiguousVariables <: LayoutStrategy end
struct ContiguousInstances <: LayoutStrategy end
struct CustomAddressing <: LayoutStrategy end


export AddressManager
export LayoutStrategy, ContiguousVariables, ContiguousInstances, CustomAddressing
export ModelAddressAssignment

### ============= End Address Manager ============= ###

### ============= Model Metadata and Traits ============= ###

"""
    ModelMetadata

"""
Base.@kwdef struct ModelMetadata
    name::Symbol
    vars::Vector{ModelVar}
    output_vars::Vector{Symbol}
    residuals::Vector{ModelResidual}
    layout::LayoutStrategy = ContiguousVariables()
end

# Bus/Node model types
abstract type NodeTrait end
struct ACNode <: NodeTrait end
struct DCNode <: NodeTrait end

# Load model types
abstract type AbstractLoad{Tv} <: AbstractModel end

export ModelMetadata
export NodeTrait, ACNode, DCNode
export AbstractLoad

### ============= End Model Metadata and Traits ============= ###

### ============= System Types ============= ###

"""
    SystemModel

A system model is a collection of models that are connected to each other.

# Fields

$TYPEDFIELDS
"""
mutable struct SystemModel{T<:NamedTuple}
    "Address manager object"
    address_manager::AddressManager

    "A named tuple of models"
    models::T

    "Array type for computational backend (Vector for CPU, CuVector for GPU)"
    array_type::Type{<:AbstractVector}

    # TODO: Boundary buses/dc nodes?

    # Inner constructor to validate component types
    function SystemModel(am::AddressManager, components::NamedTuple;
                         auto_setup::Bool=true,
                         array_type::Type{<:AbstractVector}=Vector)
        # Ensure all elements are vectors/arrays of AbstractModel subtypes
        # Check eltype since StructArray wraps the base type
        all(V -> eltype(V) <: AbstractModel, values(components)) ||
            error("All components must be AbstractModel subtypes")

        sys = new{typeof(components)}(am, components, array_type)

        # Automatically setup unless disabled
        if auto_setup
            setup!(sys)
        end

        return sys
    end
end

export SystemModel

### ========================================= ###


### ============= Data Format Types ============= ###

abstract type DataFormat end
struct PSSE <: DataFormat end
struct MATPOWER <: DataFormat end

export DataFormat, PSSE

### ========================================= ###
