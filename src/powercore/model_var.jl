
"""
    ModelVar(name::Symbol, var_type::AddressableType; kwargs...)

Constructor for internal variables
"""
function ModelVar(
    name::Symbol,
    var_type::VT;
    description::String = "",
    units::String = "",
    bounds::Union{Nothing, Tuple{Float64, Float64}} = nothing
) where {VT<:AddressableType}
    props = PropertyDict()
    
    # Add properties if provided
    !isempty(description) && (props[Description] = Description(description))
    !isempty(units) && (props[Units] = Units(units))
    !isnothing(bounds) && (props[Bounds] = Bounds(bounds[1], bounds[2]))
    
    ModelVar{Internal, VT, Nothing, Nothing, PropertyDict}(
        name,
        var_type,
        nothing,
        nothing,
        nothing,
        props
    )
end

"""
    ModelVar(name::Symbol, source_model::Type, source_name::Symbol; kwargs...)

Constructor for external variables
"""
function ModelVar(
    name::Symbol,
    var_type::VT,
    source_model::Symbol,
    source_name::Symbol,
    indexer::Symbol;
    description::String = "",
    units::String = ""
) where {VT<:AddressableType}
    props = PropertyDict()
    
    # Add properties if provided
    !isempty(description) && (props[Description] = Description(description))
    !isempty(units) && (props[Units] = Units(units))
    
    ModelVar{External, VT, Symbol, Symbol, PropertyDict}(
        name,
        var_type,
        source_model,
        source_name,
        indexer,
        props
    )
end

"""
Constructor for internal residual equations
"""
function ModelResidual(
    name::Symbol,
    eq_type::ET,
    access::RA;
    description::String = ""
) where {ET<:AddressableRes, RA<:ResAccess}
    ModelResidual{Internal, ET, RA, Nothing, Nothing, Nothing}(
        name,
        eq_type,
        access,
        nothing,
        nothing,
        nothing,
        description
    )
end

"""
Constructor for external residual equations
"""
function ModelResidual(
    name::Symbol,
    eq_type::ET,
    access::RA,
    source_model::SM,
    source_name::SR,
    indexer::IN = nothing;
    description::String = ""
) where {ET<:AddressableRes, RA<:ResAccess, SM, SR, IN}
    ModelResidual{External, ET, RA, SM, SR, IN}(
        name,
        eq_type,
        access,
        source_model,
        source_name,
        indexer,
        description
    )
end

# === Helper Functions === #
is_internal(::ModelVar{Internal, AT, SM, IN, P}) where {AT, SM, IN, P} = true
is_internal(::ModelVar{External, AT, SM, IN, P}) where {AT, SM, IN, P} = false

is_external(v::ModelVar) = !is_internal(v)

is_internal(::ModelResidual{Internal, ET, RA, SM, SR, IN}) where {ET, RA, SM, SR, IN} = true
is_internal(::ModelResidual{External, ET, RA, SM, SR, IN}) where {ET, RA, SM, SR, IN} = false

is_external(r::ModelResidual) = !is_internal(r)

# Property access helpers
function get_property(var::ModelVar, ::Type{T}) where {T<:VarProperty}
    get(var.properties, T, nothing)
end

function get_description(var::ModelVar)::Union{String, Nothing}
    desc = get_property(var, Description)
    isnothing(desc) ? nothing : desc.value
end

function get_units(var::ModelVar)::Union{String, Nothing}
    units = get_property(var, Units)
    isnothing(units) ? nothing : units.value
end

function get_bounds(var::ModelVar)::Union{Tuple{Float64, Float64}, Nothing}
    bounds = get_property(var, Bounds)
    isnothing(bounds) ? nothing : (bounds.min, bounds.max)
end

# Helper functions for metadata
function internal_vars(metadata::ModelMetadata)
    filter(is_internal, metadata.vars)
end

function external_vars(metadata::ModelMetadata)
    filter(is_external, metadata.vars)
end