
"""
    ModelVar(name::Symbol, var_type::AddressableType; kwargs...)

Constructor for internal variables
"""
function ModelVar(
    name::Symbol,
    var_type::VT;
    description::String = "",
    units::String = "",
    mass_coefficient::Union{Nothing, Symbol, MassCoefficient} = nothing
) where {VT<:AddressableType}
    props = PropertyDict()

    # Add properties if provided
    !isempty(description) && (props[Description] = Description(description))
    !isempty(units) && (props[Units] = Units(units))
    if !isnothing(mass_coefficient)
        # Accept both Symbol (backward compatible) and MassCoefficient
        props[MassCoefficient] = mass_coefficient isa MassCoefficient ? mass_coefficient : MassCoefficient(mass_coefficient)
    end

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
    ref_key::Symbol;
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
        ref_key,
        props
    )
end

"""
Constructor for external variables with composite ref_key
"""
function ModelVar(
    name::Symbol,
    var_type::VT,
    source_model::Symbol,
    source_name::Symbol,
    ref_key::Vector{Symbol};
    description::String = "",
    units::String = ""
) where {VT<:AddressableType}
    props = PropertyDict()

    !isempty(description) && (props[Description] = Description(description))
    !isempty(units) && (props[Units] = Units(units))

    ModelVar{External, VT, Symbol, Vector{Symbol}, PropertyDict}(
        name,
        var_type,
        source_model,
        source_name,
        ref_key,
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
    ref_key::RK = nothing;
    description::String = ""
) where {ET<:AddressableRes, RA<:ResAccess, SM, SR, RK}
    ModelResidual{External, ET, RA, SM, SR, RK}(
        name,
        eq_type,
        access,
        source_model,
        source_name,
        ref_key,
        description
    )
end

# === Helper Functions === #
"""
    is_internal(item::Union{ModelVar, ModelResidual}) -> Bool

Whether a model variable or residual is owned by its model ([`Internal`](@ref)), as
opposed to referencing another model's ([`External`](@ref)).
"""
is_internal(::ModelVar{Internal, AT, SM, IN, P}) where {AT, SM, IN, P} = true
is_internal(::ModelVar{External, AT, SM, IN, P}) where {AT, SM, IN, P} = false

"""
    is_external(item::Union{ModelVar, ModelResidual}) -> Bool

Whether a model variable or residual references another model ([`External`](@ref));
the negation of [`is_internal`](@ref).
"""
is_external(v::ModelVar) = !is_internal(v)

is_internal(::ModelResidual{Internal, ET, RA, SM, SR, IN}) where {ET, RA, SM, SR, IN} = true
is_internal(::ModelResidual{External, ET, RA, SM, SR, IN}) where {ET, RA, SM, SR, IN} = false

is_external(r::ModelResidual) = !is_internal(r)

# Property access helpers

"""
    get_property(var::ModelVar, ::Type{T}) where T<:VarProperty -> Union{T, Nothing}

Retrieve a property of type `T` from the variable's properties dict.

Returns `nothing` if the property is not set.
"""
function get_property(var::ModelVar, ::Type{T}) where {T<:VarProperty}
    get(var.properties, T, nothing)
end

"""
    get_description(var::ModelVar) -> Union{String, Nothing}

Get the human-readable description string for a variable.

Returns `nothing` if no `Description` property was set.
"""
function get_description(var::ModelVar)::Union{String, Nothing}
    desc = get_property(var, Description)
    isnothing(desc) ? nothing : desc.value
end

"""
    get_units(var::ModelVar) -> Union{String, Nothing}

Get the units string for a variable (e.g., "rad/s", "pu").

Returns `nothing` if no `Units` property was set.
"""
function get_units(var::ModelVar)::Union{String, Nothing}
    units = get_property(var, Units)
    isnothing(units) ? nothing : units.value
end

"""
    get_mass_coeff(var::ModelVar) -> Union{MassCoefficient, Nothing}

Get mass coefficient for a differential variable (M in M*dx/dt = f).

Returns `nothing` if not a differential variable or no coefficient set.
"""
function get_mass_coeff(var::ModelVar)::Union{MassCoefficient, Nothing}
    get_property(var, MassCoefficient)
end

# Helper functions for metadata

"""
    external_vars(metadata::ModelMetadata) -> Vector{ModelVar}

Filter metadata to return only externally-sourced variables.

These are variables gathered from other models (e.g., bus voltage).
"""
function external_vars(metadata::ModelMetadata)
    filter(is_external, metadata.vars)
end