# Base Quantity Trait System
# ============================================================================
# This module provides a trait-based system for resolving base quantities
# (voltage, power, frequency, etc.) for different model types.
#
# Models declare their relationship to base quantities via traits:
# - OwnsBase: Model stores this base quantity in a field
# - InheritsBase: Model receives this base quantity from connected models (via ExtParam)
# - NoBase: Model doesn't interact with this base quantity
#
# This eliminates hard-coded field name checks and enables extensibility for
# future stability simulation models with diverse base quantity requirements.

# ============================================================================
# Trait Type Hierarchy
# ============================================================================

"""
    BaseQuantityRole

Abstract type for base quantity trait roles.

Concrete subtypes:
- `OwnsBase`: Model stores this base quantity in a field
- `InheritsBase`: Model receives from ExtParam/terminals
- `NoBase`: Model doesn't use this base quantity
"""
abstract type BaseQuantityRole end

"""
    OwnsBase <: BaseQuantityRole

Trait indicating a model stores this base quantity directly in a field.

# Example
```julia
# Bus owns its base voltage (basekv field)
base_quantity_role(::Type{Bus}, ::Type{BaseVoltage}) = OwnsBase()
```
"""
struct OwnsBase <: BaseQuantityRole end

"""
    InheritsBase <: BaseQuantityRole

Trait indicating a model receives this base quantity from connected models
via the ExtParam system.

# Example
```julia
# PQ load inherits base voltage from connected bus
base_quantity_role(::Type{PQ}, ::Type{BaseVoltage}) = InheritsBase()
```
"""
struct InheritsBase <: BaseQuantityRole end

"""
    NoBase <: BaseQuantityRole

Trait indicating a model doesn't interact with this base quantity: it has no value on that
base. This is the default for all model-base combinations unless explicitly overridden, and
`get_base_quantity_vector` refuses it, naming the model and the quantity to convert.
"""
struct NoBase <: BaseQuantityRole end

# ============================================================================
# Core Trait API
# ============================================================================

"""
    base_quantity_role(::Type{<:AbstractModel}, ::Type{<:ReferenceType}) -> BaseQuantityRole

Determine a model's relationship to a base quantity.

Returns one of: `OwnsBase()`, `InheritsBase()`, or `NoBase()` (default).

Models override this method to declare their base quantity requirements.

# Example
```julia
# Bus owns base voltage
base_quantity_role(::Type{Bus}, ::Type{BaseVoltage}) = OwnsBase()

# PQ inherits base voltage from Bus
base_quantity_role(::Type{PQ}, ::Type{BaseVoltage}) = InheritsBase()

# PQ doesn't interact with base frequency (default)
base_quantity_role(::Type{PQ}, ::Type{BaseFrequency}) = NoBase()
```
"""
base_quantity_role(::Type{<:AbstractModel}, ::Type{<:ReferenceType}) = NoBase()

# ============================================================================
# Query Helpers
# ============================================================================

"""
    owns_base(::Type{<:AbstractModel}, ::Type{<:ReferenceType}) -> Bool

Check if a model owns a base quantity (stores it in a field).
"""
owns_base(M::Type, B::Type) = base_quantity_role(M, B) isa OwnsBase

"""
    inherits_base(::Type{<:AbstractModel}, ::Type{<:ReferenceType}) -> Bool

Check if a model inherits a base quantity from connected models.
"""
inherits_base(M::Type, B::Type) = base_quantity_role(M, B) isa InheritsBase

"""
    needs_base(::Type{<:AbstractModel}, ::Type{<:ReferenceType}) -> Bool

Check if a model needs a base quantity (owns or inherits it).
"""
needs_base(M::Type, B::Type) = !isa(base_quantity_role(M, B), NoBase)

# ============================================================================
# Base Quantity Resolution
# ============================================================================

"""
    get_base_quantity_vector(sys::SystemModel, model_name::Symbol,
                             base_type::Type{<:ReferenceType}; quantity = :values) -> Vector{Float64}

Get per-instance base quantity values for a model using the trait system.

# Arguments
- `sys`: SystemModel with solved power flow
- `model_name`: Name of the model (e.g., `:Bus`, `:PQ`)
- `base_type`: Type of base quantity (e.g., `BaseVoltage`, `BasePower`)
- `quantity`: The name of the quantity to convert on that base, for the refusal of a model
  without the base (e.g. the observation `:v`)

# Returns
Vector of base quantity values, one per model instance. Throws an `ArgumentError` naming the
model and `quantity` where the model has no such base (`NoBase`): a kV observation needs a base
voltage, which the model owns or inherits.

# Examples
```julia
# Get base voltages for PQ loads (inherited from connected buses)
base_kv = get_base_quantity_vector(sys, :PQ, BaseVoltage)

# Get base voltages for buses (owned by Bus model)
base_kv = get_base_quantity_vector(sys, :Bus, BaseVoltage)
```
"""
function get_base_quantity_vector(sys::SystemModel, model_name::Symbol,
                                   base_type::Type{<:ReferenceType}; quantity::Symbol = :values)
    model = sys.models[model_name]
    model_type = get_base_model_type(model)
    role = base_quantity_role(model_type, base_type)
    return _get_base_vector(sys, model_name, model, base_type, role, quantity)
end

# ============================================================================
# Trait Dispatch Implementations
# ============================================================================

"""
    _get_base_vector(sys, model_name, model, base_type, ::OwnsBase, quantity)

Retrieve base quantity from model's own field.
"""
function _get_base_vector(sys::SystemModel, model_name::Symbol, model,
                          base_type::Type{<:ReferenceType}, ::OwnsBase, quantity::Symbol)
    field_name = _base_field_name(base_type)
    return collect(getproperty(model, field_name))
end

"""
    _get_base_vector(sys, model_name, model, base_type, ::InheritsBase, quantity)

Retrieve base quantity from ExtParam (populated from connected models).
"""
function _get_base_vector(sys::SystemModel, model_name::Symbol, model,
                          base_type::Type{<:ReferenceType}, ::InheritsBase, quantity::Symbol)
    param_name = _ext_param_name(base_type)
    return collect(sys.external[model_name][param_name])
end

"""
    _get_base_vector(sys, model_name, model, base_type, ::NoBase, quantity)

Refuse: the model has no such base, so `quantity` cannot be converted on it.
"""
function _get_base_vector(sys::SystemModel, model_name::Symbol, model,
                          base_type::Type{<:ReferenceType}, ::NoBase, quantity::Symbol)
    throw(ArgumentError(
        "Cannot convert $model_name.$quantity on its $(nameof(base_type)): $model_name has none " *
        "(base_quantity_role($(nameof(get_base_model_type(model))), $(nameof(base_type))) is NoBase). " *
        "A kV observation needs a base voltage: declare base_quantity_role OwnsBase() or InheritsBase() for the model"))
end

# ============================================================================
# Field and Parameter Name Mapping
# ============================================================================
# These functions map ReferenceType to actual field/parameter names.
# Extend for new base quantity types.

"""
    _base_field_name(::Type{<:ReferenceType}) -> Symbol

Map ReferenceType to field name in model structs (for OwnsBase).
"""
_base_field_name(::Type{BaseVoltage}) = :basekv
# Future extensions:
# _base_field_name(::Type{BasePower}) = :mbase      # For generators with ratings
# _base_field_name(::Type{BaseInertia}) = :H        # For synchronous machines

"""
    _ext_param_name(::Type{<:ReferenceType}) -> Symbol

Map ReferenceType to parameter name in external params (for InheritsBase).
"""
_ext_param_name(::Type{BaseVoltage}) = :basekv
# Future extensions:
# _ext_param_name(::Type{BasePower}) = :sbase       # System base power
# _ext_param_name(::Type{BaseFrequency}) = :base_freq  # System base frequency
