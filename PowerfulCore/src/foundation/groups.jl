# ============================================================================
# Group System - First-Class Group Interfaces for Dynamic Models
# ============================================================================
#
# This module provides a formalized group system that makes groups first-class
# citizens in the type system. Groups define contracts (interfaces) that models
# must implement, enabling polymorphic references (e.g., TGOV1 can reference
# :generators group instead of hard-coding :GENCLS).
#
# Key concepts:
# - GroupInterface: Defines the contract for a group (required vars and residuals)
# - GroupVar: A variable that all models in the group must expose
# - GroupRes: A residual that all models in the group must expose
# - group_var_mapping: Maps group canonical names to model-specific names
#
# ============================================================================

"""
    GroupVar

Specification for a variable that all models in a group must expose.

# Fields
- `name::Symbol`: Canonical name in the group interface (e.g., :omega, :delta)
- `kind::Symbol`: Variable kind - :state or :algebraic
- `description::String`: Human-readable description

# Example
```julia
GroupVar(:omega, :state, "Rotor speed in per-unit")
GroupVar(:P_m, :algebraic, "Mechanical power input")
```
"""
struct GroupVar
    name::Symbol
    kind::Symbol  # :state or :algebraic
    description::String

    function GroupVar(name::Symbol, kind::Symbol, description::String)
        kind in (:state, :algebraic) || error("GroupVar kind must be :state or :algebraic, got :$kind")
        new(name, kind, description)
    end
end

GroupVar(name::Symbol, kind::Symbol; description::String="") = GroupVar(name, kind, description)

"""
    GroupRes

Specification for a residual that all models in a group must expose.

# Fields
- `name::Symbol`: Canonical name in the group interface (e.g., :g_P_m, :g_E_fd)
- `kind::Symbol`: Residual kind - :state or :algebraic
- `access::Symbol`: Access pattern - :private or :shared
- `description::String`: Human-readable description

# Example
```julia
GroupRes(:g_P_m, :algebraic, :shared, "Mechanical power input residual")
```
"""
struct GroupRes
    name::Symbol
    kind::Symbol    # :state or :algebraic
    access::Symbol  # :private or :shared
    description::String

    function GroupRes(name::Symbol, kind::Symbol, access::Symbol, description::String)
        kind in (:state, :algebraic) || error("GroupRes kind must be :state or :algebraic, got :$kind")
        new(name, kind, access, description)
    end
end

GroupRes(name::Symbol, kind::Symbol, access::Symbol; description::String="") =
    GroupRes(name, kind, access, description)

"""
    GroupInterface

Defines the contract that all models in a group must satisfy.

Acts as a "trait bundle" for a category of models. Models declare they implement
a group interface, and the system validates they provide the required variables and
residuals. A model initialized from a member of the group (`InitInput` of kind
`:from_group_init`) reads the initial value of one of its required variables, through
the member's `group_var_mapping`.

# Fields
- `name::Symbol`: Group name (e.g., :generators, :governors, :exciters)
- `required_vars::Vector{GroupVar}`: Variables all members must expose
- `required_residuals::Vector{GroupRes}`: Residuals all members must expose
- `optional_residuals::Vector{GroupRes}`: Residuals a member may expose; a model that delivers
  through one is routed to a member that has it (`InjectionForms`)
- `key_fields::Vector{Symbol}`: Composite key for matching (e.g., [:BUS, :ID])

# Example
```julia
GroupInterface(
    name = :generators,
    required_vars = [
        GroupVar(:omega, :state; description="Rotor speed"),
        GroupVar(:delta, :state; description="Rotor angle"),
        GroupVar(:P_m, :algebraic; description="Mechanical power input"),
    ],
    required_residuals = [
        GroupRes(:g_P_m, :algebraic, :shared; description="Mechanical power residual"),
    ],
    key_fields = [:BUS, :ID]
)
```
"""
struct GroupInterface
    name::Symbol
    required_vars::Vector{GroupVar}
    required_residuals::Vector{GroupRes}
    optional_residuals::Vector{GroupRes}
    key_fields::Vector{Symbol}
end

GroupInterface(name, required_vars, required_residuals, key_fields) =
    GroupInterface(name, required_vars, required_residuals, GroupRes[], key_fields)

# Keyword constructor for convenience
function GroupInterface(;
    name::Symbol,
    required_vars::Vector{GroupVar} = GroupVar[],
    required_residuals::Vector{GroupRes} = GroupRes[],
    optional_residuals::Vector{GroupRes} = GroupRes[],
    key_fields::Vector{Symbol} = Symbol[]
)
    GroupInterface(name, required_vars, required_residuals, optional_residuals, key_fields)
end

# ============================================================================
# Group Registry
# ============================================================================

"""
Global registry of group interfaces.

Maps group name (Symbol) to GroupInterface. Populated at module load time
with standard groups, and can be extended with custom groups.
"""
const GROUP_INTERFACES = Dict{Symbol, GroupInterface}()

"""
    register_group!(interface::GroupInterface)

Register a group interface in the global registry.

# Example
```julia
register_group!(GroupInterface(
    name = :custom_controllers,
    required_vars = [GroupVar(:output, :algebraic)],
    required_residuals = [],
    key_fields = [:BUS, :ID]
))
```
"""
function register_group!(interface::GroupInterface)
    if haskey(GROUP_INTERFACES, interface.name)
        @warn "Overwriting existing group interface" group=interface.name
    end
    GROUP_INTERFACES[interface.name] = interface
    return interface
end

"""
    get_group_interface(group::Symbol) -> GroupInterface

Retrieve a registered group interface by name.

# Throws
- `KeyError` if group is not registered
"""
function get_group_interface(group::Symbol)
    if !haskey(GROUP_INTERFACES, group)
        error("Group :$group is not registered. Available groups: $(keys(GROUP_INTERFACES))")
    end
    return GROUP_INTERFACES[group]
end

"""
    has_group_interface(group::Symbol) -> Bool

Check if a group interface is registered.
"""
has_group_interface(group::Symbol) = haskey(GROUP_INTERFACES, group)

"""
    all_groups() -> Vector{Symbol}

Return names of all registered group interfaces.
"""
all_groups() = collect(keys(GROUP_INTERFACES))

# ============================================================================
# Model Traits for Group Implementation
# ============================================================================

"""
    implements_group(::Type{M}) -> Union{Symbol, Nothing}

Declare which group interface this model implements.

Models implement this trait to declare group membership. Returns the group name
(Symbol) or `nothing` if the model doesn't implement any group.

# Examples
```julia
# GENCLS implements the :generators group
implements_group(::Type{<:GENCLS}) = :generators

# TGOV1 implements the :governors group
implements_group(::Type{<:TGOV1}) = :governors

# PQ load doesn't implement any dynamic group
implements_group(::Type{<:PQ}) = nothing
```
"""
implements_group(::Type{<:AbstractModel}) = nothing

"""
    get_controlled_group(::Type{M}) -> Union{Symbol, Nothing}

Return the group this controller controls (e.g., `:generators` for TGOV1).

Returns `nothing` for non-controller models (default).
"""
get_controlled_group(::Type{<:AbstractModel}) = nothing

"""
    group_var_mapping(::Type{M}) -> Dict{Symbol, Symbol}

Map group canonical variable names to model's actual variable names.

Keys are group names, values are model names. Example:
`Dict(:E_fd => :E_q, :omega => :omega)` maps group's `:E_fd` to GENCLS's `:E_q`.
"""
function group_var_mapping(::Type{M}) where M <: AbstractModel
    # Default: empty mapping (model doesn't implement any group)
    return Dict{Symbol, Symbol}()
end

"""
    group_res_mapping(::Type{M}) -> Dict{Symbol, Symbol}

Map group canonical residual names to model's actual residual names.

Example: `Dict(:g_E_fd => :g_E_q)` maps group's `:g_E_fd` to GENCLS's `:g_E_q`.
"""
function group_res_mapping(::Type{M}) where M <: AbstractModel
    return Dict{Symbol, Symbol}()
end

# ============================================================================
# Group Validation
# ============================================================================

"""
    validate_group_implementation(::Type{M}) -> Bool

Validate that a model correctly implements its declared group interface.

Checks that all required variables and residuals declared by
the group interface are provided by the model's metadata and mappings.

# Throws
- `AssertionError` with detailed message if validation fails

# Returns
`true` if validation passes
"""
function validate_group_implementation(::Type{M}) where M <: AbstractModel
    group = implements_group(M)

    # Models that don't implement a group pass trivially
    if group === nothing
        return true
    end

    # Get group interface
    if !has_group_interface(group)
        error("Model $M declares implements_group() = :$group but group :$group is not registered")
    end
    interface = get_group_interface(group)

    # Get model metadata and mappings
    metadata = model_metadata(M)
    var_map = group_var_mapping(M)
    res_map = group_res_mapping(M)

    # Validate key fields: model may have more fields than group (e.g., a model keyed by
    # [:BUS, :ID] for per-generator matching implementing a group keyed by [:i]).
    # Require: model has at least as many key fields as group, primary key type compatible.
    model_key = get_key_name(M)
    model_key = model_key isa Symbol ? [model_key] : collect(model_key)
    group_key = interface.key_fields
    if length(model_key) < length(group_key)
        error("Model $(nameof(M)) key fields $model_key have fewer fields than group :$group requires ($group_key)")
    elseif length(model_key) > length(group_key)
        # Allow composite keys - lookup will use primary key fallback
        @debug "Model $(nameof(M)) has composite key $model_key for group :$group (expected $group_key); primary key matching will be used."
    elseif model_key != group_key
        @warn "Model $(nameof(M)) implements group :$group with differing key field names $model_key; expected $group_key. Matching by position." maxlog=1
    end

    # Collect model variable names for lookup
    model_var_names = Set(v.name for v in metadata.vars)
    model_res_names = Set(r.name for r in metadata.residuals)

    # Validate required variables
    for gvar in interface.required_vars
        if !haskey(var_map, gvar.name)
            error("Model $M implements group :$group but missing var mapping for :$(gvar.name). " *
                  "Add `$(gvar.name) => :your_var_name` to group_var_mapping.")
        end
        model_var = var_map[gvar.name]
        if model_var ∉ model_var_names
            error("Model $M maps group var :$(gvar.name) to :$model_var but :$model_var " *
                  "not found in model_metadata vars. Available: $model_var_names")
        end
    end

    # Validate required residuals
    for gres in interface.required_residuals
        if !haskey(res_map, gres.name)
            error("Model $M implements group :$group but missing residual mapping for :$(gres.name). " *
                  "Add `$(gres.name) => :your_res_name` to group_res_mapping.")
        end
        model_res = res_map[gres.name]
        if model_res ∉ model_res_names
            error("Model $M maps group residual :$(gres.name) to :$model_res but :$model_res " *
                  "not found in model_metadata residuals. Available: $model_res_names")
        end
    end

    # An optional residual a member maps must be one it has
    for gres in interface.optional_residuals
        haskey(res_map, gres.name) || continue
        model_res = res_map[gres.name]
        if model_res ∉ model_res_names
            error("Model $M maps group residual :$(gres.name) to :$model_res but :$model_res " *
                  "not found in model_metadata residuals. Available: $model_res_names")
        end
    end

    return true
end

# ============================================================================
# Group Reference Resolution Utilities
# ============================================================================

"""
    resolve_group_var_name(group::Symbol, canonical_name::Symbol, target_type::Type) -> Symbol

Resolve a group canonical variable name to the actual variable name for a specific model type.

# Arguments
- `group::Symbol`: Group name (e.g., :generators)
- `canonical_name::Symbol`: Canonical name in group interface (e.g., :omega)
- `target_type::Type`: Target model type (e.g., GENCLS)

# Returns
Model's actual variable name (e.g., :omega for GENCLS)

# Throws
Error if mapping not found
"""
resolve_group_var_name(group::Symbol, canonical_name::Symbol, target_type::Type) =
    _resolve_group_name(group_var_mapping(target_type), "variable", group, canonical_name, target_type)

"""
    resolve_group_res_name(group::Symbol, canonical_name::Symbol, target_type::Type) -> Symbol

Resolve a group canonical residual name to the actual residual name for a specific model type.
"""
resolve_group_res_name(group::Symbol, canonical_name::Symbol, target_type::Type) =
    _resolve_group_name(group_res_mapping(target_type), "residual", group, canonical_name, target_type)

function _resolve_group_name(mapping, kind::String, group::Symbol, canonical_name::Symbol, target_type::Type)
    if !haskey(mapping, canonical_name)
        error("Model $target_type does not map group :$group $kind :$canonical_name. " *
              "Available mappings: $mapping")
    end
    return mapping[canonical_name]
end

# ============================================================================
# Helper: Check if a reference is to a group
# ============================================================================

"""
    is_group_reference(model_or_group::Symbol) -> Bool

Check if a symbol refers to a registered group or a specific model.
"""
is_group_reference(model_or_group::Symbol) = has_group_interface(model_or_group)

# ============================================================================
# Exports
# ============================================================================

export GroupVar, GroupRes, GroupInterface
export register_group!, get_group_interface, has_group_interface, all_groups
export implements_group, group_var_mapping, group_res_mapping
export validate_group_implementation
export resolve_group_var_name, resolve_group_res_name
export is_group_reference
