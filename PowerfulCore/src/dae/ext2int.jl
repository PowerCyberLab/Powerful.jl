"""
External-to-internal indexing for ExtParam retrieval.

This module handles ONLY external parameter (ExtParam) retrieval.
It is SEPARATE from ForeignVar/SharedRes which is handled by address_manager.jl.

# Scope
- **This module**: Maps external keys → internal indices for DATA retrieval (ExtParam)
- **Address manager**: Maps ForeignVar/SharedRes → addresses for COMPUTATION (separate, already exists)

# Key Functions
- `ext2int`: Build IndexOrder with index maps and foreign refs
- `build_index_map`: Create composite key → array index lookup using get_key_name()
- `build_terminal_refs`: Build foreign refs from model terminals (single source of truth for connectivity)

# See Also
- MATPOWER's ext2int.m for design inspiration
- address_manager.jl for ForeignVar/SharedRes (separate system)
"""

"""
    ext2int(models::NamedTuple) -> (NamedTuple, IndexOrder)

Build index mappings for external parameter retrieval following MATPOWER's ext2int pattern.

This function creates lookup tables to map external bus numbers (and composite keys) to
internal array indices. Called BEFORE per-unit conversion.

# Algorithm
1. For each model, build composite key → array index map using `get_key_name()`
2. For each model with ExtParam, build indexer vectors for fast lookup

# Arguments
- `models::NamedTuple`: Parsed models with external bus numbers

# Returns
- `models`: Same models (unchanged - external numbers preserved)
- `order::IndexOrder`: Index mappings for ExtParam retrieval

# Example
```julia
models, order = ext2int(models)
# order.key_indices[:Bus][(69953,)] = 67900  # Bus 69953 → array index 67900
# order.foreign_refs[:PQ][:i] = [1, 1, 2, ...]  # PQ.i → Bus array indices
# order.foreign_refs[:Branch][:i] = [5, 7, ...]  # Branch.i → from-Bus array indices
# order.foreign_refs[:Branch][:j] = [8, 9, ...]  # Branch.j → to-Bus array indices
```
"""
function ext2int(models::NamedTuple)
    # Step 1: Build model instance key → index maps using get_key_name()
    key_indices_dict = Dict{Symbol, Dict}()

    for model_name in keys(models)
        model = models[model_name]
        if isempty(model)
            key_indices_dict[model_name] = Dict()
            continue
        end

        model_type = get_base_model_type(model)
        key_fields = get_key_name(model_type)

        # Build lookup: (external_keys...) → array_index
        index_map = build_index_map(model, key_fields)
        key_indices_dict[model_name] = index_map
    end

    # Step 2: Build foreign refs from TERMINALS (unified approach)
    # This builds refs for all connectivity: ExtParams, ForeignVars, SharedRes
    foreign_refs_dict = Dict{Symbol, Dict{GroupRefKey, GroupLinkRef}}()

    for model_name in keys(models)
        model = models[model_name]
        if isempty(model)
            foreign_refs_dict[model_name] = Dict{GroupRefKey, GroupLinkRef}()
            continue
        end

        # Build refs from terminals (the single source of truth for connectivity)
        refs = build_terminal_refs(model, key_indices_dict, models)
        foreign_refs_dict[model_name] = refs
    end

    # Step 3: Create IndexOrder
    order = IndexOrder(key_indices_dict, foreign_refs_dict)

    # Models remain unchanged (external numbers preserved)
    return models, order
end

"""
    instance_key(model, fields, i) -> Tuple

Key of instance `i`: the values of `fields` (a Symbol or a collection of Symbols)
as a tuple. Single-field keys are 1-tuples, e.g. `(101,)` for Bus 101.

This is the one definition of an instance key; `key_indices`, terminal matching,
group lookups and `get_key` all use it.
"""
@inline instance_key(model, field::Symbol, i) = (getproperty(model, field)[i],)
@inline instance_key(model, fields, i) = Tuple(getproperty(model, f)[i] for f in fields)

"""
    build_index_map(model, key_fields) -> Dict

Build instance key → array index lookup for a single model (see `instance_key`).
Errors if two instances share a key.

# Examples
```julia
build_index_map(buses, :i)[(69953,)]        # → 67900
build_index_map(loads, [:i, :id])[(101, "2")]  # → 5
```
"""
function build_index_map(model, key_fields)
    index_map = Dict()
    for i in 1:length(model)
        key = instance_key(model, key_fields, i)
        haskey(index_map, key) &&
            error("Duplicate key $key for fields $key_fields at instances $(index_map[key]) and $i")
        index_map[key] = i
    end
    return index_map
end

"""
    normalize_ref_key(ref_key) -> Union{Symbol, Tuple{Vararg{Symbol}}}

Normalize a reference key specification to a canonical dictionary key format.
Symbols pass through unchanged; vectors/tuples become tuples.

This is used for both Terminal.ref_key and ModelVar.ref_key to ensure
consistent dictionary keys in foreign_refs lookups.
"""
normalize_ref_key(ref_key::Union{Symbol, Vector{Symbol}, Tuple{Vararg{Symbol}}}) =
    ref_key isa Symbol ? ref_key : Tuple(ref_key)


"""
    build_group_lookup(group_name, key_indices, models) -> Dict{Tuple, Tuple{Symbol, Int}}

Construct a lookup from instance keys to `(model_sym, index)` for all models
implementing the given group, from the pre-computed `key_indices` (which must
contain every such model).
"""
function build_group_lookup(
    group_name::Symbol,
    key_indices::AbstractDict,
    models::NamedTuple
)
    lookup = Dict{Tuple, Tuple{Symbol, Int}}()

    for (model_sym, model_data) in pairs(models)
        isempty(model_data) && continue

        model_type = get_base_model_type(model_data)
        Foundation.implements_group(model_type) == group_name || continue

        for (ext_key, idx) in key_indices[model_sym]
            lookup[ext_key] = (model_sym, idx)
        end
    end

    return lookup
end

"""
    build_terminal_refs(model, key_indices_dict, models) -> Dict{GroupRefKey, GroupLinkRef}

Build foreign refs from model TERMINALS (unified approach).

This is the single source of truth for connectivity. Terminals define how models
connect to other models or groups. This function builds GroupLinkRef for each terminal,
supporting both direct model references and group references.

# Key Format (Unified)
All references are keyed by `normalize_ref_key(terminal.ref_key)`, regardless of
whether the terminal references a group or a direct model. This provides a consistent
lookup pattern: `foreign_refs[model_name][normalize_ref_key(ref_key)]`.

# Determining Group vs Direct
- Group references produce heterogeneous GroupLinkRef (target_type is Vector{Symbol})
- Direct references produce homogeneous GroupLinkRef (all target_type entries are same)
Use `is_heterogeneous(link_ref)` to distinguish after lookup.

# Arguments
- `model`: StructArray of model instances
- `key_indices_dict`: Global key → index maps for all models
- `models`: NamedTuple of all models (for group reference resolution)

# Returns
Dictionary mapping normalized ref_key to GroupLinkRef

# Example
```julia
# PQ connects to Bus via Terminal(:Bus, :i)
# Result: Dict(:i => GroupLinkRef(:Bus, [1, 2, 3, ...]))

# GENCLS connects to bus_group via Terminal(:bus_group, :BUS)
# Result: Dict(:BUS => GroupLinkRef([:Bus, :Bus, ...], [1, 2, ...]))

# TGOV1 connects to generators via Terminal(:generators, [:BUS, :ID])
# Result: Dict((:BUS, :ID) => GroupLinkRef([:GENCLS, :GENROU, ...], [1, 1, ...]))
```
"""
function build_terminal_refs(model, key_indices_dict::Dict, models::NamedTuple)
    refs = Dict{GroupRefKey, GroupLinkRef}()

    model_type = get_base_model_type(model)
    terminals = get_terminals(model_type)

    isempty(terminals) && return refs

    n = length(model)

    for terminal in terminals
        target_ref = terminal.vertex_model
        terminal_ref_key = terminal.ref_key

        # UNIFIED: Always key by normalized ref_key
        dict_key = normalize_ref_key(terminal_ref_key)

        # Skip if already built (multiple terminals may share same key)
        haskey(refs, dict_key) && continue

        # Determine if this is a group reference
        is_group = Foundation.is_group_reference(target_ref)

        if is_group
            # Group reference - search all models implementing this group
            group_lookup = build_group_lookup(target_ref, key_indices_dict, models)

            isempty(group_lookup) &&
                error("$(model_type) terminal :$(terminal_ref_key) references group :$(target_ref), but no model in the system implements it")

            # Build heterogeneous GroupLinkRef
            target_types = Vector{Symbol}(undef, n)
            target_indices = zeros(Int, n)

            for i in 1:n
                key = instance_key(model, terminal_ref_key, i)
                if haskey(group_lookup, key)
                    target_types[i], target_indices[i] = group_lookup[key]
                else
                    error("$(model_type) terminal :$(terminal_ref_key) references key $(key) in group :$(target_ref), but no model in the group has this key")
                end
            end

            refs[dict_key] = GroupLinkRef(target_types, target_indices)
        else
            # Direct model reference
            haskey(key_indices_dict, target_ref) ||
                error("$(model_type) terminal :$(terminal_ref_key) references model :$(target_ref), which is not in the system")

            target_indices_map = key_indices_dict[target_ref]

            # Build homogeneous GroupLinkRef
            target_types = fill(target_ref, n)
            target_indices = zeros(Int, n)

            for i in 1:n
                key = instance_key(model, terminal_ref_key, i)
                if haskey(target_indices_map, key)
                    target_indices[i] = target_indices_map[key]
                else
                    error("$(model_type) terminal :$(terminal_ref_key) references key $(key) in :$(target_ref), but no such key exists")
                end
            end

            refs[dict_key] = GroupLinkRef(target_types, target_indices)
        end
    end

    return refs
end


"""
    get_index_by_key(sys::SystemModel, model::Symbol, key::Tuple) -> Int

Lookup internal array index for a model instance by its composite key.

This is a dictionary lookup: given a model's unique key (external bus number,
composite key, etc.), return the internal array index where that instance is stored.

# Arguments
- `sys`: SystemModel with IndexOrder
- `model`: Model type symbol (`:Bus`, `:PV`, `:Branch`, etc.)
- `key`: Composite key tuple (e.g., `(101,)` for Bus, `(101, "1")` for PV)

# Returns
Internal array index (Int) where the instance is stored

# Example
```julia
# Lookup Bus 101
idx = get_index_by_key(sys, :Bus, (101,))
# → 67  (Bus 101 is stored at array index 67)

# Lookup PV generator at bus 101, id "1"
idx = get_index_by_key(sys, :PV, (101, "1"))
# → 5  (This PV instance is at array index 5)

# Use to access model data
bus_data = sys.models.Bus[idx]
```

# See Also
- `get_index_by_key(sys, model, keys::Vector)` - Batch version
- `get_foreign_indices` - Get all foreign reference indices for a field
"""
function get_index_by_key(sys::SystemModel, model::Symbol, key)
    idx = find_index(sys, model, key)
    idx === nothing && error("No :$model instance with key $key")
    return idx
end

"""
    find_index(sys::SystemModel, model::Symbol, key) -> Union{Int, Nothing}

Internal array index of the `model` instance with external `key`, or `nothing`
if no such instance exists. The key is normalized with `normalize_key`. Errors if
`model` itself is not in the system.
"""
function find_index(sys::SystemModel, model::Symbol, key)
    haskey(sys.order.key_indices, model) || error("Model :$model is not in the system")
    return get(sys.order.key_indices[model], normalize_key(sys, model, key), nothing)
end

"""
    normalize_key(sys::SystemModel, model::Symbol, key) -> Tuple

Canonical form of a user-supplied instance key: a scalar becomes a 1-tuple and each
element is converted to the element type of its key field (e.g. `101 → Int32(101)`,
`"1" → String3("1")`), matching the entries of `key_indices`.
"""
function normalize_key(sys::SystemModel, model::Symbol, key)
    data = sys.models[model]
    fields = get_key_name(get_base_model_type(data))
    names = fields isa Symbol ? (fields,) : Tuple(fields)
    key_tuple = key isa Tuple ? key : (key,)
    length(key_tuple) == length(names) ||
        error("Key $key for :$model must have $(length(names)) element(s) $names")
    return Tuple(convert(eltype(getproperty(data, f)), v) for (f, v) in zip(names, key_tuple))
end

"""
    get_index_by_key(sys::SystemModel, model::Symbol, keys::AbstractVector{<:Tuple}) -> Vector{Int}

Batch lookup of internal array indices for multiple model instances.

# Arguments
- `sys`: SystemModel with IndexOrder
- `model`: Model type symbol
- `keys`: Vector of composite key tuples

# Returns
Vector of internal array indices

# Example
```julia
# Batch lookup for multiple buses
keys = [(101,), (102,), (103,)]
indices = get_index_by_key(sys, :Bus, keys)
# → [67, 68, 69]
```
"""
function get_index_by_key(sys::SystemModel, model::Symbol, keys::AbstractVector)
    return [get_index_by_key(sys, model, k) for k in keys]
end

"""
    get_field_by_key(sys::SystemModel, model::Symbol, key::Tuple, field::Symbol) -> Any

Get field value from a model instance by its composite key.

Convenience function combining `get_index_by_key` with field access.

# Arguments
- `sys`: SystemModel
- `model`: Model type symbol
- `key`: Composite key tuple
- `field`: Field name to retrieve

# Returns
Field value from the model instance

# Example
```julia
# Get setpoint from PV generator at bus 101, id "1"
pg = get_field_by_key(sys, :PV, (101, "1"), :pg)

# Equivalent to:
idx = get_index_by_key(sys, :PV, (101, "1"))
pg = sys.models.PV.pg[idx]
```
"""
function get_field_by_key(sys::SystemModel, model::Symbol, key, field::Symbol)
    idx = get_index_by_key(sys, model, key)
    return getproperty(sys.models[model], field)[idx]
end

"""
    get_field_by_key(sys::SystemModel, model::Symbol, keys::AbstractVector{<:Tuple}, field::Symbol) -> Vector

Batch retrieval of field values for multiple model instances.

# Arguments
- `sys`: SystemModel
- `model`: Model type symbol
- `keys`: Vector of composite key tuples
- `field`: Field name to retrieve

# Returns
Vector of field values

# Example
```julia
# Get setpoints from multiple PV generators
keys = [(101, "1"), (102, "1"), (103, "1")]
pg_vals = get_field_by_key(sys, :PV, keys, :pg)
# → [2.0, 1.5, 3.0]
```
"""
function get_field_by_key(sys::SystemModel, model::Symbol, keys::AbstractVector, field::Symbol)
    indices = get_index_by_key(sys, model, keys)
    field_array = getproperty(sys.models[model], field)
    return [field_array[idx] for idx in indices]
end

"""
    get_foreign_indices(sys::SystemModel, model::Symbol, field::Symbol) -> Vector{Int}

Get all foreign reference indices for a model's reference field.

This retrieves pre-computed indexer vectors that map each instance of a model
to the internal array index of the referenced model instance.

# Arguments
- `sys`: SystemModel with IndexOrder
- `model`: Model with foreign references (e.g., `:PV`, `:Branch`)
- `field`: Field holding the foreign reference (e.g., `:i`, `:j`)

# Returns
Vector of internal array indices where `result[i]` is the internal index in the
referenced model for instance `i` of the referencing model.

# Example
```julia
# Get bus indices for all PV generators
bus_indices = get_foreign_indices(sys, :PV, :i)
# → [1, 1, 3, 5, 7, ...]  (PV[1] → Bus[1], PV[2] → Bus[1], etc.)

# Get from-bus and to-bus indices for all branches
from_bus = get_foreign_indices(sys, :Branch, :i)
to_bus = get_foreign_indices(sys, :Branch, :j)
```

# See Also
- `get_index_by_key` - Single key-based lookup
- `get_foreign_link_ref` - Get full GroupLinkRef (includes target types)
"""
function get_foreign_indices(sys::SystemModel, model::Symbol, field::Union{Symbol, Tuple})
    return get_group_ref(sys.order, model, field).target_idx
end

"""
    get_foreign_link_ref(sys::SystemModel, model::Symbol, field::Symbol) -> GroupLinkRef

Get the full GroupLinkRef for a model's reference field.

Unlike `get_foreign_indices` which returns only indices, this returns the full
GroupLinkRef including target types for each instance. Useful when references
may be heterogeneous (pointing to different model types).

# Arguments
- `sys`: SystemModel with IndexOrder
- `model`: Model with foreign references (e.g., `:PV`, `:Branch`)
- `field`: Field holding the foreign reference (e.g., `:i`, `:j`)

# Returns
GroupLinkRef with target_type and target_idx vectors

# Example
```julia
# Get full reference info for generators connected to bus_group
link_ref = get_foreign_link_ref(sys, :GENCLS, :i)
# link_ref.target_type = [:Bus, :Bus, ...]  (the bus_group model of each generator)
# link_ref.target_idx = [1, 2, ...]
```
"""
function get_foreign_link_ref(sys::SystemModel, model::Symbol, field::Union{Symbol, Tuple})
    return get_group_ref(sys.order, model, field)
end

# ============================================================================
# Group-Based Control Links (heterogeneous targets)
# ============================================================================

# GroupLinkRef is now defined in Foundation.core_types and imported via DAEInfra module

# ============================================================================
# GroupLinkRef Helper Functions
# ============================================================================

"""
    gather_from_group(models::NamedTuple, link_ref::GroupLinkRef, field::Symbol) -> Vector

Gather field values from heterogeneous targets using GroupLinkRef.

For each instance, retrieves the specified field from the appropriate target model
based on the target_type and target_idx in the GroupLinkRef.

# Arguments
- `models::NamedTuple`: All models in the system
- `link_ref::GroupLinkRef`: Reference with target types and indices
- `field::Symbol`: Field name to retrieve from targets

# Returns
Vector of field values, one per source instance

# Example
```julia
# Get basekv values for PQ loads from their connected buses
link_ref = sys.order.foreign_refs[:PQ][:i]
basekv = gather_from_group(sys.models, link_ref, :basekv)
```
"""
function gather_from_group(models::NamedTuple, link_ref::GroupLinkRef, field::Symbol)
    n = length(link_ref.target_type)
    length(link_ref.target_idx) == n ||
        error("GroupLinkRef corrupted: target_type length ($n) != target_idx length ($(length(link_ref.target_idx)))")
    n == 0 && error("gather_from_group: empty GroupLinkRef for field :$field")

    T = eltype(getproperty(getfield(models, link_ref.target_type[1]), field))
    result = Vector{T}(undef, n)
    for i in 1:n
        # Direct indexing: a missing target model/field or bad index errors loudly
        result[i] = getproperty(getfield(models, link_ref.target_type[i]), field)[link_ref.target_idx[i]]
    end
    return result
end

export ext2int
export get_index_by_key, find_index, normalize_key, get_field_by_key, get_foreign_indices, get_foreign_link_ref
export normalize_ref_key
export instance_key, build_terminal_refs
export gather_from_group
