const ADDRESSABLE_TYPES = [AlgebVar, AlgebRes, StateVar, StateRes, ObservableVar]

"""
$(SIGNATURES)

Instance keys of every element of a model, as tuples (see `instance_key`):
`(101,)` for Bus 101, `(101, "1")` for a PQ load. These are the keys of
`sys.order.key_indices[model]`. Key columns are copied to the CPU first.
"""
function get_key(model_vec::StructArray{<:AbstractModel})
    fields = get_key_name(get_base_model_type(model_vec))
    names = fields isa Symbol ? (fields,) : Tuple(fields)
    cpu_cols = NamedTuple{names}(map(f -> Array(getproperty(model_vec, f)), names))
    return [instance_key(cpu_cols, names, i) for i in 1:length(model_vec)]
end

# Note: get_key_name base function is defined in Foundation (model_metadata.jl)
# Model-specific overrides remain in individual model files.

# Note: ModelAddressAssignment and ModelRetrievedAddresses constructors
# are defined in Foundation (core_types.jl)

"""
Assign addresses for model variables and residuals
"""
function assign_address!(
    am::AddressManager,
    model_name::Symbol,
    metadata::ModelMetadata,
    device_count::Int
)
    # Input sanity checks
    device_count > 0 || return

    if haskey(am.assignments, model_name)
        throw(ArgumentError("Model $model_name already has partial/complete address assignment"))
    end

    # Always use contiguous variables layout
    _assign_contiguous_addresses!(am, model_name, metadata.vars, device_count)
    _assign_contiguous_addresses!(am, model_name, metadata.residuals, device_count)

    am.assignments[model_name].is_complete = true

end

"""
    _assign_contiguous_addresses!(am, model_name, allocable, device_count)

Assign each internal variable or residual in `allocable` a contiguous block of
`device_count` addresses in its address space (`next_idx`), recording them in
`am.addresses` and the model's assignment. External items are skipped; their
addresses are retrieved from the owning model.
"""
function _assign_contiguous_addresses!(
    am::AddressManager,
    model_name::Symbol,
    allocable::Vector{T},
    device_count::Int,
) where T

    if haskey(am.assignments, model_name)
        by_type = am.assignments[model_name].by_type
    else
        am.assignments[model_name] = ModelAddressAssignment(model_name)
        by_type = am.assignments[model_name].by_type
    end

    # filter out external vars using `is_internal`
    internal_vars = filter(is_internal, allocable)

    # Assign addresses for each internal variable/residual
    for var in internal_vars
        start_idx = get!(am.next_idx, var.address_type, 1)
        range = collect(UInt, start_idx:(start_idx + device_count - 1))

        # Assign addresses
        am.addresses[(var.address_type, model_name, var.name)] = range
        by_type[var.address_type][var.name] = range

        # Update next_idx for next allocation
        am.next_idx[var.address_type] = start_idx + device_count
    end

end

"""
Gather addresses from a heterogeneous GroupLinkRef (multiple target model types).

For each instance, resolves the canonical field name to the model-specific name
and looks up the address from the appropriate target model.

# Arguments
- `am::AddressManager` - Address manager with allocated addresses
- `models::NamedTuple` - All models in the system
- `model_name::Symbol` - Model requesting addresses
- `link_ref::GroupLinkRef` - Heterogeneous link ref (varying target_type per instance)
- `address_type::AddressableType` - Type of address (StateVar, AlgebVar, etc.)
- `group_name::Symbol` - Group being referenced (e.g., :generators)
- `canonical_field::Symbol` - Canonical field name in group interface

# Returns
`Vector{UInt}` of addresses, one per instance; with `routed`, an instance it marks false (delivered
through another port, see [`injection_routes`](@ref)) has the address [`UNROUTED`](@ref)
"""
function _gather_heterogeneous_addresses(
    am::AddressManager,
    models::NamedTuple,
    model_name::Symbol,
    link_ref::GroupLinkRef,
    address_type::AddressableType,
    group_name::Symbol,
    canonical_field::Symbol;
    routed::Union{Nothing, AbstractVector{Bool}} = nothing
)
    n = length(link_ref.target_type)

    # Validate GroupLinkRef consistency
    if length(link_ref.target_idx) != n
        error("GroupLinkRef corrupted for $model_name: target_type ($n) != target_idx ($(length(link_ref.target_idx)))")
    end

    result = Vector{UInt}(undef, n)

    for i in 1:n
        # An instance not delivered through this residual has no address in it
        if !isnothing(routed) && !routed[i]
            result[i] = UNROUTED
            continue
        end
        target_model = link_ref.target_type[i]
        target_idx = link_ref.target_idx[i]

        # Validate target model exists
        if !haskey(models, target_model)
            error("Heterogeneous ref: target :$target_model not found (from $model_name[$i])")
        end

        # Resolve canonical field to model-specific field
        model_type = get_base_model_type(models[target_model])
        actual_field = _resolve_canonical_name(group_name, canonical_field, model_type, address_type)

        # Get address from actual target model
        key = (address_type, target_model, actual_field)
        if !haskey(am.addresses, key)
            error("Address not found for $target_model.$actual_field (from group :$group_name)")
        end
        src_addrs = am.addresses[key]

        # Bounds validation
        if target_idx < 1 || target_idx > length(src_addrs)
            error("Heterogeneous ref: index $target_idx out of bounds for $target_model.$actual_field (1:$(length(src_addrs))) at $model_name[$i]")
        end

        result[i] = src_addrs[target_idx]
    end
    return result
end

"""
Resolve canonical name to model-specific name for group references.
Returns source_name unchanged if source_model is not a group reference.
"""
function _resolve_canonical_name(
    source_model::Symbol,
    source_name::Symbol,
    target_model_type,
    address_type::AddressableType
)
    if !Foundation.is_group_reference(source_model)
        return source_name
    end
    is_residual = address_type isa AlgebRes || address_type isa StateRes
    return is_residual ?
        Foundation.resolve_group_res_name(source_model, source_name, target_model_type) :
        Foundation.resolve_group_var_name(source_model, source_name, target_model_type)
end

"""
    UNROUTED

The address of an instance in a residual it is not delivered through: a model that
declares an injection in one form delivers each connection through exactly one of its two
forms (see [`injection_routes`](@ref)), and the residuals of the other form have no address
for that instance. Nothing is scattered to it.
"""
const UNROUTED = typemax(UInt)

# The forms a model declares the residuals of (power and current always; voltage where it has them)
_declared_forms(forms::InjectionForms) = isempty(forms.voltage) ? (:power, :current) : (:power, :current, :voltage)

# The group and terminal every form of an injection is delivered to, and the group
# fields of each form
function _injection_ports(metadata::ModelMetadata)
    forms = metadata.injection
    model_name = metadata.name
    by_name = Dict(r.name => r for r in metadata.residuals)
    delivered = Dict(
        form => [by_name[n] for n in getfield(forms, form)] for form in _declared_forms(forms))
    first_res = first(delivered[:power])
    group, ref_key = first_res.source_model, first_res.ref_key
    for (form, residuals) in delivered, r in residuals
        (r.source_model === group && r.ref_key == ref_key) ||
            error("$model_name delivers its $form injection through $(r.name) to " *
                  "$(r.source_model) via $(r.ref_key), not to $group via $ref_key like the " *
                  "rest of its injection")
    end
    return group, ref_key, Dict(form => [r.source_name for r in rs] for (form, rs) in delivered)
end

# The form an injection of `forms` is delivered in at a bus of `target_type`: its own if the bus
# implements that form's group fields, otherwise (voltage) the current form and then the power form,
# or (power, current) the other one; `nothing` if none
function _injection_form(forms::InjectionForms, fields, target_type::Type)
    mapping = group_res_mapping(target_type)
    order = forms.natural === :voltage ? (:voltage, :current, :power) :
            forms.natural === :power ? (:power, :current) : (:current, :power)
    k = findfirst(form -> all(f -> haskey(mapping, f), fields[form]), order)
    return isnothing(k) ? nothing : order[k]
end

# The forms an injection is delivered in, for messages
_describe_ports(fields) = join(("the $form port $(fields[form])" for form in keys(fields)), ", ")

"""
    injection_form(metadata::ModelMetadata, target_type::Type) -> Symbol

The form (`:power`, `:current` or `:voltage`) a model with an `InjectionForms` delivers its injection
in at a bus of type `target_type`: the form it declares if the bus implements that form's group
residuals, otherwise the next in the order of `InjectionForms`. Throws where the bus implements none.
"""
function injection_form(metadata::ModelMetadata, target_type::Type)
    isnothing(metadata.injection) && error("$(metadata.name) declares no injection")
    group, _, fields = _injection_ports(metadata)
    form = _injection_form(metadata.injection, fields, target_type)
    isnothing(form) &&
        error("$(metadata.name) injects at $target_type, which implements none of " *
              "$(_describe_ports(fields)) of group :$group")
    return form
end

"""
    injection_routes(metadata, models, foreign_refs) -> Union{Nothing, Vector{Symbol}}

The form (`:power`, `:current` or `:voltage`) each instance of a model with an `InjectionForms`
delivers its injection in, or `nothing` for a model without one.

A connection goes through exactly one port: the one of the form the model declares if the bus
it is connected to implements that form's group residuals, otherwise the next form in the order
of `InjectionForms`, and a bus that implements none is an error that names the instance and its bus.
"""
function injection_routes(
    metadata::ModelMetadata,
    models::NamedTuple,
    foreign_refs::Union{AbstractDict, Nothing}
)
    isnothing(metadata.injection) && return nothing
    model_name = metadata.name
    group, ref_key, fields = _injection_ports(metadata)
    link_ref = get_group_ref(foreign_refs, model_name, ref_key)
    n = length(link_ref.target_type)
    routes = Vector{Symbol}(undef, n)
    for i in 1:n
        target_type = get_base_model_type(models[link_ref.target_type[i]])
        form = _injection_form(metadata.injection, fields, target_type)
        isnothing(form) &&
            error("$model_name instance $i (key $(instance_key(models[model_name], ref_key, i))) " *
                  "injects at $(link_ref.target_type[i]) (via $ref_key), which implements " *
                  "none of $(_describe_ports(fields)) of group :$group")
        routes[i] = form
    end
    return routes
end

"""
Retrieve addresses for external variables and equations for a model.

Uses unified ref_key-based lookup into `foreign_refs`. All terminal references
(both direct model references and group references) are keyed by their
`normalize_ref_key(terminal.ref_key)` value.

# Reference Types
1. **Homogeneous references**: All instances target same model type (e.g., PQ → Bus)
   - GroupLinkRef with uniform target_type
   - Fast path: direct index into source addresses
2. **Heterogeneous references**: Instances target different model types (e.g., TGOV1 → GENCLS/GENROU)
   - GroupLinkRef with varying target_type
   - Requires per-instance lookup with group name resolution

# Arguments
- `am::AddressManager` - Address manager with allocated addresses
- `metadata::ModelMetadata` - Metadata for the model requesting addresses
- `models::NamedTuple` - All models in the system
- `foreign_refs` - Group references from IndexOrder (required for external vars)
"""
function retrieve_addresses(
    am::AddressManager,
    metadata::ModelMetadata,
    models::NamedTuple;
    foreign_refs::Union{AbstractDict, Nothing} = nothing
)
    model_name = metadata.name
    retrieved = ModelRetrievedAddresses(model_name)

    # Each connection of a model with two forms of injection goes through one of them
    routes = injection_routes(metadata, models, foreign_refs)
    injected = isnothing(metadata.injection) ? Symbol[] :
        vcat(metadata.injection.power, metadata.injection.current, metadata.injection.voltage)

    # Handle both external variables and residuals
    for external in [filter(is_external, metadata.vars);
                    filter(is_external, metadata.residuals)]
        (; name, address_type, source_model, source_name, ref_key) = external

        # UNIFIED: Always lookup by ref_key
        if !isnothing(foreign_refs) && has_group_ref(foreign_refs, model_name, ref_key)
            link_ref = get_group_ref(foreign_refs, model_name, ref_key)

            if name in injected
                form = name in metadata.injection.power ? :power :
                       name in metadata.injection.current ? :current : :voltage
                retrieved.by_type[address_type][name] = _gather_heterogeneous_addresses(
                    am, models, model_name, link_ref, address_type, source_model, source_name;
                    routed = routes .=== form
                )
            elseif is_heterogeneous(link_ref)
                # Heterogeneous: multiple target model types (group reference)
                # source_model is the group name, needed for canonical name resolution
                retrieved.by_type[address_type][name] = _gather_heterogeneous_addresses(
                    am, models, model_name, link_ref, address_type, source_model, source_name
                )
            else
                # Homogeneous: all same target model type
                # Use fast path with direct index lookup from GroupLinkRef
                target_model = first(link_ref.target_type)
                model_type = get_base_model_type(models[target_model])
                actual_name = _resolve_canonical_name(source_model, source_name, model_type, address_type)

                key = (address_type, target_model, actual_name)
                haskey(am.addresses, key) || error(
                    "Address key $key not found for $model_name.$name")
                src_addrs = am.addresses[key]
                retrieved.by_type[address_type][name] = src_addrs[link_ref.target_idx]
            end
        else
            error("External variable $model_name.$name requires foreign_refs entry for ref_key :$ref_key. " *
                  "Ensure build_terminal_refs was called during system setup.")
        end
    end

    am.retrieved[model_name] = retrieved
    nothing
end

# ============================================================================
# Reverse Address Lookup
# ============================================================================

"""
    build_reverse_map!(am::AddressManager)

Build reverse address lookup map from assigned and retrieved addresses.

Populates `am.reverse_map` with mappings from DAE vector addresses to
(model_name, var_name, instance_index) tuples, organized by AddressableType.

This is called automatically during `setup_with_formulation!`. Can also be
called manually after address remapping.

# Modifies
- `am.reverse_map`: Populated with reverse lookup data

# Example
```julia
am = AddressManager()
setup_with_formulation!(am, sys, PowerInjection())
# am.reverse_map is now populated
model, var, inst = am.reverse_map.algeb_vars[UInt(42)]
```
"""
function build_reverse_map!(am::AddressManager)
    new_dict() = Dict{UInt, Tuple{Symbol, Symbol, Int}}()
    reverse_map = ReverseAddressMap(new_dict(), new_dict(), new_dict(), new_dict(), new_dict())

    # Process assignments (internal variables/residuals that this model OWNS)
    # IMPORTANT: Only process assignments, NOT retrieved!
    # Retrieved addresses are external references (including SharedRes) that
    # belong to other models. Only the owning model should be in the reverse map.
    for (model_name, assignment) in pairs(am.assignments)
        for (addr_type, vars_dict) in pairs(assignment.by_type)
            target_dict = _reverse_dict(reverse_map, addr_type)
            for (var_name, addrs) in pairs(vars_dict)
                for (instance_idx, addr) in enumerate(addrs)
                    target_dict[addr] = (model_name, var_name, instance_idx)
                end
            end
        end
    end

    am.reverse_map = reverse_map

    return nothing
end

# Reverse-map dict holding addresses of the given address type
_reverse_dict(rm::ReverseAddressMap, ::AlgebVar) = rm.algeb_vars
_reverse_dict(rm::ReverseAddressMap, ::StateVar) = rm.state_vars
_reverse_dict(rm::ReverseAddressMap, ::AlgebRes) = rm.algeb_res
_reverse_dict(rm::ReverseAddressMap, ::StateRes) = rm.state_res
_reverse_dict(rm::ReverseAddressMap, ::ObservableVar) = rm.observable_vars

"""
    lookup_address(am::AddressManager, addr_type::AddressableType, addr::UInt) -> (Symbol, Symbol, Int)

Look up metadata for a DAE vector address.

Returns `(model_name, var_name, instance_index)`. Example:
`lookup_address(am, AlgebRes(), UInt(42))` → `(:PQ, :p, 3)` (3rd PQ's active power).

Requires `build_reverse_map!(am)` to be called first.
"""
function lookup_address(
    am::AddressManager,
    addr_type::AddressableType,
    addr::UInt
)::Tuple{Symbol, Symbol, Int}
    if isnothing(am.reverse_map)
        error("Reverse address map not built. Call build_reverse_map!(am) first.")
    end

    target_dict = _reverse_dict(am.reverse_map, addr_type)

    if !haskey(target_dict, addr)
        error("Address $addr not found in reverse map for type $addr_type")
    end

    return target_dict[addr]
end

"""
    variable_name(am::AddressManager, addr_type::AddressableType, addr::UInt) -> Symbol

Return a Symbol naming the variable at the given address.

Format: `:Model_varname_instance` (e.g., `:GENCLS_omega_2`).
"""
function variable_name(
    am::AddressManager,
    addr_type::AddressableType,
    addr::UInt
)::Symbol
    model, var, inst = lookup_address(am, addr_type, addr)
    return Symbol(model, "_", var, "_", inst)
end

"""
    item_addresses(am, model_name, item) -> Vector{UInt}

Addresses of a metadata variable or residual `item` of `model_name`: the model's
own assignments when `item` is internal, its retrieved addresses otherwise.
"""
function item_addresses(am::AddressManager, model_name::Symbol, item)
    by_type = is_internal(item) ? am.assignments[model_name].by_type : am.retrieved[model_name].by_type
    return by_type[item.address_type][item.name]
end

"""
    get_addresses(am, model_name, name) -> Vector{UInt}

Get DAE addresses for a variable/residual by searching all address types.

Searches both assigned (internal) and retrieved (external) addresses.
For better performance when type is known, use the 4-argument version.
"""
function get_addresses(
    am::AddressManager,
    model_name::Symbol,
    name::Symbol
)::Vector{UInt}
    # Search in assigned addresses (internal vars/residuals)
    if haskey(am.assignments, model_name)
        for (_, vars) in am.assignments[model_name].by_type
            if haskey(vars, name)
                return vars[name]
            end
        end
    end

    # Search in retrieved addresses (external vars/residuals)
    if haskey(am.retrieved, model_name)
        for (_, vars) in am.retrieved[model_name].by_type
            if haskey(vars, name)
                return vars[name]
            end
        end
    end

    error("Variable or residual '$name' not found for model $model_name")
end

"""
    get_addresses(am, addr_type, model_name, name) -> Vector{UInt}

Get DAE addresses for a variable/residual with explicit address type (O(1) lookup).

# Arguments
- `addr_type`: AlgebVar(), StateVar(), AlgebRes(), or StateRes()
- `model_name`: Model symbol (e.g., :Bus, :PQ)
- `name`: Variable/residual name (e.g., :theta, :v)
"""
function get_addresses(
    am::AddressManager,
    addr_type::AddressableType,
    model_name::Symbol,
    name::Symbol
)::Vector{UInt}
    # Try direct lookup in main addresses dictionary
    key = (addr_type, model_name, name)
    if haskey(am.addresses, key)
        return am.addresses[key]
    end

    # Check retrieved addresses (external vars/residuals)
    if haskey(am.retrieved, model_name) &&
       haskey(am.retrieved[model_name].by_type, addr_type) &&
       haskey(am.retrieved[model_name].by_type[addr_type], name)
        return am.retrieved[model_name].by_type[addr_type][name]
    end

    error("Variable or residual '$name' of type $addr_type not found for model $model_name")
end

"""
    get_var_address_type(am::AddressManager, model_name::Symbol, var_name::Symbol) -> AddressableType

Return the address type (`StateVar()` or `AlgebVar()`) for a variable.

Example: `get_var_address_type(am, :GENROU, :delta)` → `StateVar()`.
Throws if variable not found.
"""
function get_var_address_type(
    am::AddressManager,
    model_name::Symbol,
    var_name::Symbol
)::AddressableType
    # Check assigned addresses (internal variables)
    if haskey(am.assignments, model_name)
        by_type = am.assignments[model_name].by_type
        haskey(get(by_type, StateVar(), Dict{Symbol,Vector{UInt}}()), var_name) && return StateVar()
        haskey(get(by_type, AlgebVar(), Dict{Symbol,Vector{UInt}}()), var_name) && return AlgebVar()
    end
    # Check retrieved addresses (external variables)
    if haskey(am.retrieved, model_name)
        by_type = am.retrieved[model_name].by_type
        haskey(get(by_type, StateVar(), Dict{Symbol,Vector{UInt}}()), var_name) && return StateVar()
        haskey(get(by_type, AlgebVar(), Dict{Symbol,Vector{UInt}}()), var_name) && return AlgebVar()
    end
    available = Symbol[]
    for table in (am.assignments, am.retrieved)
        haskey(table, model_name) || continue
        for addr_type in (StateVar(), AlgebVar())
            append!(available, keys(get(table[model_name].by_type, addr_type, Dict{Symbol,Vector{UInt}}())))
        end
    end
    error("Variable '$var_name' not found for model $model_name. Variables: $(join(sort!(unique!(available)), ", "))")
end


"""
    remap_addresses!(am::AddressManager, perm::Vector{Int})

Remap the algebraic addresses in-place using the permutation `perm` (e.g. an AMD
ordering of the algebraic Jacobian ∂g/∂y): new address `k` holds old address
`perm[k]`.

`perm` permutes the algebraic space only. `AlgebVar` (the `y` vector) and
`AlgebRes` (the `g` vector) are square and receive the same permutation (a
symmetric reordering of ∂g/∂y). `StateVar`, `StateRes` and `ObservableVar` are
separate address spaces (the `x`, `f` and observable vectors) and keep their
addresses; `perm` is not defined on them.
"""
function remap_addresses!(am::AddressManager, perm::Vector{Int})
    n_alg = max_address_in_space(am, AlgebVar())
    n_res = max_address_in_space(am, AlgebRes())
    n_alg == n_res || error("Cannot permute a non-square algebraic space: " *
        "$n_alg algebraic variables, $n_res algebraic residuals")
    length(perm) == n_alg || error("Address permutation has length $(length(perm)), " *
        "algebraic space has $n_alg addresses")

    # New address of old address `a` (perm[new] = old)
    inv_perm = invperm(perm)
    remap(addrs) = UInt[a == UNROUTED ? a : UInt(inv_perm[Int(a)]) for a in addrs]
    is_algebraic(addr_type) = addr_type isa AlgebVar || addr_type isa AlgebRes

    # Remap main addresses dictionary in-place
    for (key, old_addrs) in am.addresses
        is_algebraic(key[1]) && (am.addresses[key] = remap(old_addrs))
    end

    # Remap assignments (internal) and retrieved (external) addresses in-place
    for table in (am.assignments, am.retrieved), (_, entry) in table
        for (addr_type, vars_dict) in entry.by_type
            is_algebraic(addr_type) || continue
            for (var_name, old_addrs) in vars_dict
                vars_dict[var_name] = remap(old_addrs)
            end
        end
    end

    # A permutation keeps the size of the algebraic space, so next_idx is unchanged.
    # Invalidate reverse map since addresses have changed
    # User must call build_reverse_map! again if needed
    Base.setfield!(am, :reverse_map, nothing)

    return nothing
end

"""
    setup_with_formulation!(am::AddressManager, sys::SystemModel, formulation::FormulationTrait)

Assign addresses for all models based on the formulation.

Assigns internal addresses, retrieves external addresses, and applies addressing strategy.
Called automatically by problem constructors (ACPowerFlowProblem, etc.).
"""
function setup_with_formulation!(am::AddressManager, sys::SystemModel, formulation::FormulationTrait)
    # Phase 1: Assign internal addresses using formulation-specific metadata
    for (model_name, model_vec) in pairs(sys.models)
        n_instances = length(model_vec)

        # Skip models with no instances
        if n_instances == 0
            @debug "Model $model_name has 0 instances, skipping address assignment"
            continue
        end

        form_metadata = get_metadata_for_formulation(
            model_metadata(get_base_model_type(model_vec)), formulation)

        assign_address!(am, model_name, form_metadata, n_instances)
    end

    # Phase 2: Retrieve external addresses
    for (_, model_vec) in pairs(sys.models)
        # Skip if no instances
        if length(model_vec) == 0
            continue
        end

        form_metadata = get_metadata_for_formulation(
            model_metadata(get_base_model_type(model_vec)), formulation)

        # Pass group references for group reference resolution
        frgn_refs = isempty(sys.order.foreign_refs) ? nothing : sys.order.foreign_refs
        retrieve_addresses(am, form_metadata, sys.models; foreign_refs=frgn_refs)
    end

    # Phase 3: Build reverse address lookup map
    # Note: AMD/CustomPermutation reordering is applied at the solver layer
    # via apply_addressing_strategy!() which has access to jacobian_sparsity_pattern
    build_reverse_map!(am)

    return am
end
