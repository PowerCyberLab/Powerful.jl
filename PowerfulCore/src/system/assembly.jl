"""
System assembly: deriving systems by removing or adding model instances.

All paths that change which instances a system contains go through here, so the
IndexOrder and external parameters are always rebuilt the same way (`ext2int` +
`retrieve_external_params`) instead of being patched by hand.
"""

"""
    filter_instances(data::NamedTuple, keep::AbstractDict; drop_empty=false) -> NamedTuple

Keep only the instances where `keep[name]` is true, for each entry of `data`
(model StructArrays such as `sys.models`, or per-instance parameter collections
such as `sys.external` and `sys.pu_coeffs`). Entries without a mask are kept
unchanged. With `drop_empty=true`, masked entries left with no instances are removed.

# Example
```julia
keep = Dict(:PV => .!replaced)
filtered_models = filter_instances(sys.models, keep)
filtered_external = filter_instances(sys.external, keep)
```
"""
function filter_instances(data::NamedTuple, keep::AbstractDict; drop_empty::Bool = false)
    filtered = Pair{Symbol, Any}[]
    for (name, value) in pairs(data)
        if haskey(keep, name)
            mask = keep[name]
            drop_empty && !any(mask) && continue
            value = _keep_instances(value, mask)
        end
        push!(filtered, name => value)
    end
    return NamedTuple(filtered)
end

_keep_instances(arr::AbstractArray, mask) = arr[mask]
_keep_instances(nt::NamedTuple, mask) = map(v -> v[mask], nt)
_keep_instances(d::AbstractDict, mask) = typeof(d)(k => v[mask] for (k, v) in d)

"""
    filter_instances(data::ParsedSystemData, keep::AbstractDict; drop_empty=false) -> ParsedSystemData

The original physical-unit data restricted to the kept instances, so
`parsed_data.models` stays index-aligned with `sys.models`. The system base and
metadata are unchanged.
"""
filter_instances(data::ParsedSystemData, keep::AbstractDict; drop_empty::Bool = false) =
    ParsedSystemData(filter_instances(data.models, keep; drop_empty), data.system_base, data.metadata)

# Optional system data (`sys.parsed_data`, `sys.pu_coeffs`) that is absent stays absent
filter_instances(::Nothing, ::AbstractDict; drop_empty::Bool = false) = nothing

"""
    add_instances!(sys::SystemModel, additions::NamedTuple) -> SystemModel

Append model instances to `sys`, e.g. `add_instances!(sys, (Fault = faults,))`.
Instances must already be in per-unit. Models not yet in the system are added.

The system is re-indexed with `ext2int` (so group references to the added models
resolve everywhere), external parameters are retrieved for the extended models,
and the cached topology is invalidated. Errors if an added instance key already
exists, or appears twice in `additions`; `sys` is left unchanged on error.
"""
function add_instances!(sys::SystemModel, additions::NamedTuple)
    extended = Pair{Symbol, Any}[]
    for (name, new) in pairs(additions)
        isempty(new) && continue
        meta = model_metadata(get_base_model_type(new))
        declares_linked_rating(meta) && throw(ArgumentError(
            "$name is per unit on the base of the device it is linked to, which `add_instances!` does " *
            "not convert (it takes per-unit values on the system base); `derive_system` does, with the links"))
        isempty(meta.derived_params) || throw(ArgumentError(
            "$name declares derived parameters (the ratio of its rating to the system base), which come with " *
            "the per-unit conversion that `add_instances!` does not do (it takes per-unit values on the system base)"))
        key_fields = get_key_name(get_base_model_type(new))
        new_keys = [instance_key(new, key_fields, i) for i in 1:length(new)]
        allunique(new_keys) || error("Duplicate $name keys among added instances: $new_keys")

        if haskey(sys.models, name)
            for k in new_keys
                find_index(sys, name, k) === nothing || error("$name instance $k already exists")
            end
            push!(extended, name => vcat(sys.models[name], new))
        else
            push!(extended, name => new)
        end
    end
    extended_models = NamedTuple(extended)

    # Build everything before mutating, so a failure leaves sys untouched
    models = merge(sys.models, extended_models)
    _, order = ext2int(models)
    external = retrieve_external_params(extended_models, models, order; pu_coeffs = sys.pu_coeffs)

    sys.models = models
    sys.order = order
    sys.external = merge(sys.external, external)
    sys.topology = nothing
    return sys
end

"""
    derive_system(sys::SystemModel, keep::AbstractDict, new_models::NamedTuple; drop_empty=false,
                  links=ModelLinks()) -> SystemModel

A new system derived from `sys`: keep the instances selected by `keep` (see
`filter_instances`, including `drop_empty`), add `new_models`, re-index the merged
model set with `ext2int`, and retrieve external parameters for `new_models` against
it. Kept instances keep their external parameters, original parsed data and
per-unit coefficients, filtered with the same masks, and the derived system keeps
the settings of `sys`. `sys` is not modified. The external parameters of `new_models`
retrieved before the per-unit conversion (`PrePerUnit`) hold the physical value, as in `sys`
(recovered from the per-unit models through the coefficients).

`new_models` are per unit on the system base, except those that declare a rating from the
device they are linked to (`ModelRating(...; from = LinkedDevice())`): their parameters are per
unit on that device's base, and are converted here, with `links` (`ModelLinks`) naming the devices
of `sys` (which `keep` may remove) they are linked to (`convert_linked_models`). Their
coefficients join the system's `pu_coeffs`.

Shared by every derived-system builder, such as the TDS system builder (static models
plus dynamic models) and any other that assembles a system from a different mix of
models.
"""
function derive_system(sys::SystemModel, keep::AbstractDict, new_models::NamedTuple; drop_empty::Bool = false,
                       links::ModelLinks = ModelLinks())
    # The derived system owns its data: changes to it (initialization, load conversion,
    # events) must not reach the caller's model arrays. Filtering copies the models with
    # a mask; the others and the added models are copied here.
    filtered_models = filter_instances(sys.models, keep; drop_empty)
    filtered_models = NamedTuple(name => haskey(keep, name) ? models : copy(models)
                                 for (name, models) in pairs(filtered_models))
    filtered_external = filter_instances(sys.external, keep; drop_empty)
    new_models = map(copy, new_models)
    # The models given per unit on the base of the device they are linked to (the device is in
    # `sys`, which `keep` may remove it from) are converted to the system base here
    linked = convert_linked_models(sys, new_models, links)
    new_models = linked.models
    all_models = merge(filtered_models, new_models)
    _, order = ext2int(all_models)
    pu_coeffs = filter_instances(sys.pu_coeffs, keep; drop_empty)
    isempty(linked.coeffs) || (pu_coeffs = merge(something(pu_coeffs, NamedTuple()), linked.coeffs))
    new_external = retrieve_external_params(new_models, all_models, order; pu_coeffs)
    return SystemModel(all_models, merge(filtered_external, new_external), order;
                       settings = sys.settings,
                       parsed_data = filter_instances(sys.parsed_data, keep; drop_empty),
                       pu_coeffs)
end

"""
    copy(sys::SystemModel) -> SystemModel

A system that owns its model arrays, so changes to them (events, added instances)
do not reach `sys`. The index order, external parameters, parsed data, per-unit
coefficients, cached topology and settings are shared: they are replaced, never
modified in place, when a system changes (see `add_instances!`).
"""
Base.copy(sys::SystemModel) =
    SystemModel(map(copy, sys.models), sys.external, sys.order;
                settings = sys.settings, array_type = sys.array_type,
                parsed_data = sys.parsed_data, pu_coeffs = sys.pu_coeffs,
                topology = sys.topology)
