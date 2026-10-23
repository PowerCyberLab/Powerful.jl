"""
TDS Linking System - Connect dynamic models to power flow models

This module handles the linking between dynamic models and their replaced
power flow models (PV/Slack), as well as controller linking (TGOV1 → generators).
"""

# Use abstract interface to decouple from PowerfulPF (enables parallel precompilation)
using PowerfulCore: AbstractSteadyStateSolution

# ============================================================================
# Linking System - Connect dynamic models to power flow models
# ============================================================================
# Note: ModelLinks struct is defined in foundation/core_types.jl


"""
    create_model_links(sys_pf::SystemModel, dynamic_models::NamedTuple) -> ModelLinks

Build replacement links between dynamic models and power flow models.

Used for **replacement filtering**: determines which PV/Slack instances to remove
from the merged TDS system (replaced by GENCLS, GENROU, etc.).

**Note**: Controller links (TGOV1 → GENCLS) are handled via group references
built from controller pseudo-terminals and stored in `IndexOrder.foreign_refs`.

**Important**: Links are NOT required for initialization of bus-attached values.
Models can initialize V, θ via `FromSolution` using their Terminal declarations
(see `get_terminals` trait).

# Algorithm
For replacement models (GENCLS, GENROU), match to power flow models (PV, Slack)
using `get_replacement_map()` trait (see `replacement_links`).

# Returns
`ModelLinks` struct-of-arrays with (source_type, source_idx, target_type, target_idx)

# Example
```julia
links = create_model_links(sys_pf, (GENCLS=gencls_models,))
# Returns links: GENCLS[1]→PV[3], GENCLS[2]→Slack[1], etc.
```

# See Also
- `get_replacement_map`: Declares what static models a dynamic model replaces
"""
function create_model_links(sys_pf::SystemModel, dynamic_models::NamedTuple)
    links = [replacement_links(dyn_name, dyn_models, get_replacement_map(eltype(dyn_models)), sys_pf)
             for (dyn_name, dyn_models) in pairs(dynamic_models)]
    return ModelLinks(reduce(vcat, (l.source_type for l in links); init=Symbol[]),
                      reduce(vcat, (l.source_idx for l in links); init=Int[]),
                      reduce(vcat, (l.target_type for l in links); init=Symbol[]),
                      reduce(vcat, (l.target_idx for l in links); init=Int[]))
end

"""
    _as_field_value(v, T) -> T

Source value converted to a target model field's element type for matching;
string fields compare by text (e.g. a numeric DYR ID 1 matches "1").
"""
_as_field_value(v, ::Type{T}) where {T<:AbstractString} = convert(T, string(v))
_as_field_value(v, ::Type{T}) where {T} = convert(T, v)

"""
    replacement_links(name, instances, replacement_map, sys_pf) -> ModelLinks

Link each of `instances` (dynamic model instances or raw DYR records, reported as
`name`) to the static instance it replaces. Targets are tried in `replacement_map`
order; an instance links to the first target whose `target_fields` equal its
`source_fields` (converted with `_as_field_value`). Targets missing from `sys_pf` are
skipped, and each instance that links to nothing is warned about.
"""
function replacement_links(name::Symbol, instances, replacement_map, sys_pf::SystemModel)
    matches = Vector{Union{Nothing, Tuple{Symbol, Int}}}(nothing, length(instances))
    isempty(replacement_map) && return ModelLinks(Symbol[], Int[], Symbol[], Int[])

    for (target, spec) in replacement_map
        haskey(sys_pf.models, target) || continue
        target_models = sys_pf.models[target]

        # One index over the target fields, then a lookup per source instance
        target_index = build_index_map(target_models, collect(spec.target_fields))
        field_types = Tuple(eltype(getproperty(target_models, f)) for f in spec.target_fields)

        for (k, inst) in enumerate(instances)
            isnothing(matches[k]) || continue
            key = Tuple(_as_field_value(getproperty(inst, sf), T)
                        for (sf, T) in zip(spec.source_fields, field_types))
            idx = get(target_index, key, nothing)
            isnothing(idx) || (matches[k] = (target, idx))
        end
    end

    source_fields = first(values(replacement_map)).source_fields
    for (k, m) in enumerate(matches)
        isnothing(m) || continue
        values_k = map(f -> getproperty(instances[k], f), source_fields)
        @warn "No replacement target for $name[$k] with $source_fields = $values_k" maxlog=10
    end

    # Target-major order: links (and the models parsed from them) group by target
    linked = [k for target in keys(replacement_map)
                for k in eachindex(matches) if !isnothing(matches[k]) && matches[k][1] == target]
    return ModelLinks(fill(name, length(linked)), linked,
                      Symbol[matches[k][1] for k in linked], Int[matches[k][2] for k in linked])
end

"""
    build_replacement_masks(sys_pf::SystemModel, links::ModelLinks) -> Dict{Symbol, BitVector}

Build masks indicating which instances are replaced by dynamic models.

Returns a dictionary mapping model names to BitVectors where `true` indicates
the instance at that index is replaced.

# Example
```julia
links = create_model_links(sys_pf, dynamic_models)
masks = build_replacement_masks(sys_pf, links)
# masks[:PV] = BitVector([0,0,0,0,1,0,0,1,0])  # PV[5] and PV[8] replaced
```
"""
function build_replacement_masks(sys_pf::SystemModel, links::ModelLinks)
    # Initialize masks (all false = not replaced)
    masks = Dict{Symbol, BitVector}()
    for (name, models) in pairs(sys_pf.models)
        masks[name] = falses(length(models))
    end

    # Mark replaced instances (only for links that target power flow models)
    for i in 1:length(links)
        target = links.target_type[i]
        # Skip controller links (TGOV1 → GENCLS) - they don't replace PF models
        haskey(masks, target) || continue
        masks[target][links.target_idx[i]] = true
    end

    return masks
end


# ============================================================================
# System Construction - Build TDS-ready SystemModel
# ============================================================================

"""
    prepare_tds_system(sys_pf, sol_pf, dynamic_models; config) -> Tuple{SystemModel, ModelLinks}

Create TDS-ready system with dynamic models, filtering replaced static models.

# Example
```julia
# Staged API (for debugging)
sys_tds, links = prepare_tds_system(sys_pf, sol_pf, dynamic_models)
init = TDSInitialization(sys_tds, links, sol_pf)
prob = TDSProblem(init; tspan=(0.0, 10.0))

# Standard API (most users) - use setup_tds instead
init = setup_tds(sol_pf, dynamic_models)
prob = TDSProblem(init; tspan=(0.0, 10.0))
```

# Algorithm
1. Create links between dynamic and power flow models
2. Build replacement masks from links
3. `derive_system`: remove the replaced static instances, add the dynamic models
   (those per unit on the base of their linked device converted to the system base),
   re-index, and retrieve external params for the dynamic models
4. Apply ZIP conversion if enabled, to the loads of the TDS system (copies: `sys_pf`
   keeps its loads)

# Returns
- `sys_tds::SystemModel` - TDS system with filtered static + dynamic models
- `links::ModelLinks` - Mapping from dynamic to static (needed for initialization)

# Design Note
This function separates system construction from initialization:
- Construction: Parameter setup (this function)
- Initialization: State values (initialize_dynamic_models)
"""
function prepare_tds_system(
    sys_pf::SystemModel,
    sol_pf::AbstractSteadyStateSolution,
    dynamic_models::NamedTuple;
    config::SimConfig = SimConfig()
)
    # Step 1: Create links
    links = create_model_links(sys_pf, dynamic_models)

    # Step 2: Build replacement masks (compute once, use multiple times)
    masks = build_replacement_masks(sys_pf, links)

    # Step 3: Remove the replaced static instances, add the dynamic models, re-index
    # and retrieve the dynamic models' external params against the merged set
    keep = Dict(name => .!mask for (name, mask) in masks)
    sys_tds = derive_system(sys_pf, keep, dynamic_models; links)

    # Step 4: ZIP conversion of the TDS system's own loads
    tds_config = get_solver_config(config, :TDS, TDSConfig)
    tds_config.convert_pq_to_zip && apply_zip_conversion!(sys_tds, sol_pf, tds_config)

    return sys_tds, links
end

export create_model_links, prepare_tds_system
export build_replacement_masks
