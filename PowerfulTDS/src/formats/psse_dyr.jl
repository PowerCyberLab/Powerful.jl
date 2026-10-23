# PSSE DYR (dynamic model) file parser.
#
# Uses PowerDynData.jl to parse DYR files and converts the output to Powerful.jl
# dynamic model types with proper linking to power flow models.

"""
    create_model_links_for_type(model_type, dyr_data, sys_pf) -> ModelLinks

Create links between DYR instances and power flow models for a specific model type.

This is similar to `create_model_links()` but works with raw DYR data (StructVector)
before model construction, enabling us to inherit data during parsing.

# Arguments
- `model_type::Type` - Dynamic model type (e.g., GENCLS)
- `dyr_data::StructVector` - Raw DYR data from PowerDynData
- `sys_pf::SystemModel` - Power flow system

# Returns
`ModelLinks` mapping each DYR instance to its PV/Slack counterpart

Matching is done by `replacement_links` with the model type's replacement map.
"""
function create_model_links_for_type(
    model_type::Type,
    dyr_data::StructVector,
    sys_pf::SystemModel
)
    return replacement_links(nameof(model_type), dyr_data, get_replacement_map(model_type), sys_pf)
end

# ============================================================================
# Top-Level DYR Parser
# ============================================================================

"""
    parse_dynamics(dyr_file, sys_pf) -> NamedTuple

Parse PSSE DYR file and return dynamic model parameters.

Two-phase parsing:
1. **Phase 1**: Parse models that REPLACE power flow models (generators)
   - Links to PV/Slack in sys_pf via `get_replacement_map`
2. **Phase 2**: Parse models that CONTROL dynamic models (governors, exciters)
   - Links to Phase 1 results via `get_controlled_group`

# Arguments
- `dyr_file::String` - Path to PSSE DYR file
- `sys_pf::SystemModel` - Power flow system (must contain PV/Slack models)

# Returns
NamedTuple with parsed dynamic models (GENCLS, TGOV1, etc.)

# Example
```julia
dyn = parse_dynamics("case.dyr", sys)
prob = TDSProblem(setup_tds(sol, dyn); tspan=(0.0, 10.0))
```
"""
function parse_dynamics(
    dyr_file::String,
    sys_pf::SystemModel
)
    dyn_data = PowerDynData.parse_dyr(dyr_file)

    # Show parsed model summary
    if !isempty(dyn_data.models)
        sorted_models = sort(collect(dyn_data.models), by = x -> x[1])
        model_list = join(["$(name): $(length(data.data))" for (name, data) in sorted_models], ", ")
        @debug "DYR file contains: $model_list"
    end

    results = Dict{Symbol, Any}()

    # Track supported DYR model names for unsupported model detection
    supported_dyr_names = Set{String}()
    for model_type in all_models()
        dyr_name = get_dyr_model_name(model_type)
        !isnothing(dyr_name) && push!(supported_dyr_names, dyr_name)
    end

    for model_type in all_models()
        dyr_name = get_dyr_model_name(model_type)
        isnothing(dyr_name) && continue
        !haskey(dyn_data.models, dyr_name) && continue

        dyr_data = dyn_data.models[dyr_name].data
        isempty(dyr_data) && continue

        # Models with replacement maps need linking to power flow models
        if !isempty(get_replacement_map(model_type))
            links = create_model_links_for_type(model_type, dyr_data, sys_pf)
            if isempty(links)
                @warn "No $dyr_name could be linked to power flow models" maxlog=1
                continue
            end
            models = parse_dynamic_model(model_type, dyr_data, sys_pf, links)
            inherit_replaced_status!(models, links, sys_pf)
            results[nameof(model_type)] = models
        else
            # Controllers and all other models: no linking needed at parse time
            results[nameof(model_type)] = parse_dynamic_model(model_type, dyr_data, sys_pf, nothing)
        end
    end

    if isempty(results)
        @warn "No dynamic models found in DYR file: $dyr_file"
    end

    # =========================================================================
    # Report unsupported models
    # =========================================================================
    unsupported = Dict{String, Int}()
    for (model_name, model_data) in dyn_data.models
        if !(model_name in supported_dyr_names)
            unsupported[model_name] = length(model_data.data)
        end
    end

    if !isempty(unsupported)
        # Sort by count (descending) for readability
        sorted_models = sort(collect(unsupported), by = x -> -x[2])
        model_list = join(["$name ($count)" for (name, count) in sorted_models], ", ")
        @warn "Unsupported DYR models ignored: $model_list"
    end

    return NamedTuple(Pair{Symbol, Any}[k => v for (k, v) in results])
end

"""
    inherit_replaced_status!(models, links::ModelLinks, sys_pf) -> models

Take each replacing dynamic instance out of service when the static instance it
replaces (per `links`, one link per instance in order) is out of service: a unit
that is off in the power flow is off in the simulation. Errors if such a model has
no status field to express it.
"""
function inherit_replaced_status!(models, links::ModelLinks, sys_pf::SystemModel)
    length(models) == length(links) ||
        error("$(nameof(eltype(models))): $(length(models)) instances for $(length(links)) links")
    field = get_status_field(models)
    for k in eachindex(models)
        in_service(sys_pf.models[links.target_type[k]])[links.target_idx[k]] && continue
        isnothing(field) &&
            error("$(nameof(eltype(models)))[$k] replaces out-of-service " *
                  "$(links.target_type[k])[$(links.target_idx[k])] but has no status field")
        getproperty(models, field)[k] = false
    end
    return models
end

export parse_dynamics, create_model_links_for_type
