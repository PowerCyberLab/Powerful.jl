"""
TDS Initialization System - Vectorized, metadata-driven initialization

This module provides vectorized data gathering and initialization for dynamic models.
Each model declares its initialization group and dependencies, and the system
automatically resolves initialization order.
"""

# Use abstract interface to decouple from PowerfulPF (enables parallel precompilation)
using PowerfulCore: AbstractSteadyStateSolution, get_system, get_address_manager, get_state_vector
using PowerfulCore: normalize_ref_key

# ============================================================================
# Initialization Context
# ============================================================================

"""
    InitContext

Bundles all data needed for dynamic model initialization into a single context.

This eliminates parameter sprawl in `initialize_dynamic_states` and provides
a clean interface for models to access:
- Power flow solution and system
- Model links (pre-filtered per model via helpers)
- Previously initialized model states (accumulated during initialization loop)

# Fields
- `pf_sol::AbstractSteadyStateSolution`: Completed steady-state solution (power flow)
- `sys::SystemModel`: TDS system with dynamic models
- `links::ModelLinks`: All model links
- `initialized::Dict{Symbol, NamedTuple}`: Accumulated initialized states

# Usage
```julia
# Driver creates context once
ctx = InitContext(pf_sol, sys, links)

# Each model uses helper to get its data
data = gather_init_data(ctx, GENCLS)  # Returns NamedTuple from init_sources
```
"""
mutable struct InitContext
    pf_sol::AbstractSteadyStateSolution
    sys::SystemModel
    links::ModelLinks
    initialized::Dict{Symbol, NamedTuple}
end

# Constructor with empty initialized dict
InitContext(pf_sol::AbstractSteadyStateSolution, sys::SystemModel, links::ModelLinks) =
    InitContext(pf_sol, sys, links, Dict{Symbol, NamedTuple}())

"""
    gather_init_data(ctx::InitContext, ::Type{M}) -> NamedTuple

Resolve all initialization inputs for a model type using its metadata declarations.

This is the main helper for model initialization - it reads the model's `init_inputs`
metadata and resolves each input based on its `kind` field.

# Arguments
- `ctx::InitContext`: Initialization context with pf_sol, links, and accumulated states
- `::Type{M}`: Model type (e.g., GENCLS, TGOV1)

# Returns
NamedTuple with fields matching the model's `init_inputs` declarations.

# Example
```julia
function initialize_dynamic_states(model::StructArray{GENCLS{T}}, ctx::InitContext) where T
    data = gather_init_data(ctx, GENCLS)
    # data.P_gen, data.Q_gen, data.V_bus, data.theta_bus available
    # ... model-specific computation
end
```
"""
function gather_init_data(ctx::InitContext, ::Type{M}) where M
    model_name = nameof(M)
    metadata = model_metadata(M)
    return resolve_init_inputs(ctx, model_name, metadata)
end

export InitContext, gather_init_data

# ============================================================================
# Link Filtering Helpers
# ============================================================================

"""
    extract_target_info(links::ModelLinks, source_model::Symbol)
        -> (target_models::Vector{Symbol}, target_indices::Vector{Int})

Extract target model types and indices for a specific source model using vectorized mask.

# Example
```julia
target_models, target_indices = extract_target_info(links, :GENCLS)
# → target_models = [:PV, :PV, :Slack, ...]
# → target_indices = [5, 8, 3, ...]
```
"""
function extract_target_info(links::ModelLinks, source_model::Symbol)
    mask = links.source_type .== source_model
    return links.target_type[mask], links.target_idx[mask]
end


# ============================================================================
# Initialization Input Resolution
# ============================================================================

"""
    resolve_init_inputs(ctx::InitContext, model_name::Symbol, metadata::ModelMetadata) -> NamedTuple

Resolve all initialization inputs declared in model metadata.

Each input's `kind` field determines the resolution strategy:
- `:from_solution`: Fetch from DAE solution vector (V, θ)
- `:from_observation`: Compute from linked power flow model (P, Q)
- `:from_group_init`: Get from previously initialized model in group

# Arguments
- `ctx::InitContext`: Initialization context with pf_sol, sys, links, and accumulated states
- `model_name::Symbol`: Name of the model being initialized (e.g., :GENCLS)
- `metadata::ModelMetadata`: Model metadata with init_inputs declarations

# Returns
NamedTuple with one field per InitInput entry, each containing a vector of values, and one per
derived parameter of the model (`RatingRatio`, as the ratio of its rating to the system base).

# Example
```julia
metadata = model_metadata(GENCLS)
data = resolve_init_inputs(ctx, :GENCLS, metadata)

# Returns: (P_gen=[2.0, 1.5, 3.5], Q_gen=[0.5, 0.3, 1.2], V_bus=[1.05, 1.03, 1.04], ...)
data.P_gen[1]  # 2.0 - Active power for GENCLS[1] from PV[5]
```
"""
function resolve_init_inputs(
    ctx::InitContext,
    model_name::Symbol,
    metadata::ModelMetadata
)
    result_arrays = []
    result_names = Symbol[]

    for input in metadata.init_inputs
        values = _resolve_input(ctx, input, model_name)
        push!(result_names, input.name)
        push!(result_arrays, values)
    end

    # The derived parameters of the model (`RatingRatio`) are those of its external parameters
    for derived in metadata.derived_params
        push!(result_names, derived.name)
        push!(result_arrays, getproperty(getproperty(ctx.sys.external, model_name), derived.name))
    end

    return NamedTuple{tuple(result_names...)}(tuple(result_arrays...))
end

# ============================================================================
# Input Resolution (dispatches based on input.kind)
# ============================================================================

"""
    _resolve_input(ctx::InitContext, input::InitInput, model_name::Symbol) -> Vector

Resolve a single initialization input based on its kind: power flow quantities as
`Float64`, values of previously initialized models in their element type.
"""
function _resolve_input(
    ctx::InitContext,
    input::InitInput,
    model_name::Symbol
)
    if input.kind === :from_solution
        return _resolve_from_solution(ctx, input, model_name)
    elseif input.kind === :from_observation
        return _resolve_from_observation(ctx, input, model_name)
    elseif input.kind === :from_group_init
        return _resolve_from_group_init(ctx, input, model_name)
    else
        error("Unknown InitInput kind: $(input.kind). " *
              "Expected :from_solution, :from_observation, or :from_group_init")
    end
end

# ============================================================================
# Resolution Implementations
# ============================================================================

"""
    _resolve_from_observation(ctx, input, model_name) -> Vector{Float64}

Resolve input from linked power flow model (computed quantities like P, Q).

Requires replacement links - extracts target models/indices from ctx.links internally.
Handles heterogeneous linked models (e.g., some GENCLS link to PV, others to Slack),
each with their own compute_observations implementation.
"""
function _resolve_from_observation(
    ctx::InitContext,
    input::InitInput,
    model_name::Symbol
)
    # Extract links for this model (replacement links)
    target_models, target_indices = extract_target_info(ctx.links, model_name)

    if isempty(target_models)
        error("Model $model_name has no replacement links but declares FromObservation " *
              "input :$(input.name). FromObservation requires links to replaced PF models.")
    end

    pf_sol = ctx.pf_sol
    n = length(target_models)
    result = Vector{Float64}(undef, n)

    # Group by target model type for efficient batch computation
    unique_models = unique(target_models)

    for linked_model in unique_models
        # Get model type (use abstract interface)
        model = get_system(pf_sol).models[linked_model]
        model_type = get_base_model_type(model)

        # Compute observations for all instances of this model type (vectorized)
        all_obs = compute_observations(pf_sol, model_type, input.field)

        # Extract values for instances that link to this model
        mask = target_models .== linked_model
        indices_in_model = target_indices[mask]
        result[mask] = all_obs[indices_in_model]
    end

    return result
end

"""
    _resolve_from_solution(ctx, input, model_name) -> Vector{Float64}

Resolve input from DAE solution variables (V, θ).

Uses the model's Terminal declarations to find bus attachment - NO LINKS REQUIRED.
This allows standalone dynamic models (that don't replace anything) to initialize.
"""
function _resolve_from_solution(
    ctx::InitContext,
    input::InitInput,
    model_name::Symbol
)
    pf_sol = ctx.pf_sol
    am = get_address_manager(pf_sol)  # Use abstract interface
    y = get_state_vector(pf_sol)      # Use abstract interface

    # Get model and its type
    model = ctx.sys.models[model_name]
    model_type = get_base_model_type(model)
    n = length(model)

    # Get Bus terminal from model's declaration (no hardcoding!)
    bus_terminal = find_terminal(model_type, :Bus)

    if bus_terminal === nothing
        error("Model $model_name has no Bus terminal declared but uses FromSolution input. " *
              "Add `get_terminals(::Type{<:$model_name}) = (Terminal(:Bus, :field),)` to the model.")
    end

    # IMPORTANT: For FromSolution, we must use the PF solution's indices to map
    # external bus numbers to internal Bus indices. We cannot use ctx.sys foreign_refs
    # because ctx.sys may connect GENCLS.BUS to another bus model of :bus_group (other
    # indices), but we need GENCLS.BUS → Bus indices in the PF solution.
    pf_sys = get_system(pf_sol)  # Use abstract interface
    bus_key_values = getproperty(model, bus_terminal.ref_key)

    # Get state variable addresses and extract values (vectorized)
    state_addrs = get_addresses(am, :Bus, input.field)
    result = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        # Convert external bus number to internal Bus index via PF solution's indices
        bus_idx = get_index_by_key(pf_sys, :Bus, bus_key_values[i])
        result[i] = y[state_addrs[bus_idx]]
    end

    return result
end

"""
    _resolve_from_group_init(ctx, input, model_name) -> Vector

Resolve input from any model in a group (polymorphic group reference), in the element
type of the linked models' initialized values (Dual when the initialization is
differentiated with respect to their parameters).

This enables controllers (TGOV1, SEXS) to work with ANY generator model
(GENCLS, GENROU, GENSAL) without hard-coding specific model names.
"""
function _resolve_from_group_init(
    ctx::InitContext,
    input::InitInput,
    model_name::Symbol
)
    group_refs = ctx.sys.order.foreign_refs
    model = ctx.sys.models[model_name]
    base_type = get_base_model_type(model)
    controlled_group = get_controlled_group(base_type)
    controlled_group === nothing && error(
        "Model $model_name reads :$(input.field) of group :$(input.from_model) at initialization (InitInput " *
        ":$(input.name)) but declares no group it controls: declare " *
        "PowerfulCore.get_controlled_group(::Type{<:$(nameof(base_type))}) = :$(input.from_model)")
    ctrl_term = find_control_terminal(base_type)
    ctrl_term === nothing && error("Model $model_name declares controlled group $controlled_group but has no Terminal for it")
    # UNIFIED: Lookup by ref_key (the field in this model), not by group name
    ctrl_ref_key = normalize_ref_key(ctrl_term.ref_key)

    if !has_group_ref(group_refs, model_name, ctrl_ref_key)
        error("Model $model_name has no group_refs entry for ref_key :$ctrl_ref_key but declares FromGroupInit " *
              "input :$(input.name). FromGroupInit requires control links to target models.")
    end

    group_ref = get_group_ref(group_refs, model_name, ctrl_ref_key)

    group_name = input.from_model  # This is the group name (e.g., :generators)
    canonical_field = input.field   # Canonical name in group interface (e.g., :P_m)
    # The initial value of a variable the group requires (what a member offers)
    required_vars = get_group_interface(group_name).required_vars
    any(v -> v.name === canonical_field, required_vars) || error(
        "Model $model_name reads :$canonical_field of group :$group_name at initialization (InitInput " *
        ":$(input.name)), which is not a variable the group requires: $([v.name for v in required_vars])")

    target_models = group_ref.target_type
    target_indices = group_ref.target_idx

    prior_init = isempty(ctx.initialized) ? NamedTuple() : NamedTuple(ctx.initialized)

    # The initialized values of the field in each linked model
    sources = Dict{Symbol, AbstractVector}()
    for linked_model in unique(target_models)
        # Check that this model has been initialized
        if !haskey(prior_init, linked_model)
            error("Model $linked_model not found in prior_init. " *
                  "Ensure it has been initialized before $model_name (check InitInput declarations).")
        end

        model_init = prior_init[linked_model]

        # Get the model type to resolve canonical name to actual field name
        model_type = get_base_model_type(ctx.sys.models[linked_model])

        # Resolve canonical group variable to the model's own variable name
        actual_field = resolve_group_var_name(group_name, canonical_field, model_type)

        # Check that the resolved field exists
        if !haskey(model_init, actual_field)
            error("Field :$actual_field (resolved from group :$group_name canonical :$canonical_field) " *
                  "not found in initialized states for $linked_model. Available fields: $(keys(model_init))")
        end

        sources[linked_model] = model_init[actual_field]
    end

    # One value per target instance, in the promoted element type of the sources
    T = isempty(sources) ? Float64 : mapreduce(eltype, promote_type, values(sources))
    result = Vector{T}(undef, length(target_models))
    for i in eachindex(target_models)
        result[i] = sources[target_models[i]][target_indices[i]]
    end

    return result
end

# ============================================================================
# Initialization Order Resolution
# ============================================================================

"""
    InitializationOrder

Topologically sorted initialization sequence.

# Fields
- `groups::Vector{Symbol}` - Groups in initialization order
- `models_by_group::Dict{Symbol, Vector{Symbol}}` - Model types per group
"""
struct InitializationOrder
    groups::Vector{Symbol}
    models_by_group::Dict{Symbol, Vector{Symbol}}
end

"""
    compute_initialization_order(sys::SystemModel) -> InitializationOrder

Compute initialization order from model metadata using topological sort.

Dependencies are inferred automatically from `InitInput` declarations:
- `InitInput(:from_group_init, ...)` → depends on the referenced group
- `InitInput(:from_observation, ...)` or `InitInput(:from_solution, ...)` → depends on `:power_flow`

# Algorithm
1. Collect all models and the groups they implement (`implements_group`)
2. Build model → group lookup for dependency resolution
3. Infer dependencies from InitInput metadata
4. Topological sort to determine group order
5. Validate no cycles exist (Kahn's algorithm detects cycles)

# Throws
- `ErrorException` if circular dependencies detected
"""
function compute_initialization_order(sys::SystemModel)
    # Collect all model types in system
    model_types = Symbol[]
    for model_name in keys(sys.models)
        model = sys.models[model_name]
        if !isempty(model)
            push!(model_types, model_name)
        end
    end

    # Build group membership map and model → group lookup
    models_by_group = Dict{Symbol, Vector{Symbol}}()
    model_to_group = Dict{Symbol, Symbol}()

    for model_name in model_types
        model_type = get_base_model_type(sys.models[model_name])
        # Initialization group: the group the model implements (models outside any
        # group initialize as :uncategorized)
        group = something(implements_group(model_type), :uncategorized)

        if !haskey(models_by_group, group)
            models_by_group[group] = Symbol[]
        end
        push!(models_by_group[group], model_name)
        model_to_group[model_name] = group
    end

    # Infer group dependencies from InitInput metadata
    group_deps = Dict{Symbol, Set{Symbol}}()
    for (group, models) in models_by_group
        deps = Set{Symbol}()
        for model_name in models
            model_type = get_base_model_type(sys.models[model_name])
            metadata = model_metadata(model_type)

            for input in metadata.init_inputs
                if input.kind == :from_group_init
                    # Depends directly on the referenced group (polymorphic group reference)
                    target_group = input.from_model  # from_model is actually the group name
                    push!(deps, target_group)
                elseif input.kind in (:from_observation, :from_solution)
                    # Depends on power_flow for PF data
                    push!(deps, :power_flow)
                end
            end
        end
        group_deps[group] = deps
    end

    @debug "Inferred group dependencies" group_deps

    # Topological sort (Kahn's algorithm detects cycles)
    sorted_groups = _topological_sort(group_deps)

    return InitializationOrder(sorted_groups, models_by_group)
end

"""
    _topological_sort(deps::Dict{Symbol, Set{Symbol}}) -> Vector{Symbol}

Kahn's algorithm for topological sort.

# Throws
- `ErrorException` if circular dependencies detected
"""
function _topological_sort(deps::Dict{Symbol, Set{Symbol}})
    # Build in-degree map
    in_degree = Dict{Symbol, Int}()
    # Collect all unique nodes (groups that have deps + groups that are depended upon)
    all_nodes = union(keys(deps), Iterators.flatten(values(deps)))

    for node in all_nodes
        in_degree[node] = 0
    end

    for (node, node_deps) in deps
        for dep in node_deps
            in_degree[node] += 1
        end
    end

    # Find nodes with no dependencies
    queue = Symbol[]
    for (node, degree) in in_degree
        if degree == 0
            push!(queue, node)
        end
    end

    # Process queue
    result = Symbol[]
    while !isempty(queue)
        node = popfirst!(queue)
        push!(result, node)

        # Reduce in-degree for dependent nodes
        for (other, other_deps) in deps
            if node in other_deps
                in_degree[other] -= 1
                if in_degree[other] == 0
                    push!(queue, other)
                end
            end
        end
    end

    # Check for cycles
    if length(result) != length(all_nodes)
        error("Circular dependency detected in initialization groups")
    end

    return result
end

# ============================================================================
# Model-Specific Initialization (Trait)
# ============================================================================

"""
    initialize_dynamic_states(model::StructArray{<:AbstractModel}, ctx::InitContext) -> NamedTuple

Initialize ALL instances of a dynamic model at once using the initialization context.

This is a trait function that dynamic models implement. The default implementation
errors, requiring models to provide their specific initialization logic.

# Arguments
- `model::StructArray` - All instances of the dynamic model (for type dispatch)
- `ctx::InitContext` - Context with pf_sol, sys, links, and accumulated initialized states

# Returns
NamedTuple of state arrays (one vector per state variable)

# Example Implementation
```julia
function initialize_dynamic_states(model::StructArray{GENCLS{T}}, ctx::InitContext) where T
    # Use helper to gather all declared init_sources
    data = gather_init_data(ctx, GENCLS)

    # Model-specific computation using data.P_gen, data.Q_gen, etc.
    # ...

    return (delta=delta, omega=omega, ...)
end
```
"""
function initialize_dynamic_states(
    model::StructArray{<:AbstractModel},
    ctx::InitContext
)
    error("initialize_dynamic_states not implemented for $(eltype(model))")
end

# ============================================================================
# System-Level Initialization Driver
# ============================================================================

"""
    check_initialization_declared(sys::SystemModel)

Throw an `ArgumentError` naming each model of `sys` that does not declare whether it requires
initialization (`PowerfulCore.requires_initialization`): the initialization cannot tell whether
to call its `initialize_dynamic_states`, and a model with states it skips would start at zero.
"""
function check_initialization_declared(sys::SystemModel)
    refused = [String(name) for (name, model) in pairs(sys.models)
               if !declares_requires_initialization(get_base_model_type(model))]
    isempty(refused) || throw(ArgumentError(
        "The time-domain initialization cannot tell whether these models require initialization: " *
        join(refused, ", ") * " (no requires_initialization declared: a method of " *
        "PowerfulCore.requires_initialization returning true when its initialize_dynamic_states gives its " *
        "states and setpoints, false when it has no state of its own or another initializer owns them)"))
    return nothing
end

"""
    initialize_dynamic_models(
        sys::SystemModel,
        links::ModelLinks,
        pf_sol::AbstractSteadyStateSolution
    ) -> NamedTuple

Initialize the dynamic models of `sys` from the power flow solution `pf_sol`, in
dependency order using vectorized gathering.

`setup_tds` calls it for the phasor system it builds. A package that builds its own
system around the TDS dynamic models (with bus models of its own that implement
`:bus_group`) calls it to initialize them. To simulate a phasor system, use:
```julia
prob = TDSProblem(setup_tds(sol_pf, dynamic_models); tspan=(0.0, 10.0))
```

# Algorithm
1. Compute initialization order via topological sort
2. For each group in order:
   - For each model in group:
     - Call model's initialize_dynamic_states() with prior_init from previously initialized models
     - Store results for use by subsequent models
3. Build setpoints from initialized states using model metadata
4. Extract PF snapshot (v, theta)
5. Return unified init NamedTuple

# Returns
NamedTuple with three fields:
- `dyn`: Per-model initialized states (GENCLS=(delta, omega, ...), ...)
- `setpoints`: Per-model setpoints for holding equations (GENCLS=(P_m0, E_q0), ...)
- `pf`: Power flow snapshot (v=[], theta=[])
"""
function initialize_dynamic_models(
    sys::SystemModel,
    links::ModelLinks,
    pf_sol::AbstractSteadyStateSolution
)
    check_initialization_declared(sys)

    # Compute initialization order
    order = compute_initialization_order(sys)

    # Create initialization context (accumulates states as we go)
    ctx = InitContext(pf_sol, sys, links)

    for group in order.groups
        models_in_group = get(order.models_by_group, group, Symbol[])

        for model_name in models_in_group
            model = sys.models[model_name]

            if isempty(model)
                continue
            end

            # Only initialize models that require it (checked via trait)
            model_type = get_base_model_type(model)
            if !requires_initialization(model_type)
                continue
            end

            # Initialize all instances at once (vectorized)
            # Model uses gather_init_data(ctx, ModelType) to get its data
            states = initialize_dynamic_states(model, ctx)

            # Accumulate for subsequent models that depend on this one
            ctx.initialized[model_name] = states
        end
    end

    # Build dyn NamedTuple
    dyn = isempty(ctx.initialized) ? NamedTuple() : NamedTuple(ctx.initialized)

    # Build setpoints from initialized states using model metadata
    setpoints = _build_setpoints(dyn, sys)

    # Extract PF snapshot (v, theta at buses)
    pf = _extract_pf_snapshot(pf_sol)

    return (dyn = dyn, setpoints = setpoints, pf = pf)
end

"""
    _build_setpoints(dyn::NamedTuple, sys::SystemModel) -> NamedTuple

Build setpoints for holding equations from initialized dynamic states.

For each model with `initialized_params` metadata, copies source variables
to setpoint names (e.g., P_m → P_m0, E_q → E_q0).
"""
function _build_setpoints(dyn::NamedTuple, sys::SystemModel)
    setpoint_pairs = Pair{Symbol, Any}[]

    for (model_name, states) in pairs(dyn)
        model = sys.models[model_name]
        model_type = get_base_model_type(model)
        metadata = model_metadata(model_type)

        # Skip models without initialized_params
        isempty(metadata.initialized_params) && continue

        # Build setpoints for this model, of the element type of its states (Dual when
        # the initialization is differentiated with respect to model parameters)
        model_setpoints = Pair{Symbol, AbstractVector}[]
        for init_param in metadata.initialized_params
            source_vals = getfield(states, init_param.source_var)
            push!(model_setpoints, init_param.name => copy(source_vals))
        end

        push!(setpoint_pairs, model_name => NamedTuple(model_setpoints))
    end

    return isempty(setpoint_pairs) ? NamedTuple() : NamedTuple(setpoint_pairs)
end

"""
    _extract_pf_snapshot(pf_sol::AbstractSteadyStateSolution) -> NamedTuple

Extract voltage magnitude and angle snapshot from steady-state solution, with the bus
numbers (`bus`) in the order of the values and the solution's largest residual
(`residual_norm`, see `get_residual_norm`).
"""
function _extract_pf_snapshot(pf_sol::AbstractSteadyStateSolution)
    am = get_address_manager(pf_sol)  # Use abstract interface
    y = get_state_vector(pf_sol)      # Use abstract interface
    v_addrs = get_addresses(am, :Bus, :v)
    theta_addrs = get_addresses(am, :Bus, :theta)

    return (
        v = y[v_addrs],
        theta = y[theta_addrs],
        bus = copy(get_system(pf_sol).models.Bus.i),
        residual_norm = get_residual_norm(pf_sol)
    )
end

"""
    pf_snapshot_positions(pf::NamedTuple) -> Dict

Position of each bus number in the power flow snapshot `pf` (see
`_extract_pf_snapshot`): the voltage of bus `b` is `pf.v[positions[b]]`,
`pf.theta[positions[b]]`.
"""
pf_snapshot_positions(pf::NamedTuple) = Dict(b => k for (k, b) in enumerate(pf.bus))

# Export public API
export extract_target_info
export resolve_init_inputs
export InitializationOrder, compute_initialization_order
export initialize_dynamic_states  # Public trait method (model-specific initialization)
