import Base: iterate, length

# ============================================================================
# GPU Hooks (extended by a GPU package for its array types)
# ============================================================================

"""
    synchronize_if_gpu(arr::AbstractArray)

Synchronize GPU if array is a GPU array. No-op for CPU arrays.
"""
synchronize_if_gpu(arr::AbstractArray) = nothing

"""
    to_gpu(sys::SystemModel) -> SystemModel

Transfer a SystemModel (or a problem) to the GPU: a copy whose model arrays live on
the device, with the matching `array_type`; the original is not modified. A GPU package
provides the methods; without one this fallback raises an error.
"""
function to_gpu(x)
    error("No GPU support is loaded: `to_gpu` has no method for $(typeof(x)). " *
          "A GPU package provides its methods.")
end

"""
    check_gpu_support(sys::SystemModel)

Throw an `ArgumentError` naming each model of `sys` that cannot run on the GPU: a model that
declares no GPU support (`gpu_support`), and a model with instances that is
`UnsupportedOnGPU`, with its reason. A GPU package's `to_gpu` calls it before any transfer.
"""
function check_gpu_support(sys::SystemModel)
    refused = String[]
    for (name, model) in pairs(sys.models)
        M = get_base_model_type(model)
        if !declares_gpu_support(M)
            push!(refused, "$name (no GPU support declared: a method of gpu_support returning SupportedOnGPU or UnsupportedOnGPU)")
            continue
        end
        support = gpu_support(M)
        support isa UnsupportedOnGPU && !isempty(model) && push!(refused, "$name: $(support.reason)")
    end
    isempty(refused) || throw(ArgumentError("The system cannot be transferred to the GPU: " * join(refused, "; ")))
    return nothing
end

# ============================================================================
# External Parameter Retrieval
# ============================================================================

"""
    retrieve_external_params(models_to_retrieve, all_models, index_order, phase)

Retrieve external parameters for specified models at the specified phase.

External parameters are model properties that reference other models (e.g., PQ needs
`bus.ide` from Bus). Retrieval happens in two phases:
- `PrePerUnit`: Before per-unit conversion (dimensionless params like bus type, indices)
- `PostPerUnit`: After per-unit conversion (physical quantities already in per-unit)

# Arguments
- `models_to_retrieve`: NamedTuple of models to retrieve external parameters FOR
- `all_models`: NamedTuple of all models (provides lookup context for references)
- `index_order`: IndexOrder containing foreign_refs for indexing
- `phase`: RetrievalPhase type (PrePerUnit or PostPerUnit)
- `pu_coeffs`: for `PrePerUnit` retrieval from `all_models` that are already per unit (a system
  derived from a built one), the per-unit coefficients of the system (`sys.pu_coeffs`, or `nothing`
  where it has none): a parameter that reads a converted field holds the physical value, as
  one retrieved before the conversion does. `missing` (default): the values are copied as they
  are. For `PostPerUnit`, the coefficients of the models retrieved for (`convert_system_to_per_unit`,
  `derive_system`), where the derived parameters (`RatingRatio`) of their metadata are kept.

# Returns
- NamedTuple with same keys as `models_to_retrieve`, values are NamedTuples of external param arrays

# Example
```julia
# TDS case: retrieve for dynamic models, but lookup in all models (including Bus)
external = retrieve_external_params(
    dynamic_models,  # Only GENCLS
    all_models,      # GENCLS + Bus + Branch + etc.
    index_order,
    PrePerUnit
)
# Result: (GENCLS=(Bus_v=[...], Bus_a=[...]), ...)

# Normal case: retrieve for all models
external = retrieve_external_params(models, models, index_order, PrePerUnit)
```
"""
function retrieve_external_params(
    models_to_retrieve::NamedTuple,
    all_models::NamedTuple,
    index_order::IndexOrder,
    phase::Type{<:RetrievalPhase};
    pu_coeffs = missing
)
    external_pairs = Pair{Symbol, NamedTuple}[]

    for model_name in keys(models_to_retrieve)
        model = models_to_retrieve[model_name]
        base_type = get_base_model_type(model)
        metadata = model_metadata(base_type)

        # Filter ext_params by retrieval phase (metadata.ext_params is the spec)
        phase_params = filter(ep -> typeof(ep.phase) == phase, metadata.ext_params)

        # The derived parameters (`RatingRatio`) are known once the model is converted: they come
        # with the coefficients of the conversion (`PostPerUnit`)
        derived = phase === PostPerUnit ? metadata.derived_params : RatingRatio[]

        if isempty(phase_params) && isempty(derived)
            # No params for this phase - store empty NamedTuple
            push!(external_pairs, model_name => NamedTuple())
            continue
        end

        # Retrieve params matching this phase
        param_vals_dict = Dict{Symbol, Any}()

        for d in derived
            held = pu_coeffs === missing || pu_coeffs === nothing ? nothing : get(pu_coeffs, model_name, nothing)
            (held === nothing || !haskey(held, d.name)) && throw(ArgumentError(
                "$model_name declares the derived parameter $(d.name) (the ratio of its rating to the system " *
                "base), which comes with the per-unit coefficients of the model, and there are none for it"))
            param_vals_dict[d.name] = copy(held[d.name])
        end

        for ext_ref in phase_params
            # UNIFIED: Always lookup by ref_key (the field in THIS model)
            lookup_key = DAEInfra.normalize_ref_key(ext_ref.terminal.ref_key)
            link_ref = DAEInfra.get_group_ref(index_order, model_name, lookup_key)

            # Every target model must exist (group refs may span several models);
            # unmatched keys already errored when the refs were built
            for tgt in unique(link_ref.target_type)
                haskey(all_models, tgt) || error("Model $model_name requires external parameter $(ext_ref.name) " *
                                                 "from $tgt, but $tgt not found in all_models")
            end

            # Gather values per-target (supports heterogeneous group refs)
            param_values = DAEInfra.gather_from_group(all_models, link_ref, ext_ref.ref_field)
            if phase === PrePerUnit && pu_coeffs !== missing
                param_values = _physical_values(param_values, all_models, pu_coeffs, link_ref, ext_ref.ref_field,
                                                model_name, ext_ref)
            end
            param_vals_dict[ext_ref.name] = param_values
        end

        # Convert to NamedTuple
        push!(external_pairs, model_name => NamedTuple(param_vals_dict))
    end

    return NamedTuple(external_pairs)
end

"""
    merge_external(pre::NamedTuple, post::NamedTuple) -> NamedTuple

Merge external parameters from PrePerUnit and PostPerUnit phases.

# Arguments
- `pre`: External params from PrePerUnit phase
- `post`: External params from PostPerUnit phase

# Returns
- NamedTuple with merged params for each model
"""
function merge_external(pre::NamedTuple, post::NamedTuple)
    # Both should have same keys (model names)
    @assert keys(pre) == keys(post) "PrePerUnit and PostPerUnit must have same model names"

    merged_pairs = Pair{Symbol, NamedTuple}[]

    for model_name in keys(pre)
        pre_params = getfield(pre, model_name)
        post_params = getfield(post, model_name)

        # Merge the two NamedTuples
        merged = merge(pre_params, post_params)

        push!(merged_pairs, model_name => merged)
    end

    return NamedTuple(merged_pairs)
end

"""
    retrieve_external_params(models_to_retrieve, all_models, index_order) -> NamedTuple

Retrieve and merge external parameters from both the `PrePerUnit` and
`PostPerUnit` phases. Use the 4-argument form when per-unit conversion
happens between the two phases. `pu_coeffs`: see the 4-argument form (the models are per unit
already: the `PrePerUnit` parameters are recovered to the physical value).
"""
function retrieve_external_params(
    models_to_retrieve::NamedTuple,
    all_models::NamedTuple,
    index_order::IndexOrder;
    pu_coeffs = missing
)
    return merge_external(
        retrieve_external_params(models_to_retrieve, all_models, index_order, PrePerUnit; pu_coeffs),
        retrieve_external_params(models_to_retrieve, all_models, index_order, PostPerUnit; pu_coeffs),
    )
end

# The external parameters of `model_name` (its `ExtParam`s of metadata) that read `field` of
# `model`, with the group reference each reads through
function _external_readers(sys::SystemModel, model_name::Symbol, model::Symbol, field::Symbol)
    models = getfield(sys.models, model_name)
    isempty(models) && return Tuple{ExtParam, Any}[]
    # (the system holds the external parameters of a model, retrieved from its declarations: a model
    # it holds none for reads no field, and its metadata, built at each call, is not asked for)
    held = get(sys.external, model_name, nothing)
    held !== nothing && isempty(held) && return Tuple{ExtParam, Any}[]
    found = Tuple{ExtParam, Any}[]
    for ep in model_metadata(get_base_model_type(models)).ext_params
        ep.ref_field === field || continue
        link_ref = DAEInfra.get_group_ref(sys.order, model_name,
                                          DAEInfra.normalize_ref_key(ep.terminal.ref_key))
        model in link_ref.target_type && push!(found, (ep, link_ref))
    end
    return found
end

"""
    external_dependents(sys, model::Symbol, field::Symbol) -> Vector{Symbol}

The models of `sys` with an external parameter (`ExtParam`) that is a copy of `field` of
`model` (an event that changes the field changes what their kernels read, and what their
automata's values read: `refresh_external!`).
"""
function external_dependents(sys::SystemModel, model::Symbol, field::Symbol)
    sys.order === nothing && return Symbol[]
    return [name for name in keys(sys.models) if !isempty(_external_readers(sys, name, model, field))]
end

"""
    refresh_external!(sys, model::Symbol, field::Symbol)

The external parameters retrieved from `field` of `model` are copies, taken when the
system was built: after an event changed the field they are copied again, in place (the
arrays the models' workspaces hold), as they were taken. A parameter retrieved after the
per-unit conversion (`PostPerUnit`) is a copy of the per-unit field. One retrieved before
it (`PrePerUnit`) holds the physical value, which the model now holds per unit: the physical
value is recovered through the per-unit coefficient of the field (`sys.pu_coeffs`, which
converted it), as it is at the build of the system and of the systems derived from it. A
field that no model declares as a parameter to convert (an index, a base voltage) was not
converted, and is copied as it is. An `ArgumentError` names a field that is converted
and has no coefficient in `sys`, whose physical value cannot be recovered.
"""
function refresh_external!(sys::SystemModel, model::Symbol, field::Symbol)
    sys.order === nothing && return nothing
    for name in keys(sys.models), (ep, link_ref) in _external_readers(sys, name, model, field)
        values = gather_from_group(sys.models, link_ref, field)
        ep.phase isa PrePerUnit &&
            (values = _physical_values(values, sys.models, sys.pu_coeffs, link_ref, field, name, ep))
        copyto!(getfield(getfield(sys.external, name), ep.name), values)
    end
    return nothing
end

# The physical values of the `field` of the targets of `link_ref`, from their per-unit `values`
# (per-unit = physical × coefficient): `values` as they are where no target converts it
_physical_values(values, models, pu_coeffs, link_ref, field::Symbol, reader::Symbol, ep::ExtParam) =
    (coefficients = _pu_coefficients(models, pu_coeffs, link_ref, field, reader, ep);
     coefficients === nothing ? values : values ./ coefficients)

# The per-unit coefficients (per-unit = physical × coefficient) of `field` of the targets of
# `link_ref`, in order; `nothing` where no target converts it. An `ArgumentError` where one does
# and there is no coefficient for it.
function _pu_coefficients(models, pu_coeffs, link_ref, field::Symbol, reader::Symbol, ep::ExtParam)
    n = length(link_ref.target_type)
    coefficients = ones(n)
    converted = false
    # Whether a target model converts `field` is asked once per run of instances of one model
    # (building a model's metadata is not free)
    last_target, last_converts = nothing, false
    for i in 1:n
        target = link_ref.target_type[i]
        if target !== last_target
            last_converts = any(p -> p.name === field, model_metadata(get_base_model_type(getfield(models, target))).params)
            last_target = target
        end
        last_converts || continue
        converted = true
        coefficient = pu_coeffs === nothing || !haskey(pu_coeffs, target) ||
                      !haskey(getfield(pu_coeffs, target), field) ? nothing :
                      getfield(pu_coeffs, target)[field]
        (coefficient === nothing || link_ref.target_idx[i] > length(coefficient)) && throw(ArgumentError(
            "$reader reads $field of $target as its external parameter $(ep.name), retrieved before the " *
            "per-unit conversion (PrePerUnit), but the system has no per-unit coefficient for " *
            "$target.$field to recover the physical value from the per-unit field"))
        coefficients[i] = coefficient[link_ref.target_idx[i]]
    end
    return converted ? coefficients : nothing
end

export retrieve_external_params, merge_external, external_dependents, refresh_external!

# Helper type for iteration
struct ComponentIterator
    system::SystemModel
    n::Int  # Current iteration number
end

Base.length(sys::SystemModel) = length(sys.models)

# Iteration over models
function Base.iterate(iter::ComponentIterator)
    models_tuple = values(iter.system.models)
    if iter.n > length(models_tuple)
        return nothing
    end

    return (models_tuple[iter.n], ComponentIterator(iter.system, iter.n + 1))
end

# Make SystemModel iterable
Base.iterate(sys::SystemModel) =
    iterate(ComponentIterator(sys, 1))

Base.iterate(sys::SystemModel, state::ComponentIterator) =
    iterate(state)


# ============================================================================
# System Residual Workspace
# ============================================================================

# ============================================================================
# System-Level Workspace
# ============================================================================

"""
    check_formulation_support(sys::SystemModel, formulation::FormulationTrait)

Throw an `ArgumentError` naming every model of `sys` that is not defined in
`formulation` (see [`supports_formulation`](@ref)). Problem constructors call it before
building anything, so an unsupported formulation fails at once instead of evaluating a
model with variables that do not exist.
"""
function check_formulation_support(sys::SystemModel, formulation::FormulationTrait)
    F = typeof(formulation)
    unsupported = [name for (name, model_vec) in pairs(sys.models)
                   if !supports_formulation(get_base_model_type(model_vec), F)]
    isempty(unsupported) || throw(ArgumentError(
        "The $(nameof(F)) formulation is not supported by $(join(unsupported, ", ")) " *
        "(no variables and residuals are defined for it)"))
    return nothing
end

"""
    OwnRows

Rows held for the out-of-service instances of a model: `algebraic[k][i]` is the address
of the k-th own algebraic residual of instance `i`, which is also the address of its
k-th own algebraic variable; `state[k][i]` the address of the k-th own zero-mass state
of instance `i` (its residual row and variable), and `zero_mass[k][i]` whether that
instance has zero mass.
"""
struct OwnRows
    algebraic::Vector{Vector{UInt}}
    state::Vector{Vector{UInt}}
    zero_mass::Vector{BitVector}
end


"""
    SystemBuffers{T, F}

System-level workspace with pre-allocated `ResidualBuffers` for every model type.

All memory is allocated during construction; the hot path
`evaluate_system_residuals!` performs zero allocation.

# Type Parameters
- `T`: Element type (Float64, or a ForwardDiff Dual type)
- `F`: Formulation trait type

# Fields
- `formulation`: Formulation trait
- `evaluators`: One closure per model, `(f, g, x, y, t) -> nothing`, in evaluation order
- `evaluators_with_recorder`: The same with an `InjectionRecorder` argument
- `evaluators_accumulate`: Accumulate-only closures for `evaluate_generator_contributions!`
- `workspaces_by_name`: Model name → `ResidualBuffers` (cold path, debugging, recorder)
- `models_by_name`: Model name → the model arrays its evaluator reads: the system's,
  except the columns given a copy of element type `T` (`typed_fields`)
- `metadata_by_name`: Model name → formulation-specific `ModelMetadata`
- `automaton_models`: The models that declare automata (their workspaces hold an
  `AutomatonState` of them), in evaluation order
- `own_models`: Models with own algebraic variables (variables without a source model)
  or own states whose mass is zero for some instance
- `own_rows`: Per entry of `own_models`, the rows it holds (see `OwnRows`); built on
  first use (`nothing` until then)
- `zero_mass`: Per entry of `own_models`, the own states with a mass coefficient that
  is zero for some instance, as `name => mask` pairs (`mask[i]`: instance `i` has zero
  mass); the mass matrix is fixed during a simulation, so they are found once
- `out_of_service`: Per entry of `own_models`, whether the model had an
  out-of-service instance in the last evaluation (set by its evaluator closure)
- `hold`: Values at which the own algebraic variables of out-of-service instances are
  held (`NaN` until set with [`set_hold!`](@ref)), shared by the workspaces of every
  element type built for the same solve
- `hold_x`: The same for their own zero-mass states

# Out-of-service instances

An instance that is out of service (`false` in the status dict) contributes nothing,
but its own algebraic variables are still entries of `y`. After all models are
evaluated, [`hold_out_of_service!`](@ref) replaces each of their residual rows with
`y[a] - hold[a] = 0`, so they keep a nonsingular equation: a bus de-energized by a
contingency, or the controller inputs of a tripped machine, stay at their held values.
Its own states with zero mass are algebraic too and are held the same way
(`x[a] - hold_x[a] = 0`); states with mass keep `M·ẋ = 0`, i.e. their values.

# Example
```julia
prob = ACPowerFlowProblem(sys)
sys_ws = SystemBuffers(Float64, prob.address_manager, prob.sys, prob.formulation)
status = build_status_dict(prob.sys)

# Hot path (zero allocation); a power flow's state is the operating point at t = 0
evaluate_system_residuals!(sys_ws, dae, status, 0.0)
```
"""
# ============================================================================
# DESIGN: Vector{Function} for Zero-Allocation Hot Path
# ============================================================================
#
# Problem: NamedTuple type parameters embed all model types, causing compilation
# explosion when model count grows (10 type params → exponential specializations).
#
# Solution: Replace typed tuple iteration with Vector{Function} closures.
# Each closure captures concrete types internally but presents type-erased
# Function interface externally → only 2 type parameters.
#
# Hot path: iterate Vector{Function}, each closure is fully specialized inside.
# Cold path: Dict{Symbol, Any} for debugging/InjectionRecorder access.
#
# Cost: One dynamic dispatch per model per evaluation (~nanoseconds).
# Savings: 50+ seconds compilation time reduction.
# ============================================================================
struct SystemBuffers{T<:Number, F<:FormulationTrait}
    formulation::F                              # Formulation trait for type dispatch
    evaluators::Vector{Function}                # Hot path: (f,g,x,y,t) -> nothing
    evaluators_with_recorder::Vector{Function}  # Recorder path: (f,g,x,y,t,recorder) -> nothing
    evaluators_accumulate::Vector{Function}     # Cold path: for evaluate_generator_contributions!
    workspaces_by_name::Dict{Symbol, Any}       # Cold path: for InjectionRecorder/debugging
    models_by_name::Dict{Symbol, Any}           # Cold path: model arrays read by the evaluators
    metadata_by_name::Dict{Symbol, Any}         # Cold path: for InjectionRecorder/debugging
    automaton_models::Vector{Symbol}            # Models that declare automata
    own_models::Vector{Symbol}                  # Models with own algebraic or zero-mass rows
    own_rows::Vector{Union{Nothing, OwnRows}}   # Per own model, built on first use
    zero_mass::Vector{Vector{Pair{Symbol, BitVector}}}  # Per own model, own zero-mass states
    out_of_service::Vector{Bool}                # Per own model, set by the evaluators
    hold::Vector{Float64}                       # Held values of out-of-service own variables
    hold_x::Vector{Float64}                     # Held values of their own zero-mass states
    foreign_refresh::Dict{Symbol, Base.RefValue{Any}}  # name => a (x,y,t)->nothing to run before this
                                                 # model's own calc_residuals!, or nothing (the default);
                                                 # set after the fact by PowerfulTDS where this model
                                                 # gathers a governed state of another model (`Governed.foreign`)
end

# ============================================================================
# Constructor Helpers
# ============================================================================

"""
    _compute_evaluation_order(sys, am, formulation) -> Vector{Symbol}

Compute model evaluation order by `evaluation_pass`, filtering to models that will have workspaces.

This is called BEFORE building workspaces so that NamedTuples are constructed in the
correct order for zero-allocation iteration via `values()`.
"""
function _compute_evaluation_order(sys::SystemModel, am::AddressManager, formulation::FormulationTrait)
    model_pairs = Pair{Symbol, Int}[]

    for (model_name, model_vec) in pairs(sys.models)
        # Skip models with no address assignments
        if !haskey(am.assignments, model_name)
            continue
        end

        form_metadata = get_metadata_for_formulation(
            model_metadata(get_base_model_type(model_vec)), formulation)

        # Skip models with no variables and no residuals
        if isempty(form_metadata.vars) && isempty(form_metadata.residuals)
            continue
        end

        push!(model_pairs, model_name => form_metadata.evaluation_pass)
    end

    # Sort by evaluation_pass
    sort!(model_pairs, by=last)

    return Symbol[first(p) for p in model_pairs]
end

# ============================================================================
# Constructor
# ============================================================================

"""
    SystemBuffers(::Type{T}, am::AddressManager, sys::SystemModel, formulation::FormulationTrait) where T

Construct system workspace with all model workspaces pre-allocated.

Iterates through all models in the system, creates a ResidualBuffers
for each using formulation-specific metadata. This function performs ALL
memory allocation - subsequent residual evaluations are zero-allocation.

The array type (CPU vs GPU) is inferred from the model data.
Use `to_gpu(sys)` to transfer the system to GPU before creating workspaces.

# Arguments
- `T`: Float type for calculations
- `am`: AddressManager with formulation-specific addresses
- `sys`: SystemModel containing model vectors; the system-wide settings come from `sys.settings`
- `formulation`: PowerInjection() or CurrentInjection()
- `hold`, `hold_x`: Held values of out-of-service own variables and zero-mass states
  (optional); pass those of the Float64 workspace to share them with a workspace of
  another element type
- `automata`: Model name → `AutomatonState` to share (optional; `nothing`: new ones); pass
  those of the Float64 workspace (`automaton_states`): a workspace of another element type
  shares their modes and discrete quantities and holds its memory and values in its own
  element type (`AutomatonState(state, T)`)
- `typed_fields`: Model name → fields whose columns the evaluators read from a copy
  of element type `T` instead of the system's (optional), so that values of type `T`
  (e.g. parameters as Dual numbers) can be written there; the other columns are shared

# Returns
- `SystemBuffers`: Fully allocated workspace ready for calculations

# Example
```julia
# Create problem which sets up infrastructure
prob = ACPowerFlowProblem(sys, CurrentInjection())
# SystemBuffers created internally
```
"""
function SystemBuffers(
    ::Type{T},
    am::AddressManager,
    sys::SystemModel,
    formulation::FormulationTrait;
    hold::Vector{Float64}=fill(NaN, Int(max_address_in_space(am, AlgebVar()))),
    hold_x::Vector{Float64}=fill(NaN, Int(max_address_in_space(am, StateVar()))),
    typed_fields::AbstractDict{Symbol, <:AbstractVector{Symbol}}=Dict{Symbol, Vector{Symbol}}(),
    automata::Union{Nothing, AbstractDict{Symbol, <:AutomatonState}}=nothing
) where T<:Number
    # Compute sorted evaluation order from sys.models
    sorted_names = _compute_evaluation_order(sys, am, formulation)
    own_models = Symbol[]
    own_rows = Union{Nothing, OwnRows}[]
    zero_mass = Vector{Pair{Symbol, BitVector}}[]
    out_of_service = Bool[]

    # Build closures and cold-path Dicts for each model
    evaluators = Function[]
    evaluators_with_recorder = Function[]
    evaluators_accumulate = Function[]
    workspaces_by_name = Dict{Symbol, Any}()
    models_by_name = Dict{Symbol, Any}()
    metadata_by_name = Dict{Symbol, Any}()
    automaton_models = Symbol[]
    foreign_refresh = Dict{Symbol, Base.RefValue{Any}}()

    for model_name in sorted_names
        # Note: _compute_evaluation_order already filtered for address assignments and vars/residuals
        model_vec = haskey(typed_fields, model_name) ?
            _with_typed_columns(getfield(sys.models, model_name), typed_fields[model_name], T) :
            getfield(sys.models, model_name)

        # Create formulation-specific metadata (cache to avoid hot-loop allocations)
        form_metadata = get_metadata_for_formulation(
            model_metadata(get_base_model_type(model_vec)), formulation)

        # Create workspace for this model with system's array type
        ws = ResidualBuffers(T, form_metadata, am, model_vec, sys.array_type;
                             system=sys.settings,
                             automata=automata === nothing ? nothing : get(automata, model_name, nothing))
        has_automata(form_metadata) && push!(automaton_models, model_name)

        # Store in Dicts for cold path access (InjectionRecorder, debugging)
        workspaces_by_name[model_name] = ws
        models_by_name[model_name] = model_vec
        metadata_by_name[model_name] = form_metadata
        # Index of this model in own_models (0: nothing to hold); its evaluators record
        # whether it has out-of-service instances for hold_out_of_service!
        k_own = 0
        zero_mass_states = _own_zero_mass_states(form_metadata, model_vec)
        if any(_is_own_algebraic, form_metadata.vars) || !isempty(zero_mass_states)
            push!(own_models, model_name); push!(own_rows, nothing)
            push!(zero_mass, zero_mass_states); push!(out_of_service, false)
            k_own = length(own_models)
        end

        # Capture concrete types in closures - these are fully specialized internally
        # but present type-erased Function interface externally
        model = model_vec
        ext = getfield(sys.external, model_name)
        meta = form_metadata
        form = formulation

        # A box PowerfulTDS may later point at a (x,y,t)->nothing to run before this model's own
        # calc_residuals!, where this model gathers a governed state of another model as a foreign
        # reference (`Governed.foreign`; see `_gather_and_calc!`'s own docstring). `nothing` (the
        # default, and every model today) costs one pointer load and a comparison.
        refresh_ref = Ref{Any}(nothing)
        foreign_refresh[model_name] = refresh_ref

        # Hot path closure - captures ALL with concrete types
        # Status is passed as Dict{Symbol, Vector{Bool}}, lookup active_mask by model_name
        #
        # `refresh_ref[]` is `Ref{Any}`-typed (it must be: it is set, if ever, to a closure of a type
        # only known once PowerfulTDS builds `Governed`, long after this function returns). Branching on
        # it HERE, not inside `_gather_and_calc!`, matters: the (overwhelmingly common) `=== nothing`
        # branch calls `_gather_and_calc!` with no `foreign_refresh!` keyword at all, so that call keeps
        # its existing, fully concrete specialization (the type-narrowing Julia already does after an
        # `=== nothing` check does not, by itself, make passing the still-`Any`-typed box as a keyword
        # argument free) -- passing `refresh_ref[]` itself as a keyword on every call, as an earlier draft
        # did, made every model's own hot path `Any`-typed and allocating, confirmed by
        # `test_tds_jacobian_closure.jl`'s allocation budget (180800 vs <= 10400). Every model takes this
        # branch today (`Governed.foreign` is empty everywhere); only a future cross-model foreign reader
        # would take the other one.
        let mn = model_name, k = k_own, refresh_ref = refresh_ref  # Capture model_name for Dict lookup
            push!(evaluators, (f, g, x, y, t, status) -> begin
                active = status[mn]  # Dict lookup: ~2ns, single concrete type
                rf = refresh_ref[]
                if rf === nothing
                    _gather_and_calc!(ws, model, ext, f, g, x, y, meta, form, t, active)
                else
                    _gather_and_calc!(ws, model, ext, f, g, x, y, meta, form, t, active; foreign_refresh! = rf)
                end
                scatter_residuals!(ws, f, g, active)
                k > 0 && (@inbounds out_of_service[k] = !_all_in_service(active))
            end)
        end

        # Recorder path closure - also captures model_name for recording
        let mn = model_name, k = k_own, refresh_ref = refresh_ref  # Capture model_name in let block for closure
            push!(evaluators_with_recorder, (f, g, x, y, t, status, recorder) -> begin
                active = status[mn]
                rf = refresh_ref[]
                if rf === nothing
                    _gather_and_calc!(ws, model, ext, f, g, x, y, meta, form, t, active)
                else
                    _gather_and_calc!(ws, model, ext, f, g, x, y, meta, form, t, active; foreign_refresh! = rf)
                end
                record_model_contributions!(recorder, mn, ws, t, active)
                scatter_residuals!(ws, f, g, active)
                k > 0 && (@inbounds out_of_service[k] = !_all_in_service(active))
            end)
        end

        # Accumulate path closure - for evaluate_generator_contributions!
        let mn = model_name, refresh_ref = refresh_ref
            push!(evaluators_accumulate, (f, g, x, y, t, status) -> begin
                active = status[mn]
                rf = refresh_ref[]
                if rf === nothing
                    _gather_and_calc!(ws, model, ext, f, g, x, y, meta, form, t, active)
                else
                    _gather_and_calc!(ws, model, ext, f, g, x, y, meta, form, t, active; foreign_refresh! = rf)
                end
                scatter_residuals_accumulate!(ws, f, g, active)  # Uses += instead of =
            end)
        end
    end

    return SystemBuffers{T, typeof(formulation)}(
        formulation,
        evaluators,
        evaluators_with_recorder,
        evaluators_accumulate,
        workspaces_by_name,
        models_by_name,
        metadata_by_name,
        automaton_models,
        own_models,
        own_rows,
        zero_mass,
        out_of_service,
        hold,
        hold_x,
        foreign_refresh
    )
end

"""
    _with_typed_columns(model::StructArray, fields, ::Type{T}) -> StructArray

`model` with the columns `fields` copied to element type `T`; the other columns are the
same arrays. The element type of the array stays that of `model`, which dispatch on the
model type needs, so rows of a copy are not accessed as structs.
"""
function _with_typed_columns(model::StructArray{M}, fields::AbstractVector{Symbol}, ::Type{T}) where {M, T}
    cols = StructArrays.components(model)
    for field in fields
        hasproperty(cols, field) || error("Model $(M) has no field :$field")
    end
    typed = map(keys(cols), values(cols)) do name, col
        name in fields ? T.(col) : col
    end
    return StructArray{M}(NamedTuple{keys(cols)}(typed))
end

"""
    with_float_columns(model::StructArray, ::Type{T}) -> StructArray

`model` with each of its floating-point columns (its real parameters) in a copy of element
type `T`, its other columns shared: for evaluating a model's kernels with parameters of
another type (the trace of the variables its automaton values read, PowerfulTDS
`trace_value_reads`, which makes them tracers).
"""
function with_float_columns(model::StructArray, ::Type{T}) where T
    cols = StructArrays.components(model)
    return _with_typed_columns(model, [name for (name, col) in pairs(cols) if eltype(col) <: AbstractFloat], T)
end

"""
    automaton_states(sys_ws::SystemBuffers) -> Dict{Symbol, AutomatonState}

The `AutomatonState` of each model of `sys_ws` with automata, to share with a workspace
of another element type (`SystemBuffers` keyword `automata`).
"""
automaton_states(sys_ws::SystemBuffers) =
    Dict{Symbol, AutomatonState}(name => sys_ws.workspaces_by_name[name].automata
                                 for name in sys_ws.automaton_models)

"""
    sync_automaton_memory!(dst::SystemBuffers, src::SystemBuffers)

Copy the memory of the automata of every model of `src` (the Float64 workspace, where
the simulation sets it) into the memory of `dst` (a workspace of another element type
built beside it, which holds its own, see `AutomatonState(state, T)`), converting the
element type: values without derivatives.
"""
function sync_automaton_memory!(dst::SystemBuffers, src::SystemBuffers)
    for name in src.automaton_models
        _copy_memory!(dst.workspaces_by_name[name], src.workspaces_by_name[name])
    end
    return nothing
end

"""
    automaton_memory_synced(dst::SystemBuffers, src::SystemBuffers) -> Bool

Whether the memory of the automata of every model of `dst` equals that of `src` (the Float64
workspace) as `sync_automaton_memory!` would leave it: the same values (`NaN` before a start
counted equal), without derivatives. The check of the rates' test hook, which asserts that nothing
leaves the workspaces apart.
"""
function automaton_memory_synced(dst::SystemBuffers, src::SystemBuffers)
    return all(name -> _memory_synced(dst.workspaces_by_name[name], src.workspaces_by_name[name]), src.automaton_models)
end

_memory_synced(dst::ResidualBuffers, src::ResidualBuffers) =
    _memory_synced(dst.automata.memory, src.automata.memory)
_memory_synced(dst::NamedTuple{K}, src::NamedTuple{K}) where K = all(map(_column_synced, values(dst), values(src)))
# (`isequal` of a Dual and a number: equal values and no partials; `NaN` equals `NaN`)
_column_synced(dst::AbstractVector, src::AbstractVector) = all(i -> isequal(dst[i], src[i]), eachindex(src))

# The memory of the automata of a model's workspace, copied from `src`'s into `dst`'s (a function
# barrier: the workspaces are found by name, so their types are known at run time only)
_copy_memory!(dst::ResidualBuffers, src::ResidualBuffers) = _copy_memory!(dst.automata.memory, src.automata.memory)

# The columns of a model's memory, unrolled
_copy_memory!(dst::NamedTuple{K}, src::NamedTuple{K}) where K = (map(_copy_column!, values(dst), values(src)); nothing)
_copy_column!(dst::AbstractVector, src::AbstractVector) = (dst .= src; nothing)

# Whether every instance is in service. On the CPU a branch-free AND that vectorizes
# (the short-circuiting `all` checks one element at a time); `all` reduces on the GPU.
@inline function _all_in_service(active::Vector{Bool})
    m = true
    @inbounds @simd for x in active
        m &= x
    end
    return m
end
_all_in_service(active::AbstractVector{Bool}) = all(active)

# An own algebraic variable or residual: algebraic, without a source model
_is_own_algebraic(item) = (item.address_type isa AlgebVar || item.address_type isa AlgebRes) &&
                          item.source_model === nothing

# An own state: differential, without a source model
_is_own_state(item) = (item.address_type isa StateVar || item.address_type isa StateRes) &&
                      item.source_model === nothing

# Own states whose mass coefficient is zero for some instance, as name => mask pairs
function _own_zero_mass_states(metadata, model_vec)
    states = Pair{Symbol, BitVector}[]
    for var in metadata.vars
        _is_own_state(var) || continue
        mass_coeff = get_mass_coeff(var)
        isnothing(mass_coeff) && continue
        mask = BitVector(iszero.(Array(eval_mass_coeff(mass_coeff, model_vec))))
        any(mask) && push!(states, var.name => mask)
    end
    return states
end

"""
    _own_rows(ws, metadata, model_name, zero_mass) -> OwnRows

The rows held for the out-of-service instances of a model (see `OwnRows`). Errors unless
its own algebraic residuals share addresses with its own algebraic variables, in the
same order, and each own zero-mass state has an own residual row at its address.
"""
function _own_rows(ws, metadata, model_name::Symbol, zero_mass)
    algebraic = _own_algebraic_rows(ws, metadata, model_name)
    own_state_rows = [getfield(ws.res_addrs, r.name) for r in metadata.residuals if _is_own_state(r)]
    state = Vector{Vector{UInt}}()
    for (name, _) in zero_mass
        addrs = getfield(ws.var_addrs, name)
        addrs in own_state_rows || error("$model_name has out-of-service instances, but its " *
            "zero-mass state $name has no own residual row at its address, so it cannot be held")
        push!(state, addrs)
    end
    return OwnRows(algebraic, state, [mask for (_, mask) in zero_mass])
end

"""
    _own_algebraic_rows(ws, metadata, model_name) -> Vector{Vector{UInt}}

Addresses of the own algebraic residuals of a model, per residual and instance.
Errors unless they are the addresses of its own algebraic variables, in the same
order, which is what lets `hold_out_of_service!` hold those variables.
"""
function _own_algebraic_rows(ws, metadata, model_name::Symbol)
    own_var_addrs = [getfield(ws.var_addrs, v.name) for v in metadata.vars if _is_own_algebraic(v)]
    own_res_addrs = [getfield(ws.res_addrs, r.name) for r in metadata.residuals if _is_own_algebraic(r)]
    own_res_addrs == own_var_addrs || error("$model_name has out-of-service instances, but " *
        "its own algebraic residuals do not share addresses with its own algebraic " *
        "variables, so they cannot be held")
    return Vector{Vector{UInt}}(own_res_addrs)
end

"""
    set_hold!(sys_ws::SystemBuffers, y::AbstractVector)
    set_hold!(sys_ws::SystemBuffers, x::AbstractVector, y::AbstractVector)

Hold the own algebraic variables of out-of-service instances at their values in `y`
(e.g. the initial guess of a solve), and their own zero-mass states at their values in
`x`; see [`hold_out_of_service!`](@ref).
"""
set_hold!(sys_ws::SystemBuffers, y::AbstractVector) = (copyto!(sys_ws.hold, y); nothing)
function set_hold!(sys_ws::SystemBuffers, x::AbstractVector, y::AbstractVector)
    copyto!(sys_ws.hold_x, x)
    copyto!(sys_ws.hold, y)
    return nothing
end

"""
    hold_out_of_service!(sys_ws, f, g, x, y, status)

Replace every own algebraic residual row `a` of an out-of-service instance by
`g[a] = y[a] - hold[a]`, and the row of each own zero-mass state by
`f[a] = x[a] - hold_x[a]`. Called by `evaluate_system_residuals!` after all models are
evaluated (whose evaluator closures record in `sys_ws.out_of_service` which models
have out-of-service instances, so the pass is free when every instance is in service),
so the row keeps an equation (its Jacobian row is the identity) although the instance
contributes nothing. Errors if an out-of-service instance belongs to a model whose own
rows cannot be held, or if the held values were never set.
"""
function hold_out_of_service!(sys_ws::SystemBuffers, f::AbstractVector, g::AbstractVector,
                              x::AbstractVector, y::AbstractVector,
                              status::Dict{Symbol, <:AbstractVector{Bool}})
    @inbounds for k in eachindex(sys_ws.out_of_service)
        sys_ws.out_of_service[k] || continue
        model_name = sys_ws.own_models[k]
        rows = sys_ws.own_rows[k]
        if rows === nothing  # first use: cold path
            rows = _own_rows(sys_ws.workspaces_by_name[model_name],
                             sys_ws.metadata_by_name[model_name], model_name, sys_ws.zero_mass[k])
            sys_ws.own_rows[k] = rows
        end
        active = status[model_name]
        _hold_rows!(g, y, sys_ws.hold, rows.algebraic, active, model_name)
        for j in eachindex(rows.state)
            _hold_rows!(f, x, sys_ws.hold_x, (rows.state[j],), active, model_name, rows.zero_mass[j])
        end
    end
    return nothing
end

# Scalar loop over the (few) out-of-service instances; GPU arrays reject scalar
# indexing, so out-of-service instances are a CPU feature. `selected[i]`: instance `i`'s
# row is held when it is out of service.
function _hold_rows!(r::AbstractVector, v::AbstractVector, hold::Vector{Float64}, rows, active,
                     model_name::Symbol, selected = nothing)
    @inbounds for addrs in rows, i in eachindex(addrs)
        active[i] && continue
        selected === nothing || selected[i] || continue
        a = addrs[i]
        h = hold[a]
        isnan(h) && error("$model_name instance $i is out of service, but the values at which " *
                          "its own variables are held were never set (set_hold!)")
        r[a] = v[a] - h
    end
    return nothing
end


# ============================================================================
# System Residual Evaluation (Hot Path - Zero Allocation)
# ============================================================================

"""
    evaluate_system_residuals!(sys_ws::SystemBuffers{T}, dae::DAE{T}) where T

Evaluate residuals for all models in the system - ZERO ALLOCATION.

This is the hot path called in every Newton-Raphson iteration. Orchestrates
the gather-compute-scatter pattern across all model types:

1. Zero out global residual vector
2. For each model:
   - Gather variables from DAE into local workspace
   - Zero local residual buffers
   - Calculate model-specific residuals
   - Scatter results back to global DAE

# Arguments
- `sys_ws`: System workspace with pre-allocated buffers
- `dae`: DAE system with global variable and residual vectors

# Performance Notes
- Zero allocation after initial workspace construction
- Type-stable dispatch on each model's calc_residuals!
- Models are evaluated sequentially (can be parallelized if SharedRes uses atomics)
- Compiler inlines and optimizes the entire call chain

# Example
```julia
# In Newton-Raphson loop
for iter in 1:maxiters
    evaluate_system_residuals!(sys_ws, dae)  # ← Zero allocation here

    if norm(dae.g, Inf) < tol
        break
    end

    calculate_jacobian!(J, sys_ws, dae)
    Δy = J \\ dae.g
    dae.y .-= Δy
end
```
"""

"""
    evaluate_system_residuals!(sys_ws, dae, status, t)

DAE-based API for residual evaluation (convenience wrapper), at time `t`, the time of the
state in `dae`.

Extracts vectors from DAE and delegates to 7-arg implementation.
Handles both power flow (empty f,x) and TDS (non-empty f,x).

# Arguments
- `sys_ws::SystemBuffers{T}`: Pre-allocated workspace
- `dae::DAE`: DAE container with f, g, x, y vectors
- `status::Dict{Symbol, Vector{Bool}}`: Active status for each model
- `t`: Time of the state (no default: a model's residuals can depend on it directly)
"""
function evaluate_system_residuals!(
    sys_ws::SystemBuffers{T},
    dae::DAE{AT},
    status::Dict{Symbol, <:AbstractVector{Bool}},
    t::Real
) where {T, AT}
    evaluate_system_residuals!(sys_ws, dae.f, dae.g, dae.x, dae.y, status, t)
    return nothing
end

"""
    evaluate_system_residuals!(sys_ws, f, g, x, y, status, t; recorder=nothing)

DAE-based API for residual evaluation with differential states, at time `t`.

This interface supports time-domain simulation where models have both
differential states (x) and algebraic variables (y), producing both
differential residuals (f) and algebraic residuals (g).

# Arguments
- `sys_ws::SystemBuffers`: Pre-allocated workspace
- `f::AbstractVector`: Differential residual vector (dx/dt contributions)
- `g::AbstractVector`: Algebraic residual vector
- `x::AbstractVector`: Differential state vector
- `y::AbstractVector`: Algebraic variable vector
- `status::Dict{Symbol, Vector{Bool}}`: Active status for each model (required)
- `t`: Time of the state `(x, y)` (required: a model's residuals can depend on it directly)
- `recorder`: Optional InjectionRecorder for tracking model contributions

# Implementation
For each model:
- If model has differential states: gather from both x and y, scatter to both f and g
- If model is algebraic-only: gather from y, scatter to g (same as power flow)

# Example
```julia
# TDS problem with dynamic models
prob = TDSProblem(setup_tds(sol_pf, dynamic_models))
adcache = create_adcache(prob.address_manager, sys, prob.formulation, ForwardDiffBackend())

# Evaluate DAE residuals
x = zeros(n_diff)  # Differential states
y = zeros(n_alg)   # Algebraic variables
f = zeros(n_diff)  # Differential residuals
g = zeros(n_alg)   # Algebraic residuals
status = build_status_dict(sys)

adcache(f, g, x, y, status, t)  # Uses this method internally
```
"""
function evaluate_system_residuals!(
    sys_ws::SystemBuffers,
    f::AbstractVector,
    g::AbstractVector,
    x::AbstractVector,
    y::AbstractVector,
    status::Dict{Symbol, <:AbstractVector{Bool}},
    t::Real;
    recorder = nothing  # Optional InjectionRecorder for tracking model contributions
)
    # Zero out global residual vectors
    fill!(f, zero(eltype(f)))
    fill!(g, zero(eltype(g)))

    # Evaluate each model using Vector{Function} iteration
    # Each closure is fully specialized internally, but calling through
    # Vector{Function} requires dynamic dispatch (~16 bytes per call).
    # For 6 models: 96 bytes total. This is the intentional trade-off for
    # 50+ second compile time reduction.
    if recorder === nothing
        # Fast path: iterate evaluator closures
        @inbounds for eval_fn in sys_ws.evaluators
            eval_fn(f, g, x, y, t, status)
        end
    else
        # Recording path: iterate recorder closures
        @inbounds for eval_fn in sys_ws.evaluators_with_recorder
            eval_fn(f, g, x, y, t, status, recorder)
        end
    end

    # Own rows of out-of-service instances keep an equation
    hold_out_of_service!(sys_ws, f, g, x, y, status)

    # Synchronize GPU kernels if needed
    synchronize_if_gpu(g)

    nothing
end

"""
    evaluate_generator_contributions!(sys_ws::SystemBuffers{T}, y::AbstractVector{T}, status::Dict{Symbol, Vector{Bool}}) where T

Evaluate generator contributions (reactive power injections) after power flow converges.

This cold-path function uses `scatter_residuals_accumulate!` instead of `scatter_residuals!`,
which skips `:set` mode and allows generator control equations (PV voltage control, Slack
power balance) to accumulate their contributions instead of replacing bus power balance
equations. This enables extraction of actual reactive power injections from generators.

# Implementation
Uses `_eval_models_accumulate!` which calls `scatter_residuals_accumulate!` for each model.
This treats `SharedRes(:set)` as `:add` during scatter, allowing control equations to
contribute to power balance instead of overwriting it.

# Physical Interpretation
After power flow converges, residuals g = 0. With accumulate-only scatter:
- For PV generators: g[bus_Q] = actual Q injection to maintain V = Vs
- For Slack generators: g[bus_P] and g[bus_Q] = actual P and Q injections

# Arguments
- `sys_ws`: Pre-allocated system workspace
- `y`: Converged solution vector (from power flow)
- `status`: Dict{Symbol, Vector{Bool}} with active status for each model

# Returns
- `g`: Residual vector containing generator contributions at bus injection positions

# Example
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())
ws = setup_powerflow_newton(prob)
status = build_status_dict(sys)
g_gen = evaluate_generator_contributions!(get_float_ws(ws.adcache), sol.dae.y, status)
```

# See Also
- `evaluate_system_residuals!`: Base residual evaluation function
- `scatter_residuals_accumulate!`: Cold-path scatter that skips `:set` mode
- `get_generator_residuals`: Lazy accessor for cached generator contributions in PowerFlowSolution
"""
function evaluate_generator_contributions!(
    sys_ws::SystemBuffers{T},
    y::AbstractVector{T},
    status::Dict{Symbol, <:AbstractVector{Bool}}
) where T
    # Allocate temporary residual and empty differential vectors (PF has no differential states)
    g = zeros(T, length(y))
    f = Vector{T}(undef, 0)
    x = Vector{T}(undef, 0)

    # Cold path - iterate accumulate closures (allocation OK in cold path)
    # Uses scatter_residuals_accumulate! which skips :set mode
    for eval_fn in sys_ws.evaluators_accumulate
        eval_fn(f, g, x, y, zero(T), status)
    end

    return g
end


# ============================================================================
# Evaluate Model Residuals
# ============================================================================

"""
    evaluate_model_residuals!(ws::ResidualBuffers{T}, model_vec::AbstractModel,
                              external::NamedTuple, dae::DAE{T}, metadata::ModelMetadata;
                              force_accumulate::Bool=false) where T

Evaluate residuals for a single model type - zero allocation.

Generic function that works for any model by dispatching to model-specific
calc_residuals! implementations.

# Steps
1. Gather variables from DAE.y into workspace
2. Gather accumulated residuals for SharedRes(:set) models (if needed)
3. Calculate residuals (model-specific dispatch)
4. Scatter residuals to DAE.g

# Arguments
- `ws`: Pre-allocated workspace for this model
- `model_vec`: Model vector with data
- `external`: External parameters for this model (topology-derived references)
- `g, y` or `dae`: Residual and variable vectors (or DAE container)
- `metadata`: Model metadata for gather/scatter operations
- `formulation`: Formulation trait (PowerInjection or CurrentInjection)
"""

# Gather a model's inputs and compute its residuals into the workspace (no scatter).
#
# Where `metadata.automata` governs a state (`governs_states`), the governed state's own column of
# `ws.var_vals` is substituted with its relation's current value before `calc_residuals!` reads it
# (`enforce_governed_affine!`, round 19a Phase D): the residual a Newton iterate or a Jacobian/tgrad
# perturbation reads must see the relation, not the value the integration drifted to between the
# discrete writes that land it there exactly (`_enter_governed!`). Unlike the values path
# (`evaluate_automaton_values!`), and unlike this path's own earlier (Phase C) shape, this needs no
# `calc_automaton_values!` kernel call at all: `calc_residuals!` reads `ws.var_vals` directly, never
# `ws.automata.values`, so a governed state's relation can be read as its pre-extracted affine form (a
# dot product over the already-gathered variables, `AffineRelation`/`refresh_affine_governed!`) instead
# of a second kernel invocation -- the relation is, by `check_governed`'s own admissibility proof,
# exactly affine in those variables with parameter-dependent coefficients, computed once per parameter
# set (at setup, and whenever an event retraces the model), not per residual evaluation. A model without a
# governed state pays only one `Bool` field read, `metadata.governs_any_state` (cached at `ModelMetadata`
# construction, `Foundation/core_types.jl`): `governs_states(metadata.automata)` itself is NOT free to call
# here -- `automata`'s declared field type, `Tuple{Vararg{Automaton}}`, is not concrete, so calling it on a
# `metadata::ModelMetadata` VALUE is a dynamic dispatch on every call, regardless of which model (measured
# and fixed in round 19a Phase D's second pass: the private round 19a design/timing notes
# record the dispatch found via `@code_warntype`, the fix, and a
# paired, simultaneous-process re-measurement of this exact path for two non-governed models -- GENROU and
# PQ -- that, once the dynamic dispatch is gone, shows no reproducible regression against `main` left: the
# gap an earlier single-trial measurement showed was noise, not demonstrated by a repeated, paired trial).
# A `@generated`, model-TYPE-dispatched trait (to fold the branch away entirely at compile time for a model
# whose type never governs a state, not just make the runtime check cheap) was tried and reverted: its
# generator calls `model_metadata(M)`, which is defined per model in packages that load AFTER PowerfulCore
# (PowerfulPF, PowerfulTDS, ...), and a `@generated` function's generator only sees methods that exist at
# ITS OWN definition time -- confirmed by a `MethodError` reproduced directly, including long after every
# package has finished loading, not merely during precompilation.
#
# `foreign_refresh!`, where not `nothing` (set by PowerfulTDS after an `AutomatonRuntime`'s `Governed`
# finds this model gathers *another* model's governed state as a foreign reference,
# `SystemBuffers.foreign_refresh`): run after `gather_variables!`, before `calc_residuals!`, and
# overwrites this model's own `ws.var_vals` columns that mirror the foreign governed state with its
# owner's current relation (never the raw, possibly drifted value `x`/`y` holds at that address) --
# the same substitution `enforce_governed_affine!` makes for a model's own governed state, carried to
# every other model that reads it, since each model's own `gather_variables!` call would otherwise read
# the raw address and overwrite whatever an earlier pass wrote. No model in this codebase does this today
# (`Governed.foreign` is empty for every installed model); the hook costs nothing where it is `nothing`.
@inline function _gather_and_calc!(ws, model_vec, external, f, g, x, y, metadata, formulation, t, active_mask;
                                   foreign_refresh! = nothing)
    gather_variables!(ws, x, y)
    foreign_refresh! === nothing || foreign_refresh!(x, y, t)
    metadata.governs_any_state && enforce_governed_affine!(ws, metadata)
    gather_residuals!(ws, f, g)
    calc_residuals!(ws, model_vec, external, formulation, t; active_mask=active_mask)
    nothing
end

# Unified evaluate_model_residuals! - handles both PF (empty f,x) and TDS (non-empty f,x).
# Vectors may hold Duals while the workspace is Float64 (ForwardDiff Jacobians).
@inline function evaluate_model_residuals!(
    ws::ResidualBuffers{T},
    model_vec::StructArray{M},
    external::NamedTuple,
    f::AbstractVector,
    g::AbstractVector,
    x::AbstractVector,
    y::AbstractVector,
    metadata::ModelMetadata,
    formulation::FormulationTrait,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, M<:AbstractModel}
    _gather_and_calc!(ws, model_vec, external, f, g, x, y, metadata, formulation, t, active_mask)
    scatter_residuals!(ws, f, g, active_mask)
    nothing
end

# DAE-based evaluate_model_residuals! - convenience wrapper
@inline function evaluate_model_residuals!(
    ws::ResidualBuffers{T},
    model_vec::StructArray{M},
    external::NamedTuple,
    dae::DAE,
    metadata::ModelMetadata,
    formulation::FormulationTrait,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, M<:AbstractModel}
    evaluate_model_residuals!(ws, model_vec, external, dae.f, dae.g, dae.x, dae.y, metadata, formulation, t; active_mask=active_mask)
end


# ============================================================================
# Display Methods
# ============================================================================

"""
    Base.show(io::IO, sys::SystemModel)

Display component counts and array backend.
"""
function Base.show(io::IO, sys::SystemModel)
    AT = sys.array_type
    backend_name = AT <: Vector ? "Vector" : string(nameof(AT))

    if get(io, :compact, false)
        print(io, "SystemModel{$backend_name}(...)")
    else
        println(io, "SystemModel{$backend_name}:")
        println(io, "  Components:")
        for (name, data) in pairs(sys.models)
            println(io, "    $name: $(length(data)) elements")
        end
    end
end

# ============================================================================
# Status Dict Helpers (for parallel ACCC)
# ============================================================================

"""
    build_status_dict(sys::SystemModel) -> Dict{Symbol, Vector{Bool}}

Build a dictionary mapping model names to their in-service status: a **copy** of
each model's status column (see `get_status_field`), or all `true` for models
without one.

The dict belongs to whoever builds it (a solver workspace, a contingency): changing
it never touches `sys`, and changing `sys` afterwards does not reach it. Events
that switch devices during a simulation update the workspace's dict (see
`apply!(sys, status, mod)`).

# Example
```julia
status = build_status_dict(sys)
status[:Branch][5] = false  # Trip branch 5 for this solve only
```
"""
function build_status_dict(sys::SystemModel)
    AT = sys.array_type
    status = Dict{Symbol, AT{Bool}}()
    for (name, model_vec) in pairs(sys.models)
        sf = get_status_field(model_vec)
        if sf !== nothing
            # Copy (Vector or CuVector), converted to Bool if stored as Int
            status[name] = copy(getproperty(model_vec, sf))
        else
            # Models without status field are always active
            # Create array of correct type filled with true
            arr = AT{Bool}(undef, length(model_vec))
            fill!(arr, true)
            status[name] = arr
        end
    end
    return status
end

"""
    copy_status_dict(status::Dict{Symbol, Vector{Bool}}) -> Dict{Symbol, Vector{Bool}}

Copy a status dictionary, e.g. one per thread or contingency from a shared base.

# Example
```julia
base_status = build_status_dict(sys)

Threads.@threads for i in 1:n_contingencies
    status = copy_status_dict(base_status)
    status[:Branch][contingency_idx[i]] = false
    # ... solve with thread-local status
end
```
"""
function copy_status_dict(status::Dict{Symbol, V}) where V<:AbstractVector{Bool}
    Dict{Symbol, V}(k => copy(v) for (k, v) in status)
end

# ============================================================================
# Exports
# ============================================================================

export SystemBuffers
export evaluate_system_residuals!, evaluate_generator_contributions!
export build_status_dict, copy_status_dict, check_formulation_support
export set_hold!, hold_out_of_service!
export automaton_states, sync_automaton_memory!, automaton_memory_synced, with_float_columns
