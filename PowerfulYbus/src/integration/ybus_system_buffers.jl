# ============================================================================
# Y-Bus System Buffers
# ============================================================================
#
# Metadata-driven workspace for Y-bus TDS evaluation using closure pattern.
# Similar to SystemBuffers but separates generator and controller evaluation.
#
# ============================================================================

using PowerfulCore: SystemModel, AddressManager, ResidualBuffers, FormulationTrait
using PowerfulCore: PowerInjection, model_metadata, get_base_model_type, SimConfig
using PowerfulCore: gather_variables!, gather_residuals!, calc_residuals!, scatter_residuals!
using PowerfulCore: StateVar, StateRes, PrivateRes, SharedRes
using PowerfulCore: build_res_scatter_tuple, scatter_res_add!, scatter_res_set!
using PowerfulCore: build_status_dict

"""
    YBusSystemBuffers{T, F}

System-level workspace for Y-bus TDS evaluation with metadata-driven discovery.

Uses closure pattern (like SystemBuffers) to avoid compilation explosion.
Separates generators (need network currents) from controllers (standard evaluation).

# Type Parameters
- `T`: Float type (Float64, Float32, etc.)
- `F`: FormulationTrait type (PowerInjection, etc.)

# Fields
- `formulation::F`: Formulation trait for type dispatch
- `generator_evaluators::Vector{Function}`: Y-bus generators: (f,g,x,y,t,net_sol) -> nothing
- `internal_voltage_gatherers::Vector{Function}`: Per machine model: (E,x,y) -> nothing,
  setting the internal voltages `E[k]` of its in-service instances (network generator
  order `model_to_gen_idx`)
- `controller_evaluators::Vector{Function}`: Controllers: (f,g,x,y,t,V_bus) -> nothing
- `workspaces_by_name::Dict{Symbol, Any}`: For debugging access
- `metadata_by_name::Dict{Symbol, Any}`: For debugging access
- `models_by_name::Dict{Symbol, Any}`: Model arrays for gather_internal_voltages!
- `externals_by_name::Dict{Symbol, Any}`: External parameters of each model (e.g. `Mratio`, its
  base over the system base, for the current sources)
- `model_to_gen_idx::Dict{Symbol, Vector{Int}}`: Maps model_name → [model_idx → gen_buses_idx]
- `automaton_models::Vector{Symbol}`: The models that declare automata (their workspaces
  hold an `AutomatonState`, run by the automaton runtime of PowerfulTDS)
- `status`: Model name → in-service mask, read by the evaluators (and the automata)

# Design
- **Zero hard-coding**: Uses `network_partition_role` trait for discovery
- **Metadata-driven**: Uses `model_metadata()` for structure
- **Type-erased closures**: Each closure fully specialized internally
- **Allocation-free hot path**: All memory pre-allocated during construction
- **One source of bus voltage**: controller variables that are bus voltages
  (`bus_voltage_inputs`) are set from the network's bus voltages

# Example
```julia
# Create Y-bus system buffers (computes gen index mapping internally)
sys_ws = YBusSystemBuffers(Float64, am, sys, PowerInjection(), tds_prob.init)

# Evaluation (with voltage reconstruction)
network_sol = solve_network!(cache, sys)
V_bus = cache.V_bus
evaluate_ybus_residuals!(sys_ws, dae, network_sol, V_bus, t)
```
"""
struct YBusSystemBuffers{T<:Number, F<:FormulationTrait}
    formulation::F
    generator_evaluators::Vector{Function}       # (f,g,x,y,t,net_sol) -> nothing
    internal_voltage_gatherers::Vector{Function} # (E,x,y) -> nothing
    controller_evaluators::Vector{Function}      # (f,g,x,y,t,V_bus) -> nothing
    workspaces_by_name::Dict{Symbol, Any}        # Cold path: for debugging
    metadata_by_name::Dict{Symbol, Any}          # Cold path: for debugging
    models_by_name::Dict{Symbol, Any}            # Cold path: model arrays for gather_internal_voltages!
    externals_by_name::Dict{Symbol, Any}         # Cold path: external parameters (current sources' Mratio)
    model_to_gen_idx::Dict{Symbol, Vector{Int}}  # model_name → [model_idx → gen_buses_idx]
    automaton_models::Vector{Symbol}             # Models that declare automata
    status::Dict{Symbol, <:AbstractVector{Bool}} # In-service masks of the evaluators
end

"""
    YBusSystemBuffers(::Type{T}, am::AddressManager, sys::SystemModel,
                      formulation::FormulationTrait, tds_init; automata=nothing) where T

Construct Y-bus system workspace with metadata-driven discovery.

The setpoints of every workspace (E_q0, P_m0, ...) are set from the TDS
initialization `tds_init` (`sync_ybus_setpoints!`), so buffers of any number type
built from the same initialization evaluate the same residuals.

**Zero hard-coding**: Uses `network_partition_role` trait to automatically discover
which models are generators (MachineWithInternalVoltage) vs controllers (NotInYBusFormulation).

# Arguments
- `T`: Float type for calculations
- `am`: AddressManager with formulation-specific addresses
- `sys`: SystemModel containing model vectors
- `formulation`: PowerInjection() (or others when implemented)
- `tds_init`: TDS initialization (`tds_prob.init`) holding the setpoint values
- `automata`: Model name → `AutomatonState` to share (optional; `nothing`: new ones), as
  `automaton_states` of another `YBusSystemBuffers` returns them: the workspaces of
  every number type, and the trapezoidal solver's and the DAE residual's, read the
  discrete state the automaton runtime sets in one

# Returns
- `YBusSystemBuffers`: Fully allocated workspace ready for Y-bus evaluation

# Discovery Process
0. Refuse a model without a declared role or with an unsupported one (`check_ybus_support`)
1. Compute generator index mapping (model_idx → gen_buses_idx) from sys.models
2. Iterate all models in `sys.models`
3. Check `network_partition_role(model_type)`
4. If `MachineWithInternalVoltage`: create generator evaluator with `calc_ybus_residuals!`
5. If `NotInYBusFormulation`: create controller evaluator with `calc_residuals!`
6. If `NonConformingLoad` with state variables (e.g., REGCA1): create controller evaluator
   (current injection still handled by NC solver, but state residuals need evaluation)
7. Others (FixedNetworkElement, StaticGenerator, pure algebraic NonConformingLoad): skip
   (handled by the Y matrix and the network solve)
"""
function YBusSystemBuffers(
    ::Type{T},
    am::AddressManager,
    sys::SystemModel,
    formulation::FormulationTrait,
    tds_init;
    automata::Union{Nothing, AbstractDict{Symbol, <:AutomatonState}}=nothing
) where T<:Number
    # Every model declares its role, and none with instances is unsupported: a model is
    # refused by name before any role is read
    check_ybus_support(sys)

    # Compute generator index mapping FIRST (needed for closures)
    # This maps model_name → [model_idx → gen_buses_idx]
    model_to_gen_idx = compute_gen_index_mapping(sys)

    # Build status dict for active mask lookup in closures
    status = build_status_dict(sys)

    # Closure vectors
    generator_evaluators = Function[]
    internal_voltage_gatherers = Function[]
    controller_evaluators = Function[]
    workspaces_by_name = Dict{Symbol, Any}()
    metadata_by_name = Dict{Symbol, Any}()
    models_by_name = Dict{Symbol, Any}()
    externals_by_name = Dict{Symbol, Any}()
    automaton_models = Symbol[]

    # Metadata-driven discovery - NO hard-coding!
    @debug "YBusSystemBuffers: Processing models" model_names=collect(keys(sys.models)) am_models=collect(keys(am.assignments))
    for (model_name, model_vec) in pairs(sys.models)
        # Skip models with no address assignments
        if !haskey(am.assignments, model_name)
            @debug "  Skipping model (no address assignments)" model_name
            continue
        end

        # Get model type and role
        base_type = get_base_model_type(model_vec)
        base_metadata = model_metadata(base_type)
        role = network_partition_role(base_type)
        @debug "  Processing model" model_name role

        form_metadata = get_metadata_for_formulation(base_metadata, formulation)

        # Skip models with no variables and no residuals
        (isempty(form_metadata.vars) && isempty(form_metadata.residuals)) && continue

        # Create workspace
        ws = ResidualBuffers(T, form_metadata, am, model_vec, sys.array_type;
                             system=sys.settings,
                             automata=automata === nothing ? nothing : get(automata, model_name, nothing))
        if has_automata(form_metadata)
            # The automaton runtime evaluates a model's residuals with `calc_residuals!`:
            # a machine's Y-bus residuals (`calc_ybus_residuals!`) are not those
            role isa MachineWithInternalVoltage && throw(ArgumentError(
                "The Y-bus simulation (PowerfulYbus) does not realize the automata of $model_name " *
                "($(join((a.name for a in form_metadata.automata), ", "))), a machine with an internal " *
                "voltage: its residuals in the Y-bus formulation are calc_ybus_residuals!, which the " *
                "automaton runtime does not evaluate"))
            push!(automaton_models, model_name)
        end

        # Store for debugging and external access
        workspaces_by_name[model_name] = ws
        metadata_by_name[model_name] = form_metadata
        models_by_name[model_name] = model_vec
        externals_by_name[model_name] = getfield(sys.external, model_name)

        # Capture for closures
        model = model_vec
        ext = getfield(sys.external, model_name)
        meta = form_metadata
        form = formulation

        # Dispatch based on network partition role
        if role isa MachineWithInternalVoltage
            # Generator: uses calc_ybus_residuals! with network solution
            # Capture the mapping vector directly (computed at construction time)
            gen_idx_vec = model_to_gen_idx[model_name]

            # Pre-compute filtered scatter tuple that excludes bus injections (p_inj, q_inj)
            # Bus power balance is handled by network solve (I = Y*E) in Y-bus mode
            internal_scatter_tup = build_res_scatter_tuple(
                ws.res_addrs, ws.res_vals, form_metadata;
                skip = is_bus_injection
            )
            # Capture model_name for status dict lookup
            mn = model_name

            push!(generator_evaluators, (f, g, x, y, t, network_sol) -> begin
                active_mask = status[mn]  # Status dict lookup: ~2ns
                gather_variables!(ws, x, y)
                gather_residuals!(ws, f, g)
                # Captured Vector{Int} - O(1) access, no Dict lookup
                calc_ybus_residuals!(ws, model, ext, network_sol, form, gen_idx_vec; active_mask=active_mask)
                # Use pre-computed filtered tuple for zero-allocation scatter
                scatter_res_add!(internal_scatter_tup, f, g, active_mask)
                scatter_res_set!(internal_scatter_tup, f, g, active_mask)
            end)
            push!(internal_voltage_gatherers, _internal_voltage_gatherer(ws, model, meta, gen_idx_vec))

        elseif role isa NotInYBusFormulation
            # Controller: uses standard calc_residuals!, with its bus voltage variables
            # (e.g. SEXS.v) from the network's bus voltages
            push!(controller_evaluators, _controller_evaluator(ws, model, ext, meta, form, status, model_name, sys))

        elseif role isa NonConformingLoad
            # NC loads: current injection handled by NC solver in solve_network!
            # BUT: if NC load has state variables (e.g., REGCA1), we need to evaluate
            # its state residuals here. Pure algebraic NC loads (PQ) can be skipped.
            has_states = has_state_vars(form_metadata)

            if has_states
                # NC load with dynamics (e.g., REGCA1): evaluate state residuals
                # Current injection is still handled by NC solver separately
                @debug "  Adding NC load with dynamics to controller_evaluators" model_name
                push!(controller_evaluators, _controller_evaluator(ws, model, ext, meta, form, status, model_name, sys))
            end
            # else: pure NC load (PQ) - fully handled by NC solver, skip

        # else: FixedNetworkElement (in Y) or StaticGenerator (network solve
        # constraint, no variables of its own) - skip
        end
    end

    sys_ws = YBusSystemBuffers{T, typeof(formulation)}(
        formulation,
        generator_evaluators,
        internal_voltage_gatherers,
        controller_evaluators,
        workspaces_by_name,
        metadata_by_name,
        models_by_name,
        externals_by_name,
        model_to_gen_idx,
        automaton_models,
        status
    )
    sync_ybus_setpoints!(sys_ws, tds_init)
    return sys_ws
end

"""
    automaton_states(sys_ws::YBusSystemBuffers) -> Dict{Symbol, AutomatonState}

The `AutomatonState` of each model of `sys_ws` with automata, to share with another
`YBusSystemBuffers` (keyword `automata`).
"""
automaton_states(sys_ws::YBusSystemBuffers) =
    Dict{Symbol, AutomatonState}(name => sys_ws.workspaces_by_name[name].automata
                                 for name in sys_ws.automaton_models)

"""
    sync_automaton_memory!(dst::YBusSystemBuffers, src::YBusSystemBuffers)

Copy the memory of the automata of every model of `src` (the Float64 workspace, where
the simulation sets it) into the memory of `dst` (a workspace of another number type
sharing its automaton state, which holds memory of its own), as values without
derivatives.
"""
function sync_automaton_memory!(dst::YBusSystemBuffers, src::YBusSystemBuffers)
    for name in src.automaton_models
        dst_memory = dst.workspaces_by_name[name].automata.memory
        for (m, column) in pairs(src.workspaces_by_name[name].automata.memory)
            getfield(dst_memory, m) .= column
        end
    end
    return nothing
end

"""
    automaton_memory_synced(dst::YBusSystemBuffers, src::YBusSystemBuffers) -> Bool

Whether the memory of the automata of `dst` equals that of `src` as `sync_automaton_memory!`
leaves it (see the method for `SystemBuffers`).
"""
function automaton_memory_synced(dst::YBusSystemBuffers, src::YBusSystemBuffers)
    return all(src.automaton_models) do name
        all(map((d, s) -> all(i -> isequal(d[i], s[i]), eachindex(s)),
                values(dst.workspaces_by_name[name].automata.memory), values(src.workspaces_by_name[name].automata.memory)))
    end
end

"""
    copy_setpoints!(dst::YBusSystemBuffers, src::YBusSystemBuffers)

Copy the setpoints of every model of `src` into those of `dst` (a workspace of another number
type, built from the initialization, whose setpoints are those it had when `src`'s were
last changed by an event: `set_setpoint!`).
"""
function copy_setpoints!(dst::YBusSystemBuffers, src::YBusSystemBuffers)
    for (name, ws) in pairs(src.workspaces_by_name)
        dst_setpoints = dst.workspaces_by_name[name].setpoints
        for (setpoint, values) in pairs(ws.setpoints)
            copyto!(getproperty(dst_setpoints, setpoint), values)
        end
    end
    return nothing
end

"""
    setpoint_array(sys_ws::YBusSystemBuffers, model::Symbol, setpoint::Symbol) -> AbstractVector

The array of `setpoint` of `model` in `sys_ws`; an `ArgumentError` names the model or the
setpoint that does not exist.
"""
function setpoint_array(sys_ws::YBusSystemBuffers, model::Symbol, setpoint::Symbol)
    haskey(sys_ws.workspaces_by_name, model) || throw(ArgumentError(
        "SetpointChange of $model.$setpoint: $model has no equations in this problem"))
    setpoints = sys_ws.workspaces_by_name[model].setpoints
    hasproperty(setpoints, setpoint) || throw(ArgumentError(
        "SetpointChange of $model.$setpoint: $model has no setpoint $setpoint. Setpoints: " *
        join(keys(setpoints), ", ")))
    return getproperty(setpoints, setpoint)
end

# Closure setting the internal voltages E[k] of the in-service instances of one machine
# model from its gathered variables (k = gen_idx_vec[i], 0 when out of service)
function _internal_voltage_gatherer(ws, model, meta, gen_idx_vec::Vector{Int})
    instances = [(k, i) for (i, k) in enumerate(gen_idx_vec) if k > 0]
    return (E, x, y) -> begin
        gather_variables!(ws, x, y)
        @inbounds for (k, i) in instances
            E[k] = get_internal_voltage(ws, model, i)
        end
    end
end

# Closure evaluating one model's residuals with the standard gather/calc/scatter path;
# its bus voltage variables are set from the bus voltages V_bus after the gather
function _controller_evaluator(ws, model, ext, meta, form, status, model_name::Symbol, sys)
    bus_inputs = bus_voltage_inputs(ws, meta, sys, model_name)
    return (f, g, x, y, t, V_bus) -> begin
        active_mask = status[model_name]
        gather_variables!(ws, x, y)
        set_bus_voltages!(bus_inputs, V_bus)
        gather_residuals!(ws, f, g)
        calc_residuals!(ws, model, ext, form, t; active_mask=active_mask)
        scatter_residuals!(ws, f, g, active_mask)
    end
end

"""
    compute_gen_index_mapping(sys::SystemModel) -> Dict{Symbol, Vector{Int}}

Compute generator index mapping from model array indices to gen_buses indices.

The gen_buses ordering is determined by sorting the internal bus indices of the
in-service generators (`online_machines`). This matches the ordering used by
`FullYMatrix.gen_buses` and the NetworkSolution generator arrays.

# Returns
Dict mapping model_name → Vector{Int} where vec[model_idx] = gen_buses_idx
"""
function compute_gen_index_mapping(sys::SystemModel)
    # Every machine model gets an entry; out-of-service instances map to 0
    model_to_gen_idx = Dict{Symbol, Vector{Int}}()
    for (model_name, model_vec) in pairs(sys.models)
        network_partition_role(eltype(model_vec)) isa MachineWithInternalVoltage || continue
        model_to_gen_idx[model_name] = zeros(Int, length(model_vec))
    end

    # gen_buses ordering: in-service machines sorted by internal bus index
    machines = sort(online_machines(sys); by = m -> m[3])
    for (gen_buses_idx, (model_name, model_idx, _)) in enumerate(machines)
        model_to_gen_idx[model_name][model_idx] = gen_buses_idx
    end

    return model_to_gen_idx
end

"""
    evaluate_ybus_residuals!(sys_ws::YBusSystemBuffers, dae, network_sol, V_bus, t)

Evaluate all residuals using Y-bus network solution - ZERO ALLOCATION.

# Algorithm
1. Evaluate generators with their terminal voltages from the network solution
   (calc_ybus_residuals!)
2. Evaluate controllers with standard evaluation (calc_residuals!), their bus
   voltage variables set from `V_bus`

# Arguments
- `sys_ws`: Y-bus system workspace
- `dae`: DAE with (f,g,x,y) vectors (or just (g,y) for power flow)
- `network_sol`: NetworkSolution with generator currents
- `V_bus`: Bus voltages (complex) of the network solution
- `t`: Current time

# Performance
- Zero allocation after initial construction
- Type-stable dispatch via closures
- Each closure is fully specialized internally
"""
function evaluate_ybus_residuals!(
    sys_ws::YBusSystemBuffers{T},
    dae,
    network_sol::NetworkSolution{T},
    V_bus::AbstractVector{Complex{T}},
    t::Real=0.0
) where T
    # Extract vectors (supports both full DAE and algebraic-only)
    f = hasfield(typeof(dae), :f) ? dae.f : T[]
    g = dae.g
    x = hasfield(typeof(dae), :x) ? dae.x : T[]
    y = dae.y

    # Zero residuals
    # Use eltype() instead of T to support mixed Float64/Dual during autodiff
    fill!(f, zero(eltype(f)))
    fill!(g, zero(eltype(g)))

    # Evaluate generators with network solution
    for eval_fn in sys_ws.generator_evaluators
        eval_fn(f, g, x, y, t, network_sol)
    end

    # Evaluate controllers with standard path (bus voltages from V_bus)
    for eval_fn in sys_ws.controller_evaluators
        eval_fn(f, g, x, y, t, V_bus)
    end

    nothing
end

"""
    sync_ybus_setpoints!(sys_ws::YBusSystemBuffers, tds_init)

Synchronize setpoints in Y-bus workspaces from TDS initialization data.

After TDS initialization computes setpoints (E_q0, P_m0, etc.), this function
copies them from `tds_init.dyn` into the Y-bus workspace setpoint arrays.

# Arguments
- `sys_ws`: YBusSystemBuffers with workspaces containing setpoint arrays
- `tds_init`: TDS initialization data (from `tds_prob.init`)

# Algorithm
For each model workspace:
1. Get `initialized_params` from metadata (defines var → setpoint mapping)
2. Read initial variable values from `tds_init.dyn[model_name]`
3. Copy to workspace setpoint arrays

Called by the `YBusSystemBuffers` constructor.
"""
function sync_ybus_setpoints!(sys_ws::YBusSystemBuffers, tds_init)
    # Iterate through all workspaces
    for (model_name, ws) in pairs(sys_ws.workspaces_by_name)
        # Skip if no initialization data for this model
        haskey(tds_init.dyn, model_name) || continue

        # Get initialization data for this model
        init_data = tds_init.dyn[model_name]

        # Get metadata to find initialized_params mapping
        meta = sys_ws.metadata_by_name[model_name]

        # Copy each initialized parameter from init_data to setpoints
        for init_param in meta.initialized_params
            setpoint_name = init_param.name        # e.g., :E_q0
            var_name = init_param.source_var       # e.g., :E_q

            # Read initial value from init_data
            if haskey(init_data, var_name)
                init_values = getproperty(init_data, var_name)

                # Write to workspace setpoints
                if haskey(ws.setpoints, setpoint_name)
                    setpoint_array = getproperty(ws.setpoints, setpoint_name)
                    copyto!(setpoint_array, init_values)
                end
            end
        end
    end

    nothing
end

# ============================================================================
# Y-Bus Specific Helpers
# ============================================================================

"""
    is_bus_injection(res_meta) -> Bool

Check if a residual is a bus power/current injection.

Bus injections (p_inj, q_inj) scatter to Bus power balance residuals, which are
handled implicitly by the network solve (I = Y*E) in Y-bus formulation.

Used as a predicate for `build_res_scatter_tuple` to filter out bus injections
at construction time, enabling zero-allocation scatter in the hot path.

# Returns
- `true` if `source_model === :bus_group` (network injection to bus)
- `false` otherwise (internal residual or control signal to other models)

# Examples
- `GENCLS.p_inj` → `source_model = :bus_group` → true (skip in Y-bus)
- `TGOV1.tm` → `source_model = :generators` → false (keep - governor control)
- `GENCLS.g_P_m` → `source_model = nothing` → false (keep - internal)

# Usage
```julia
# Pre-compute filtered scatter tuple at construction time
internal_scatter_tup = build_res_scatter_tuple(
    ws.res_addrs, ws.res_vals, metadata;
    skip = is_bus_injection
)

# Hot path: zero-allocation scatter using pre-computed tuple
scatter_res_add!(internal_scatter_tup, f, g, active_mask)
scatter_res_set!(internal_scatter_tup, f, g, active_mask)
```
"""
@inline is_bus_injection(res_meta) = res_meta.source_model === :bus_group

export YBusSystemBuffers, evaluate_ybus_residuals!, copy_setpoints!, setpoint_array
export compute_gen_index_mapping
export is_bus_injection
