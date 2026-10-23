"""
    InjectionRecorder - Track per-model contributions to SharedRes addresses

This module provides a way to record which models contribute what values to
shared residual equations (like bus power/current balance) during TDS simulation.

# Usage
```julia
# 1. Define what to monitor (residuals p, q of buses 1 and 2)
monitor = InjectionMonitor(:Bus => [:p => [(1,), (2,)], :q => [(1,), (2,)]])

# 2. Create recorder (builds address mappings)
recorder = InjectionRecorder(monitor, sys_ws, am)

# 3. Pass recorder to evaluate_system_residuals! during TDS
evaluate_system_residuals!(sys_ws, f, g, x, y, status, t; recorder=recorder)

# 4. Retrieve results
table = contributions_table(recorder, :Bus, :p, (1,))
```
"""

# Import SharedRes from Foundation (parent module)
using ..Foundation: SharedRes, AddressableRes, is_shared_res

# get_addresses is already available within DAEInfra (no import needed)

# ============================================================================
# Types (Fully Parameterized for Type Stability)
# ============================================================================

"""
    InjectionMonitor{K}

Specifies which SharedRes addresses to monitor.

# Type Parameters
- `K`: External key type (e.g., `Int` for bus numbers)

# Example
```julia
# Monitor the active and reactive power balance of buses 1, 2, 3 (K = Tuple{Int})
monitor = InjectionMonitor(
    :Bus => [:p => [(1,), (2,), (3,)], :q => [(1,), (2,), (3,)]]
)
```
"""
struct InjectionMonitor{K}
    targets::Dict{Symbol, Vector{Pair{Symbol, Vector{K}}}}
end

# Convenience constructor - infers K from first key
function InjectionMonitor(pairs::Pair{Symbol}...)
    # Infer key type from first pair's first var's first key
    first_pair = first(pairs)
    first_var_spec = first(first_pair.second)
    K = eltype(first_var_spec.second)

    targets = Dict{Symbol, Vector{Pair{Symbol, Vector{K}}}}()
    for (model, var_specs) in pairs
        targets[model] = [Symbol(v.first) => collect(K, v.second) for v in var_specs]
    end
    return InjectionMonitor{K}(targets)
end

"""
    ContributorInfo

Information about a model that contributes to a monitored address.
"""
struct ContributorInfo
    model_name::Symbol
    res_name::Symbol
    local_indices::Vector{Int}
end

"""
    AddressMapping{K}

Pre-computed mapping for efficient recording.

# Type Parameters
- `K`: External key type
"""
struct AddressMapping{K}
    contributors::Dict{Tuple{Symbol, Symbol, K}, Vector{ContributorInfo}}
    dae_addresses::Dict{Tuple{Symbol, Symbol, K}, UInt64}
end

"""
    InjectionRecorder{T, K}

Records per-model contributions to monitored SharedRes addresses.

Uses a two-phase recording approach:
1. `record_model_contributions!` caches values in `pending` (called many times per step)
2. `commit_contributions!` pushes pending to timeseries (called once per accepted step)

# Type Parameters
- `T`: Numeric type (Float64, etc.)
- `K`: External key type (Int for bus numbers, etc.)

# Fields
- `monitor`: What addresses are being monitored
- `address_map`: Pre-computed contributor mappings
- `pending`: Cached values from latest evaluation (not yet committed)
- `contributions`: Time series storage (target => model => values)
- `timestamps`: Time points recorded (one per accepted step)
"""
mutable struct InjectionRecorder{T, K}
    monitor::InjectionMonitor{K}
    address_map::AddressMapping{K}
    pending::Dict{Tuple{Symbol, Symbol, K}, Dict{Symbol, T}}  # Cached values
    contributions::Dict{Tuple{Symbol, Symbol, K}, Dict{Symbol, Vector{T}}}
    timestamps::Vector{T}
end

# ============================================================================
# Setup Functions
# ============================================================================

"""
    InjectionRecorder(monitor::InjectionMonitor{K}, sys_ws, sys, am) where K

Advanced constructor, for a caller that already has a `SystemBuffers` and does not want
one built from a `TDSProblem`. Most callers use [`InjectionRecorder`](@ref)`(monitor,
prob::TDSProblem)` instead.

# Extended help
- `monitor`: InjectionMonitor specifying which addresses to track
- `sys_ws`: SystemBuffers (provides model workspaces and metadata)
- `sys`: SystemModel (provides index order for key lookups)
- `am`: AddressManager (provides DAE address mappings)
"""
function InjectionRecorder(
    monitor::InjectionMonitor{K},
    sys_ws,  # SystemBuffers
    sys,     # SystemModel (for index lookups)
    am       # AddressManager
) where K
    T = Float64

    # Build address mapping
    address_map = build_address_mapping(monitor, sys_ws, sys, am)

    # Initialize pending cache and contributions storage
    pending = Dict{Tuple{Symbol, Symbol, K}, Dict{Symbol, T}}()
    contributions = Dict{Tuple{Symbol, Symbol, K}, Dict{Symbol, Vector{T}}}()

    for key in keys(address_map.contributors)
        pending[key] = Dict{Symbol, T}()
        contributions[key] = Dict{Symbol, Vector{T}}()
        for contrib in address_map.contributors[key]
            pending[key][contrib.model_name] = zero(T)
            contributions[key][contrib.model_name] = T[]
        end
    end

    return InjectionRecorder{T, K}(
        monitor,
        address_map,
        pending,
        contributions,
        T[]
    )
end

# NOTE: The InjectionRecorder(monitor, prob::TDSProblem) convenience constructor
# is defined in PowerfulCore.jl after SystemBuffers is available, to avoid
# circular dependencies between DAEInfra and System modules.

"""
    build_address_mapping(monitor, sys_ws, sys, am) -> AddressMapping{K}
"""
function build_address_mapping(
    monitor::InjectionMonitor{K},
    sys_ws,
    sys,  # SystemModel for index lookups
    am
) where K
    contributors = Dict{Tuple{Symbol, Symbol, K}, Vector{ContributorInfo}}()
    dae_addresses = Dict{Tuple{Symbol, Symbol, K}, UInt64}()

    for (target_model, var_specs) in monitor.targets
        for (target_var, ext_keys) in var_specs
            target_type = monitored_residual_type(sys_ws, target_model, target_var)
            target_addrs = get_addresses(am, target_type, target_model, target_var)

            for ext_key in ext_keys
                int_idx = get_index_by_key(sys, target_model, ext_key)

                dae_addr = target_addrs[int_idx]
                key = (target_model, target_var, ext_key)
                dae_addresses[key] = dae_addr
                contributors[key] = ContributorInfo[]

                find_contributors!(contributors[key], target_type, dae_addr, sys_ws)
            end
        end
    end

    return AddressMapping{K}(contributors, dae_addresses)
end

"""
    monitored_residual_type(sys_ws, model, residual) -> AddressableRes

Address space (`AlgebRes()` or `StateRes()`) of the monitored residual `residual`
of `model`. Errors if `model` has no such residual.
"""
function monitored_residual_type(sys_ws, model::Symbol, residual::Symbol)
    haskey(sys_ws.metadata_by_name, model) ||
        error("Cannot monitor $model: no such model in the system")
    residuals = sys_ws.metadata_by_name[model].residuals
    idx = findfirst(r -> r.name === residual, residuals)
    idx === nothing && error("Cannot monitor $model.$residual: $model has no residual " *
        "$residual. Residuals: $(join([r.name for r in residuals], ", "))")
    return residuals[idx].address_type
end

"""
    find_contributors!(result, target_type, dae_addr, sys_ws)

Append to `result` every `SharedRes(:add)` residual (of any model) that adds into
address `dae_addr` of the address space `target_type`. Contributors are matched by
address, so residuals that reach the equation through a group reference (e.g. a
generator's injection through `:bus_group`) are found too.
"""
function find_contributors!(
    result::Vector{ContributorInfo},
    target_type::AddressableRes,
    dae_addr::UInt64,
    sys_ws
)
    for (model_name, ws) in sys_ws.workspaces_by_name
        metadata = sys_ws.metadata_by_name[model_name]

        for res in metadata.residuals
            is_shared_res(:add)(res) || continue
            res.address_type === target_type || continue  # separate address spaces

            res_addrs = getfield(ws.res_addrs, res.name)
            local_idxs = findall(==(dae_addr), res_addrs)

            if !isempty(local_idxs)
                push!(result, ContributorInfo(model_name, res.name, local_idxs))
            end
        end
    end
end



# ============================================================================
# Recording Functions (Hot Path - Must Be Type Stable)
# ============================================================================

"""
    record_model_contributions!(recorder, model_name, ws, t, active_mask)

Cache contributions from a model's workspace into the pending buffer.

This function is called at every `evaluate_system_residuals!` invocation
(including Newton iterations). Values are cached in `pending` and will be
committed to the timeseries when `commit_contributions!` is called at
accepted steps.

# Arguments
- `recorder`: InjectionRecorder to store contributions
- `model_name`: Symbol identifying the model
- `ws`: ResidualBuffers workspace
- `t`: Current time
- `active_mask`: Vector{Bool} indicating active status (passed from caller)
"""
@inline function record_model_contributions!(
    recorder::InjectionRecorder{T, K},
    model_name::Symbol,
    ws,
    t::Real,
    active_mask::AbstractVector{Bool}
) where {T, K}
    # Cache contributions in pending buffer (overwrites on each call). A model's
    # contribution is the sum over all of its residuals that add into the address.
    for (key, contrib_list) in recorder.address_map.contributors
        found = false
        total = zero(T)
        for contrib in contrib_list
            contrib.model_name == model_name || continue
            found = true

            res_vals = getfield(ws.res_vals, contrib.res_name)
            @inbounds for idx in contrib.local_indices
                if active_mask[idx]
                    total += res_vals[idx]
                end
            end
        end

        # Cache in pending (will be committed at accepted step)
        found && (recorder.pending[key][model_name] = total)
    end
    nothing
end

"""
    commit_contributions!(recorder, t)

Commit cached pending values to the timeseries at time `t`.

This should be called once per accepted step (via a callback), not during
every residual evaluation. This ensures the timeseries contains exactly
one entry per accepted step with the converged values.
"""
function commit_contributions!(recorder::InjectionRecorder{T, K}, t::Real) where {T, K}
    # Push timestamp
    push!(recorder.timestamps, T(t))

    # Commit all pending values to timeseries
    for (key, model_dict) in recorder.pending
        for (model_name, value) in model_dict
            push!(recorder.contributions[key][model_name], value)
        end
    end
    nothing
end

# ============================================================================
# Retrieval Functions
# ============================================================================

"""
    get_contributions(recorder, target_model, target_var, ext_key) -> Dict{Symbol, Vector{T}}

Retrieve all model contributions to a monitored address.

Returns a Dict mapping model names to their contribution time series.
Throws if address is not monitored.
"""
function get_contributions(
    recorder::InjectionRecorder{T, K},
    target_model::Symbol,
    target_var::Symbol,
    ext_key::K
) where {T, K}
    key = (target_model, target_var, ext_key)
    if !haskey(recorder.contributions, key)
        error("Address ($target_model, $target_var, $ext_key) not monitored")
    end
    return recorder.contributions[key]
end

"""
    get_timeseries(recorder, target_model, target_var, ext_key, model_name) -> (timestamps, values)

Retrieve a specific model's contribution timeseries to a monitored address.

Returns tuple of (timestamps, values) for the specified contributing model.
Throws if address not monitored or model doesn't contribute to it.
"""
function get_timeseries(
    recorder::InjectionRecorder{T, K},
    target_model::Symbol,
    target_var::Symbol,
    ext_key::K,
    model_name::Symbol
) where {T, K}
    contribs = get_contributions(recorder, target_model, target_var, ext_key)
    if !haskey(contribs, model_name)
        error("Model $model_name does not contribute to ($target_model, $target_var, $ext_key)")
    end
    return (recorder.timestamps, contribs[model_name])
end

"""
    contributions_table(recorder, target_model, target_var, ext_key)

Returns NamedTuple with :time, one field per model, and :total
"""
function contributions_table(
    recorder::InjectionRecorder{T, K},
    target_model::Symbol,
    target_var::Symbol,
    ext_key::K
) where {T, K}
    contribs = get_contributions(recorder, target_model, target_var, ext_key)

    n = length(recorder.timestamps)
    total = zeros(T, n)

    cols = Pair{Symbol, Vector{T}}[]
    push!(cols, :time => copy(recorder.timestamps))

    for (model_name, values) in contribs
        push!(cols, model_name => copy(values))
        total .+= values
    end

    push!(cols, :total => total)
    return NamedTuple(cols)
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    list_monitored_addresses(recorder::InjectionRecorder) -> Vector

Return the monitored targets of `recorder` as `(model, variable, external_key)`
tuples, for example `(:Bus, :p, (1,))`.
"""
function list_monitored_addresses(recorder::InjectionRecorder)
    return collect(keys(recorder.address_map.dae_addresses))
end

"""
    list_contributors(recorder::InjectionRecorder, target_model::Symbol,
                      target_var::Symbol, ext_key) -> Vector{Symbol}

Return the names of the models that contribute to the monitored residual
`target_var` of instance `ext_key` of `target_model`. Errors if that target is not
monitored.
"""
function list_contributors(
    recorder::InjectionRecorder{T, K},
    target_model::Symbol,
    target_var::Symbol,
    ext_key::K
) where {T, K}
    key = (target_model, target_var, ext_key)
    if !haskey(recorder.address_map.contributors, key)
        error("Address ($target_model, $target_var, $ext_key) not monitored")
    end
    return unique!([c.model_name for c in recorder.address_map.contributors[key]])
end

