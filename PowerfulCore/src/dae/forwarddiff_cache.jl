"""
ForwardDiff-based autodiff cache with Dual workspace caching.

This module provides the ForwardDiff backend implementation:
- ForwardDiffCache: Caches SystemBuffers{Dual} for zero-allocation after warmup
- get_tds_dual_buffers: TDS-specific Dual buffer caching
"""

export ForwardDiffCache, get_tds_dual_buffers, get_typed_workspace!

# ============================================================================
# ForwardDiffCache
# ============================================================================

"""
    ForwardDiffCache{F} <: AbstractADCache{F}

ForwardDiff-based autodiff cache with Dual workspace caching.

Pre-allocates Float64 workspace and caches Dual workspaces on demand.
This avoids the overhead of creating new workspaces on every Jacobian
evaluation (typically ~1ms per allocation for large systems).

**Zero allocation**: After warmup (first call per unique Dual type)

# Type Parameters
- `F`: FormulationTrait type (PowerInjection or CurrentInjection)

# Fields
- `core::ADCacheCore{F}`: Shared state (sys, formulation, address_manager, float_ws)
- `dual_workspaces::Dict{DataType, SystemBuffers}`: Cached Dual workspaces
- `tds_dual_buffers::Dict{DataType, NamedTuple}`: TDS-specific Dual buffers
- `empty_vecs::Dict{DataType, Tuple{Vector, Vector}}`: Cached empty (f, x) vectors per Dual type
- `parameter_fields::Dict{Symbol, Vector{Symbol}}`: Model fields that every Dual workspace
  reads from its own copy of its element type (see `register_parameter_fields!`)

# Example
```julia
# Create via factory function (recommended)
adcache = create_adcache(am, sys, formulation, ForwardDiffBackend())

# Use as callable
adcache(g, y)        # Power flow
adcache(f, g, x, y)  # TDS
```
"""
mutable struct ForwardDiffCache{F<:FormulationTrait} <: AbstractADCache{F}
    core::ADCacheCore{F}
    dual_workspaces::Dict{DataType, SystemBuffers}
    tds_dual_buffers::Dict{DataType, NamedTuple}
    empty_vecs::Dict{DataType, Tuple{Vector, Vector}}
    parameter_fields::Dict{Symbol, Vector{Symbol}}
end

# ============================================================================
# Factory Function
# ============================================================================

"""
    create_adcache(am, sys, formulation, ::ForwardDiffBackend)

Create ForwardDiffCache for ForwardDiff backend.
"""
function create_adcache(am::AddressManager, sys, formulation::F,
                        ::ForwardDiffBackend) where {F<:FormulationTrait}
    core = ADCacheCore(am, sys, formulation)
    ForwardDiffCache{F}(
        core,
        Dict{DataType, SystemBuffers}(),
        Dict{DataType, NamedTuple}(),
        Dict{DataType, Tuple{Vector, Vector}}(),
        Dict{Symbol, Vector{Symbol}}()
    )
end

"""
    get_backend(::ForwardDiffCache)

Returns ForwardDiffBackend() for ForwardDiffCache.
"""
get_backend(::ForwardDiffCache) = ForwardDiffBackend()

# ============================================================================
# Callable Interfaces
# ============================================================================

"""
    (cache::ForwardDiffCache)(g::AbstractVector, y::AbstractVector, status, t)

Evaluate algebraic residuals (power flow interface) at time `t`, the time of the state `y`.

Automatically dispatches to Float64 or Dual number path based on input type.
"""
function (cache::ForwardDiffCache)(g::AbstractVector, y::AbstractVector,
                                   status::Dict{Symbol, <:AbstractVector{Bool}}, t::Real)
    T = eltype(g)
    # Get or create cached empty vectors for this type (allocates once per unique T)
    (f, x) = get!(cache.empty_vecs, T) do
        (Vector{T}(undef, 0), Vector{T}(undef, 0))
    end
    # Delegate to unified 6-arg interface
    cache(f, g, x, y, status, t)
    return nothing
end

"""
    (cache::ForwardDiffCache)(f, g, x, y, status, t; recorder=nothing)

Evaluate DAE residuals (both differential and algebraic) at time `t`, the time of the state
`(x, y)`: a model's residuals can depend on it directly, not only through the state.

Automatically dispatches to Float64 or Dual number path based on input type.
Caches Dual workspaces on first use to avoid reallocation overhead.
"""
function (cache::ForwardDiffCache)(
    f::AbstractVector,
    g::AbstractVector,
    x::AbstractVector,
    y::AbstractVector,
    status::Dict{Symbol, <:AbstractVector{Bool}},
    t::Real;
    recorder = nothing
)
    # Determine workspace type from algebraic output g
    T = eltype(g)

    if T == Float64
        # Float64 path: Use cached workspace
        evaluate_system_residuals!(cache.core.float_ws, f, g, x, y, status, t; recorder=recorder)
    else
        # Dual number path (no recording - used in Jacobian computation)
        evaluate_system_residuals!(get_typed_workspace!(cache, T), f, g, x, y, status, t)
    end

    return nothing
end

# ============================================================================
# Typed Workspaces
# ============================================================================

"""
    get_typed_workspace!(cache::ForwardDiffCache, ::Type{T}) -> SystemBuffers

Workspace for element type `T`. For `Float64`, the cache's Float64 workspace. For
Dual types, a workspace cached per `T`, created on first use with every setpoint
copied from the Float64 workspace, and the modes and discrete quantities of its automata
shared with their memory copied (see
`AutomatonState(state, T)`), so Dual evaluation computes the same residual function
(setpoints such as time constants or references are not left at zero). Its
evaluators read the registered parameter fields (`register_parameter_fields!`) from
copies of element type `T`, initialized from the system.
"""
get_typed_workspace!(cache::ForwardDiffCache, ::Type{Float64}) = cache.core.float_ws

function get_typed_workspace!(cache::ForwardDiffCache, ::Type{T}) where T
    get!(() -> _typed_workspace(T, cache.core, cache.parameter_fields), cache.dual_workspaces, T)
end

# The workspace of element type `T` beside the Float64 workspace of `core`
function _typed_workspace(::Type{T}, core::ADCacheCore, typed_fields::Dict{Symbol, Vector{Symbol}}) where T
    float_ws = core.float_ws
    ws = SystemBuffers(T, core.address_manager, core.sys, core.formulation;
                       hold = float_ws.hold, hold_x = float_ws.hold_x, typed_fields,
                       automata = automaton_states(float_ws))
    _copy_setpoints!(ws, float_ws)
    return ws
end

"""
    sync_automaton_memory!(cache::ForwardDiffCache)

Copy the memory of the automata from the Float64 workspace, where the simulation sets
it, into every Dual workspace of `cache` (which hold their own, see
`AutomatonState(state, T)`), as values without derivatives. The simulation calls it
where a start or a transition changes memory; a Dual workspace built later starts from
the Float64 memory.
"""
function sync_automaton_memory!(cache::ForwardDiffCache)
    for ws in values(cache.dual_workspaces)
        sync_automaton_memory!(ws, cache.core.float_ws)
    end
    return nothing
end

"""
    _copy_setpoints!(dst::SystemBuffers, src::SystemBuffers)

Copy every setpoint array from `src` into `dst` (converting element type).
Both must be built from the same address manager and system.
"""
function _copy_setpoints!(dst::SystemBuffers, src::SystemBuffers)
    for (model_name, src_ws) in src.workspaces_by_name
        dst_setpoints = dst.workspaces_by_name[model_name].setpoints
        for (name, src_arr) in pairs(src_ws.setpoints)
            getfield(dst_setpoints, name) .= src_arr
        end
    end
    nothing
end

"""
    register_parameter_fields!(cache::ForwardDiffCache, model_name::Symbol, fields)

Give the fields `fields` of model `model_name` a copy of element type `T` in every Dual
workspace `T` of `cache`, so that parameter values of that type (e.g. Dual numbers of a
parameter Jacobian) can be written there without changing the system. Dual workspaces
built before are dropped and rebuilt on first use.

The copies take their values from the system when a workspace is built; write the
fields through `scatter_parameters!`, which updates the copies too.
"""
function register_parameter_fields!(cache::ForwardDiffCache, model_name::Symbol,
                                    fields::AbstractVector{Symbol})
    registered = get!(cache.parameter_fields, model_name, Symbol[])
    new_fields = setdiff(fields, registered)
    isempty(new_fields) && return nothing
    append!(registered, new_fields)
    invalidate_dual_workspaces!(cache)
    return nothing
end

# ============================================================================
# TDS Dual Buffers
# ============================================================================

"""
    get_tds_dual_buffers(cache::ForwardDiffCache, ::Type{T}, n_diff::Int, n_alg::Int)

Get or create TDS-specific Dual buffers for type T.

Returns a NamedTuple with (x_buf, y_buf, f_buf, g_buf) of type Vector{T}.
Buffers are cached for reuse across ForwardDiff evaluations.
"""
function get_tds_dual_buffers(cache::ForwardDiffCache, ::Type{T}, n_diff::Int, n_alg::Int) where T
    get!(cache.tds_dual_buffers, T) do
        (
            x_buf = Vector{T}(undef, n_diff),
            y_buf = Vector{T}(undef, n_alg),
            f_buf = Vector{T}(undef, n_diff),
            g_buf = Vector{T}(undef, n_alg)
        )
    end
end

# ============================================================================
# Dual Workspace Invalidation
# ============================================================================

"""
    invalidate_dual_workspaces!(adcache::ForwardDiffCache)

Invalidate cached Dual workspaces when Float64 setpoints change.

Called after Float64 workspace setpoints are updated (e.g. by a parameter system).
Next Dual evaluation will create fresh workspace from updated Float64 values.
"""
function invalidate_dual_workspaces!(adcache::ForwardDiffCache)
    empty!(adcache.dual_workspaces)
    nothing
end
