# ============================================================================
# Parameters: model fields and setpoints written by element type
# ============================================================================
#
# A parameter vector p maps, slice by slice, to a field of the instances of one model:
# a model field (e.g. the exciter gain SEXS.K), a setpoint of the model's workspace
# (e.g. the governor reference TGOV1.pref0) or a memory of its automata (e.g. the time a
# timer was armed, which moves with the parameters that move that time). Float64 values are written where the
# Float64 evaluation reads them, the system's model arrays and the Float64 workspace;
# values of another element type (Dual numbers of a Jacobian with respect to p) are
# written into the workspace of that type, whose evaluators read registered model
# fields from their own typed copy (`register_parameter_fields!`). The system is never
# rebuilt, so the Float64 evaluation keeps its concrete field types.
# ============================================================================

"""
    ParameterSlice

The positions `range` of a parameter vector, mapped to the field `field` of the
instances `instances` of model `model`: a model field (`source = :model`, e.g. the
exciter gain `SEXS.K`), a setpoint of the model's workspace (`source = :setpoint`,
e.g. the governor reference `TGOV1.pref0`) or a memory of the model's automata
(`source = :memory`, see `AutomatonState`).

Build it with `ParameterSlice(cache, model, field, instances, first; source)`, which
checks the model, field and instances against `cache`.
"""
struct ParameterSlice{I<:AbstractVector{Int}}
    model::Symbol
    field::Symbol
    source::Symbol
    instances::I
    range::UnitRange{Int}
end

"""
    ParameterSlice(cache, model, field, instances, first; source=:model) -> ParameterSlice

The slice of a parameter vector starting at position `first` that holds `field` of the
`instances` of `model` in `cache`. Errors if the model has no equations in the cache,
the field or setpoint does not exist or does not hold real numbers, or an instance is
out of range.
"""
function ParameterSlice(cache::AbstractADCache, model::Symbol, field::Symbol,
                        instances::AbstractVector{<:Integer}, first::Int;
                        source::Symbol = :model)
    source in (:model, :setpoint, :memory) ||
        throw(ArgumentError("Invalid source :$source. Must be :model, :setpoint or :memory"))
    float_ws = get_float_ws(cache)
    haskey(float_ws.workspaces_by_name, model) ||
        error("No workspace found for model :$model: it has no equations in this problem")
    values = _parameter_array(float_ws, model, field, source)
    eltype(values) <: AbstractFloat ||
        error(":$model.$field holds $(eltype(values)), not real numbers")
    idx = collect(Int, instances)
    all(i -> 1 <= i <= length(values), idx) ||
        throw(ArgumentError("Instances $instances out of range for :$model.$field " *
                            "($(length(values)) instances)"))
    array_type = get_system(cache).array_type
    device_idx = is_gpu_array_type(array_type) ? array_type{Int}(idx) : idx
    return ParameterSlice(model, field, source, device_idx, first:(first + length(idx) - 1))
end

"""
    _parameter_array(ws::SystemBuffers, model, field, source) -> AbstractVector

The array that the evaluators of `ws` read `field` of `model` from: a column of the
model arrays (`:model`), a setpoint of the model's workspace (`:setpoint`) or a memory
of its automata (`:memory`).
"""
function _parameter_array(ws::SystemBuffers, model::Symbol, field::Symbol, source::Symbol)
    if source === :model
        arrays = ws.models_by_name[model]
        hasproperty(arrays, field) || error("Model :$model has no field :$field")
        return getproperty(arrays, field)
    elseif source === :memory
        memory = ws.workspaces_by_name[model].automata.memory
        hasproperty(memory, field) ||
            error("Memory :$field not found for :$model. Available memory: $(keys(memory))")
        return getfield(memory, field)
    else
        setpoints = ws.workspaces_by_name[model].setpoints
        hasproperty(setpoints, field) ||
            error("Setpoint :$field not found for :$model. Available setpoints: $(keys(setpoints))")
        return getfield(setpoints, field)
    end
end

_parameter_array(ws::SystemBuffers, s::ParameterSlice) = _parameter_array(ws, s.model, s.field, s.source)

"""
    parameter_array(cache::ForwardDiffCache, slice::ParameterSlice, ::Type{T}) -> AbstractVector

The array holding the values of `slice` for evaluations with element type `T`: for
`Float64`, the system's model field or the setpoint of the Float64 workspace; for
another `T`, the typed copy in the workspace of type `T` (built on first use). A model
field has a typed copy only once registered (`register_parameters!`).
"""
parameter_array(cache::ForwardDiffCache, s::ParameterSlice, ::Type{T}) where T =
    _parameter_array(get_typed_workspace!(cache, T), s)

"""
    register_parameters!(cache::ForwardDiffCache, slices)

Register the model fields of `slices` (`register_parameter_fields!`), so that values of
any element type can be written to them with `scatter_parameters!`. Setpoints and memory
need no registration: every workspace has its own.
"""
function register_parameters!(cache::ForwardDiffCache, slices::AbstractVector{<:ParameterSlice})
    for s in slices
        s.source === :model && register_parameter_fields!(cache, s.model, [s.field])
    end
    return nothing
end

"""
    gather_parameters(cache, slices) -> AbstractVector{Float64}

The current values of the parameter vector described by `slices`, from the Float64
evaluation (the system's model fields and the Float64 workspace's setpoints).
"""
function gather_parameters(cache::AbstractADCache, slices::AbstractVector{<:ParameterSlice})
    float_ws = get_float_ws(cache)
    isempty(slices) && return Float64[]
    n = sum(s -> length(s.range), slices)
    p = similar(_parameter_array(float_ws, first(slices)), Float64, n)
    for s in slices
        gather!(view(p, s.range), _parameter_array(float_ws, s), s.instances)
    end
    return p
end

"""
    scatter_parameters!(cache::ForwardDiffCache, slices, p::AbstractVector{T})

Write the parameter vector `p` into the arrays that evaluations with element type `T`
read (`parameter_array`). Float64 values also update the typed copies in the Dual
workspaces already built, so every workspace evaluates the same residual function.
"""
function scatter_parameters!(cache::ForwardDiffCache, slices::AbstractVector{<:ParameterSlice},
                             p::AbstractVector{T}) where T
    for s in slices
        scatter_unmasked!(parameter_array(cache, s, T), view(p, s.range), s.instances)
    end
    if T === Float64
        for ws in values(cache.dual_workspaces), s in slices
            scatter_unmasked!(_parameter_array(ws, s), view(p, s.range), s.instances)
        end
    end
    return nothing
end
