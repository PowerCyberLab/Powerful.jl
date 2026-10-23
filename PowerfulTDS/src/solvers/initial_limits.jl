# ============================================================================
# Initial Values Outside Their Limits
# ============================================================================
#
# The bounds a model declares as `VariableLimit`s (what its automata enforce) are
# checked against the initialized states before a simulation starts.

export enforce_initial_limits!


# The values of a bound of a `VariableLimit` for the instances of `models`, whose
# initialized variables are `initial`: nothing, the field it names, or what its function
# computes (one number, or one per instance) from the model arrays, or from them and the
# initialized variables where it takes both
_bound_values(::Nothing, models, initial) = nothing
_bound_values(field::Symbol, models, initial) = getproperty(models, field)
_bound_values(f::Function, models, initial) = applicable(f, models, initial) ? f(models, initial) : f(models)
_at(values::Number, k) = values
_at(values::AbstractVector, k) = values[k]

# Whether a `VariableLimit` applies at instance `k`: where its `when` (nothing: everywhere) holds
_applies(::Nothing, models, initial, k) = true
_applies(f::Function, models, initial, k) = _at(f(models, initial), k)

# A value as the messages show it: the value of a Dual (the sensitivity runs seed the limits with the parameters)
_shown(x::Real) = round(ForwardDiff.value(x), sigdigits = 6)

"""
    enforce_initial_limits!(sys, dyn, adjust::Bool) -> Vector{NamedTuple}

Check the initial value of every limited variable (from `dyn`, the initialized states)
against its bounds, for the in-service instances of `sys`: those the model declares as
`VariableLimit`s (`variable_limits`).

A value outside its bounds is an operating point the limits do not allow: clamping it
would move the system away from the power flow solution. With `adjust`, the violated
bound is moved to the initial value, so the simulation starts in steady state with the
device at its limit, and a warning lists every adjustment. Otherwise an
`ArgumentError` lists the violations.

Returns the violations as `(model, instance, variable, value, bound, limit)` named
tuples, with `limit` the bound before any adjustment.

    enforce_initial_limits!(initialize, sys, adjust::Bool) -> result

`initialize()` computes the initialization (an object with `dyn`, the initialized states).
Where a bound was moved, the start is computed again from the data as moved: what was
computed from a bound (a reference clamped to the limits that the start then violated) is
that of the new limits. The limits then hold at the start, which is checked.
"""
function enforce_initial_limits!(initialize::Function, sys::SystemModel, adjust::Bool)
    result = initialize()
    isempty(enforce_initial_limits!(sys, result.dyn, adjust)) && return result
    result = initialize()
    enforce_initial_limits!(sys, result.dyn, false)
    return result
end

function enforce_initial_limits!(sys::SystemModel, dyn::NamedTuple, adjust::Bool)
    violations = NamedTuple[]
    for (model_name, models) in pairs(sys.models)
        isempty(models) && continue
        declared = variable_limits(model_metadata(get_base_model_type(models)))
        isempty(declared) && continue
        status_field = get_status_field(models)
        in_service = isnothing(status_field) ? trues(length(models)) : getproperty(models, status_field)
        instances = get_key(models)
        for limit in declared
            x0 = dyn[model_name][limit.variable]
            initial = dyn[model_name]
            lower, upper = _bound_values(limit.lower, models, initial), _bound_values(limit.upper, models, initial)
            for k in eachindex(x0)
                in_service[k] || continue
                _applies(limit.when, models, initial, k) || continue
                bound, limits = upper !== nothing && x0[k] > _at(upper, k) ? (limit.upper, upper) :
                                lower !== nothing && x0[k] < _at(lower, k) ? (limit.lower, lower) : (nothing, nothing)
                isnothing(bound) && continue
                push!(violations, (model = model_name, instance = instances[k], variable = limit.variable,
                                   value = x0[k], bound = bound isa Symbol ? bound : :computed,
                                   limit = _at(limits, k)))
                adjust || continue
                bound isa Symbol || throw(ArgumentError(
                    "$(model_name) ($(join(string.(instances[k]), ", "))): $(limit.variable) = " *
                    "$(_shown(x0[k])) is outside its limit $(_shown(_at(limits, k))), " *
                    "which the model computes from its parameters and cannot be moved to it"))
                limits[k] = x0[k]
            end
        end
    end
    isempty(violations) && return violations

    lines = ["  $(v.model) ($(join(string.(v.instance), ", "))): $(v.variable) = $(_shown(v.value)) is outside " *
             "$(v.bound) = $(_shown(v.limit))" for v in violations]
    adjust || throw(ArgumentError(
        "Initial values outside their limits (TDSConfig(adjust_limits = true) moves the " *
        "limits to them):\n" * join(lines, "\n")))
    @warn "Limits moved to the initial values outside them, so the simulation starts in " *
          "steady state with these devices at their limits:\n" * join(lines, "\n")
    return violations
end
