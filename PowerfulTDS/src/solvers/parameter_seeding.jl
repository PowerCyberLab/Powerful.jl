"""
Parameters of Dynamic Models Seeded with Dual Numbers

The shared basis of the sensitivities of the TDS system to parameters of its dynamic
models (controller gains, time constants, inertias, ...): the eigenvalue sensitivities
(`eigen_sensitivity.jl`) and other packages' own sensitivity integrations.

A parameter moves the operating point: the power flow is held (dispatch and bus
voltages), and the initialization computes the states and setpoints of the dynamic
models from it and from the parameters. Running the initialization with the parameters
as Dual numbers gives the operating point and its derivatives, and the model fields and
setpoints that move with the parameters, which are written as Dual values into the
workspaces that evaluate the residuals.
"""

using ForwardDiff: ForwardDiff
using StructArrays: StructArray, components
using PowerfulCore: ParameterSlice, get_float_ws

# Tag of the Dual numbers carrying the parameter derivatives
struct ParameterSensitivityTag end

# Widths of the Dual numbers: the number of parameters differentiated in one pass. A width
# per parameter count would compile the Dual initialization and the nested Jacobian anew
# for every count. The narrower width serves the few parameters of a typical tuning; the
# wider one keeps up to its count in one pass (each further pass initializes again).
const SENSITIVITY_WIDTHS = (8, 16)

"""
    _parameter_list(dyn_models, params) -> Vector{Tuple{Symbol, Symbol, Int}}

The parameters `params`, `(model, field, instances)` tuples, as one
`(model, field, instance)` per parameter. Errors if there is none, one is given twice,
or one is not a real field of an instance of one of `dyn_models`.
"""
function _parameter_list(dyn_models::NamedTuple, params)
    param_list = [(m, f, Int(i)) for (m, f, instances) in params for i in instances]
    isempty(param_list) && throw(ArgumentError("No parameters given"))
    allunique(param_list) || throw(ArgumentError("A parameter is given twice: $param_list"))
    _check_parameters(dyn_models, param_list)
    return param_list
end

"""
    _sensitivity_dual_type(np) -> Type{<:ForwardDiff.Dual}

The Dual type for `np` parameters: the narrower of `SENSITIVITY_WIDTHS` if it takes them
all, the wider otherwise (in passes of its width).
"""
function _sensitivity_dual_type(np::Int)
    narrow, wide = SENSITIVITY_WIDTHS
    return ForwardDiff.Dual{ParameterSensitivityTag, Float64, np <= narrow ? narrow : wide}
end

"""
    _seeded_operating_point(sol_pf, dyn_models, params, ::Type{D}, prob, adapter, cache; config)
        -> NamedTuple

The operating point of `prob` from the initialization with the parameters `params`
seeded (partial `j` along parameter `j`, the other partials of `D` zero): the state
`u_d = [x; y]` and the mass coefficients `M_d` of the differential states as Dual
numbers, and the model fields and setpoints that move with the parameters (`slices` of
`cache` and their Dual `values`).
"""
function _seeded_operating_point(sol_pf, dyn_models, params, ::Type{D}, prob, adapter, cache;
                                 config) where D
    init_d = setup_tds(sol_pf, _seed_parameters(dyn_models, params, D); config)
    _check_external_parameters(init_d.sys, D)

    am = prob.address_manager
    x_d, y_d = zeros(D, adapter.n_diff), zeros(D, adapter.n_alg)
    init_tds_buffers!(x_d, y_d, am, init_d.sys,
                      TDSInitData(dyn = init_d.dyn, setpoints = init_d.setpoints, pf = init_d.pf))
    u_d = [x_d; y_d]
    isapprox(ForwardDiff.value.(u_d), pack_state(adapter); rtol = 1e-12, atol = 1e-12) ||
        error("The initialization with Dual parameters does not reproduce the operating point")

    M_d = diag(build_tds_mass_matrix(init_d.sys, am, x_d, y_d))[1:adapter.n_diff]
    slices, values = _moving_quantities(cache, init_d.sys, init_d.setpoints, D)
    return (; u_d, M_d, slices, values)
end

"""
    _check_parameters(dyn_models, param_list)

Error unless each parameter `(model, field, instance)` of `param_list` is a real field of
an instance of one of the models `dyn_models`.
"""
function _check_parameters(dyn_models::NamedTuple, param_list)
    for (m, f, i) in param_list
        haskey(dyn_models, m) ||
            throw(ArgumentError("Model :$m is not one of the dynamic models $(keys(dyn_models))"))
        model = dyn_models[m]
        hasproperty(model, f) || error("Model :$m has no field :$f")
        eltype(getproperty(model, f)) === Float64 || error(":$m.$f holds $(eltype(getproperty(model, f))), not real numbers")
        1 <= i <= length(model) ||
            throw(ArgumentError("Instance $i out of range for :$m ($(length(model)) instances)"))
    end
    return nothing
end

"""
    _seed_parameters(dyn_models, param_list, ::Type{D}) -> NamedTuple

`dyn_models` rebuilt with element type `D`, their real fields converted and the
parameter `j` of `param_list` seeded with the `j`-th unit partial (`param_list` has at
most as many parameters as `D` has partials). Every model is converted, so that each
initialization receives the values of the models before it in its own element type.
"""
function _seed_parameters(dyn_models::NamedTuple, param_list, ::Type{D}) where D
    width = ForwardDiff.npartials(D)
    length(param_list) <= width || throw(ArgumentError("More parameters than partials of $D"))

    seeded = map(keys(dyn_models)) do m
        model = dyn_models[m]
        cols = map(c -> eltype(c) === Float64 ? D.(c) : c, components(model))
        for (j, (pm, f, i)) in enumerate(param_list)
            pm === m || continue
            col = getproperty(cols, f)
            col[i] = D(ForwardDiff.value(col[i]), ForwardDiff.Partials(ntuple(k -> Float64(k == j), width)))
        end
        StructArray{get_base_model_type(model){D}}(cols)
    end
    return NamedTuple{keys(dyn_models)}(Tuple(seeded))
end

# Whether any entry of a vector carries a nonzero derivative
_moves(col) = any(x -> any(!iszero, ForwardDiff.partials(x)), col)

"""
    _moving_quantities(cache, sys_d, setpoints_d, ::Type{D}) -> (slices, values)

The model fields of the Dual system `sys_d` and the setpoints `setpoints_d` that carry
derivatives, as parameter slices of `cache` (every instance) and their Dual values in
slice order.
"""
function _moving_quantities(cache, sys_d::SystemModel, setpoints_d::NamedTuple, ::Type{D}) where D
    slices = ParameterSlice[]
    values = D[]
    workspaces = get_float_ws(cache).workspaces_by_name
    function add!(model_name, field, col, source)
        eltype(col) === D && _moves(col) || return
        slice = ParameterSlice(cache, model_name, field, 1:length(col), length(values) + 1; source)
        push!(slices, slice)
        append!(values, col)
    end
    for (model_name, model) in pairs(sys_d.models)
        haskey(workspaces, model_name) || continue
        for field in propertynames(model)
            add!(model_name, field, getproperty(model, field), :model)
        end
    end
    for (model_name, setpoints) in pairs(setpoints_d), (name, col) in pairs(setpoints)
        add!(model_name, name, col, :setpoint)
    end
    return slices, values
end

"""
    _check_external_parameters(sys, ::Type{D})

Error if an external parameter of a model (gathered from another model when the system
is built) carries a derivative: the residual evaluation reads those values from the
Float64 system, so the derivative would be lost.
"""
function _check_external_parameters(sys::SystemModel, ::Type{D}) where D
    for (model_name, ext) in pairs(sys.external), (name, col) in pairs(ext)
        col isa AbstractVector && eltype(col) === D && _moves(col) &&
            error("Parameter enters :$model_name through its external parameter :$name, " *
                  "which the sensitivities do not differentiate")
    end
    return nothing
end
