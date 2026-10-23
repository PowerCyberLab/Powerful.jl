"""
    variable_accessors.jl

Generic accessor functions for reading/writing variables and residuals in SystemModel
using external keys (bus numbers, load IDs, etc.) or internal array indices.

These functions combine:
1. IndexOrder: external key → internal array index
2. AddressManager: model/variable → DAE addresses
3. DAE: addresses → actual values
"""

#==============================================================================
# Variable Accessors (DAE.y)
==============================================================================#

"""
    variable_addresses(am, model, var) -> Vector{UInt}

Addresses of variable `var` of `model` in the variable vector `y`. Residual names are
rejected: residual addresses index `g`, not `y`.
"""
function variable_addresses(am::AddressManager, model::Symbol, var::Symbol)
    return get_addresses(am, get_var_address_type(am, model, var), model, var)
end

"""
    get_var(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol, ext_key) -> T

Get variable value using external key (e.g., bus number).

Generic function that works for any variable in any model by combining:
1. IndexOrder: external key → internal array index
2. AddressManager: model/variable → DAE addresses
3. State vector: addresses → actual values

# Arguments
- `prob::ACPowerFlowProblem`: Power flow problem (contains sys and address_manager)
- `y::AbstractVector`: State vector (from solution, workspace, or custom)
- `model::Symbol`: Model name (`:Bus`, `:PQ`, `:PV`, etc.)
- `var::Symbol`: Variable name (`:v`, `:theta`, `:p`, etc.)
- `ext_key`: External key for lookup (bus number, or tuple for composite keys)

# Returns
- Variable value (type matches `y` element type)

# Examples
```julia
# After solving - query solution
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, RobustNR())
v = get_var(prob, sol.dae.y, :Bus, :v, (4001,))

# During iteration - query workspace state
workspace = setup_powerflow_newton(prob)
v_iter = get_var(prob, workspace.cache.iter_dae.y, :Bus, :v, (4001,))

# Query trial point during line search
v_trial = get_var(prob, workspace.y_trial, :Bus, :v, (4001,))

# Bus voltage angle
theta = get_var(prob, sol.dae.y, :Bus, :theta, (4001,))
```

Only variables (entries of `y`) are accessible. Residual names such as a load's `:p`
or a branch's `:p_i` raise an error; injections and flows are reported by
[`collect_model_table`](@ref).

# Errors
- If `prob.sys.order` is `nothing`: Error (IndexOrder required)
- If model or external key not found: error from `get_index_by_key`
- If model/variable not found: Error from AddressManager
"""
function get_var(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol, ext_key)
    sys = prob.sys
    am = prob.address_manager

    # Validate preconditions
    isnothing(sys.order) && error(
        "SystemModel has no IndexOrder - cannot map external keys. " *
        "Use get_var_by_index() for direct internal index access."
    )

    # Step 1: External key → internal index
    internal_idx = get_index_by_key(sys, model, ext_key)

    # Step 2: Get addresses for this variable
    addrs = variable_addresses(am, model, var)

    # Step 3: Return value from state vector
    return y[addrs[internal_idx]]
end

"""
    get_var_by_index(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol, idx::Int) -> T

Get variable value using internal array index (1-based indexing).

Faster than `get_var()` since it skips external→internal mapping.
Use when you already have internal indices.

# Arguments
- `prob::ACPowerFlowProblem`: Power flow problem
- `y::AbstractVector`: State vector
- `model::Symbol`: Model name
- `var::Symbol`: Variable name
- `idx::Int`: Internal array index (1-based)

# Example
```julia
# Direct access by internal index (no external key lookup)
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, RobustNR())
v = get_var_by_index(prob, sol.dae.y, :Bus, :v, 1)  # First bus in internal ordering
```
"""
function get_var_by_index(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol, idx::Int)
    addrs = variable_addresses(prob.address_manager, model, var)
    return y[addrs[idx]]
end

"""
    get_all_vars(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol) -> Vector

Get all values for a variable across all instances (in internal order).

Returns a view into the state vector, not a copy.

# Examples
```julia
# Get all bus voltages from solution
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, RobustNR())
voltages = get_all_vars(prob, sol.dae.y, :Bus, :v)

# Get all bus voltage angles
angles = get_all_vars(prob, sol.dae.y, :Bus, :theta)
```
"""
function get_all_vars(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol)
    addrs = variable_addresses(prob.address_manager, model, var)
    return y[addrs]
end

"""
    get_vars_dict(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol) -> Dict

Get all values as a dictionary keyed by external keys.

Useful for exporting results with semantic bus numbers.

# Example
```julia
# Get bus voltages indexed by external bus number
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, RobustNR())
voltages = get_vars_dict(prob, sol.dae.y, :Bus, :v)
# Dict(1 => 1.06, 2 => 1.045, 4001 => 1.02, ...)
```
"""
function get_vars_dict(prob::ACPowerFlowProblem, y::AbstractVector, model::Symbol, var::Symbol)
    sys = prob.sys
    isnothing(sys.order) && error("SystemModel has no IndexOrder")

    result = Dict()
    addrs = variable_addresses(prob.address_manager, model, var)

    for (ext_key, internal_idx) in sys.order.key_indices[model]
        result[ext_key] = y[addrs[internal_idx]]
    end

    return result
end

# String convenience methods for Python interop (strings don't auto-convert to
# Symbols in juliacall): convert model and variable names, then dispatch.
for f in (:get_var, :get_var_by_index, :get_all_vars, :get_vars_dict)
    @eval begin
        $f(sol::Union{PowerFlowSolution, MultiIslandSolution}, model::AbstractString, var::AbstractString, args...) =
            $f(sol, Symbol(model), Symbol(var), args...)
        $f(prob::ACPowerFlowProblem, y::AbstractVector, model::AbstractString, var::AbstractString, args...) =
            $f(prob, y, Symbol(model), Symbol(var), args...)
    end
end


#==============================================================================
# Bus Voltages (formulation-neutral)
==============================================================================#

"""
    bus_voltages(sol::PowerFlowSolution) -> (v, theta)

Voltage magnitude (pu) and angle (rad) of every bus, in internal order, whatever the
formulation of the solution: `PowerInjection` solves them directly,
`CurrentInjection` solves the rectangular components (`v = |v_r + j·v_i|`), and
`DCPowerFlow` solves the angles, with the magnitudes at the case voltages `vm` that
its load and shunt models use.

Observations and reports read bus voltages through this function, so they are
defined for every formulation.
"""
bus_voltages(sol::PowerFlowSolution) = _bus_voltages(sol, sol.prob.formulation)

function _bus_voltages(sol::PowerFlowSolution, ::PowerInjection)
    am = sol.prob.address_manager
    return (sol.dae.y[get_addresses(am, AlgebVar(), :Bus, :v)],
            sol.dae.y[get_addresses(am, AlgebVar(), :Bus, :theta)])
end

function _bus_voltages(sol::PowerFlowSolution, ::CurrentInjection)
    am = sol.prob.address_manager
    v_r = sol.dae.y[get_addresses(am, AlgebVar(), :Bus, :v_r)]
    v_i = sol.dae.y[get_addresses(am, AlgebVar(), :Bus, :v_i)]
    return (hypot.(v_r, v_i), atan.(v_i, v_r))
end

function _bus_voltages(sol::PowerFlowSolution, ::DCPowerFlow)
    theta = sol.dae.y[get_addresses(sol.prob.address_manager, AlgebVar(), :Bus, :theta)]
    return (collect(Float64, sol.prob.sys.models.Bus.vm), theta)
end

"""
    terminal_voltages(sol::PowerFlowSolution, model::Symbol, terminal::Symbol) -> (v, theta)

Voltage magnitude and angle of the bus at terminal field `terminal` (e.g. `:i`) of
every instance of `model`, from [`bus_voltages`](@ref).
"""
function terminal_voltages(sol::PowerFlowSolution, model::Symbol, terminal::Symbol)
    v, theta = bus_voltages(sol)
    idx = get_foreign_indices(sol.prob.sys, model, terminal)
    return (v[idx], theta[idx])
end

"""Generator models whose units share the generator power of their bus."""
const GENERATOR_MODELS = (:PV, :Slack)

"""
    generator_power(sol::PowerFlowSolution, model::Symbol) -> (p, q)

Active and reactive power (pu) that each generator of `model` (`:PV` or `:Slack`)
supplies. The power flow solves the total generator power at each bus: the balance of
every other contribution to the bus power equations (`get_generator_residuals`). In
`CurrentInjection` the bus equations are current balances, and the slack bus rows
belong to the slack equations, so the power balances of the `PowerInjection` form are
evaluated at the solved voltages instead.

The units at a bus share its total:
- A PV unit supplies its scheduled active power `pg`. Slack units share the active
  power balance of their bus in proportion to their MVA base `mbase`, as units with
  equal per-unit droop do.
- All units at a bus share its reactive power in proportion to `rmpct`, the percent of
  the total Mvar each unit supplies.

An out-of-service unit supplies nothing. `DCPowerFlow` has no reactive power: `q` is
`nothing`.
"""
function generator_power(sol::PowerFlowSolution, model::Symbol)
    model in GENERATOR_MODELS ||
        error("generator_power: :$model is not a generator model $(GENERATOR_MODELS)")
    sys = sol.prob.sys
    units = sys.models[model]
    total = _generator_power(sol, get_foreign_indices(sys, model, :i), sol.prob.formulation)
    p = model === :PV ? scheduled_power(units) :
        total.p .* _bus_shares(sys, model, :mbase, (:Slack,))
    q = isnothing(total.q) ? nothing : total.q .* _bus_shares(sys, model, :rmpct, GENERATOR_MODELS)
    return (p = p, q = q)
end

"""
    scheduled_power(units) -> Vector

Scheduled active power (pu) of generator units: `pg` in service, zero out of service.
"""
scheduled_power(units) = units.status .* units.pg

"""
    _bus_shares(sys, model, weight, models) -> Vector{Float64}

Share of each unit of `model` in a power its bus supplies: its `weight` field over the
sum of the weights of the in-service units of `models` at the bus. A unit alone in
service at its bus supplies all of it, and an out-of-service unit none.
"""
function _bus_shares(sys, model::Symbol, weight::Symbol, models)
    n_bus = length(sys.models.Bus)
    total = zeros(n_bus)
    n_on = zeros(Int, n_bus)
    for m in models
        haskey(sys.models, m) || continue
        units = sys.models[m]
        w = getproperty(units, weight)
        for (k, b) in enumerate(get_foreign_indices(sys, m, :i))
            units.status[k] || continue
            total[b] += w[k]
            n_on[b] += 1
        end
    end

    units = sys.models[model]
    w = getproperty(units, weight)
    buses = get_foreign_indices(sys, model, :i)
    return map(eachindex(buses)) do k
        b = buses[k]
        units.status[k] || return 0.0
        n_on[b] == 1 && return 1.0
        total[b] > 0 || error("The $(n_on[b]) generators in service at bus " *
                              "$(sys.models.Bus.i[b]) have no $weight to share their power by")
        return w[k] / total[b]
    end
end

# Total generator power at each of `buses` (internal bus indices), per formulation
function _generator_power(sol, buses, ::PowerInjection)
    return _power_balance_generators(get_generator_residuals(sol), sol.prob.address_manager, buses)
end

function _generator_power(sol, buses, ::CurrentInjection)
    sys = sol.prob.sys
    polar = ACPowerFlowProblem(sys, PowerInjection(); auto_split = false)
    am = polar.address_manager
    v, theta = bus_voltages(sol)
    y = zeros(Int(max_address_in_space(am, AlgebVar())))
    length(y) == 2 * length(v) ||
        error("generator_power: the PowerInjection variables of this system are not only the bus voltages")
    y[get_addresses(am, AlgebVar(), :Bus, :v)] .= v
    y[get_addresses(am, AlgebVar(), :Bus, :theta)] .= theta
    g = evaluate_generator_contributions!(SystemBuffers(Float64, am, sys, PowerInjection()), y,
                                          build_status_dict(sys))
    return _power_balance_generators(g, am, buses)
end

function _generator_power(sol, buses, ::DCPowerFlow)
    g = get_generator_residuals(sol)
    # Negate: the bus residuals have the sign opposite to the generator output
    return (p = -g[get_addresses(sol.prob.address_manager, AlgebRes(), :Bus, :p)[buses]], q = nothing)
end

# Generator output from accumulated power balances (negated: the bus residuals have the
# sign opposite to the generator output)
_power_balance_generators(g, am, buses) =
    (p = -g[get_addresses(am, AlgebRes(), :Bus, :p)[buses]],
     q = -g[get_addresses(am, AlgebRes(), :Bus, :q)[buses]])
