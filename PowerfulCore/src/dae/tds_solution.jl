# ============================================================================
# TDS Solution - Unified solution type for time-domain simulation results
#
# This provides a unified interface for both:
# - TDS (DiffEq ODESolution) - preserves interpolation and plotting
# - Ybus (trapezoidal) - stores time vector and state matrix directly
#
# Key design: "Worse is Better" - different internal representations,
# consolidated accessors dispatch internally. Zero data duplication.
# ============================================================================

using ..Foundation: SystemModel, AddressManager, StateVar, AlgebVar

"""
    TDSSolution{SOL, T, SYS, AM}

Unified solution type for time-domain simulation results.

Works for both TDS (DiffEq ODESolution) and Ybus (trapezoidal) solvers.
Provides ergonomic `get_series(sol, :Model, :var, key)` API without needing
the problem object.

# Type Parameters
- `SOL`: ODESolution type (for TDS) or `Nothing` (for Ybus)
- `T`: Numeric type (typically Float64)
- `SYS`: SystemModel type
- `AM`: AddressManager type

# Fields
- `sol`: ODESolution (TDS) or `nothing` (Ybus)
- `t`: Time vector (Ybus only) or `nothing` (TDS uses sol.t)
- `u`: State matrix (Ybus only) or `nothing` (TDS uses sol.u)
- `retcode`: Solution status (:Success, :MaxIters, etc.)
- `n_diff`: Number of differential states
- `n_alg`: Number of algebraic states
- `sys`: SystemModel reference (for external key lookups)
- `address_manager`: AddressManager reference (for address lookups)

# Example
```julia
# Get time series for bus voltage
t, v = get_series(sol, :Bus, :v, 5001)

# Interpolation (TDS only)
state_at_0_5 = sol(0.5)

# Access time vector
t = sol.t

# Index into timesteps
u_at_step_10 = sol[10]
```
"""
struct TDSSolution{SOL, T, SYS<:SystemModel, AM<:AddressManager}
    sol::SOL                              # ODESolution (TDS) or nothing (Ybus)
    t::Union{Nothing, Vector{T}}          # Ybus: vector, TDS: nothing (use sol.t)
    u::Union{Nothing, Matrix{T}}          # Ybus: matrix, TDS: nothing (use sol.u)
    retcode::Symbol                       # :Success, :MaxIters, etc.
    n_diff::Int                           # Differential state count
    n_alg::Int                            # Algebraic state count
    sys::SYS                              # SystemModel reference
    address_manager::AM                   # For get_addresses()
end

# Constructor for TDS (ODESolution) - infers T from sol.t
function TDSSolution(
    sol::SOL,
    ::Nothing,
    ::Nothing,
    retcode::Symbol,
    n_diff::Int,
    n_alg::Int,
    sys::SYS,
    address_manager::AM
) where {SOL, SYS<:SystemModel, AM<:AddressManager}
    # Infer T from sol.t element type
    T = eltype(sol.t)
    return TDSSolution{SOL, T, SYS, AM}(sol, nothing, nothing, retcode, n_diff, n_alg, sys, address_manager)
end

# ============================================================================
# Consolidated Access Layer (Single Point of Truth)
#
# All internal access to time/state data goes through these inline functions.
# This is the ONLY place where TDS vs Ybus conditionals appear.
# ============================================================================

"""
    _get_t(sol::TDSSolution) -> AbstractVector

Get the time vector from a TDSSolution.
Uses `getfield` to avoid infinite recursion with `getproperty`.
"""
@inline function _get_t(sol::TDSSolution)
    sol_inner = getfield(sol, :sol)
    isnothing(sol_inner) ? getfield(sol, :t) : sol_inner.t
end

"""
    _get_u_at_step(sol::TDSSolution, step::Int) -> AbstractVector

Get the state vector at a specific timestep.
Uses `getfield` to avoid infinite recursion with `getproperty`.
"""
@inline function _get_u_at_step(sol::TDSSolution, step::Int)
    sol_inner = getfield(sol, :sol)
    isnothing(sol_inner) ? view(getfield(sol, :u), :, step) : sol_inner.u[step]
end

"""
    _get_value(sol::TDSSolution, addr::Int, step::Int) -> Number

Get a single value at (address, timestep).
Uses `getfield` to avoid infinite recursion with `getproperty`.
"""
@inline function _get_value(sol::TDSSolution, addr::Int, step::Int)
    sol_inner = getfield(sol, :sol)
    isnothing(sol_inner) ? getfield(sol, :u)[addr, step] : sol_inner.u[step][addr]
end

"""
    _n_steps(sol::TDSSolution) -> Int

Get the number of timesteps in the solution.
"""
@inline _n_steps(sol::TDSSolution) = length(_get_t(sol))

"""
    _has_interpolation(sol::TDSSolution) -> Bool

Check if interpolation is available (TDS solutions only).
"""
@inline _has_interpolation(sol::TDSSolution) = !isnothing(sol.sol)

# ============================================================================
# Solution Interface (from solution_interface.jl pattern)
# ============================================================================

"""
    get_system(sol::TDSSolution) -> SystemModel

Get the system model from a TDS solution.
"""
get_system(sol::TDSSolution) = sol.sys

"""
    get_address_manager(sol::TDSSolution) -> AddressManager

Get the address manager from a TDS solution.
"""
get_address_manager(sol::TDSSolution) = sol.address_manager

"""
    get_state_vector(sol::TDSSolution) -> AbstractVector

Get the final state vector from a TDS solution.
"""
get_state_vector(sol::TDSSolution) = _get_u_at_step(sol, _n_steps(sol))

"""
    is_converged(sol::TDSSolution) -> Bool

Check if the TDS solution converged successfully.
"""
is_converged(sol::TDSSolution) = sol.retcode === :Success

# ============================================================================
# Interface Forwarding (for DiffEq compatibility)
# ============================================================================

"""
    (sol::TDSSolution)(t)

Interpolate the solution at time `t`. Only available for TDS solutions.
"""
function (sol::TDSSolution)(t)
    _has_interpolation(sol) || error("Interpolation not available for Ybus solutions")
    return sol.sol(t)
end

"""
    Base.getindex(sol::TDSSolution, i::Int)

Get the state vector at timestep `i`.
"""
Base.getindex(sol::TDSSolution, i::Int) = _get_u_at_step(sol, i)

"""
    Base.getproperty(sol::TDSSolution, s::Symbol)

Property forwarding for common fields.
- `sol.t` returns the time vector
- `sol.u` returns the state data (vector-of-vectors for TDS, matrix for Ybus)
- `sol.retcode` is the stored status Symbol (`:Success`, `:MaxIters`, ...), as for
  `PowerFlowSolution`; the DiffEq `ReturnCode` of a TDS run is `sol.sol.retcode`
- `sol.stats` is the integrator statistics of a TDS run (`naccept`, `nreject`, `nf`,
  ...); Ybus solutions come from a fixed-step integrator and have none
"""
function Base.getproperty(sol::TDSSolution, s::Symbol)
    s === :t && return _get_t(sol)
    if s === :stats
        sol_inner = getfield(sol, :sol)
        isnothing(sol_inner) && error("A Ybus TDSSolution has no integrator statistics; " *
                                      "it was integrated with fixed steps (see `sol.t`)")
        return sol_inner.stats
    end
    # Forward .u to sol.sol.u for TDS (ODESolution), or return stored matrix for Ybus
    if s === :u
        sol_inner = getfield(sol, :sol)
        return isnothing(sol_inner) ? getfield(sol, :u) : sol_inner.u
    end
    return getfield(sol, s)
end

"""
    Base.length(sol::TDSSolution) -> Int

Return the number of timesteps in the solution.
"""
Base.length(sol::TDSSolution) = _n_steps(sol)

"""
    Base.lastindex(sol::TDSSolution) -> Int

Return the last valid index (for `sol[end]` syntax).
"""
Base.lastindex(sol::TDSSolution) = _n_steps(sol)

# ============================================================================
# Main API: get_series
# ============================================================================

"""
    get_series(sol::TDSSolution, model::Symbol, var::Symbol, ext_key) -> (t, values)

Extract a time series for a specific variable from the solution.

This is the primary API for extracting simulation results. It handles:
1. External key → internal index conversion
2. Variable type detection (StateVar vs AlgebVar)
3. Address calculation with proper offset
4. Value extraction from either TDS or Ybus storage format

# Arguments
- `sol::TDSSolution`: The simulation solution
- `model::Symbol`: Model name (e.g., `:Bus`, `:GENROU`)
- `var::Symbol`: Variable name (e.g., `:v`, `:delta`, `:omega`)
- `ext_key`: External key (e.g., `5001` or `(5001,)` for buses)

# Returns
- `(t, values)`: Tuple of time vector and value vector

# Example
```julia
# Get bus voltage
t, v = get_series(sol, :Bus, :v, 5001)

# Get generator rotor angle (key is a tuple for generators)
t, delta = get_series(sol, :GENROU, :delta, (101, "1"))

# Plot results
using Plots
plot(t, v, label="Bus 5001 voltage")
```
"""
function get_series(sol::TDSSolution, model::Symbol, var::Symbol, ext_key)
    return (_get_t(sol), _series_values(sol, state_address(sol, model, var, ext_key)))
end

"""
    state_address(sol::TDSSolution, model::Symbol, var::Symbol, ext_key) -> Int

Position of the variable `var` of the instance `ext_key` (external key tuple) of `model`
in the combined state vector `[x; y]` of `sol`.
"""
function state_address(sol::TDSSolution, model::Symbol, var::Symbol, ext_key)
    addrs = _series_addresses(sol, model, var)
    # External → Internal (key_indices covers static and dynamic models)
    return addrs[get_index_by_key(sol.sys, model, ext_key)]
end

# Absolute address of `var` for every instance of `model` in the combined [x; y] vector
function _series_addresses(sol::TDSSolution, model::Symbol, var::Symbol)
    haskey(sol.sys.models, model) ||
        error("Model :$model not found in system. Available: $(keys(sol.sys.models))")

    am = sol.address_manager
    addr_type = get_var_address_type(am, model, var)
    # StateVar addresses are 1:n_diff, AlgebVar addresses need offset
    offset = addr_type isa StateVar ? 0 : sol.n_diff
    return [Int(offset + a) for a in get_addresses(am, addr_type, model, var)]
end

# Values at one absolute address over all stored steps
function _series_values(sol::TDSSolution, addr::Int)
    n = _n_steps(sol)
    values = Vector{eltype(_get_t(sol))}(undef, n)
    @inbounds for i in 1:n
        values[i] = _get_value(sol, addr, i)
    end
    return values
end

"""
    _format_external_key(key) -> String

Format an instance key tuple for use in variable names (elements joined with underscores).
"""
_format_external_key(key::Tuple) = join(key, "_")

"""
    get_series(sol::TDSSolution, model::Symbol, var::Symbol) -> NamedTuple

Extract time series for ALL instances of a variable from the solution.

Returns a NamedTuple with `t` (time vector) and one column per instance,
named `{Model}_{var}_{external_key}`.

# Arguments
- `sol::TDSSolution`: The simulation solution
- `model::Symbol`: Model name (e.g., `:Bus`, `:GENROU`)
- `var::Symbol`: Variable name (e.g., `:v`, `:delta`, `:omega`)

# Returns
NamedTuple with fields:
- `t`: Time vector
- `{Model}_{var}_{key}`: One field per instance

# Example
```julia
# Get all bus voltages
data = get_series(sol, :Bus, :v)
# Returns (t=..., Bus_v_1=..., Bus_v_2=..., Bus_v_3=..., ...)

# Plot all voltages
using Plots
plot(data.t, [data[k] for k in keys(data) if k != :t])

# Convert to DataFrame
using DataFrames
df = DataFrame(data)
```

# See Also
- `get_series(sol, model, var, ext_key)` - Extract single variable
"""
function get_series(sol::TDSSolution, model::Symbol, var::Symbol)
    addrs = _series_addresses(sol, model, var)
    # External keys come directly from the model (works for both static and dynamic models)
    ext_keys = get_key(sol.sys.models[model])

    t = _get_t(sol)
    data_pairs = Pair{Symbol, Vector{eltype(t)}}[:t => collect(t)]
    for (int_idx, addr) in enumerate(addrs)
        col_name = Symbol(model, "_", var, "_", _format_external_key(ext_keys[int_idx]))
        push!(data_pairs, col_name => _series_values(sol, addr))
    end

    return NamedTuple(data_pairs)
end

# ============================================================================
# Show method for nice printing
# ============================================================================

function Base.show(io::IO, sol::TDSSolution{SOL, T}) where {SOL, T}
    solver_type = isnothing(sol.sol) ? "Ybus" : "TDS"
    n = _n_steps(sol)
    t = _get_t(sol)
    tspan = n > 0 ? "($(t[1]), $(t[end]))" : "(empty)"
    print(io, "TDSSolution{$solver_type, $T}: $n timesteps, tspan=$tspan, retcode=$(sol.retcode)")
end

function Base.show(io::IO, ::MIME"text/plain", sol::TDSSolution{SOL, T}) where {SOL, T}
    solver_type = isnothing(sol.sol) ? "Ybus (trapezoidal)" : "TDS (DiffEq)"
    n = _n_steps(sol)
    t = _get_t(sol)
    tspan = n > 0 ? "($(t[1]), $(t[end]))" : "(empty)"
    println(io, "TDSSolution{$T}")
    println(io, "  Solver: $solver_type")
    println(io, "  Timesteps: $n")
    println(io, "  Time span: $tspan")
    println(io, "  States: $(sol.n_diff) differential, $(sol.n_alg) algebraic")
    println(io, "  Return code: $(sol.retcode)")
    if _has_interpolation(sol)
        println(io, "  Interpolation: available")
    end
end
