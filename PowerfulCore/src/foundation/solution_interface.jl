# ============================================================================
# Abstract Steady-State Solution Interface
#
# This interface decouples downstream packages (e.g. PowerfulTDS) from
# the concrete PowerFlowSolution type in PowerfulPF, enabling parallel
# precompilation.
#
# Any steady-state solver solution (power flow, initialization, etc.) can
# implement this interface to be used for dynamic simulation initialization.
# ============================================================================

"""
    AbstractSteadyStateSolution

Abstract type for steady-state solution results.

Packages that need steady-state solutions for initialization (e.g., TDS and any
simulation mode built on it) should depend on this interface rather than concrete
solution types.

# Required Interface Methods
Implementations must provide:
- `get_system(sol)` → `SystemModel`
- `get_address_manager(sol)` → `AddressManager`
- `get_state_vector(sol)` → `AbstractVector{<:Real}`
- `get_residual_norm(sol)` → `Real`

# Optional Interface Methods
- `is_converged(sol)` → `Bool` (default: `true`)

# Example Implementation
```julia
struct PowerFlowSolution <: AbstractSteadyStateSolution
    prob::ACPowerFlowProblem
    dae::DAE
    converged::Bool
    # ...
end

get_system(sol::PowerFlowSolution) = sol.prob.sys
get_address_manager(sol::PowerFlowSolution) = sol.prob.address_manager
get_state_vector(sol::PowerFlowSolution) = sol.dae.y
is_converged(sol::PowerFlowSolution) = sol.converged
```

# Usage in TDS
```julia
function setup_tds(sol::AbstractSteadyStateSolution, dyn_models::NamedTuple; config)
    sys = get_system(sol)
    am = get_address_manager(sol)
    y = get_state_vector(sol)
    # ... use sys, am, y for initialization
end
```
"""
abstract type AbstractSteadyStateSolution end

"""
    get_system(sol::AbstractSteadyStateSolution) -> SystemModel

Get the system model from a steady-state solution.

This provides access to:
- `sys.models` - Model collections
- `sys.order.key_indices` - External key to internal index mapping
- `sys.config` - System configuration
"""
function get_system end

"""
    get_address_manager(sol::AbstractSteadyStateSolution) -> AddressManager

Get the address manager from a steady-state solution.

This provides:
- Variable address lookups via `get_addresses(am, model, var)`
- Address space dimensions
"""
function get_address_manager end

"""
    get_state_vector(sol::AbstractSteadyStateSolution) -> AbstractVector{<:Real}

Get the solution state vector (algebraic variables).

For power flow solutions, this is the `y` vector containing:
- Bus voltage magnitudes (v)
- Bus voltage angles (theta)
- Other algebraic variables per formulation
"""
function get_state_vector end

"""
    is_converged(sol::AbstractSteadyStateSolution) -> Bool

Check if the solution converged.

Default implementation returns `true`. Override for solutions that track convergence.
"""
is_converged(::AbstractSteadyStateSolution) = true

"""
    get_residual_norm(sol::AbstractSteadyStateSolution) -> Real

The largest residual of the solution's equations (infinity norm) when it stopped: how
closely the solution meets each of them, e.g. a slack bus's `v - vs`.
"""
function get_residual_norm end

# Export interface
export AbstractSteadyStateSolution
export get_system, get_address_manager, get_state_vector, is_converged, get_residual_norm
