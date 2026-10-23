# Solution types for power flow analysis
# This file is included early to resolve circular dependencies with model files

export PowerFlowSolution, MultiIslandSolution

"""
    PowerFlowSolution{T, PROB, DAE_T}

Result from power flow calculation.

Contains:
- Solution state (converged, iterations, residual)
- Reference to problem (for system access)
- Reference to DAE (for variable access)
- Lazy-computed observations (e.g., generator residuals)

# Fields
- `converged::Bool`: Whether the solver converged
- `residual_norm::T`: Final residual norm (||g||∞)
- `iterations::Int`: Number of iterations taken
- `retcode::Any`: Return code from solver
- `stats::Any`: Solver statistics (timing, etc.)
- `prob::PROB`: Reference to the problem (ACPowerFlowProblem)
- `dae::DAE_T`: Reference to workspace DAE (enables sol.dae.y access)

# Usage
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

if sol.converged
    # Access solution variables
    v = get_var(sol, :Bus, :v, 1)

    # Access raw state vector
    y = sol.dae.y
end
```

# Lazy Evaluation
Generator residuals are lazily computed and cached on first access. This avoids
unnecessary computation when these values aren't needed.

# DAE Reference
The `dae` field stores a reference to the workspace DAE from the solver. This enables
efficient access to the solution state vector (sol.dae.y) for TDS initialization.
Since it's a reference, memory overhead is minimal (~8 bytes).
"""
struct PowerFlowSolution{T, PROB, DAE_T} <: AbstractSteadyStateSolution
    converged::Bool
    residual_norm::T
    iterations::Int
    retcode::Any
    stats::Any
    prob::PROB  # ACPowerFlowProblem (contains sys, address_manager, formulation)
    dae::DAE_T  # Reference to workspace DAE (enables sol.dae.y access for state vector)
    _generator_residuals::Base.RefValue{Union{Nothing, Vector{T}}}  # Lazy cache
end

"""
    Base.show(io::IO, sol::PowerFlowSolution)

Display convergence status and key metrics.
"""
function Base.show(io::IO, sol::PowerFlowSolution{T}) where {T}
    if get(io, :compact, false)
        print(io, "PowerFlowSolution{$T}(...)")
    else
        println(io, "PowerFlowSolution{$T}:")
        println(io, "  Converged: $(sol.converged)")
        println(io, "  Iterations: $(sol.iterations)")
        println(io, "  Residual: $(round(sol.residual_norm, sigdigits=3))")
        println(io, "  Variables: $(length(sol.dae.y))")
    end
end

# Accessor for address_manager (convenience)
Base.getproperty(sol::PowerFlowSolution, s::Symbol) = s === :address_manager ? sol.prob.address_manager : getfield(sol, s)

# ============================================================================
# AbstractSteadyStateSolution Interface Implementation
# ============================================================================
# This allows PowerfulTDS to depend on PowerfulCore instead of PowerfulPF

import PowerfulCore: get_system, get_address_manager, get_state_vector, is_converged, get_residual_norm

"""
    get_system(sol::PowerFlowSolution) -> SystemModel

Get the system model from a power flow solution.
"""
get_system(sol::PowerFlowSolution) = sol.prob.sys

"""
    get_address_manager(sol::PowerFlowSolution) -> AddressManager

Get the address manager from a power flow solution.
"""
get_address_manager(sol::PowerFlowSolution) = sol.prob.address_manager

"""
    get_state_vector(sol::PowerFlowSolution) -> Vector

Get the solution state vector (algebraic variables y).
"""
get_state_vector(sol::PowerFlowSolution) = sol.dae.y

"""
    is_converged(sol::PowerFlowSolution) -> Bool

Check if the power flow solution converged.
"""
is_converged(sol::PowerFlowSolution) = sol.converged

"""
    get_residual_norm(sol::PowerFlowSolution) -> Real

The largest power flow mismatch when the solve stopped (`sol.residual_norm`).
"""
get_residual_norm(sol::PowerFlowSolution) = sol.residual_norm

"""
    MultiIslandSolution{T}

Solution from power flow calculation on a multi-island system.

Transparent wrapper that makes multiple islands look like a single solution.
Users can query variables using the same API as PowerFlowSolution.

# Fields
- `islands::Vector{SystemModel}`: Individual island systems
- `solutions::Vector{PowerFlowSolution{T}}`: Per-island solutions
- `converged::Bool`: True if ALL islands converged
- `iterations::Int`: Maximum iterations across all islands
- `residual_norm::T`: Maximum residual norm across all islands
- `retcode::Any`: Aggregated return code
- `stats::Any`: Aggregated statistics

# Example
```julia
# Automatic island handling (recommended)
prob = ACPowerFlowProblem(sys)  # Auto-detects islands
sol = solve(prob, RobustNR())

# Access works the same regardless of islands
v = get_var(sol, :Bus, :v, (4001,))  # Finds correct island automatically
```
"""
struct MultiIslandSolution{T}
    islands::Vector{SystemModel}
    solutions::Vector  # Vector of PowerFlowSolution (type inferred from construction)
    converged::Bool
    iterations::Int
    residual_norm::T
    retcode::Any
    stats::Any
end

"""
    Base.show(io::IO, sol::MultiIslandSolution)

Display convergence status for multi-island solution.
"""
function Base.show(io::IO, sol::MultiIslandSolution{T}) where {T}
    if get(io, :compact, false)
        print(io, "MultiIslandSolution{$T}(...)")
    else
        n_converged = count(s -> s.converged, sol.solutions)
        total_vars = sum(s -> length(s.dae.y), sol.solutions)

        println(io, "MultiIslandSolution{$T}:")
        println(io, "  Islands: $(length(sol.islands))")
        println(io, "  Converged: $n_converged/$(length(sol.islands))")
        println(io, "  Max iterations: $(sol.iterations)")
        println(io, "  Max residual: $(round(sol.residual_norm, sigdigits=3))")
        println(io, "  Total variables: $total_vars")
    end
end

# ============================================================================
# Multi-Island Solution Accessors
# ============================================================================

"""
    get_var(sol::MultiIslandSolution, model::Symbol, var::Symbol, ext_key) -> T

Get variable value from multi-island solution using external key.

Transparently searches through islands to find which contains the key.

# Examples
```julia
# After solving multi-island system
prob = ACPowerFlowProblem(sys)  # auto-detects islands
sol = solve(prob, RobustNR())

# Access variables without knowing island structure
v = get_var(sol, :Bus, :v, (4001,))
theta = get_var(sol, :Bus, :theta, (4001,))
```
"""
function get_var(sol::MultiIslandSolution, model::Symbol, var::Symbol, ext_key)
    # Linear search through islands to find which contains this key
    for (i, island) in enumerate(sol.islands)
        # Check if island has IndexOrder
        if isnothing(island.order)
            error("Island $i has no IndexOrder - cannot use get_var with external keys")
        end

        # Check if model exists in this island
        if !haskey(island.order.key_indices, model)
            continue  # Model not in this island, try next
        end

        # Check if key exists in this model
        if find_index(island, model, ext_key) !== nothing
            # Found the island - access from corresponding solution
            return get_var(sol.solutions[i], model, var, ext_key)
        end
    end

    error("External key $ext_key not found in model :$model in any island")
end

"""
    get_all_vars(sol::MultiIslandSolution, model::Symbol, var::Symbol) -> Vector

Get all values for a variable across all islands (concatenated in island order).

# Example
```julia
sol = solve(ACPowerFlowProblem(sys), RobustNR())
voltages = get_all_vars(sol, :Bus, :v)  # All bus voltages from all islands
```
"""
function get_all_vars(sol::MultiIslandSolution, model::Symbol, var::Symbol)
    result = []
    for (i, island) in enumerate(sol.islands)
        if haskey(island.models, model)
            island_vars = get_all_vars(sol.solutions[i], model, var)
            append!(result, island_vars)
        end
    end
    return result
end

"""
    get_vars_dict(sol::MultiIslandSolution, model::Symbol, var::Symbol) -> Dict

Get all values as a dictionary keyed by external keys (merged from all islands).

# Example
```julia
sol = solve(ACPowerFlowProblem(sys), RobustNR())
voltages = get_vars_dict(sol, :Bus, :v)
# Dict(1 => 1.06, 2 => 1.045, 4001 => 1.02, ...)  # Includes all islands
```
"""
function get_vars_dict(sol::MultiIslandSolution, model::Symbol, var::Symbol)
    result = Dict()
    for (i, island) in enumerate(sol.islands)
        if haskey(island.order.key_indices, model)
            island_dict = get_vars_dict(sol.solutions[i], model, var)
            merge!(result, island_dict)
        end
    end
    return result
end

