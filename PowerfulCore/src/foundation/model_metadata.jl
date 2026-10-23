"""
Interface for model metadata.

All models must implement model_metadata(::Type{ModelType}) to describe
their variables and residuals.
"""

# ============================================================================
# Type Utilities
# ============================================================================

"""
    get_base_model_type(model_vec::StructArray{T}) -> Type

Extract base model type from a StructArray (e.g., `StructArray{PQ{Float64}}` → `PQ`).

Returns the unparameterized type for metadata lookup.
"""
function get_base_model_type(model_vec::StructArray{T}) where T <: AbstractModel
    # For StructArray, the eltype IS the base type
    # Return the unparameterized base type (e.g., Bus instead of Bus{Float64})
    return Base.unwrap_unionall(T).name.wrapper
end

# ============================================================================
# Model Metadata Interface
# ============================================================================

"""
    model_metadata(::Type{M}) where M<:AbstractModel -> ModelMetadata

Get metadata describing a model's variables, residuals, and structure.

Each model type must implement this trait to define its DAE interface.
Returns `ModelMetadata` with `vars`, `residuals`, and `observable_vars` fields.
"""
function model_metadata(model::Type)
    error("No `model_metadata` function defined for model $model. " *
          "Each model must implement model_metadata(::Type{$model}) to describe its structure.")
end

"""
    calc_residuals!(ws::ResidualBuffers, model_vec::AbstractModel)

Calculate residuals for a model. Each model type must implement this trait.

- `ws`: Pre-allocated workspace with gathered variables in `ws.var_vals`
- `model_vec`: Model vector instance

Implementations read from `ws.var_vals` and write results to `ws.res_vals`.
"""
function calc_residuals! end

# ============================================================================
# Observation Computation Interface
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{ModelT}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all instances of model type `ModelT` in the solution.

This is the **primary API** for models to provide derived outputs (power flows, generations, etc.)
that are computed from the solved state rather than being direct state variables.

# Arguments
- `sol::PowerFlowSolution`: Complete power flow solution
  - `sol.dae.y`: State vector (voltages, angles)
  - `sol.prob.sys`: System model (parameters, topology)
  - `sol.prob.address_manager`: Variable address mapping
  - Residuals available via `get_generator_residuals(sol)` for power balance extraction
- `::Type{ModelT}`: Model type (e.g., Branch, PV, Slack, PQ)
- `var::Symbol`: Observation variable name (defined in model's observations metadata)

# Returns
- `Vector{Float64}`: Values for all instances, length = number of model instances

# Implementation Guidelines
- **Vectorized**: Compute all instances at once for efficiency (use broadcasting, SIMD, etc.)
- **Complete access**: Solution object provides everything needed (state, parameters, residuals)
- **Internal helpers allowed**: Models can use scalar helpers internally if needed
- **Error handling**: Throw descriptive error for unknown variable names

# Example Implementations

```julia
# Simple case: direct computation from state and parameters
function compute_observations(sol::PowerFlowSolution, ::Type{Branch}, var::Symbol)
    branches = sol.prob.sys.models.Branch
    am = sol.prob.address_manager

    # Gather voltages for all branches (vectorized)
    theta_i = sol.dae.y[get_addresses(am, :Branch, :theta_i)]
    theta_j = sol.dae.y[get_addresses(am, :Branch, :theta_j)]
    v_i = sol.dae.y[get_addresses(am, :Branch, :v_i)]
    v_j = sol.dae.y[get_addresses(am, :Branch, :v_j)]

    # Compute power flows for all branches
    if var == :p_from
        return @. v_i^2 * branches.gi - v_i * v_j * (branches.gi * cos(theta_i - theta_j) + branches.bi * sin(theta_i - theta_j))
    # ... other variables
    end
end

# Complex case: extract from residuals
function compute_observations(sol::PowerFlowSolution, ::Type{Slack}, var::Symbol)
    slack = sol.prob.sys.models.Slack

    if var == :q_gen
        # Extract reactive power from bus power balance residuals
        g_gen = get_generator_residuals(sol)
        bus_indices = slack.i
        q_addrs = get_addresses(sol.prob.address_manager, :Bus, :q)
        return g_gen[q_addrs[bus_indices]]  # Vectorized gather
    end
end
```

# Generic Fallback
This method provides a generic fallback that automatically handles:
1. **State variables**: Extract from solution vector (e.g., :theta, :v for Bus)
2. **Parameters**: Return from model struct (e.g., :pg for PV)

Models only need custom implementations for computed observations (e.g., power flows).

# See Also
- [`model_metadata`](@ref): Defines observations available for each model
- [`collect_model_table`](@ref PowerfulPF.collect_model_table): Uses this API for report generation
"""
function compute_observations end

# Generic fallback implementation in src/solvers/powerflow.jl
# (after PowerFlowSolution type is defined)

"""
    get_status_field(model_type::Type{<:AbstractModel}) -> Union{Symbol, Nothing}

Get the status field name (`:status`) for a model type.

Returns `nothing` if the model has no status field (e.g., Bus).
Every model with an in-service flag names it `status` (Bool).
"""
function get_status_field(model_type::Type{<:AbstractModel})
    metadata = model_metadata(model_type)
    return metadata.status_field
end

"""
    get_status_field(model_vec::StructArray) -> Union{Symbol, Nothing}

Convenience wrapper that extracts base type and calls `get_status_field(Type)`.
"""
function get_status_field(model_vec::StructArray{T}) where {T<:AbstractModel}
    base_type = get_base_model_type(model_vec)
    return get_status_field(base_type)
end

# ============================================================================
# Metadata Query API
# ============================================================================

"""
    has_state_vars(metadata::ModelMetadata) -> Bool

Check if model has any differential state variables (`StateVar`).

Returns `true` for dynamic models with time derivatives (e.g., GENCLS delta, omega),
`false` for static power flow models (e.g., PQ).
"""
has_state_vars(metadata::ModelMetadata) = any(is_state_var, metadata.vars)

"""
    get_terminals(::Type{ModelType}) -> Tuple{Terminal, ...}

Declare which vertices this model connects to in the network graph.

Terminal count determines behavior:
- 0 terminals: Topological model (Area, Zone) - no spatial location
- 1 terminal: Single-attachment device (PQ, PV) - island membership only
- 2+ terminals: Multi-terminal device (Branch) - creates edges

Default returns empty tuple (no terminals).
"""
function get_terminals(::Type{M}) where M <: AbstractModel
    return ()  # Default: no terminals (non-connecting device)
end

# ============================================================================
# Terminal Query API
# ============================================================================

"""
    find_terminal(model_type::Type, vertex_model::Symbol) -> Terminal or nothing

Find the first terminal connecting to a specific vertex type (e.g., `:Bus`).

- `model_type`: Model type (e.g., `GENCLS`, `Branch`, `PQ`)
- `vertex_model`: Vertex type to find (e.g., `:Bus`, `:DCBus`)

Returns `nothing` if no terminal connects to that vertex type.
"""
function find_terminal(model_type::Type, vertex_model::Symbol)
    terminals = get_terminals(model_type)
    for t in terminals
        if t.vertex_model == vertex_model
            return t
        end
        # Also match if looking for :Bus and terminal uses :bus_group
        # (:bus_group is a polymorphic reference that resolves to any model implementing
        # it, such as Bus)
        if t.vertex_model == :bus_group && vertex_model == :Bus
            return t
        end
    end
    return nothing
end

# ============================================================================
# Key Field Trait
# ============================================================================

"""
    get_key_name(::Type{M}) where M <: AbstractModel -> Symbol or Vector{Symbol}

Trait that specifies the field(s) that uniquely identify instances of a model type.

Key fields are used for:
1. **Indexing**: Building external ID → internal index mappings
2. **Lookup**: Finding specific instances by external identifier
3. **Group validation**: Verifying model key fields match group interface

# Returns
- `Symbol`: Single key field (e.g., `:i` for Bus)
- `Vector{Symbol}`: Composite key fields (e.g., `[:i, :id]` for PQ, `[:i, :j, :ckt]` for Branch)

# Default Implementation
Throws `ArgumentError` - all models must implement this trait.

# Examples
```julia
# Single key field (bus number)
get_key_name(::Type{Bus}) = :i

# Composite key (bus + load ID)
get_key_name(::Type{<:PQ}) = [:i, :id]

# Three-field composite key (from bus, to bus, circuit)
get_key_name(::Type{Branch}) = [:i, :j, :ckt]
```

# See Also
- [`model_metadata`](@ref): Model structure metadata
- [`validate_group_implementation`](@ref): Uses key fields for group validation
"""
function get_key_name(::Type{M}) where M
    throw(ArgumentError("No unique index name defined for model type $M"))
end

export get_base_model_type
export model_metadata, calc_residuals!, compute_observations, get_status_field
export has_state_vars, get_terminals
export find_terminal
export get_key_name
