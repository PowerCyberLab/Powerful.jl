"""
TDS initialization validation and diagnostics.
"""

using PrettyTables
using Printf: @sprintf

# ============================================================================
# Validation Result
# ============================================================================

"""
    ResidualValidationResult

Residual validation check results.

# Fields
- `max_diff_residual::Float64` - Maximum absolute differential residual
- `max_alg_residual::Float64` - Maximum absolute algebraic residual
- `non_zero_diff::Vector{Tuple{UInt, Float64, Symbol, Symbol}}` - (address, value, model, equation)
- `non_zero_alg::Vector{Tuple{UInt, Float64, Symbol, Symbol}}` - (address, value, model, equation)
- `is_valid::Bool` - True if all residuals below tolerance
"""
struct ResidualValidationResult
    max_diff_residual::Float64
    max_alg_residual::Float64
    non_zero_diff::Vector{Tuple{UInt, Float64, Symbol, Symbol}}
    non_zero_alg::Vector{Tuple{UInt, Float64, Symbol, Symbol}}
    is_valid::Bool
end

# ============================================================================
# Validation Function
# ============================================================================

"""
    validate_tds_initialization(adapter::TDSAdapter, am::AddressManager, t; atol=1e-6, max_show=20)
    validate_tds_initialization(prob::TDSProblem; atol=1e-6, max_show=20)

Validate TDS initialization by checking that the DAE starts in steady state.

Evaluates the residuals at the initial time and checks the steady-state conditions
`g = 0` and `f = M·rate`: the rates are zero at an equilibrium, and those of the
sinusoidal steady state otherwise (`initial_rates(prob.init.steady_state)`).

# Methods
- `(adapter, am, t)` - Use existing adapter, whose initial state is the state at `t` (the
  start of its problem, `prob.tspan[1]`; called by `solve(; validate=true)`)
- `(prob)` - Standalone validation (creates adapter internally)

# Arguments
- `adapter::TDSAdapter` - TDS adapter with initialized buffers
- `am::AddressManager` - Address manager for residual lookup
- `prob::TDSProblem` - Initialized TDS problem
- `atol::Float64` - Absolute tolerance for zero check (default: 1e-6)
- `max_show::Int` - Maximum violations to include in result (default: 20)

# Returns
`ResidualValidationResult` with diagnostics

# Example
```julia
prob = TDSProblem(setup_tds(sol_pf, dynamic_models))
result = validate_tds_initialization(prob)

if !result.is_valid
    println("Max residual: ", result.max_alg_residual)
end
```
"""
function validate_tds_initialization(
    adapter::TDSAdapter,
    am::AddressManager,
    t::Real;
    atol::Float64 = 1e-6,
    max_show::Int = 20
)
    # Evaluate at the initial time `t`: in steady state g = 0 and f = M·rate, where the rates
    # are zero at an equilibrium and those of the sinusoidal steady state otherwise
    u0 = vcat(adapter.x_buf, adapter.y_buf)
    du = similar(u0)
    evaluate!(du, adapter, u0, t, ODEInterface())

    mass = diag(adapter.mass_matrix)
    f0 = du[1:adapter.n_diff] .- mass[1:adapter.n_diff] .* adapter.rate_buf
    g0 = du[adapter.n_diff+1:end]

    # Find non-zero residuals
    non_zero_diff = Tuple{UInt, Float64, Symbol, Symbol}[]
    non_zero_alg = Tuple{UInt, Float64, Symbol, Symbol}[]

    max_diff = maximum(abs, f0; init=0.0)
    max_alg = maximum(abs, g0; init=0.0)

    # Collect differential residual violations
    for (i, val) in enumerate(f0)
        if abs(val) > atol && length(non_zero_diff) < max_show
            model, var, _ = lookup_address(am, StateRes(), UInt(i))
            push!(non_zero_diff, (UInt(i), val, model, var))
        end
    end

    # Collect algebraic residual violations
    for (i, val) in enumerate(g0)
        if abs(val) > atol && length(non_zero_alg) < max_show
            model, var, _ = lookup_address(am, AlgebRes(), UInt(i))
            push!(non_zero_alg, (UInt(i), val, model, var))
        end
    end

    is_valid = (max_diff <= atol) && (max_alg <= atol)

    return ResidualValidationResult(max_diff, max_alg, non_zero_diff, non_zero_alg, is_valid)
end

# Convenience method - creates adapter internally
function validate_tds_initialization(prob::TDSProblem; atol::Float64 = 1e-6, max_show::Int = 20)
    adcache = create_adcache(prob)
    adapter = TDSAdapter(adcache, prob)
    return validate_tds_initialization(adapter, prob.address_manager, prob.tspan[1]; atol, max_show)
end

# ============================================================================
# Reporting Functions
# ============================================================================

"""
    print_validation_report(result::ResidualValidationResult; io::IO=stdout)

Print formatted TDS initialization validation report.

Displays maximum residuals and tables of non-zero residuals (violations).
Useful for diagnosing initialization problems before running simulation.

# Arguments
- `result::ResidualValidationResult`: Validation result from `validate_tds_initialization`
- `io::IO=stdout`: Output stream (default: stdout)

# Output
- If valid: Shows "✓ All residuals within tolerance" with max values
- If invalid: Shows "✗ Non-zero residuals detected!" with tables listing
  violating equations (address, model, equation name, residual value)

# Example
```julia
prob = TDSProblem(init; tspan=(0.0, 10.0))
result = validate_tds_initialization(prob)

if !result.is_valid
    print_validation_report(result)
    # ======================================================================
    # TDS Initialization Validation Report
    # ======================================================================
    # ✗ Non-zero residuals detected!
    #   Max differential residual: 1.23e-02
    #   Max algebraic residual:    5.67e-03
    # ...
end
```

# See Also
- [`validate_tds_initialization`](@ref): Generate validation result
- [`write_validation_report`](@ref): Write report to file
"""
function print_validation_report(result::ResidualValidationResult; io::IO = stdout)
    println(io, "\n" * "="^70)
    println(io, "TDS Initialization Validation Report")
    println(io, "="^70)

    if result.is_valid
        println(io, "✓ All residuals within tolerance")
        println(io, "  Max differential residual: ", @sprintf("%.2e", result.max_diff_residual))
        println(io, "  Max algebraic residual:    ", @sprintf("%.2e", result.max_alg_residual))
        return
    end

    println(io, "✗ Non-zero residuals detected!")
    println(io, "  Max differential residual: ", @sprintf("%.2e", result.max_diff_residual))
    println(io, "  Max algebraic residual:    ", @sprintf("%.2e", result.max_alg_residual))
    println(io)

    # Print differential residual violations
    if !isempty(result.non_zero_diff)
        println(io, "Differential Residuals (f != 0):")
        println(io, "-" * "-"^69)

        data = []
        for (addr, val, model, eq) in result.non_zero_diff
            push!(data, [addr, string(model), string(eq), @sprintf("%.6e", val)])
        end

        data_matrix = permutedims(hcat(data...))
        column_labels = ["Address", "Model", "Equation", "Residual"]
        pretty_table(io, data_matrix; column_labels = column_labels)
        println(io)
    end

    # Print algebraic residual violations
    if !isempty(result.non_zero_alg)
        println(io, "Algebraic Residuals (g != 0):")
        println(io, "-" * "-"^69)

        data = []
        for (addr, val, model, eq) in result.non_zero_alg
            push!(data, [addr, string(model), string(eq), @sprintf("%.6e", val)])
        end

        data_matrix = permutedims(hcat(data...))
        column_labels = ["Address", "Model", "Equation", "Residual"]
        pretty_table(io, data_matrix; column_labels = column_labels)
        println(io)
    end

    println(io, "="^70)
end

"""
    write_validation_report(result::ResidualValidationResult, filename::String)

Write validation report to file.
"""
function write_validation_report(result::ResidualValidationResult, filename::String)
    open(filename, "w") do io
        print_validation_report(result; io=io)
    end
    println("Validation report written to $filename")
end

"""
    _build_filter_mask(sys, keys, model_name, n_instances) -> BitVector

Build boolean mask for filtering model instances based on provided keys.
Errors if a requested key does not exist.

# Arguments
- `sys` - SystemModel holding the model
- `keys` - Optional Dict mapping model names to lists of keys to show
- `model_name` - Symbol name of the model
- `n_instances` - Total number of instances

# Returns
Boolean vector where `true` indicates instance should be shown.
"""
function _build_filter_mask(sys, keys, model_name, n_instances)
    if keys === nothing || !haskey(keys, model_name)
        # No filter - show all instances
        return trues(n_instances)
    end

    mask = falses(n_instances)
    for k in keys[model_name]
        mask[get_index_by_key(sys, model_name, k)] = true
    end
    return mask
end

"""
    print_initialization_summary(prob::TDSProblem; io::IO=stdout, keys=nothing, vars=nothing, include_observable_vars=false)
    print_initialization_summary(initial_states::NamedTuple, sys::SystemModel; io::IO=stdout, keys=nothing, vars=nothing)

Print summary tables of initialized dynamic model states.

For each dynamic model, displays a table with:
- Instance identification columns (from `get_key_name`, splatted if composite)
- State variable columns (delta, omega, E_q, etc.)
- Output variables (P_e, Q_e) if `include_observable_vars=true` (TDSProblem only)

# Arguments
- `prob::TDSProblem` or `initial_states, sys` - Source of initialization data
- `io::IO=stdout` - Output stream
- `keys=nothing` - Optional filter for specific instances
  - If `nothing`, shows all instances
  - If provided, only shows instances matching the specified keys
  - Format: `Dict(:ModelName => [key1, key2, ...])` where keys match `get_key_name` format
- `vars=nothing` - Optional filter for state variables to display
  - If `nothing`, shows all state variables
  - If provided, only shows specified variables
  - Format: `Vector{Symbol}` (e.g., `[:delta, :omega]`)
- `include_observable_vars=false` - Include output variables (P_e, Q_e) by evaluating residuals (TDSProblem only)

# Examples
```julia
# Show all instances
prob = TDSProblem(setup_tds(sol_pf, dynamic_models); tspan=(0.0, 10.0))
print_initialization_summary(prob)

# Show only specific GENCLS instances by composite key (BUS, ID)
print_initialization_summary(prob; keys = Dict(:GENCLS => [(1, "1"), (2, "1")]))

# Show only specific state variables
print_initialization_summary(prob; vars = [:delta, :omega])

# Include output variables (P_e, Q_e for GENCLS)
print_initialization_summary(prob; include_observable_vars = true)

# Combine filters: specific instances and specific variables
print_initialization_summary(prob;
    keys = Dict(:GENCLS => [(1, "1"), (2, "1")]),
    vars = [:delta, :omega])

# Direct from initial_states (advanced)
print_initialization_summary(initial_states, sys)
```
"""
function print_initialization_summary end

function _print_initialization_summary_impl(
    initial_states::NamedTuple,
    sys::SystemModel,
    model_workspaces;
    io::IO = stdout,
    keys = nothing,
    vars = nothing
)
    if isempty(initial_states)
        println(io, "No dynamic models initialized.")
        return
    end

    println(io, "\n" * "="^70)
    println(io, "TDS Initialization Summary")
    println(io, "="^70)

    for (model_name, states) in pairs(initial_states)
        n_instances = length(first(states))

        if n_instances == 0
            continue
        end

        # Get model type and key names
        model = sys.models[model_name]
        model_type = get_base_model_type(model)
        key_name = get_key_name(model_type)
        metadata = model_metadata(model_type)

        # Filter instances
        row_mask = _build_filter_mask(sys, keys, model_name, n_instances)
        n_filtered = count(row_mask)
        if n_filtered == 0
            continue
        end

        println(io, "\n$model_name: $n_filtered instances")
        println(io, "-"^70)

        # Build column names and data
        col_names = Symbol[]
        col_data = Vector[]

        # Add key columns (splat if composite)
        if key_name isa Vector  # Composite key: [:BUS, :ID]
            for field in key_name
                push!(col_names, field)
                push!(col_data, getproperty(model, field)[row_mask])
            end
        else  # Single key: :i
            push!(col_names, key_name)
            push!(col_data, getproperty(model, key_name)[row_mask])
        end

        # Add state variable columns (filter by vars if provided)
        for (var_name, var_values) in pairs(states)
            # Skip if vars filter is active and this variable is not in the list
            if vars !== nothing && !(var_name in vars)
                continue
            end

            push!(col_names, var_name)
            push!(col_data, var_values[row_mask])
        end

        # Add output variables if workspace provided and model has them
        if model_workspaces !== nothing && haskey(model_workspaces, model_name)
            ws = model_workspaces[model_name]
            for output_var in metadata.observable_vars
                # Skip if vars filter is active and this variable is not in the list
                if vars !== nothing && !(output_var in vars)
                    continue
                end

                # Observable variables are stored with the residuals (ResidualBuffers)
                push!(col_names, output_var)
                push!(col_data, getfield(ws.res_vals, output_var)[row_mask])
            end
        end

        # Build data matrix for PrettyTables
        data_matrix = hcat([Vector(col) for col in col_data]...)

        # Count how many columns are keys vs states
        n_key_cols = length(key_name isa Vector ? key_name : [key_name])

        # Format state columns as floats (skip key columns)
        state_col_indices = collect((n_key_cols + 1):length(col_names))

        # Print table with float formatting for state columns only
        formatters = isempty(state_col_indices) ? nothing : [PrettyTables.fmt__printf("%.6f", state_col_indices)]
        pretty_table(io, data_matrix;
                    column_labels = col_names,
                    formatters = formatters,
                    display_size = (-1, -1))  # Show all rows and columns

        println(io)
    end

    println(io, "="^70)
end

# Public API - from initial_states (no output vars)
print_initialization_summary(initial_states::NamedTuple, sys::SystemModel; io::IO = stdout, keys = nothing, vars = nothing) =
    _print_initialization_summary_impl(initial_states, sys, nothing; io = io, keys = keys, vars = vars)

# Public API - from TDSProblem
function print_initialization_summary(
    prob::TDSProblem;
    io::IO = stdout,
    keys = nothing,
    vars = nothing,
    include_observable_vars::Bool = false
)
    if !include_observable_vars
        # Fast path: no workspace needed
        return _print_initialization_summary_impl(prob.init.dyn, prob.sys, nothing; io = io, keys = keys, vars = vars)
    end

    adcache = create_adcache(prob)
    adapter = TDSAdapter(adcache, prob)

    # Evaluate to populate output vars in workspaces, at the time of the initial state
    u0 = vcat(adapter.x_buf, adapter.y_buf)
    du = similar(u0)
    evaluate!(du, adapter, u0, prob.tspan[1], ODEInterface())

    # Print with output vars (access workspaces from ADCache)
    _print_initialization_summary_impl(prob.init.dyn, prob.sys, get_float_ws(adcache).workspaces_by_name; io = io, keys = keys, vars = vars)
end

# Export public API
export validate_tds_initialization
export print_validation_report, write_validation_report
export ResidualValidationResult
export print_initialization_summary
