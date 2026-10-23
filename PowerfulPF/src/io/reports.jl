"""
Power flow report generation with configurable units and formatting.

This module provides a flexible system for generating reports from power flow solutions.
Models declare their variables via ModelMetadata, and this module uses that metadata
to automatically collect and format data with configurable unit conversions.

# Key Features
- Unit conversions based on existing PUType system
- Configurable formatting and column selection
- Works with SystemModel and PowerFlowSolution
- Leverages Tables.jl interface for export flexibility

# Observation API
Models compute derived quantities (observations) via `compute_observations`:
- Defined in `src/foundation/model_metadata.jl` (foundation layer)
- Signature: `compute_observations(sol::PowerFlowSolution, ::Type{ModelT}, var::Symbol) -> Vector{Float64}`
- Always vectorized: returns values for ALL instances at once
- Full solution access: models can read state, residuals, parameters, etc.
- See [`compute_observations`](@ref) for implementation details
"""

using Statistics

# ============================================================================
# Configuration Types
# ============================================================================

"""
    UnitPreferences

Configuration for unit conversions in reports.

# Fields
- `power::Symbol`: Power units - `:physical` (MW/MVAR) or `:pu`
- `voltage::Symbol`: Voltage units - `:pu` (default) or `:kv`
- `angle::Symbol`: Angle units - `:deg` (default) or `:rad`
- `impedance::Symbol`: Impedance units - `:pu` or `:ohm`
- `base_mva::Float64`: System base power in MVA (default: 100.0)

# Examples
```julia
# Default: physical power, pu voltage, degrees
UnitPreferences()

# All physical units
UnitPreferences(voltage=:kv, power=:physical)

# All per-unit
UnitPreferences(power=:pu, voltage=:pu, angle=:rad)
```
"""
Base.@kwdef struct UnitPreferences
    power::Symbol = :physical      # :physical (MW/MVAR) or :pu
    voltage::Symbol = :pu           # :pu or :kv
    angle::Symbol = :deg            # :deg or :rad
    impedance::Symbol = :pu         # :pu or :ohm
    base_mva::Float64 = 100.0
end

"""
    FormattingPreferences

Configuration for value formatting in reports.

# Fields
- `float_precision::Int`: Default decimal places (default: 4)
"""
Base.@kwdef struct FormattingPreferences
    float_precision::Int = 4
end

"""
    ColumnPreferences

Configuration for column selection and naming.

# Fields
- `use_descriptions::Bool`: Use ModelParam.description vs. field name
- `include_units::Bool`: Append "(MW)", "(pu)" to column names
- `column_style::Symbol`: Naming style - `:compact` (p_from), `:short` (p_from (MW)), `:verbose` (full description)
- `include_id::Bool`: Include instance ID column
- `include_vars::Union{Nothing, Vector{Symbol}}`: If set, only include these variables
- `exclude_vars::Vector{Symbol}`: Exclude these variables
- `include_params::Bool`: Include model parameters
- `exclude_params::Vector{Symbol}`: Exclude these parameters
"""
Base.@kwdef struct ColumnPreferences
    use_descriptions::Bool = true
    include_units::Bool = true
    column_style::Symbol = :short  # :compact, :short, :verbose
    include_id::Bool = true
    include_vars::Union{Nothing, Vector{Symbol}} = nothing
    exclude_vars::Vector{Symbol} = Symbol[]
    include_params::Bool = true
    exclude_params::Vector{Symbol} = Symbol[]
end

"""
    ReportConfig

Complete configuration for power flow report generation.

# Fields
- `units::UnitPreferences`: Unit conversion preferences
- `formatting::FormattingPreferences`: Value formatting preferences
- `columns::ColumnPreferences`: Column selection and naming

# Examples
```julia
# Default configuration
config = ReportConfig()

# Customize units only
config = ReportConfig(units=UnitPreferences(voltage=:kv, angle=:rad))

# Customize multiple aspects
config = ReportConfig(
    units = UnitPreferences(voltage=:kv, power=:pu),
    formatting = FormattingPreferences(float_precision=6),
    columns = ColumnPreferences(exclude_vars=[:internal_state])
)

# Use preset configurations
config = AllPerUnit
config = AllPhysical
```
"""
Base.@kwdef struct ReportConfig
    units::UnitPreferences = UnitPreferences()
    formatting::FormattingPreferences = FormattingPreferences()
    columns::ColumnPreferences = ColumnPreferences()
end

# Allow partial updates via keyword arguments
function ReportConfig(base::ReportConfig;
                     units=base.units,
                     formatting=base.formatting,
                     columns=base.columns)
    ReportConfig(units, formatting, columns)
end

# ============================================================================
# Preset Configurations
# ============================================================================

"""Default configuration: MW, pu voltage, degrees (excludes vm/va initial values)"""
const Default = ReportConfig(
    columns = ColumnPreferences(
        exclude_params = [:vm, :va]  # Exclude initial values - show only computed solution
    )
)

"""All per-unit: power, voltage in pu; angles in radians"""
const AllPerUnit = ReportConfig(
    units = UnitPreferences(
        power = :pu,
        voltage = :pu,
        angle = :rad
    )
)

"""All physical units: MW, kV, degrees"""
const AllPhysical = ReportConfig(
    units = UnitPreferences(
        power = :physical,
        voltage = :kv,
        angle = :deg
    )
)

"""Minimal report: only variables, no parameters"""
const Minimal = ReportConfig(
    columns = ColumnPreferences(
        include_params = false,
        use_descriptions = false
    )
)

"""Compact format: fewer decimal places"""
const Compact = ReportConfig(
    formatting = FormattingPreferences(
        float_precision = 2
    )
)

# ============================================================================
# Unit Conversion Functions
# ============================================================================

"""
    to_report_units(pu_type::PUType, value, config::ReportConfig[, base_kv]) -> converted_value

Convert a per-unit value to report units based on configuration and PUType.

For observations (computed values with no stored coefficients):
- Most types only need config (ActivePower, ReactivePower, AngleDegrees, NoPU)
- VoltageNominal needs base_kv for kV conversion

# Arguments
- `pu_type`: The PUType from observation metadata
- `value`: The value to convert (in per-unit or internal representation)
- `config`: Report configuration specifying desired units
- `base_kv`: Base voltage in kV (only for VoltageNominal)

# Examples
```julia
config = ReportConfig(units=UnitPreferences(power=:physical, base_mva=100.0))
to_report_units(ActivePower(), 1.5, config)  # Returns 150.0 MW

config = ReportConfig(units=UnitPreferences(voltage=:kv))
to_report_units(VoltageNominal(), 1.03, config, 230.0)  # Returns 236.9 kV
```
"""
function to_report_units end

# Active power: pu → MW (observations only)
function to_report_units(::ActivePower, pu_val, config::ReportConfig)
    config.units.power == :physical ? pu_val * config.units.base_mva : pu_val
end

# Reactive power: pu → MVAR (observations only)
function to_report_units(::ReactivePower, pu_val, config::ReportConfig)
    config.units.power == :physical ? pu_val * config.units.base_mva : pu_val
end

# Voltage: pu → kV (observations only, needs base_kv)
function to_report_units(::VoltageNominal, pu_val, config::ReportConfig, base_kv)
    config.units.voltage == :kv ? pu_val * base_kv : pu_val
end

# Angle: radians → degrees (observations only)
function to_report_units(::AngleDegrees, rad_val, config::ReportConfig)
    config.units.angle == :deg ? rad2deg(rad_val) : rad_val
end

# No conversion needed
to_report_units(::NoPU, val, config::ReportConfig) = val

# Default fallback (for any future observation types)
to_report_units(::PUType, val, config::ReportConfig) = val

# ============================================================================
# Column Naming
# ============================================================================

"""
    _capitalize_var_name(s::String) -> String

Capitalize variable names for formal presentation.
Examples: "p_from" → "P_from", "theta" → "Theta", "basekv" → "Basekv"
"""
function _capitalize_var_name(s::String)
    isempty(s) && return s
    return uppercase(s[1:1]) * s[2:end]
end

"""
    get_column_name(var_name::Symbol, param_meta, config::ReportConfig) -> String

Get column name based on configuration preferences.
Uses metadata description if available and appends units suffix if requested.
"""
function get_column_name(var_name::Symbol,
                        param_meta::Union{ModelParam, ModelObservable, Nothing},
                        config::ReportConfig)
    # Handle column_style
    if config.columns.column_style == :compact
        # Compact: just variable name, no units, capitalized
        return _capitalize_var_name(string(var_name))
    end

    # Determine base name
    base_name = if config.columns.column_style == :verbose &&
                   config.columns.use_descriptions && !isnothing(param_meta) &&
                   !isempty(param_meta.description)
        # Verbose: use full description
        param_meta.description
    else
        # Short: use capitalized variable name
        _capitalize_var_name(string(var_name))
    end

    # Add units suffix if requested and metadata available
    if !config.columns.include_units || isnothing(param_meta)
        return base_name
    end

    # Determine units suffix based on PUType and config
    suffix = if param_meta.pu_type isa ActivePower
        config.units.power == :physical ? " (MW)" : " (pu)"
    elseif param_meta.pu_type isa ReactivePower
        config.units.power == :physical ? " (MVAR)" : " (pu)"
    elseif param_meta.pu_type isa PowerNominal
        # Active or reactive power: the unit the parameter declares (MW or MVAR)
        isempty(param_meta.units) && error("Power parameter $var_name declares no units (MW or MVAR)")
        config.units.power == :physical ? " (" * param_meta.units * ")" : " (pu)"
    elseif param_meta.pu_type isa VoltageNominal
        config.units.voltage == :kv ? " (kV)" : " (pu)"
    elseif param_meta.pu_type isa AngleDegrees
        config.units.angle == :deg ? " (deg)" : " (rad)"
    elseif param_meta.pu_type isa ImpedanceNominal
        config.units.impedance == :ohm ? " (Ω)" : " (pu)"
    elseif param_meta.pu_type isa AdmittanceNominal
        config.units.impedance == :ohm ? " (S)" : " (pu)"
    elseif param_meta.pu_type isa CurrentNominal
        config.units.power == :physical ? " (A)" : " (pu)"
    elseif !isempty(param_meta.units)
        " (" * param_meta.units * ")"
    else
        ""
    end

    return base_name * suffix
end

# ============================================================================
# Discovery API - Output Variables
# ============================================================================

"""
    get_available_outputs(model_type::Type{<:AbstractModel}) -> Vector{Symbol}

Get all available output variables for a model (DAE variables + observations).
These are the variables that can be included in reports via `collect_model_table`.

# Example
```julia
julia> get_available_outputs(Branch)
6-element Vector{Symbol}:
 :p_from
 :q_from
 :p_to
 :q_to
 :p_loss
 :q_loss

# Use in report generation
config = ReportConfig(columns=ColumnPreferences(include_vars=[:p_from, :p_loss]))
report = collect_model_table(sys, :Branch; config=config)
```
"""
function get_available_outputs(model_type::Type{<:AbstractModel})
    metadata = model_metadata(model_type)
    outputs = Symbol[]

    # DAE output variables
    append!(outputs, metadata.observable_vars)

    # Computed observations (skip if already in observable_vars - those are metadata-only)
    for obs in metadata.observations
        obs.name in metadata.observable_vars && continue
        push!(outputs, obs.name)
    end

    return outputs
end

"""
    get_default_outputs(model_type::Type{<:AbstractModel}) -> Vector{Symbol}

Get the default output variables for a model (those shown in reports by default).
This is a subset of `get_available_outputs` - observations with `default_output=false`
are excluded unless explicitly requested via `include_vars`.

# Example
```julia
julia> get_default_outputs(Branch)
4-element Vector{Symbol}:
 :p_from
 :q_from
 :p_to
 :q_to

# Losses are available but not shown by default
julia> get_available_outputs(Branch)
6-element Vector{Symbol}:
 :p_from
 :q_from
 :p_to
 :q_to
 :p_loss  # Available but not default
 :q_loss  # Available but not default
```
"""
function get_default_outputs(model_type::Type{<:AbstractModel})
    metadata = model_metadata(model_type)
    defaults = Symbol[]

    # DAE output variables (always included by default)
    append!(defaults, metadata.observable_vars)

    # Observations marked as default_output=true (skip if already in observable_vars - those are metadata-only)
    for obs in metadata.observations
        obs.name in metadata.observable_vars && continue
        if obs.default_output
            push!(defaults, obs.name)
        end
    end

    return defaults
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
    _find_param_for_var(metadata::ModelMetadata, var_name::Symbol) -> Union{ModelParam, Nothing}

Find parameter metadata for a variable name.
Returns nothing if no matching parameter is found.
"""
function _find_param_for_var(metadata::ModelMetadata, var_name::Symbol)
    for param in metadata.params
        param.name == var_name && return param
    end
    return nothing
end

"""
    _find_observation_for_var(metadata::ModelMetadata, var_name::Symbol) -> Union{ModelObservable, Nothing}

Find observation metadata for a variable name.
Returns nothing if no matching observation is found.
"""
function _find_observation_for_var(metadata::ModelMetadata, var_name::Symbol)
    for obs in metadata.observations
        obs.name == var_name && return obs
    end
    return nothing
end

# Base quantity resolution now handled by trait system in base_quantity_traits.jl
# See get_base_quantity_vector(sys, model_name, BaseVoltage)

# ============================================================================
# Main Collection Function (SystemModel)
# ============================================================================

"""
    observation_defined(formulation, pu_type) -> Bool

Whether a quantity of type `pu_type` is solved in `formulation`. `DCPowerFlow` has no
reactive power and does not solve voltage magnitudes, so reports leave those columns
out.
"""
observation_defined(::FormulationTrait, ::PUType) = true
observation_defined(::DCPowerFlow, ::Union{ReactivePower, VoltageNominal}) = false

"""
    collect_model_table(sol::PowerFlowSolution, model_name::Symbol;
                       config::ReportConfig=Default) -> Union{NamedTuple, Nothing}

Collect model data with full configuration control.
Returns a Tables.jl-compatible NamedTuple of columns.

# Arguments
- `sol`: PowerFlowSolution from solve_powerflow
- `model_name`: Model to collect (e.g., `:Bus`, `:PV`, `:Branch`)
- `config`: Report configuration (default: `Default`)

# Returns
NamedTuple with columns for variables and parameters, or `nothing` if model not found.

# Examples
```julia
# Solve and collect data
prob = ACPowerFlowProblem(sys)
sol = solve(prob, RobustNR())
bus_data = collect_model_table(sol, :Bus)

# Custom units
config = ReportConfig(units=UnitPreferences(voltage=:kv))
bus_data = collect_model_table(sol, :Bus, config=config)

# Use preset
bus_data = collect_model_table(sol, :Bus, config=AllPhysical)

# Export to DataFrame
using DataFrames
df = DataFrame(collect_model_table(sol, :Bus))

# Export to CSV
using CSV
CSV.write("buses.csv", collect_model_table(sol, :Bus))
```
"""

"""
    _convert_to_report_units(vals::Vector, param_name::Symbol, pu_type::PUType,
                             sys::SystemModel, model_name::Symbol, config::ReportConfig, is_observation::Bool)

Convert per-unit values to report units.
- **Parameters**: Use stored coefficients (exact inverse: `physical = pu / coeff`)
- **Observations**: Use PUType-based conversion (computed values have no coefficients)
"""
function _convert_to_report_units(vals::Vector, param_name::Symbol, pu_type::PUType,
                                  sys::SystemModel, model_name::Symbol,
                                  config::ReportConfig, is_observation::Bool)
    if is_observation
        # Observations are computed values - no coefficients exist
        # Use PUType-based conversion
        if pu_type isa VoltageNominal
            # VoltageNominal needs per-instance base_kv for kV conversion
            # Use trait system to resolve base voltage (no hard-coded field checks)
            base_kv_vec = get_base_quantity_vector(sys, model_name, BaseVoltage; quantity = param_name)
            return [to_report_units(pu_type, vals[i], config, base_kv_vec[i]) for i in eachindex(vals)]
        else
            # Other observation types only need config
            return [to_report_units(pu_type, v, config) for v in vals]
        end
    else
        # Parameters - must have coefficients in modern workflow
        if isnothing(sys.pu_coeffs) || !haskey(sys.pu_coeffs, model_name) ||
           !haskey(sys.pu_coeffs[model_name], param_name)
            error("Missing coefficients for parameter :$param_name in model :$model_name. " *
                  "System must be created with per-unit conversion enabled.")
        end

        # Use coefficient-based conversion (exact inverse of forward conversion)
        coeffs = sys.pu_coeffs[model_name][param_name]
        return [vals[i] / coeffs[i] for i in eachindex(vals)]
    end
end

function collect_model_table(sol::PowerFlowSolution, model_name::Symbol;
                             config::ReportConfig=Default)
    sys = sol.prob.sys
    y = sol.dae.y

    !haskey(sys.models, model_name) && return nothing

    model = sys.models[model_name]
    base_type = get_base_model_type(model)
    metadata = model_metadata(base_type)

    n = length(model)
    n == 0 && return nothing

    # Build columns in order: keys → observations → params
    col_names = Symbol[]
    col_data = Vector[]

    # Instance IDs (if requested) - unpack composite keys
    if config.columns.include_id
        keys = get_key(model)
        key_name = get_key_name(base_type)
        key_fields = key_name isa Symbol ? (key_name,) : key_name

        # One column per key field (keys are tuples), converted to String
        for (idx, field) in enumerate(key_fields)
            push!(col_names, Symbol(_capitalize_var_name(string(field))))
            push!(col_data, [String(string(k[idx])) for k in keys])
        end
    end

    # Validate user input upfront
    available = get_available_outputs(base_type)

    if !isnothing(config.columns.include_vars)
        # Validate include_vars
        requested = Set{Symbol}(config.columns.include_vars)
        invalid = setdiff(requested, available)
        if !isempty(invalid)
            error("Invalid variables in include_vars for :$model_name: $(sort(collect(invalid))). " *
                  "Available: $(sort(available))")
        end
    end

    if !isempty(config.columns.exclude_vars)
        # Validate exclude_vars
        invalid = setdiff(config.columns.exclude_vars, available)
        if !isempty(invalid)
            error("Invalid variables in exclude_vars for :$model_name: $(sort(collect(invalid))). " *
                  "Available: $(sort(available))")
        end
    end

    # Determine which variables to include based on column preferences
    vars_to_include = if !isnothing(config.columns.include_vars)
        # User explicitly specified variables - already validated above
        Set{Symbol}(config.columns.include_vars)
    else
        # No explicit specification - use defaults
        default_vars = Set{Symbol}()

        # Include all DAE output variables (always shown by default)
        union!(default_vars, metadata.observable_vars)

        # Include observations marked as default_output=true
        for obs in metadata.observations
            if obs.default_output
                push!(default_vars, obs.name)
            end
        end

        # Apply exclusions (already validated above)
        filter(v -> !(v in config.columns.exclude_vars), default_vars)
    end

    # Process DAE variables from observable_vars (in metadata order)
    formulation = sol.prob.formulation
    for var_name in metadata.observable_vars
        var_name in vars_to_include || continue

        # Find metadata for unit conversion (observations for metadata)
        obs_meta = _find_observation_for_var(metadata, var_name)
        param_meta = _find_param_for_var(metadata, var_name)
        var_meta = !isnothing(obs_meta) ? obs_meta : param_meta
        !isnothing(var_meta) && !observation_defined(formulation, var_meta.pu_type) && continue

        # Values in any formulation (e.g. bus voltages from the rectangular components)
        vals = compute_observations(sol, base_type, var_name)

        if !isnothing(var_meta)
            # DAE output variables are observations (e.g., Bus.theta, Bus.v)
            vals = _convert_to_report_units(vals, var_name, var_meta.pu_type, sys, model_name, config, true)
            col_name = Symbol(get_column_name(var_name, var_meta, config))
            push!(col_names, col_name)
            push!(col_data, vals)
        else
            @warn "No metadata found for DAE variable :$var_name in model $model_name - skipping"
        end
    end

    # Process computed observations (in metadata order)
    for obs in metadata.observations
        obs.name in vars_to_include || continue
        obs.name in metadata.observable_vars && continue  # Skip if already processed as DAE variable
        observation_defined(formulation, obs.pu_type) || continue

        # Computed observation - call model's compute_observations (vectorized)
        model_type = get_base_model_type(model)

        # Compute all observations at once (vectorized, no iteration)
        vals = compute_observations(sol, model_type, obs.name)

        # Convert to report units (observations are computed values)
        vals = _convert_to_report_units(vals, obs.name, obs.pu_type, sys, model_name, config, true)
        col_name = Symbol(get_column_name(obs.name, obs, config))
        push!(col_names, col_name)
        push!(col_data, vals)
    end

    # Parameters (if requested)
    if config.columns.include_params
        params_to_include = filter(p -> !(p.name in config.columns.exclude_params),
                                   metadata.params)

        for param in params_to_include
            # Every metadata parameter is a model field
            vals = collect(getproperty(model, param.name))

            # Convert units using coefficients
            vals = _convert_to_report_units(vals, param.name, param.pu_type, sys, model_name, config, false)

            # Get column name
            col_name = Symbol(get_column_name(param.name, param, config))
            push!(col_names, col_name)
            push!(col_data, vals)
        end
    end

    return NamedTuple(zip(col_names, col_data))
end

"""
    collect_model_table(sol::MultiIslandSolution, model_name::Symbol; config=Default)

Collect model data from a multi-island solution, merging results from all islands.

Returns a unified table containing data from all islands. Instances are included
from all islands based on their external keys.

**Note**: `MultiIslandSolution` is defined in `src/solvers/solution_types.jl`.
This function handles the I/O-specific concern of merging data tables across islands.

# Arguments
- `sol::MultiIslandSolution`: Multi-island power flow solution
- `model_name::Symbol`: Name of model type (e.g., `:Bus`, `:Branch`)
- `config::ReportConfig`: Configuration for units, formatting, and columns

# Example
```julia
sol = solve(ACPowerFlowProblem(sys), RobustNR())  # May return MultiIslandSolution
bus_data = collect_model_table(sol, :Bus)
df = DataFrame(bus_data)
```
"""
function collect_model_table(sol::MultiIslandSolution, model_name::Symbol;
                             config::ReportConfig=Default)
    # Collect data from each island
    island_tables = NamedTuple[]

    for (island_sys, island_sol) in zip(sol.islands, sol.solutions)
        # Check if this island has the requested model
        if haskey(island_sys.models, model_name)
            # Collect data from this island's solution
            table = collect_model_table(island_sol, model_name; config=config)
            if !isnothing(table)
                push!(island_tables, table)
            end
        end
    end

    # If no islands have this model, return nothing
    isempty(island_tables) && return nothing

    # Merge all island tables
    # All tables should have same column names (from same model type)
    col_names = keys(island_tables[1])
    merged_cols = []

    for col_name in col_names
        # Concatenate column data from all islands
        col_data = vcat([table[col_name] for table in island_tables]...)
        push!(merged_cols, col_data)
    end

    return NamedTuple(zip(col_names, merged_cols))
end

# ============================================================================
# Pretty Printing with PrettyTables.jl
# ============================================================================

using PrettyTables

"""
    format_table_pretty(io::IO, table::NamedTuple; config=Default, kwargs...)

Helper function to format a NamedTuple table with PrettyTables.
Applies rounding based on config.formatting.float_precision.

Used internally by both print_model_table() and format_text().
"""
function format_table_pretty(io::IO, table::NamedTuple; config::ReportConfig=Default, backend::Symbol=:text, kwargs...)
    isnothing(table) && return nothing

    # Apply formatting (rounding)
    if config.formatting.float_precision > 0
        prec = config.formatting.float_precision
        rounded_data = map(collect(pairs(table))) do (key, vals)
            if eltype(vals) <: AbstractFloat
                key => round.(vals, digits=prec)
            else
                key => vals
            end
        end
        table = NamedTuple(rounded_data)
    end

    # Pretty print - display_size only valid for text backend
    if backend == :text
        PrettyTables.pretty_table(io, table;
            backend = backend,
            display_size = (-1, -1),
            kwargs...)
    else
        PrettyTables.pretty_table(io, table;
            backend = backend,
            kwargs...)
    end

    return nothing
end

"""
    print_model_table(sol, model_name::Symbol; config=Default, kwargs...)

Pretty-print model data as formatted table using Pretty Tables.jl.

This function activates the `config.formatting` preferences (float_precision, etc.)
which are otherwise unused by `collect_model_table()`.

# Arguments
- `sol`: PowerFlowSolution or MultiIslandSolution
- `model_name::Symbol`: Model to print (e.g., `:Bus`, `:Branch`)
- `config::ReportConfig`: Report configuration
- `kwargs...`: Passed to PrettyTables.jl

# Examples
```julia
# Terminal output (default)
print_model_table(sol, :Bus)

# Customized precision using config
config = ReportConfig(formatting = FormattingPreferences(float_precision = 2))
print_model_table(sol, :Bus; config=config)

# Compact format
print_model_table(sol, :Branch; config=Compact)
```
"""
function print_model_table(sol::Union{PowerFlowSolution, MultiIslandSolution},
                          model_name::Symbol;
                          config::ReportConfig=Default,
                          kwargs...)
    # Collect data using existing API
    data = collect_model_table(sol, model_name; config=config)
    isnothing(data) && return nothing

    # Delegate to formatting helper
    format_table_pretty(stdout, data; config=config, kwargs...)

    return nothing
end


# ============================================================================
# Exports
# ============================================================================

export ReportConfig, UnitPreferences, FormattingPreferences, ColumnPreferences
export Default, AllPerUnit, AllPhysical, Minimal, Compact
export collect_model_table, print_model_table
