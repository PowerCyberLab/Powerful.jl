# ============================================================================
# Configuration System
# ============================================================================
#
# Configuration for Powerful.jl:
# - System-wide settings (frequency base, MVA base) are owned by the SystemModel
#   (`sys.settings`, a SystemConfig defined in Foundation)
# - SimConfig holds how the system is simulated: solver-specific options (TDS,
#   power flow, etc.)
#
# Config is validated against Julia struct definitions (no separate schema).
#
# NOTE: SimConfig is defined in config_types.jl (included early in PowerfulCore
# for use in DAEInfra).
# ============================================================================

using TOML

# ============================================================================
# Config Access API (Type-Stable)
# ============================================================================

"""
    get_solver_config(config::SimConfig, solver_name::Symbol, ::Type{T}) where T

Get config for a solver. Fully type-stable via NamedTuple field access.

# Arguments
- `config`: SimConfig instance
- `solver_name`: Solver identifier (e.g., :TDS)
- `T`: Config type (e.g., TDSConfig)

# Returns
- Config struct instance, fully typed
- Default config if not configured

# Example
```julia
tds_config = get_solver_config(config, :TDS, TDSConfig)
```
"""
@inline function get_solver_config(config::SimConfig, solver_name::Symbol, ::Type{T}) where T
    if hasproperty(config.solver, solver_name)
        return getproperty(config.solver, solver_name)
    else
        return T()
    end
end

# ============================================================================
# TOML Loading with Validation
# ============================================================================

"""
    load_config(path::String) -> SimConfig

Load configuration from a TOML file with validation.

Validates all keys against Julia struct definitions and provides helpful
error messages with suggestions for typos.

# Arguments
- `path`: Path to TOML configuration file

# Returns
- `SimConfig` with validated settings

# Errors
- Unknown top-level sections
- Unknown keys within sections (with typo suggestions)
- Type conversion errors

# Example
```julia
config = load_config("simulation.toml")
```
"""
function load_config(path::String)
    toml = TOML.parsefile(path)
    load_config_from_dict(toml, path)
end

"""
    load_config_from_string(toml_string::String) -> SimConfig

Load simulation configuration from a TOML string.

Parses TOML content and returns a validated `SimConfig`. Useful for testing
and programmatic configuration without file I/O.

# Arguments
- `toml_string::String`: TOML-formatted configuration string

# Returns
- `SimConfig`: Validated configuration with solver settings

# Errors
- Unknown section names or keys (with typo suggestions)
- Type conversion errors

# Example
```julia
config = load_config_from_string(\"\"\"
    [solver.TDS]
    abstol = 1e-8
    convert_pq_to_zip = true
    p_zip = [0.3, 0.3, 0.4]
\"\"\")

# Access settings
config.solver.TDS.abstol  # 1e-8
```

# See Also
- [`load_config`](@ref): Load from file
- [`SimConfig`](@ref): Configuration struct
"""
function load_config_from_string(toml_string::String)
    toml = TOML.parse(toml_string)
    load_config_from_dict(toml, "<string>")
end

"""
Internal: Load config from parsed TOML dict.
"""
function load_config_from_dict(toml::Dict, source::String)
    # Validate top-level sections
    haskey(toml, "system") &&
        error("Section [system] in $source: system settings belong to the SystemModel, " *
              "e.g. SystemModel(case; frequency_base = 50.0)")
    valid_sections = Set(["solver"])
    for key in keys(toml)
        if key ∉ valid_sections
            error("Unknown config section '$key' in $source. Valid sections: $(join(sort(collect(valid_sections)), ", "))")
        end
    end

    # Parse solver settings
    solver = parse_solver_section(get(toml, "solver", Dict()), source)

    SimConfig(solver)
end

"""
    parse_section(::Type{T}, dict::Dict, section_name::String, source::String) where T

Parse a TOML section into a struct with validation.

Validates all keys against struct field names and suggests corrections for typos.
"""
function parse_section(::Type{T}, dict::Dict, section_name::String, source::String) where T
    isempty(dict) && return T()

    valid_keys = Set(string.(fieldnames(T)))
    provided_keys = Set(keys(dict))

    # Check for unknown keys
    unknown = setdiff(provided_keys, valid_keys)
    if !isempty(unknown)
        suggestions = [find_closest_match(k, valid_keys) for k in unknown]
        error_lines = ["  '$k' - did you mean '$s'?" for (k, s) in zip(unknown, suggestions)]
        error("""
            Unknown settings in [$section_name] in $source:
            $(join(error_lines, "\n"))

            Valid keys: $(join(sort(collect(valid_keys)), ", "))
            """)
    end

    # Build kwargs for constructor
    kwargs = Dict{Symbol, Any}()
    for name in fieldnames(T)
        key = string(name)
        if haskey(dict, key)
            field_type = fieldtype(T, name)
            try
                kwargs[name] = convert_toml_value(field_type, dict[key])
            catch e
                error("Cannot convert value '$(dict[key])' to $field_type for '$key' in [$section_name] in $source: " *
                      sprint(showerror, e))
            end
        end
    end

    T(; kwargs...)
end

"""
    convert_toml_value(::Type{F}, value) -> F

Convert a parsed TOML value to the config field type `F`.

TOML has no tuple type: arrays parse as `Vector`. A tuple-typed field (e.g.
`NTuple{3,Float64}`) therefore takes a TOML array, which must have exactly as many
elements as the tuple; each element is converted to its tuple slot's type.
"""
convert_toml_value(::Type{F}, value) where {F} = convert(F, value)

function convert_toml_value(::Type{F}, value::AbstractVector) where {F<:Tuple}
    if isconcretetype(F)
        n = fieldcount(F)
        length(value) == n ||
            throw(ArgumentError("expected an array of $n elements for $F, got $(length(value))"))
    end
    return convert(F, Tuple(value))
end

"""
    parse_solver_section(dict::Dict, source::String) -> NamedTuple

Parse the [solver] section of TOML config.

Solver config types are those registered with `register_solver_config!`.
Returns a NamedTuple for type stability.
"""
function parse_solver_section(dict::Dict, source::String)
    isempty(dict) && return (;)

    registry = SOLVER_CONFIG_REGISTRY
    keys_vec = Symbol[]
    vals_vec = Any[]

    for (solver_name, config_dict) in dict
        sym = Symbol(solver_name)

        if !haskey(registry, sym)
            valid = sort(collect(keys(registry)))
            if isempty(valid)
                error("Unknown solver '$solver_name' in [solver] in $source. No solvers have config defined.")
            else
                suggestion = find_closest_match(solver_name, Set(string.(valid)))
                error("""
                    Unknown solver '$solver_name' in [solver] in $source.
                    Did you mean '$suggestion'?

                    Solvers with config: $(join(valid, ", "))
                    """)
            end
        end

        T = registry[sym]
        push!(keys_vec, sym)
        push!(vals_vec, parse_section(T, config_dict, "solver.$solver_name", source))
    end

    return NamedTuple{Tuple(keys_vec)}(Tuple(vals_vec))
end

# ============================================================================
# Config Registry
# ============================================================================

# Solver config registry (populated by solver modules)
const SOLVER_CONFIG_REGISTRY = Dict{Symbol, Type}()

"""
    register_solver_config!(name::Symbol, ::Type{T}) where T

Register a solver config type.

Called by solver modules to make their config available for TOML loading.

# Example
```julia
# In tds.jl
register_solver_config!(:TDS, TDSConfig)
```
"""
function register_solver_config!(name::Symbol, ::Type{T}) where T
    SOLVER_CONFIG_REGISTRY[name] = T
end

# ============================================================================
# Typo Detection Utilities
# ============================================================================

"""
    find_closest_match(typo::String, valid::Set{String}) -> String

Find the closest matching string using Levenshtein distance.
"""
function find_closest_match(typo::String, valid::Set{String})
    isempty(valid) && return typo

    best_match = first(valid)
    best_distance = levenshtein_distance(typo, best_match)

    for candidate in valid
        dist = levenshtein_distance(typo, candidate)
        if dist < best_distance
            best_distance = dist
            best_match = candidate
        end
    end

    best_match
end

"""
    levenshtein_distance(s1::String, s2::String) -> Int

Compute the Levenshtein (edit) distance between two strings.
"""
function levenshtein_distance(s1::String, s2::String)
    m, n = length(s1), length(s2)

    # Early termination for empty strings
    m == 0 && return n
    n == 0 && return m

    # Use two rows instead of full matrix for memory efficiency
    prev_row = collect(0:n)
    curr_row = Vector{Int}(undef, n + 1)

    for i in 1:m
        curr_row[1] = i
        for j in 1:n
            cost = s1[i] == s2[j] ? 0 : 1
            curr_row[j + 1] = min(
                prev_row[j + 1] + 1,      # Deletion
                curr_row[j] + 1,          # Insertion
                prev_row[j] + cost        # Substitution
            )
        end
        prev_row, curr_row = curr_row, prev_row
    end

    prev_row[n + 1]
end

# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, config::SimConfig)
    print(io, "SimConfig(solver=$(length(keys(config.solver))) solvers)")
end

function Base.show(io::IO, ::MIME"text/plain", config::SimConfig)
    println(io, "SimConfig:")

    solver_keys = keys(config.solver)
    if !isempty(solver_keys)
        println(io, "  Solver:")
        for name in solver_keys
            settings = getproperty(config.solver, name)
            println(io, "    $name: $settings")
        end
    end
end

# ============================================================================
# Exports
# ============================================================================

export load_config, load_config_from_string
export get_solver_config
export register_solver_config!
