# ============================================================================
# Configuration Types
# ============================================================================
#
# Core configuration structs for Powerful.jl, separated for early loading.
# These types have no dependencies on traits or registry.
# ============================================================================

# ============================================================================
# Simulation Config Container
# ============================================================================

"""
    SimConfig{S}

Configuration of how a system is simulated: the solver settings.
System-wide settings (frequency and MVA base) are owned by the `SystemModel`
(`sys.settings`, see `SystemConfig`), not by `SimConfig`.

Type-stable via a NamedTuple type parameter for the solver configs.

# Type Parameters
- `S<:NamedTuple`: Solver configs, e.g., `(TDS = TDSConfig(),)`

# Fields
- `solver::S`: Solver-name -> settings struct (fully type-stable)

# Construction
```julia
# Direct construction
config = SimConfig(solver = (TDS = TDSConfig(convert_pq_to_zip = true),))

# From TOML file
config = load_config("simulation.toml")

# Default (all defaults)
config = SimConfig()
```

# TOML Format
```toml
[solver.TDS]
convert_pq_to_zip = true
```
"""
struct SimConfig{S<:NamedTuple}
    solver::S
end

# Constructor with keyword arguments (also serves as default)
SimConfig(; solver::NamedTuple = (;)) = SimConfig(solver)

# ============================================================================
# Exports
# ============================================================================

export SimConfig
