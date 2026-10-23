# ============================================================================
# System Settings
# ============================================================================

"""
    SystemConfig{T}

System-wide settings of a power system model, owned by its `SystemModel`
(`sys.settings`) so that every analysis of the system uses the same values. Solver
options live in `SimConfig` instead.

# Type Parameter
- `T<:AbstractFloat`: Float type matching simulation (Float64, Float32)

# Fields
- `frequency_base::T`: System frequency base in Hz (default: 60.0)
- `mva_base::T`: System MVA base (default: 100.0)
- `omega_base::T`: Angular frequency = 2pi * frequency_base (computed)
- `constant_power_vmin::T`: Voltage (pu) below which constant-power load (PQ `pl`,
  `ql`) is reduced smoothly to zero at zero voltage (default: 0.7, the default of
  PSS/E PQBRAK and of PowerWorld's minimum voltage for constant power loads); see
  `constant_power_factor`

# Example
```julia
sys = SystemModel(case)                                 # settings from the case (base MVA, frequency)
sys = SystemModel(case; frequency_base = 50.0)          # override
sys = SystemModel(case; constant_power_vmin = 0.8)      # constant-power breakpoint
sys.settings.omega_base
```
"""
struct SystemConfig{T<:AbstractFloat}
    frequency_base::T
    mva_base::T
    omega_base::T
    constant_power_vmin::T
end

# Parameterized constructor with defaults
function SystemConfig{T}(;
    frequency_base::T = T(60.0),
    mva_base::T = T(100.0),
    omega_base::T = T(2 * pi) * frequency_base,
    constant_power_vmin::T = T(0.7)
) where T<:AbstractFloat
    frequency_base > 0 || error("frequency_base must be positive, got $frequency_base")
    mva_base > 0 || error("mva_base must be positive, got $mva_base")
    0 < constant_power_vmin <= 2 ||
        error("constant_power_vmin must be in (0, 2] pu, got $constant_power_vmin")
    SystemConfig{T}(frequency_base, mva_base, omega_base, constant_power_vmin)
end

# Default to Float64
SystemConfig(; kwargs...) = SystemConfig{Float64}(; kwargs...)
