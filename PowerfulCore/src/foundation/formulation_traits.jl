"""
    Formulation traits for power flow equations.

Power flow problems can be formulated using different equation types.
The formulation is a property of the Problem, not the System.
"""

# ============================================================================
# Formulation Trait Hierarchy
# ============================================================================

"""
    FormulationTrait

Abstract base for power flow equation formulations.

Formulations determine the mathematical representation of network equations:
- `PowerInjection`: Power balance equations S = P + jQ (default, nonlinear)
- `CurrentInjection`: Current balance equations I = Ir + j*Ii (nonlinear)
- `DCPowerFlow`: Active power balance P = 0 only (linear approximation)

The formulation affects:
1. What equations are solved (residuals)
2. What variables are used (polar vs rectangular vs angles-only)
3. How kernels compute contributions (nonlinear vs linear)
4. Solver selection (Newton vs direct linear solve)
"""
abstract type FormulationTrait end

"""
    PowerInjection <: FormulationTrait

Power injection formulation (default, existing implementation).

# Variables
- Polar coordinates: (θ, V) - angle and magnitude

# Residuals
- Power balance at each bus: P = 0, Q = 0
- Branch contributes: P_i, Q_i, P_j, Q_j (power flows)
- Load contributes: P, Q (consumption)
- Generator contributes: P (injection), V_setpoint (control)

# Characteristics
- Most common formulation in power flow
- Well-understood convergence properties
- Trigonometric functions in Jacobian
"""
struct PowerInjection <: FormulationTrait end

"""
    CurrentInjection <: FormulationTrait

Current injection formulation (alternative, better numerics).

# Variables
- Rectangular coordinates: (Vr, Vi) - real and imaginary voltage

# Residuals
- Current balance at each bus: Ir = 0, Ii = 0
- Branch contributes: I_r_i, I_i_i, I_r_j, I_i_j (current flows)
- Load contributes: I_r, I_i (consumption as current)
- Generator contributes: I_r (from P), V_mag (control)

# Characteristics
- Better numerical properties for ill-conditioned systems
- No trigonometric functions in Jacobian (algebraic only)
- Natural for rectangular voltage coordinates
- I = Y*V is linear in V (rectangular)

# Benefits
- Better convergence near voltage collapse
- Simpler derivatives (no sin/cos)
- Can leverage Y-bus sparsity more directly
"""
struct CurrentInjection <: FormulationTrait end

"""
    DCPowerFlow <: FormulationTrait

DC power flow formulation (linearized approximation).

# Variables
- Voltage angles only: θ - no voltage magnitudes
- All voltage magnitudes assumed at 1.0 pu

# Residuals
- Active power balance only: P = 0
- Reactive power ignored (Q equations removed)
- Branch contributes: P_i, P_j (active power flows only)
- Load contributes: P (active consumption only)
- Generator contributes: P (active injection only)

# Characteristics
- Linear equations: P = B'θ where B' is the DC admittance matrix
- Fast solution via direct linear solve
- No convergence issues (linear system always solvable if connected)
- Commonly used for contingency analysis and economic dispatch

# Assumptions
- All voltage magnitudes = 1.0 pu
- Small angle differences: sin(θ_i - θ_j) ≈ θ_i - θ_j
- Branch resistance negligible: r << x
- Reactive power negligible: P_ij ≈ (θ_i - θ_j) / x

# Benefits
- Extremely fast (linear solve, no iteration needed)
- Guaranteed convergence for connected networks
- Good approximation for transmission networks with high X/R ratio
- Industry standard for planning and market applications

# Limitations
- Ignores voltage magnitude variations
- Ignores reactive power flows
- Less accurate for distribution networks (low X/R)
- Cannot capture voltage collapse phenomena
"""
struct DCPowerFlow <: FormulationTrait end

"""
    DEFAULT_FORMULATION

The formulation used when none is specified: [`PowerInjection`](@ref)`()`.
"""
const DEFAULT_FORMULATION = PowerInjection()

export FormulationTrait, PowerInjection, CurrentInjection, DCPowerFlow, DEFAULT_FORMULATION

# ============================================================================
# Formulation Support Query
# ============================================================================

"""
    supports_formulation(::Type{M}, ::Type{F}) -> Bool where {M<:AbstractModel, F<:FormulationTrait}

Whether model type `M` is defined in formulation `F`.

By default this follows the model metadata, the definition of the model's variables
and residuals: the base metadata is the `PowerInjection` form, and another formulation
is supported when the model registers variables or residuals for it
(`vars_by_formulation` / `residuals_by_formulation`). A model without variables or
residuals (e.g. an area or zone) takes part in every formulation. Problem constructors
check every model of the system with
[`check_formulation_support`](@ref PowerfulCore.System.check_formulation_support).

# Example
```julia
supports_formulation(Branch, CurrentInjection)  # true: Branch registers (v_r, v_i) metadata
supports_formulation(GENROU, CurrentInjection)  # false: defined in the PowerInjection form only
```
"""
function supports_formulation(::Type{M}, ::Type{F}) where {M<:AbstractModel, F<:FormulationTrait}
    metadata = model_metadata(M)
    isempty(metadata.vars) && isempty(metadata.residuals) && return true
    F <: PowerInjection && return true
    return haskey(metadata.vars_by_formulation, F) || haskey(metadata.residuals_by_formulation, F)
end

export supports_formulation
