export PQ, PQVec

"""
    $TYPEDEF

Static PQ load data records.

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct PQ{Tv} <: AbstractModel

    "Bus number or extended bus name enclosed in single quotes"
    i::Int32

    """
    One- or two-character uppercase non-blank alphanumeric load identifier used to distinguish
    multiple loads at bus i. For buses with a single load, recommended to use '1'
    """
    id::String3

    "Load status: true for in-service, false for out-of-service"
    status::Bool

    "Area number assigned to this load (1 to maximum areas)"
    area::Int16

    "Zone number assigned to this load (1 to maximum zones)"
    zone::Int16

    "Active power component of constant MVA load in MW"
    pl::Tv

    "Reactive power component of constant MVA load in Mvar"
    ql::Tv

    "Active power component of constant current load in MW at 1.0 pu voltage"
    ip::Tv

    "Reactive power component of constant current load in Mvar at 1.0 pu voltage"
    iq::Tv

    "Active power component of constant admittance load in MW at 1.0 pu voltage"
    yp::Tv

    """
    Reactive power component of constant admittance load in Mvar at 1.0 pu voltage.
    Positive for inductive (consuming) load, like `ql` and `iq`. The RAW-file `YQ`
    field uses the opposite sign (negative for inductive) and is negated at parse time.
    """
    yq::Tv

    "Owner number assigned to this load (1 to maximum owners)"
    owner::Int16

    "Load scaling flag: true for scalable load, false for fixed load"
    scale::Union{Bool,Missing}

    "Interruptible load flag: true for interruptible load, false for non-interruptible"
    intrpt::Union{Bool,Missing}
end

@register_model PQ

# Type aliases for StructArray-based vectorized types
"""
    PQVec{T}

Alias for `StructArray{PQ{T}}`, the struct-of-arrays container holding all
[`PQ`](@ref) instances of a system.
"""
const PQVec{T} = StructArray{PQ{T}}

get_key_name(::Type{T}) where T<:PQ = [:i, :id]

# Base quantity traits - PQ inherits base voltage from connected bus
base_quantity_role(::Type{<:PQ}, ::Type{BaseVoltage}) = InheritsBase()

# Connectivity declaration - single source of truth for topology AND ExtParams
const PQ_TERMINALS = (Terminal(:Bus, :i),)  # Single-attachment device

get_terminals(::Type{<:PQ}) = PQ_TERMINALS

const PQ_VARS = [
    ModelVar(
        :theta,
        AlgebVar(),
        :Bus,
        :theta,
        :i;
        description="Bus voltage phase angle",
        units="rad",
    ),
    ModelVar(
        :v,
        AlgebVar(),
        :Bus,
        :v,
        :i;
        description="Voltage magnitude",
        units="pu",
    ),

]

const PQ_observable_vars = []

# CurrentInjection formulation - rectangular coordinates
const PQ_VARS_CURRENT = [
    ModelVar(
        :v_r,
        AlgebVar(),
        :Bus,
        :v_r,
        :i;
        description="Voltage real part",
        units="pu",
    ),
    ModelVar(
        :v_i,
        AlgebVar(),
        :Bus,
        :v_i,
        :i;
        description="Voltage imaginary part",
        units="pu",
    ),
]

# DC Power Flow formulation - angle only
const PQ_VARS_DC = [
    ModelVar(
        :theta,
        AlgebVar(),
        :Bus,
        :theta,
        :i;
        description="Bus voltage phase angle",
        units="rad",
    ),
]

const PQ_RESIDUALS = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :p,
        :i,
        description="Active power injection from load",
    ),
    ModelResidual(
        :q,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :q,
        :i,
        description="Reactive power injection from load",
    )
]

# CurrentInjection formulation residuals
const PQ_RESIDUALS_CURRENT = [
    ModelResidual(
        :i_r,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :i_r,
        :i,
        description="Real current injection from load",
    ),
    ModelResidual(
        :i_i,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :i_i,
        :i,
        description="Imaginary current injection from load",
    )
]

# DC Power Flow formulation residuals - active power only
const PQ_RESIDUALS_DC = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :p,
        :i,
        description="Active power injection from load (DC approximation)",
    )
]

# PQ parameters for per-unit conversion (only fields that need conversion)
const PQ_PARAMS = [
    ModelParam(:pl, PowerNominal, description="Constant power P", units="MW"),
    ModelParam(:ql, PowerNominal, description="Constant power Q", units="MVAR"),
    ModelParam(:ip, PowerNominal, description="Constant current P", units="MW"),
    ModelParam(:iq, PowerNominal, description="Constant current Q", units="MVAR"),
    ModelParam(:yp, PowerNominal, description="Constant admittance P", units="MW"),
    ModelParam(:yq, PowerNominal, description="Constant admittance Q", units="MVAR"),
]

# PQ observations (computed from ZIP model and voltage)
const PQ_OBSERVATIONS = [
    ModelObservable(:p_consumed, ActivePower, units="pu", description="Active power consumed", default_output=true),
    ModelObservable(:q_consumed, ReactivePower, units="pu", description="Reactive power consumed", default_output=true),
]

# PQ has no device ratings
const PQ_RATINGS = ModelRating[]

# PQ references bus for base voltage - uses PQ_TERMINALS
const PQ_EXT_PARAMS = [
    ExtParam(:basekv, PQ_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
    ExtParam(:ide, PQ_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type (1=PQ, 2=PV, 3=Slack)"),
    ExtParam(:vm, PQ_TERMINALS[1], :vm, Dimensionless(), PostPerUnit(),
             required=true, description="Case bus voltage magnitude (pu), the operating voltage of the DC model")
]

function model_metadata(::Type{PQ})
    metadata = ModelMetadata(
        name = :PQ,
        vars = PQ_VARS,
        observable_vars = PQ_observable_vars,
        residuals = PQ_RESIDUALS,
        params = PQ_PARAMS,
        observations = PQ_OBSERVATIONS,
        ratings = PQ_RATINGS,
        ext_params = PQ_EXT_PARAMS,
        status_field = :status,
    )

    # Register CurrentInjection formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = PQ_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = PQ_RESIDUALS_CURRENT

    # Register DC Power Flow formulation overrides
    metadata.vars_by_formulation[DCPowerFlow] = PQ_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = PQ_RESIDUALS_DC

    return metadata
end

# Declare formulation support
supports_formulation(::Type{PQ}, ::Type{PowerInjection}) = true
supports_formulation(::Type{PQ}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{PQ}, ::Type{DCPowerFlow}) = true

### === Model Format Adapters ===

supports_format(::Type{PQ}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{PQ}, ::Type{PSSE}) = (:loads,)

function parse_model(model::Type{PQ}, raw, ::FormatSupport{PSSE})
    loads = raw.loads

    # Create StructArray directly from named tuple of vectors
    return StructArray{PQ{Float64}}((
        i = loads.i,
        id = loads.id,
        status = loads.status,
        area = loads.area,
        zone = loads.zone,
        pl = loads.pl,
        ql = loads.ql,
        ip = loads.ip,
        iq = loads.iq,
        yp = loads.yp,
        yq = .-loads.yq,  # RAW YQ is negative for inductive load; store consumption-positive
        owner = loads.owner,
        scale = loads.scale,
        intrpt = loads.intrpt,
    ))
end


### === Residual Calculation === ###


# ============================================================================
# ZIP Load Physics (Reusable for Residuals and Reports)
# ============================================================================

"""
    constant_power_factor(v, vmin) -> f

Fraction of the constant-power load (`pl`, `ql`) drawn at voltage magnitude `v` (pu).
At or above the breakpoint `vmin` (the system's `constant_power_vmin`) the load is
constant power, f = 1; below it the load falls off as a raised cosine,

    f(v) = (1 - cos(π v / vmin)) / 2,

so that f(0) = 0, f(vmin) = 1 and f'(vmin) = 0: the load and its voltage derivative are
continuous at the breakpoint and the load current stays finite as v → 0. This is the
low-voltage characteristic documented for PowerWorld's minimum per-unit voltage for
constant power loads (cosine fall-off to zero at zero voltage, continuous derivative,
default 0.7 pu); PSS/E applies a smooth characteristic below PQBRAK (default 0.7 pu).
Evaluated as the equal `sin(π v / (2 vmin))²`, which has no cancellation at small `v`.
Branchless (`ifelse`) for SIMD and GPU kernels.
"""
@inline function constant_power_factor(v, vmin)
    low = sin(π * v / (2 * vmin))^2
    return ifelse(v < vmin, low, one(low))
end

"""
    constant_power_current_factor(v_mag, vmin) -> f(v_mag) / v_mag²

Factor of the constant-power load current, `I = factor · conj(S) · V`. Above the
breakpoint it is 1/|V|² (constant power); below it f/|V|² with `f` from
`constant_power_factor`, written as `(π/(2 vmin))² sinc(|V|/(2 vmin))²`, which is
exact and finite at |V| = 0 (limit π²/(4 vmin²)).
"""
@inline function constant_power_current_factor(v_mag, vmin)
    low = (π / (2 * vmin))^2 * sinc(v_mag / (2 * vmin))^2
    return ifelse(v_mag < vmin, low, oftype(low, inv(v_mag * v_mag)))
end

"""
    constant_power_current_factor_derivative(v_mag, vmin) -> dg/dv_mag

Derivative of `g = constant_power_current_factor(v_mag, vmin)` with respect to the
voltage magnitude: -2/|V|³ at or above `vmin`; below it
`(π/(2 vmin))² · sinc(x) · cosc(x) / vmin` with `x = |V|/(2 vmin)` (`cosc` is the
derivative of `sinc`), which is finite at |V| = 0.
"""
@inline function constant_power_current_factor_derivative(v_mag, vmin)
    x = v_mag / (2 * vmin)
    low = (π / (2 * vmin))^2 * sinc(x) * cosc(x) / vmin
    return ifelse(v_mag < vmin, low, oftype(low, -2 / (v_mag * v_mag * v_mag)))
end

export constant_power_factor, constant_power_current_factor, constant_power_current_factor_derivative

"""
    calc_zip_load(v, pl, ql, ip, iq, yp, yq, vmin) -> (p=..., q=...)

Calculate ZIP load model power consumption.

# ZIP Load Model
- P = pl·f(v) + ip*v + yp*v²  (constant power + current + impedance)
- Q = ql·f(v) + iq*v + yq*v²

where `f = constant_power_factor(v, vmin)` is 1 at or above the constant-power
breakpoint `vmin` and falls off smoothly to zero below it.

# Arguments
- `v`: Voltage magnitude (pu)
- `pl, ql`: Constant power components (MW, MVAR)
- `ip, iq`: Constant current components (MW@1pu, MVAR@1pu)
- `yp, yq`: Constant impedance components (MW@1pu, MVAR@1pu)
- `vmin`: Constant-power breakpoint (pu), `sys.settings.constant_power_vmin`

# Returns
Named tuple `(p=..., q=...)` with power consumption in pu.

# Note
This is a pure physics function compatible with CPU, GPU, and reporting.
"""
@inline function calc_zip_load(v, pl, ql, ip, iq, yp, yq, vmin)
    v2 = v * v
    f = constant_power_factor(v, vmin)
    p = pl * f + ip * v + yp * v2
    q = ql * f + iq * v + yq * v2
    return (p=p, q=q)
end

"""
    calc_load_current(v_r, v_i, pl, ql, ip, iq, yp, yq, vmin) -> (i_r=..., i_i=...)

Calculate ZIP load current injection from rectangular voltage.

Converts ZIP load model from power form to current form:
- **Constant Power**: I_p = S*/V* = (P - jQ) / (Vr - jVi)
- **Constant Current**: I_c = (Ip - jIq) * (V/|V|)
- **Constant Admittance**: I_y = Y*V = (Yp - jYq) * (Vr + jVi)

# Arguments
- `v_r, v_i`: Real and imaginary voltage [pu]
- `pl, ql`: Constant power components [pu]
- `ip, iq`: Constant current components [pu @ 1pu]
- `yp, yq`: Constant admittance components [pu @ 1pu]
- `vmin`: Constant-power breakpoint (pu), `sys.settings.constant_power_vmin`

# Returns
Named tuple `(i_r=..., i_i=...)` with current injection in pu.

# Error Handling
Throws error if |V| < 1e-12 (voltage collapse condition).

# Sign Convention
Positive current = consumption (current flowing INTO the load, OUT of the bus).
This matches S = V * conj(I) with load convention.

# Physics Details

## Constant Power Component
S_p = f(|V|)·(P + jQ), with `f = constant_power_factor(|V|, vmin)`
I_p = S_p*/V* = f(|V|)·(P - jQ)(Vr + jVi) / |V|²
The factor f/|V|² comes from `constant_power_current_factor`.

## Constant Current Component
At 1 pu voltage, load draws current with magnitude based on ip, iq.
At actual voltage V, current magnitude scales with 1, but direction follows V:
I_c = (ip - j*iq) * (V / |V|) = (ip - j*iq) * (Vr + j*Vi) / |V|

## Constant Admittance Component
Y = yp - j*yq (load admittance)
I_y = Y * V = (yp - j*yq) * (Vr + j*Vi)
"""
@inline function calc_load_current(v_r, v_i, pl, ql, ip, iq, yp, yq, vmin)
    # Voltage magnitude squared
    v_mag_sq = v_r * v_r + v_i * v_i

    # Error on voltage collapse (the constant-current direction V/|V| is undefined at 0)
    v_mag_sq < 1e-12 && error("Voltage collapse detected: |V| ≈ 0 at load")

    v_mag = sqrt(v_mag_sq)
    inv_v = one(v_r) / v_mag

    # Constant power component: I_p = f(|V|)(P - jQ)(Vr + jVi) / |V|²
    # (P - jQ)(Vr + jVi) = P*Vr + Q*Vi + j(P*Vi - Q*Vr)
    cp = constant_power_current_factor(v_mag, vmin)
    i_p_r = (pl * v_r + ql * v_i) * cp
    i_p_i = (pl * v_i - ql * v_r) * cp

    # Constant current component: I_c = (ip - j*iq) * (Vr + j*Vi) / |V|
    # (ip - j*iq)(Vr + j*Vi) = ip*Vr + iq*Vi + j(ip*Vi - iq*Vr)
    i_c_r = (ip * v_r + iq * v_i) * inv_v
    i_c_i = (ip * v_i - iq * v_r) * inv_v

    # Constant admittance: I_y = (yp - j*yq) * (Vr + j*Vi)
    # (yp - j*yq)(Vr + j*Vi) = yp*Vr + yq*Vi + j(yp*Vi - yq*Vr)
    i_y_r = yp * v_r + yq * v_i
    i_y_i = yp * v_i - yq * v_r

    # Total current (all components add)
    i_r = i_p_r + i_c_r + i_y_r
    i_i = i_p_i + i_c_i + i_y_i

    return (i_r=i_r, i_i=i_i)
end

"""
    calc_load_current_jacobian(v_r, v_i, pl, ql, ip, iq, yp, yq, vmin)
        -> (rr=∂i_r/∂v_r, ri=∂i_r/∂v_i, ir=∂i_i/∂v_r, ii=∂i_i/∂v_i)

Exact derivative of [`calc_load_current`](@ref) with respect to the rectangular
voltage (derived with SymPy). With m = |V|, g = `constant_power_current_factor(m, vmin)`,
g' = `constant_power_current_factor_derivative(m, vmin)`,
a_r + j a_i = (pl - j ql)(v_r + j v_i) and c_r + j c_i = (ip - j iq)(v_r + j v_i):

    ∂i_r/∂v_r =  pl g + g' v_r a_r / m + ip / m - v_r c_r / m³ + yp
    ∂i_r/∂v_i =  ql g + g' v_i a_r / m + iq / m - v_i c_r / m³ + yq
    ∂i_i/∂v_r = -ql g + g' v_r a_i / m - iq / m - v_r c_i / m³ - yq
    ∂i_i/∂v_i =  pl g + g' v_i a_i / m + ip / m - v_i c_i / m³ + yp

The current is not holomorphic in V (it depends on |V|), so the derivative is a real
2×2 block rather than a complex scalar. Same sign convention as `calc_load_current`
(consumption). Errors if |V| < 1e-6, like `calc_load_current`.
"""
@inline function calc_load_current_jacobian(v_r, v_i, pl, ql, ip, iq, yp, yq, vmin)
    v_mag_sq = v_r * v_r + v_i * v_i
    v_mag_sq < 1e-12 && error("Voltage collapse detected: |V| ≈ 0 at load")

    v_mag = sqrt(v_mag_sq)
    inv_v = one(v_r) / v_mag
    inv_v3 = inv_v * inv_v * inv_v

    g = constant_power_current_factor(v_mag, vmin)
    dg_m = constant_power_current_factor_derivative(v_mag, vmin) * inv_v  # g'/m

    a_r = pl * v_r + ql * v_i
    a_i = pl * v_i - ql * v_r
    c_r = ip * v_r + iq * v_i
    c_i = ip * v_i - iq * v_r

    rr = pl * g + dg_m * v_r * a_r + ip * inv_v - v_r * c_r * inv_v3 + yp
    ri = ql * g + dg_m * v_i * a_r + iq * inv_v - v_i * c_r * inv_v3 + yq
    ir = -ql * g + dg_m * v_r * a_i - iq * inv_v - v_r * c_i * inv_v3 - yq
    ii = pl * g + dg_m * v_i * a_i + ip * inv_v - v_i * c_i * inv_v3 + yp
    return (rr=rr, ri=ri, ir=ir, ii=ii)
end

# ============================================================================
# Reporting Interface - Vectorized Observation Computation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:PQ}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all PQ loads in the solution (vectorized).

# Supported Variables
- `:p_consumed`: Total active power consumed (pu) accounting for ZIP load model
- `:q_consumed`: Total reactive power consumed (pu) accounting for ZIP load model

# Implementation Notes
Reuses `calc_zip_load` physics function to calculate total consumption, with the
system's constant-power breakpoint:
- P_total = pl·f(|V|) + ip*|V| + yp*|V|²
- Q_total = ql·f(|V|) + iq*|V| + yq*|V|²

Voltages are extracted from solution via AddressManager (same pattern as old scalar API).
"""
function compute_observations(sol::PowerFlowSolution, ::Type{<:PQ}, var::Symbol)
    pq = sol.prob.sys.models.PQ
    vmin = sol.prob.sys.settings.constant_power_vmin
    var in (:p_consumed, :q_consumed) || error("Unknown PQ output variable: $var")

    # Terminal bus voltages (any formulation); out-of-service loads consume nothing
    v_vals, _ = terminal_voltages(sol, :PQ, :i)
    n = length(pq)
    result = Vector{Float64}(undef, n)

    if sol.prob.formulation isa DCPowerFlow
        var == :p_consumed || error("PQ observation :$var is not defined in the DCPowerFlow formulation")
        @inbounds for i in 1:n
            result[i] = pq.status[i] * dc_load_power(pq.pl[i], pq.ip[i], pq.yp[i], v_vals[i], vmin)
        end
        return result
    end

    @inbounds for i in 1:n
        power = calc_zip_load(v_vals[i], pq.pl[i], pq.ql[i],
                              pq.ip[i], pq.iq[i], pq.yp[i], pq.yq[i], vmin)
        result[i] = pq.status[i] * (var == :p_consumed ? power.p : power.q)
    end

    return result
end

# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    pq_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for PQ load residuals using ZIP load model.

This is a **pure math kernel** with no addressing logic, optimized for:
- CPU vectorization (`@simd`)
- GPU execution (no branching in hot loop)
- Easy testing and verification

Dispatches to CPU or GPU implementation based on array type parameter.

# ZIP Load Model
- P_load = pl·f(v) + ip*v + yp*v²
- Q_load = ql·f(v) + iq*v + yq*v²

with `f = constant_power_factor(v, params.vmin)` (constant power above the breakpoint).

# Masking
- Inactive loads (status==false): contribute 0
- Bus type filtering handled by scatter SharedRes mode (:add vs :set)

# Arguments
- `outputs::NamedTuple{(:p, :q)}`: Preallocated output buffers
- `params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)}`: Load parameters and the
  system's constant-power breakpoint (scalar)
- `vars::NamedTuple{(:v,)}`: Voltage magnitude from workspace
- `active::AbstractVector{Bool}`: Mask for in-service loads
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)

# Performance Notes
- Branchless computation using mask multiplication
- `@simd` friendly (no conditionals in loop)
- Zero allocation
- Automatic CPU/GPU dispatch

# Example
```julia
# CPU dispatch
pq_residual_kernel!(
    (p = p_res, q = q_res),
    (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq, vmin = 0.7),
    (v = voltage_vector,),
    active_mask,
    Vector{Float64}
)

# GPU dispatch (automatically uses CUDA kernel)
pq_residual_kernel!(
    (p = p_res_gpu, q = q_res_gpu),
    (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq, vmin = 0.7),
    (v = voltage_vector_gpu,),
    active_mask_gpu,
    CuVector{Float64}
)
```
"""
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    # Pure vectorized computation - no addressing logic
    @inbounds @simd for i in eachindex(outputs.p)
        # Branchless masking - only active status
        # Bus type filtering handled by scatter mode (:add vs :set)
        mask = active[i]

        # Calculate load using extracted physics function
        power = calc_zip_load(vars.v[i], params.pl[i], params.ql[i],
                              params.ip[i], params.iq[i], params.yp[i], params.yq[i], params.vmin)

        # Sign convention: load consumption = negative injection into bus
        outputs.p[i] = mask * (-power.p)
        outputs.q[i] = mask * (-power.q)
    end

    nothing
end


# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, pq::StructArray{<:PQ}, external::NamedTuple) where {T, AT}

Calculate PQ load residual contributions using ZIP load model.

Implements the composite load model with constant MVA (P), constant current (I),
and constant impedance (Z) components:
- P_load = pl·f(v) + ip*v + yp*v²
- Q_load = ql·f(v) + iq*v + yq*v²

where `f = constant_power_factor(v, ws.system.constant_power_vmin)`.

Residual convention: positive values represent load consumption contributions.
Only active (in-service) loads contribute to residuals.
Bus type filtering handled by scatter SharedRes mode (:add vs :set).

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `pq`: PQ load model vector with numerical fields
- `_external`: External parameters (unused after refactoring)
- `active_mask`: Boolean vector indicating which instances are in-service

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `active_mask` for out-of-service loads
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: `ws.var_vals.v` (voltage magnitudes)
- Writes: `ws.res_vals.p`, `ws.res_vals.q` (power residuals)
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pq::StructArray{<:PQ},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract StructArray fields to NamedTuples
    # AT parameter dispatches to CPU (Vector) or GPU (CuArray) implementation
    pq_residual_kernel!(
        (p = ws.res_vals.p, q = ws.res_vals.q),              # outputs: NamedTuple of arrays
        (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq,
         vmin = ws.system.constant_power_vmin),              # params: StructArray fields + breakpoint
        (v = ws.var_vals.v,),                                # vars: NamedTuple of arrays
        active_mask,                                          # mask: array
        AT                                                     # dispatch: type parameter
    )
end

# CurrentInjection formulation - uses rectangular voltages
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pq::StructArray{<:PQ},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and call current injection kernel
    pq_current_kernel!(
        (i_r = ws.res_vals.i_r, i_i = ws.res_vals.i_i),     # outputs: current residuals
        (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq,
         vmin = ws.system.constant_power_vmin),              # params: ZIP load parameters + breakpoint
        (v_r = ws.var_vals.v_r, v_i = ws.var_vals.v_i),     # vars: rectangular voltages
        active_mask,                                          # mask: array
        AT                                                     # dispatch: type parameter
    )
end

# ============================================================================
# Current Injection Kernel (CPU)
# ============================================================================

"""
    pq_current_kernel!(outputs, params, vars, active, ::Type{<:Vector})

Compute PQ load current injections for CurrentInjection formulation.

Uses rectangular voltage coordinates (V_r, V_i) and computes current consumption
based on ZIP load model. Sign convention: Load WITHDRAWS current from bus (negative injection).
"""
@inline function pq_current_kernel!(
    outputs::NamedTuple{(:i_r, :i_i)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)},
    vars::NamedTuple{(:v_r, :v_i)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.i_r)
        mask = active[k]

        # Get current from helper function
        current = calc_load_current(
            vars.v_r[k], vars.v_i[k],
            params.pl[k], params.ql[k],
            params.ip[k], params.iq[k],
            params.yp[k], params.yq[k], params.vmin
        )

        # Sign convention: load consumption = negative current injection into bus
        outputs.i_r[k] = mask * (-current.i_r)
        outputs.i_i[k] = mask * (-current.i_i)
    end

    nothing
end

# ============================================================================
# DC Power Flow Kernel (CPU)
# ============================================================================

"""
    pq_dc_kernel!(outputs, params, active, ::Type{<:Vector})

Compute PQ load active power consumption for DC Power Flow formulation.

In DC power flow the load draws its power at the case bus voltage `vm`:
`pl·f(vm) + ip·vm + yp·vm²` (see `dc_load_power`).

Sign convention: Load CONSUMES power (negative contribution to bus balance).
"""
@inline function pq_dc_kernel!(
    outputs::NamedTuple{(:p,)},
    params::NamedTuple{(:pl, :ip, :yp, :vm, :vmin)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.p)
        mask = active[k]

        # Sign convention: load consumption = negative injection into bus
        outputs.p[k] = mask * (-dc_load_power(params.pl[k], params.ip[k], params.yp[k], params.vm[k], params.vmin))
    end

    nothing
end

# DC Power Flow calc_residuals! dispatcher
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pq::StructArray{<:PQ},
    external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and call DC kernel
    pq_dc_kernel!(
        (p = ws.res_vals.p,),           # outputs: active power only
        (pl = pq.pl, ip = pq.ip, yp = pq.yp, vm = external.vm,
         vmin = ws.system.constant_power_vmin),  # ZIP active power at case voltage
        active_mask,                     # mask: array
        AT                               # dispatch: type parameter
    )
end

