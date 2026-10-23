export Fault, FaultVec

"""
    $TYPEDEF

Three-phase fault model for time-domain simulation.

Represents a 3-phase-to-ground fault as a shunt admittance to ground.
The fault absorbs reactive power (large negative susceptance) to cause
voltage depression at the faulted bus.

# Physics
A 3-phase fault with impedance Zf = Rf + jXf has admittance:
- Yf = 1/Zf = Gf + jBf
- Gf = Rf / (Rf² + Xf²)
- Bf = -Xf / (Rf² + Xf²)

Power consumed by fault (ANDES convention, positive = consumption):
- P = Gf x V²
- Q = -Bf x V² (negative Bf → positive Q → absorbs reactive power)

For a **bolted fault** (Zf → 0), the fault path is predominantly inductive,
so we use Gf ≈ 0 and Bf << 0 (large negative susceptance).

# Fields

$TYPEDFIELDS

# Usage
```julia
# Create fault programmatically (not from RAW file)
fault = Fault(i=10, id="F1", rf=0.0, xf=0.0001)  # Near-bolted fault at bus 10

# Add to system
push!(sys.models.Fault, fault)

# Apply/clear via events
events = [
    FieldModification(1.0, :Fault, (Int32(10), String3("F1")), :status, :set, 1.0),
    FieldModification(1.1, :Fault, (Int32(10), String3("F1")), :status, :set, 0.0),
]
```
"""
@kwdef mutable struct Fault{Tv} <: AbstractModel
    "Bus number where fault occurs"
    i::Int32

    "Fault identifier (e.g., \"F1\", \"1\")"
    id::String3

    "Fault resistance in per-unit on system base"
    rf::Tv

    "Fault reactance in per-unit on system base (small positive → large negative Bf)"
    xf::Tv

    "Fault status: false = cleared (off), true = applied (on)"
    status::Bool
end

@register_model Fault

# Type alias for StructArray-based vectorized types
"""
    FaultVec{T}

Alias for `StructArray{Fault{T}}`, the struct-of-arrays container holding all
[`Fault`](@ref) instances of a system.
"""
const FaultVec{T} = StructArray{Fault{T}}

get_key_name(::Type{T}) where T<:Fault = [:i, :id]

# Base quantity traits - Fault inherits base voltage from connected bus
base_quantity_role(::Type{<:Fault}, ::Type{BaseVoltage}) = InheritsBase()

# Connectivity declaration - connects via bus_group (works with any bus model, e.g. Bus)
const FAULT_TERMINALS = (Terminal(:bus_group, :i),)

get_terminals(::Type{<:Fault}) = FAULT_TERMINALS

# Variables: voltage from bus_group (e.g. Bus)
const FAULT_VARS = [
    ModelVar(
        :theta,
        AlgebVar(),
        :bus_group,
        :theta,
        :i;
        description="Bus voltage phase angle",
        units="rad",
    ),
    ModelVar(
        :v,
        AlgebVar(),
        :bus_group,
        :v,
        :i;
        description="Bus voltage magnitude",
        units="pu",
    ),
]

const FAULT_VARS_CURRENT = [
    ModelVar(
        :v_r,
        AlgebVar(),
        :bus_group,
        :v_r,
        :i;
        description="Bus voltage real component",
        units="pu",
    ),
    ModelVar(
        :v_i,
        AlgebVar(),
        :bus_group,
        :v_i,
        :i;
        description="Bus voltage imaginary component",
        units="pu",
    ),
]

# Residuals: power injection to bus_group (e.g. Bus)
const FAULT_RESIDUALS = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :bus_group,
        :p,
        :i,
        description="Active power consumed by fault",
    ),
    ModelResidual(
        :q,
        AlgebRes(),
        SharedRes(),
        :bus_group,
        :q,
        :i,
        description="Reactive power consumed by fault",
    )
]

const FAULT_RESIDUALS_CURRENT = [
    ModelResidual(
        :i_r,
        AlgebRes(),
        SharedRes(:add),
        :bus_group,
        :i_r,
        :i,
        description="Real current injection from fault",
    ),
    ModelResidual(
        :i_i,
        AlgebRes(),
        SharedRes(:add),
        :bus_group,
        :i_i,
        :i,
        description="Imaginary current injection from fault",
    )
]

# No per-unit conversion needed - fault parameters are already in per-unit
const FAULT_PARAMS = ModelParam[]

# Observations
const FAULT_OBSERVATIONS = [
    ModelObservable(:p_consumed, ActivePower, units="pu", description="Active power consumed by fault", default_output=true),
    ModelObservable(:q_consumed, ReactivePower, units="pu", description="Reactive power consumed by fault", default_output=true),
]

# No device ratings
const FAULT_RATINGS = ModelRating[]

# External parameters from connected bus
const FAULT_EXT_PARAMS = [
    ExtParam(:basekv, FAULT_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
    ExtParam(:ide, FAULT_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type (1=PQ, 2=PV, 3=Slack)")
]

function model_metadata(::Type{Fault})
    metadata = ModelMetadata(
        name = :Fault,
        vars = FAULT_VARS,
        observable_vars = [],
        residuals = FAULT_RESIDUALS,
        params = FAULT_PARAMS,
        observations = FAULT_OBSERVATIONS,
        ratings = FAULT_RATINGS,
        ext_params = FAULT_EXT_PARAMS,
        status_field = :status,
    )

    # Register CurrentInjection formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = FAULT_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = FAULT_RESIDUALS_CURRENT

    return metadata
end

# Declare formulation support
# No state of its own: its current follows from the bus voltage it reads
requires_initialization(::Type{<:Fault}) = false
supports_formulation(::Type{Fault}, ::Type{PowerInjection}) = true
supports_formulation(::Type{Fault}, ::Type{CurrentInjection}) = true

# ============================================================================
# Fault Constructors - User-Friendly API
# ============================================================================

"""
    create_fault(bus::Integer; id="F1", rf=0.0, xf=1e-4, status=false) -> Fault

Create a Fault instance with computed admittance values.

# Arguments
- `bus`: Bus number where fault occurs
- `id`: Fault identifier (default: "F1")
- `rf`: Fault resistance in per-unit (default: 0.0)
- `xf`: Fault reactance in per-unit (default: 1e-4 for near-bolted fault)
- `status`: Initial status, false=cleared, true=applied (default: false)

# Physics
For near-bolted fault (rf≈0, xf→0⁺):
- Gf ≈ 0
- Bf = -1/xf (large negative, absorbs reactive power)

# Example
```julia
# Bolted fault at bus 10, initially cleared
fault = create_fault(10)

# Fault with impedance at bus 5
fault = create_fault(5; rf=0.01, xf=0.05, id="F2")
```
"""
function create_fault(bus::Integer; id::AbstractString="F1", rf::Real=0.0, xf::Real=1e-4, status::Bool=false)
    T = Float64
    rf >= 0 && xf >= 0 || error("Fault impedance must be non-negative: rf=$rf, xf=$xf")

    return Fault{T}(
        i = Int32(bus),
        id = String3(id),
        rf = T(rf),
        xf = T(xf),
        status = status,
    )
end

"""
    fault_admittance(rf, xf) -> (gf, bf)

Fault admittance `gf + j·bf = 1 / (rf + j·xf)`. A bolted fault (`rf² + xf² < 1e-12`)
is represented by `gf = 0`, `bf = -1e4`.

The residual kernels evaluate it from `rf` and `xf` at every call, so an event that
changes the fault impedance takes effect. Branch-free, for SIMD kernels.
"""
@inline function fault_admittance(rf, xf)
    z_sq = rf^2 + xf^2
    bolted = z_sq < oftype(z_sq, 1e-12)
    inv_z_sq = inv(ifelse(bolted, one(z_sq), z_sq))
    gf = ifelse(bolted, zero(z_sq), rf * inv_z_sq)
    bf = ifelse(bolted, oftype(z_sq, -1e4), -xf * inv_z_sq)  # negative for inductive xf > 0
    return (gf, bf)
end

# ============================================================================
# Fault Physics (Reusable for Residuals and Reports)
# ============================================================================

"""
    calc_fault_power(v, gf, bf) -> (p=..., q=...)

Calculate the power consumed by an applied fault. Whether the fault is applied is the
in-service mask of the caller (built from the `status` field and updated by events).

# Physics
- P = gf × v²  (resistive losses, typically ≈ 0)
- Q = (-bf) × v²  (reactive power absorbed)

For bolted fault (bf << 0): Q >> 0 (absorbs large reactive power, causes voltage drop)

# Arguments
- `v`: Voltage magnitude (pu)
- `gf`: Fault conductance (pu)
- `bf`: Fault susceptance (pu, negative for inductive)

# Returns
Named tuple `(p=..., q=...)` with power in pu (ANDES convention: positive = consumption).
"""
@inline function calc_fault_power(v, gf, bf)
    v2 = v * v
    p = gf * v2
    q = (-bf) * v2  # -bf because bf < 0 for inductive fault
    return (p=p, q=q)
end

"""
    calc_fault_current(v_r, v_i, gf, bf) -> (i_r=..., i_i=...)

Calculate the current drawn by an applied fault to ground. Whether the fault is applied
is the in-service mask of the caller.

# Current Injection Model
I = Yf × V where Yf = gf + j×bf

- I_r = gf×v_r - bf×v_i
- I_i = gf×v_i + bf×v_r

# Arguments
- `v_r`: Voltage real component (pu)
- `v_i`: Voltage imaginary component (pu)
- `gf`: Fault conductance (pu)
- `bf`: Fault susceptance (pu)

# Returns
Named tuple `(i_r=..., i_i=...)` with current in pu.
"""
# A fault is a switched shunt admittance Yf = gf + j*bf
@inline calc_fault_current(v_r, v_i, gf, bf) = PowerfulPF.calc_shunt_current(v_r, v_i, gf, bf)

# ============================================================================
# Reporting Interface
# ============================================================================

function compute_observations(sol::PowerFlowSolution, ::Type{<:Fault}, var::Symbol)
    fault = sol.prob.sys.models.Fault

    # Terminal bus voltages (any formulation)
    v_vals, _ = PowerfulPF.terminal_voltages(sol, :Fault, :i)

    n = length(fault)
    result = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        gf, bf = fault_admittance(fault.rf[i], fault.xf[i])
        power = calc_fault_power(v_vals[i], gf, bf)

        # A cleared fault consumes nothing
        if var == :p_consumed
            result[i] = fault.status[i] * power.p
        elseif var == :q_consumed
            result[i] = fault.status[i] * power.q
        else
            error("Unknown Fault output variable: $var")
        end
    end

    return result
end

# ============================================================================
# Pure Computation Kernels
# ============================================================================

"""
    fault_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for fault residuals (PowerInjection formulation).

Fault consumes power proportional to voltage squared, with the admittance
`gf + j·bf` from [`fault_admittance`](@ref):
- P = gf × v²
- Q = (-bf) × v²

The in-service mask `active` (the fault status, updated by events) switches it.

For bolted fault (bf << 0): absorbs large reactive power, causes voltage drop.
"""
@inline function fault_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:rf, :xf)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        mask = active[i]

        # Calculate power using physics function
        gf, bf = fault_admittance(params.rf[i], params.xf[i])
        power = calc_fault_power(vars.v[i], gf, bf)

        # Sign convention: positive = injection into bus
        # Fault consumes power, so injection is negative (same as Shunt)
        outputs.p[i] = mask * (-power.p)
        outputs.q[i] = mask * (-power.q)
    end

    nothing
end

"""
    fault_current_kernel!(outputs, params, external, vars, active, ::Type{AT})

Pure computation kernel for fault current residuals (CurrentInjection formulation).

Fault draws current I = Yf × V while the in-service mask `active` (the fault status,
updated by events) is set.
"""
@inline function fault_current_kernel!(
    outputs::NamedTuple{(:i_r, :i_i)},
    params::NamedTuple{(:rf, :xf)},
    external::NamedTuple,
    vars::NamedTuple{(:v_r, :v_i)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.i_r)
        # Mask out slack buses
        mask = active[i] & (external.ide[i] != BUS_SLACK)

        # Calculate current using physics function
        gf, bf = fault_admittance(params.rf[i], params.xf[i])
        current = calc_fault_current(vars.v_r[i], vars.v_i[i], gf, bf)

        # Sign convention: Y*V gives current INTO fault (to ground), negate for current INTO bus
        outputs.i_r[i] = mask * (-current.i_r)
        outputs.i_i[i] = mask * (-current.i_i)
    end

    nothing
end

# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================


@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    fault::StructArray{<:Fault},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    fault_residual_kernel!(
        (p = ws.res_vals.p, q = ws.res_vals.q),
        (rf = fault.rf, xf = fault.xf),
        (v = ws.var_vals.v,),
        active_mask,
        AT
    )
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    fault::StructArray{<:Fault},
    external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    fault_current_kernel!(
        (i_r = ws.res_vals.i_r, i_i = ws.res_vals.i_i),
        (rf = fault.rf, xf = fault.xf),
        external,
        (v_r = ws.var_vals.v_r, v_i = ws.var_vals.v_i),
        active_mask,
        AT
    )
end

"""
    add_fault!(sys::SystemModel, fault::Fault) -> SystemModel

Add a Fault instance to `sys` (see `add_instances!`). Errors if a fault with the
same (bus, id) exists or the bus is not in the system.
"""
add_fault!(sys::SystemModel, fault::Fault) = add_instances!(sys, (Fault = StructArray([fault]),))

