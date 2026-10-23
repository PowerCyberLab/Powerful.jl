export PV, PVVec

"""
    $TYPEDEF

PV generator data records.

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct PV{Tv} <: AbstractModel
    "Bus number or extended bus name enclosed in single quotes"
    i::Int32

    """
    One or two character uppercase non-blank alphanumeric machine identifier.
    Used to distinguish multiple machines at bus i.
    Recommended to use '1' for buses with single machine.
    """
    id::String3

    "Generator active power output in MW. Default: 0.0"
    pg::Tv

    """
    Generator reactive power output in Mvar.
    Only needed if case is to be treated as solved. Default: 0.0
    """
    qg::Tv

    """
    Maximum generator reactive power output in Mvar.
    Must equal fixed Mvar output for non-regulating generators. Default: 9999.0
    """
    qt::Tv

    """
    Minimum generator reactive power output in Mvar.
    Must equal fixed Mvar output for non-regulating generators. Default: -9999.0
    """
    qb::Tv

    "Regulated voltage setpoint in per unit. Default: 1.0"
    vs::Tv

    """
    Remote bus number to regulate voltage.
    Must be type 1 or 2 bus. Zero means generator regulates own voltage.
    Must be zero for type 3 (swing) bus. Default: 0
    """
    ireg::Int32

    """
    Total MVA base of units represented by this machine.
    Required for switching studies, fault analysis and dynamic simulation.
    Default: system base MVA
    """
    mbase::Tv

    """
    Real part of machine impedance ZSORCE in pu on MBASE.
    Used for switching studies, fault analysis and dynamic simulation.
    For dynamics: set to unsaturated subtransient/transient impedance.
    For short-circuit: use saturated values. Default: 0.0
    """
    zr::Tv

    "Imaginary part of machine impedance ZSORCE. See zr. Default: 1.0"
    zx::Tv

    """
    Real part of step-up transformer impedance in pu on MBASE.
    Set to zero if transformer explicitly modeled. Default: 0.0
    """
    rt::Tv

    "Imaginary part of step-up transformer impedance. See rt. Default: 0.0"
    xt::Tv

    """
    Step-up transformer off-nominal turns ratio in pu.
    Only used if transformer impedance nonzero. Default: 1.0
    """
    gtap::Tv

    """
    Machine status: true = in service. Default: true
    """
    status::Bool

    """
    Percent of total Mvar for voltage control at regulated bus.
    Required if multiple devices control voltage at bus IREG.
    Must be positive. Default: 100.0
    """
    rmpct::Tv

    "Maximum generator active power output in MW. Default: 9999.0"
    pt::Tv

    "Minimum generator active power output in MW. Default: -9999.0"
    pb::Tv

    "Primary owner number (1 to max owners). Default: owner of bus i"
    o1::Int16

    "Fraction of ownership for owner 1 (>0). Normalized to sum to 1.0. Default: 1.0"
    f1::Union{Missing,Tv}

    "Second owner number. Default: 0"
    o2::Union{Missing,Int16}

    "Ownership fraction for owner 2. See f1"
    f2::Union{Missing,Tv}

    "Third owner number. Default: 0"
    o3::Union{Missing,Int16}

    "Ownership fraction for owner 3. See f1"
    f3::Union{Missing,Tv}

    "Fourth owner number. Default: 0"
    o4::Union{Missing,Int16}

    "Ownership fraction for owner 4. See f1"
    f4::Union{Missing,Tv}

    """
    Wind machine control mode:
    0 = not wind machine
    1 = wind machine with QT/QB limits
    2 = wind machine with limits from output and wpf
    3 = fixed reactive power from output and wpf
    Default: 0
    """
    wmod::Union{Missing,Int8}

    """
    Power factor for reactive power calculation when wmod is 2 or 3.
    Sign determines reactive power direction. Default: 1.0
    """
    wpf::Union{Missing,Tv}
end

@register_model PV

# Type aliases for StructArray-based vectorized types
"""
    PVVec{T}

Alias for `StructArray{PV{T}}`, the struct-of-arrays container holding all
[`PV`](@ref) instances of a system.
"""
const PVVec{T} = StructArray{PV{T}}

get_key_name(::Type{T}) where T<:PV = [:i, :id]

# Base quantity traits - PV inherits base voltage from connected bus
base_quantity_role(::Type{<:PV}, ::Type{BaseVoltage}) = InheritsBase()

# Connectivity declaration - single source of truth for topology AND ExtParams
const PV_TERMINALS = (Terminal(:Bus, :i),)  # Single-attachment device

get_terminals(::Type{<:PV}) = PV_TERMINALS

const PV_VARS = [
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
        description="Bus voltage magnitude",
        units="pu",
    )
]

const PV_observable_vars = []

# CurrentInjection formulation - rectangular coordinates
const PV_VARS_CURRENT = [
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
    )
]

# DC Power Flow formulation - angle only
const PV_VARS_DC = [
    ModelVar(
        :theta,
        AlgebVar(),
        :Bus,
        :theta,
        :i;
        description="Bus voltage phase angle",
        units="rad",
    )
]

const PV_RESIDUALS = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :p,
        :i,
        description="Active power injection from generator",
    ),
    ModelResidual(
        :v_ctrl,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: overwrites Bus Q with voltage control equation
        :Bus,
        :q,  # PV voltage control REPLACES the Q equation at PV buses
        :i,
        description="Voltage magnitude control equation (replaces Q balance)",
    )
]

# CurrentInjection formulation residuals
const PV_RESIDUALS_CURRENT = [
    ModelResidual(
        :i_r,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: overwrites Bus I_r equation
        :Bus,
        :i_r,
        :i,
        description="Real current/power balance or V_r setpoint (for slack)",
    ),
    ModelResidual(
        :i_i_or_vmag,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: overwrites Bus I_i equation
        :Bus,
        :i_i,
        :i,
        description="Voltage magnitude control or V_i setpoint",
    )
]

# DC Power Flow formulation residuals - active power only
const PV_RESIDUALS_DC = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :p,
        :i,
        description="Active power injection from generator (DC approximation)",
    )
]

# PV parameters for per-unit conversion (only fields that need conversion)
# NOTE: PSS/E stores PV.pg in MW (nominal), not pu on mbase
const PV_PARAMS = [
    ModelParam(:pg, PowerNominal, description="Active power output", units="MW"),
    ModelParam(:qg, PowerNominal, description="Reactive power output", units="MVAR"),
]

# PV observations (outputs from solved power flow)
const PV_OBSERVATIONS = [
    ModelObservable(:p_gen, ActivePower, units="pu", description="Active power generated", default_output=true),
    ModelObservable(:q_gen, ReactivePower, units="pu", description="Reactive power generated", default_output=true),
]

# PV has device power rating
const PV_RATINGS = [
    ModelRating(:mbase, NominalPower, description="Generator MVA base")
]

# PV references bus for base voltage - uses PV_TERMINALS
const PV_EXT_PARAMS = [
    ExtParam(:basekv, PV_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage")
]

function model_metadata(::Type{PV})
    metadata = ModelMetadata(
        name = :PV,
        vars = PV_VARS,
        observable_vars = PV_observable_vars,
        residuals = PV_RESIDUALS,
        params = PV_PARAMS,
        observations = PV_OBSERVATIONS,
        ratings = PV_RATINGS,
        ext_params = PV_EXT_PARAMS,
        status_field = :status,
        # PV uses evaluation_pass=1 to run AFTER branches/transformers in CurrentInjection
        # (reads accumulated currents and replaces with power balance equations)
        evaluation_pass = 1,
    )

    # Register CurrentInjection formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = PV_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = PV_RESIDUALS_CURRENT

    # Register DC Power Flow formulation overrides
    metadata.vars_by_formulation[DCPowerFlow] = PV_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = PV_RESIDUALS_DC

    return metadata
end

# Declare formulation support
supports_formulation(::Type{PV}, ::Type{PowerInjection}) = true
supports_formulation(::Type{PV}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{PV}, ::Type{DCPowerFlow}) = true

### === Model Format Adapters ===

supports_format(::Type{PV}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{PV}, ::Type{PSSE}) = (:generators,)

function parse_model(::Type{PV}, raw, ::FormatSupport{PSSE})
    generators = raw.generators

    # Create StructArray directly from named tuple of vectors
    return StructArray{PV{Float64}}((
        i = generators.i,
        id = generators.id,
        pg = generators.pg,
        qg = generators.qg,
        qt = generators.qt,
        qb = generators.qb,
        vs = generators.vs,
        ireg = generators.ireg,
        mbase = generators.mbase,
        zr = generators.zr,
        zx = generators.zx,
        rt = generators.rt,
        xt = generators.xt,
        gtap = generators.gtap,
        status = generators.stat,
        rmpct = generators.rmpct,
        pt = generators.pt,
        pb = generators.pb,
        o1 = generators.o1,
        f1 = generators.f1,
        o2 = generators.o2,
        f2 = generators.f2,
        o3 = generators.o3,
        f3 = generators.f3,
        o4 = generators.o4,
        f4 = generators.f4,
        wmod = generators.wmod,
        wpf = generators.wpf,
    ))
end

### === Residual Calculation === ###


# ============================================================================
# Reporting Interface - Vectorized Observation Computation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:PV}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all PV generators in the solution (vectorized).

# Supported Variables
- `:p_gen`: Active power generation (pu) - the scheduled `pg` of in-service units
- `:q_gen`: Reactive power generation (pu) - each unit's share of the generator
  reactive power at its bus (see [`generator_power`](@ref))
"""
function compute_observations(sol::PowerFlowSolution, ::Type{<:PV}, var::Symbol)
    if var == :p_gen
        return scheduled_power(sol.prob.sys.models.PV)
    elseif var == :q_gen
        q = generator_power(sol, :PV).q
        isnothing(q) && error("PV observation :q_gen is not defined in the DCPowerFlow formulation")
        return q
    else
        error("Unknown PV output variable: $var")
    end
end

# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    pv_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for PV generator residuals.

PV generators have two control equations:
1. **Active power injection**: P_inj = -Pg (negative = generation)
2. **Voltage control**: V - Vs = 0 (replaces Q balance equation)

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:p, :v_ctrl)}`: Preallocated output buffers
- `params::NamedTuple{(:pg, :vs)}`: Generator parameters
- `vars::NamedTuple{(:theta, :v)}`: Bus voltage variables from workspace
- `active::AbstractVector{Bool}`: Mask for in-service generators
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)
"""
@inline function pv_residual_kernel!(
    outputs::NamedTuple{(:p, :v_ctrl)},
    params::NamedTuple{(:pg, :vs)},
    vars::NamedTuple{(:theta, :v)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        mask = active[i]

        # PV equations: power injection (positive = into bus) and voltage control
        outputs.p[i] = mask * params.pg[i]
        outputs.v_ctrl[i] = mask * (vars.v[i] - params.vs[i])
    end

    nothing
end

# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, pv::StructArray{<:PV}) where {T, AT}

Calculate PV generator residual contributions.

PV generators have two control equations:
1. **Active power injection**: Injects specified Pg into bus P balance
2. **Voltage control**: Controls voltage to setpoint Vs (replaces Q balance)

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Physics

For a PV bus, the standard bus equations are modified:
- **P equation**: ∑P_inj - ∑P_load - P_branch = 0  (generator contributes +Pg)
- **Q equation REPLACED by**: V - Vs = 0  (voltage control)

The generator's reactive power Q is NOT specified - it adjusts automatically
to maintain voltage within limits [Qb, Qt]. If limits are violated, the bus
converts to PQ type (Q fixed, V free) - handled externally.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `pv`: PV generator model vector with numerical fields

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `active_mask` for out-of-service generators
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: ws.var_vals.v (current voltage)
- Writes: ws.res_vals.p (P injection), ws.res_vals.v_ctrl (V control)
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pv::StructArray{<:PV},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    pv_residual_kernel!(
        # Outputs: residual buffers
        (p = ws.res_vals.p, v_ctrl = ws.res_vals.v_ctrl),

        # Parameters: generator settings (no slack handling needed)
        (pg = pv.pg, vs = pv.vs),

        # Variables: voltage angle and magnitude from workspace
        (theta = ws.var_vals.theta, v = ws.var_vals.v),

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

# CurrentInjection formulation - uses rectangular voltages
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pv::StructArray{<:PV},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    pv_current_kernel!(
        # Outputs: current/voltage residual buffers
        (i_r = ws.res_vals.i_r, i_i_or_vmag = ws.res_vals.i_i_or_vmag),

        # Parameters: generator settings (no slack handling needed)
        (pg = pv.pg, vs = pv.vs),

        # Variables: rectangular voltage components from workspace
        (v_r = ws.var_vals.v_r, v_i = ws.var_vals.v_i),

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

# ============================================================================
# Current Injection Kernel (CPU)
# ============================================================================

"""
    pv_current_kernel!(outputs, params, vars, active, ::Type{<:Vector})

Compute PV generator equations for CurrentInjection formulation.

PV buses have two equations:
  - Equation 1 (i_r): P balance - convert accumulated current to power and compare with -Pg
    (V_r * I_r + V_i * I_i) + Pg = 0  (current times voltage = power)
  - Equation 2 (i_i_or_vmag): V magnitude constraint
    V_r^2 + V_i^2 - V_s^2 = 0

NOTE: The current residuals (i_r, i_i) are accumulated from branches and loads
before this function is called. This function reads those accumulated values to compute
the power balance.
"""
@inline function pv_current_kernel!(
    outputs::NamedTuple{(:i_r, :i_i_or_vmag)},
    params::NamedTuple{(:pg, :vs)},
    vars::NamedTuple{(:v_r, :v_i)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.i_r)
        mask = active[k]

        # PV bus: P balance + V magnitude control
        # P balance: accumulated current times voltage should equal generation
        # NOTE: outputs.i_r[k] and outputs.i_i_or_vmag[k] contain accumulated current from branches/loads
        # Sign convention: positive = injection into bus, so total injection = 0
        p_from_network = vars.v_r[k] * outputs.i_r[k] + vars.v_i[k] * outputs.i_i_or_vmag[k]
        pg_pu = params.pg[k]
        pv_res_p = p_from_network + pg_pu  # P balance: sum of injections = 0

        # V magnitude constraint: |V|^2 - |V_s|^2 = 0
        v_mag_sq = vars.v_r[k] * vars.v_r[k] + vars.v_i[k] * vars.v_i[k]
        v_setpoint_sq = params.vs[k] * params.vs[k]
        pv_res_vmag = v_mag_sq - v_setpoint_sq

        # Write masked results (setter mode: completely replace bus residuals)
        outputs.i_r[k] = mask * pv_res_p
        outputs.i_i_or_vmag[k] = mask * pv_res_vmag
    end

    nothing
end

# ============================================================================
# DC Power Flow Kernel (CPU)
# ============================================================================

"""
    pv_dc_kernel!(outputs, params, vars, active, ::Type{<:Vector})

Compute PV generator equations for DC Power Flow formulation.

PV buses inject active power:
  - Active power injection: P = -Pg (generation)

In DC power flow:
- Voltage magnitudes assumed 1.0 pu (no V control equation)
- Reactive power ignored (no Q equations)
- Only active power balance matters

Sign convention: Generator INJECTS power (negative contribution to bus balance).
"""
@inline function pv_dc_kernel!(
    outputs::NamedTuple{(:p,)},
    params::NamedTuple{(:pg,)},
    vars::NamedTuple{(:theta,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.p)
        mask = active[k]

        # PV bus: Active power injection only
        # Sign convention: positive = injection into bus
        outputs.p[k] = mask * params.pg[k]
    end

    nothing
end

# DC Power Flow calc_residuals! dispatcher
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pv::StructArray{<:PV},
    _external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to DC kernel
    pv_dc_kernel!(
        # Outputs: active power residual only
        (p = ws.res_vals.p,),

        # Parameters: active power setpoint (no slack handling needed)
        (pg = pv.pg,),

        # Variables: angle only (voltage magnitude assumed 1.0 pu)
        (theta = ws.var_vals.theta,),

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end
