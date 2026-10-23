export Slack, SlackVec

"""
    $TYPEDEF

Slack (swing/reference) bus generator data records.

Slack buses provide the voltage and angle reference for the power system.
They share the same data structure as PV generators but have different
equation semantics:
- PV: Power injection + voltage control
- Slack: Angle reference + voltage reference

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct Slack{Tv} <: AbstractModel
    "Bus number or extended bus name enclosed in single quotes"
    i::Int32

    """
    One or two character uppercase non-blank alphanumeric machine identifier.
    Used to distinguish multiple machines at bus i.
    Recommended to use '1' for buses with single machine.
    """
    id::String3

    "Generator active power output in MW (measured, not setpoint). Default: 0.0"
    pg::Tv

    """
    Generator reactive power output in Mvar (measured, not setpoint).
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

@register_model Slack

# Type aliases for StructArray-based vectorized types
"""
    SlackVec{T}

Alias for `StructArray{Slack{T}}`, the struct-of-arrays container holding all
[`Slack`](@ref) instances of a system.
"""
const SlackVec{T} = StructArray{Slack{T}}

get_key_name(::Type{T}) where T<:Slack = [:i, :id]

# Base quantity traits - Slack inherits base voltage from connected bus
base_quantity_role(::Type{<:Slack}, ::Type{BaseVoltage}) = InheritsBase()

# Connectivity declaration - single source of truth for topology AND ExtParams
const SLACK_TERMINALS = (Terminal(:Bus, :i),)  # Single-attachment device

get_terminals(::Type{<:Slack}) = SLACK_TERMINALS

# ============================================================================
# Metadata - PowerInjection Formulation
# ============================================================================

const SLACK_VARS = [
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

# Slack uses SharedRes(:set) for BOTH equations - it's the reference bus
const SLACK_RESIDUALS = [
    ModelResidual(
        :theta_ref,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: θ = θ_ref (angle reference)
        :Bus,
        :p,  # Slack angle reference REPLACES the P equation at slack buses
        :i,
        description="Voltage angle reference equation (replaces P balance)",
    ),
    ModelResidual(
        :v_ref,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: V = V_ref (voltage reference)
        :Bus,
        :q,  # Slack voltage reference REPLACES the Q equation at slack buses
        :i,
        description="Voltage magnitude reference equation (replaces Q balance)",
    )
]

# CurrentInjection formulation - rectangular coordinates
const SLACK_VARS_CURRENT = [
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

const SLACK_RESIDUALS_CURRENT = [
    ModelResidual(
        :v_r_ref,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: V_r = V_s*cos(θ_s)
        :Bus,
        :i_r,
        :i,
        description="Voltage real part reference (replaces I_r balance)",
    ),
    ModelResidual(
        :v_i_ref,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: V_i = V_s*sin(θ_s)
        :Bus,
        :i_i,
        :i,
        description="Voltage imaginary part reference (replaces I_i balance)",
    )
]

# DC Power Flow formulation - angle only
const SLACK_VARS_DC = [
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

const SLACK_RESIDUALS_DC = [
    ModelResidual(
        :theta_ref,
        AlgebRes(),
        SharedRes(:set),  # Setter mode: θ = θ_ref
        :Bus,
        :p,
        :i,
        description="Voltage angle reference equation (replaces P balance)",
    )
]

# Slack parameters for per-unit conversion (same as PV)
const SLACK_PARAMS = [
    ModelParam(:pg, PowerNominal, description="Active power output (measured)", units="MW"),
    ModelParam(:qg, PowerNominal, description="Reactive power output (measured)", units="MVAR"),
]

# Slack observations (outputs from solved power flow)
const SLACK_OBSERVATIONS = [
    ModelObservable(:p_gen, ActivePower, units="pu", description="Active power generated", default_output=true),
    ModelObservable(:q_gen, ReactivePower, units="pu", description="Reactive power generated", default_output=true),
]

# Slack has device power rating
const SLACK_RATINGS = [
    ModelRating(:mbase, NominalPower, description="Generator MVA base")
]

# Slack references bus for base voltage and gets bus angle setpoint
const SLACK_EXT_PARAMS = [
    ExtParam(:basekv, SLACK_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
    ExtParam(:va, SLACK_TERMINALS[1], :va, Dimensionless(), PostPerUnit(),
             required=true, description="Bus voltage angle setpoint in radians")
]

function model_metadata(::Type{Slack})
    metadata = ModelMetadata(
        name = :Slack,
        vars = SLACK_VARS,
        observable_vars = [],
        residuals = SLACK_RESIDUALS,
        params = SLACK_PARAMS,
        observations = SLACK_OBSERVATIONS,
        ratings = SLACK_RATINGS,
        ext_params = SLACK_EXT_PARAMS,
        status_field = :status,
        evaluation_pass = 1000,  # Slack runs in last pass (like PV)
    )

    # Register CurrentInjection formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = SLACK_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = SLACK_RESIDUALS_CURRENT

    # Register DC Power Flow formulation overrides
    metadata.vars_by_formulation[DCPowerFlow] = SLACK_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = SLACK_RESIDUALS_DC

    return metadata
end

# Declare formulation support
supports_formulation(::Type{Slack}, ::Type{PowerInjection}) = true
supports_formulation(::Type{Slack}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{Slack}, ::Type{DCPowerFlow}) = true

# ============================================================================
# Model Format Adapters
# ============================================================================

supports_format(::Type{Slack}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Slack}, ::Type{PSSE}) = (:generators,)

function parse_model(::Type{Slack}, raw, ::FormatSupport{PSSE})
    generators = raw.generators

    # Filter for slack generators only (ide == 3)
    # Note: We need to find which generators are at slack buses
    # This requires bus information to be available
    # For now, we'll return an empty array and let SystemModel do the filtering

    # Return empty StructArray with correct type
    return StructArray{Slack{Float64}}((
        i = Int32[],
        id = String3[],
        pg = Float64[],
        qg = Float64[],
        qt = Float64[],
        qb = Float64[],
        vs = Float64[],
        ireg = Int32[],
        mbase = Float64[],
        zr = Float64[],
        zx = Float64[],
        rt = Float64[],
        xt = Float64[],
        gtap = Float64[],
        status = Bool[],
        rmpct = Float64[],
        pt = Float64[],
        pb = Float64[],
        o1 = Int16[],
        f1 = Union{Missing,Float64}[],
        o2 = Union{Missing,Int16}[],
        f2 = Union{Missing,Float64}[],
        o3 = Union{Missing,Int16}[],
        f3 = Union{Missing,Float64}[],
        o4 = Union{Missing,Int16}[],
        f4 = Union{Missing,Float64}[],
        wmod = Union{Missing,Int8}[],
        wpf = Union{Missing,Float64}[],
    ))
end

# ============================================================================
# Residual Calculation - Reuse PV Kernels with Simplified Logic
# ============================================================================


"""
    slack_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for Slack bus reference equations.

Slack buses provide the voltage and angle reference:
1. **Angle reference**: θ - θ_ref = 0
2. **Voltage reference**: V - V_ref = 0

These equations REPLACE the bus P and Q balance equations (via SharedRes(:set)).

# Arguments
- `outputs::NamedTuple{(:theta_ref, :v_ref)}`: Preallocated output buffers
- `params::NamedTuple{(:vs, :va)}`: Voltage and angle setpoints
- `vars::NamedTuple{(:theta, :v)}`: Bus voltage variables from workspace
- `active::AbstractVector{Bool}`: Mask for in-service generators
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)
"""
@inline function slack_residual_kernel!(
    outputs::NamedTuple{(:theta_ref, :v_ref)},
    params::NamedTuple{(:vs, :va)},
    vars::NamedTuple{(:theta, :v)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.theta_ref)
        mask = active[i]

        # Slack equations: enforce setpoints
        outputs.theta_ref[i] = mask * (vars.theta[i] - params.va[i])
        outputs.v_ref[i] = mask * (vars.v[i] - params.vs[i])
    end

    nothing
end

# ============================================================================
# Data Adapter for PowerInjection Formulation
# ============================================================================

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    slack::StructArray{<:Slack},
    external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    slack_residual_kernel!(
        # Outputs: residual buffers
        (theta_ref = ws.res_vals.theta_ref, v_ref = ws.res_vals.v_ref),

        # Parameters: voltage and angle setpoints
        (vs = slack.vs, va = external.va),

        # Variables: voltage angle and magnitude from workspace
        (theta = ws.var_vals.theta, v = ws.var_vals.v),

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

# ============================================================================
# CurrentInjection Formulation
# ============================================================================

@inline function slack_current_kernel!(
    outputs::NamedTuple{(:v_r_ref, :v_i_ref)},
    params::NamedTuple{(:vs, :va)},
    vars::NamedTuple{(:v_r, :v_i)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.v_r_ref)
        mask = active[k]

        # Compute rectangular voltage setpoints
        v_r_setpoint = params.vs[k] * cos(params.va[k])
        v_i_setpoint = params.vs[k] * sin(params.va[k])

        # Slack equations: V_r = V_r_setpoint, V_i = V_i_setpoint
        outputs.v_r_ref[k] = mask * (vars.v_r[k] - v_r_setpoint)
        outputs.v_i_ref[k] = mask * (vars.v_i[k] - v_i_setpoint)
    end

    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    slack::StructArray{<:Slack},
    external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    slack_current_kernel!(
        (v_r_ref = ws.res_vals.v_r_ref, v_i_ref = ws.res_vals.v_i_ref),
        (vs = slack.vs, va = external.va),
        (v_r = ws.var_vals.v_r, v_i = ws.var_vals.v_i),
        active_mask,
        AT
    )
end

# ============================================================================
# DCPowerFlow Formulation
# ============================================================================

@inline function slack_dc_kernel!(
    outputs::NamedTuple{(:theta_ref,)},
    params::NamedTuple{(:va,)},
    vars::NamedTuple{(:theta,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.theta_ref)
        mask = active[k]

        # Slack equation: θ = θ_ref
        outputs.theta_ref[k] = mask * (vars.theta[k] - params.va[k])
    end

    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    slack::StructArray{<:Slack},
    external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to DC kernel
    slack_dc_kernel!(
        (theta_ref = ws.res_vals.theta_ref,),
        (va = external.va,),
        (theta = ws.var_vals.theta,),
        active_mask,
        AT
    )
end

# ============================================================================
# Reporting Interface - Vectorized Observation Computation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:Slack}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all Slack generators in the solution (vectorized).

# Supported Variables
- `:p_gen`: Active power generation (pu) - each unit's share of the active power balance
  at its bus
- `:q_gen`: Reactive power generation (pu) - each unit's share of the generator
  reactive power at its bus

See [`generator_power`](@ref) for how units at a bus share its power.
"""
function compute_observations(sol::PowerFlowSolution, ::Type{<:Slack}, var::Symbol)
    if var == :p_gen
        return generator_power(sol, :Slack).p
    elseif var == :q_gen
        q = generator_power(sol, :Slack).q
        isnothing(q) && error("Slack observation :q_gen is not defined in the DCPowerFlow formulation")
        return q
    else
        error("Unknown Slack output variable: $var")
    end
end
