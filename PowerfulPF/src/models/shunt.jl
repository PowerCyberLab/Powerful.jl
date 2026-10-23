export Shunt, ShuntVec, AbstractShunt

"""
    AbstractShunt <: AbstractModel

A device that puts a fixed admittance `g + j·b` to ground at one bus: the fixed shunt
([`Shunt`](@ref)) and the switched shunt held at its initial admittance
([`SwitchedShunt`](@ref)). They share the variables, residuals, kernels and observations
written here. Each model gives its admittance parameters through `shunt_admittance_params`,
and the admittance of its instance `i` is `shunt_admittance(params, i)`. The network
builders (Y-bus, fast-decoupled B″, DC injections, and others) find these
models with `fixed_shunt_models`.
"""
abstract type AbstractShunt <: AbstractModel end

"""
    $TYPEDEF

Fixed shunt compensator data records.

Fixed shunts are voltage-dependent reactive power devices used for voltage support
and reactive power compensation in power systems. They typically represent:
- Capacitor banks (negative BL, provides reactive support)
- Shunt reactors (positive BL, absorbs reactive power)

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct Shunt{Tv} <: AbstractShunt

    "Bus number or extended bus name enclosed in single quotes"
    i::Int32

    """
    One- or two-character uppercase non-blank alphanumeric shunt identifier.
    Used to distinguish multiple shunts at bus i. Recommended to use '1' for single shunt.
    """
    id::String3

    "Shunt status: true = in service"
    status::Bool

    """
    Active component of shunt admittance in MW at 1.0 pu voltage.
    Typically zero for fixed shunts. Positive = resistive load.
    """
    gl::Tv

    """
    Reactive component of shunt admittance in Mvar at 1.0 pu voltage.
    PSS/E convention:
    Positive = capacitive (provides reactive support)
    Negative = inductive (absorbs reactive power)
    """
    bl::Tv
end

@register_model Shunt

# Type aliases for StructArray-based vectorized types
"""
    ShuntVec{T}

Alias for `StructArray{Shunt{T}}`, the struct-of-arrays container holding all
[`Shunt`](@ref) instances of a system.
"""
const ShuntVec{T} = StructArray{Shunt{T}}

get_key_name(::Type{T}) where T<:Shunt = [:i, :id]

# Base quantity traits - a shunt inherits base voltage from connected bus
base_quantity_role(::Type{<:AbstractShunt}, ::Type{BaseVoltage}) = InheritsBase()

"""
    shunt_admittance_params(model::StructArray{<:AbstractShunt}) -> NamedTuple

The admittance parameters of the instances of a shunt model, as the shunt kernels read
them with `shunt_admittance`: `(gl, bl)` for a `Shunt`, `(binit,)` for a `SwitchedShunt`.
"""
shunt_admittance_params(shunt::StructArray{<:AbstractShunt}) = (gl = shunt.gl, bl = shunt.bl)

"""
    shunt_admittance(params, i) -> (g, b)

Admittance `g + j·b` (pu) of instance `i` of a shunt model, from its admittance
parameters (`shunt_admittance_params`): `(gl[i], bl[i])` of a fixed shunt, and
`(0, binit[i])` of a switched shunt, which is a pure susceptance. The method is chosen
by the parameters' names at compile time, so one kernel loop serves every shunt model
on the CPU and the GPU.
"""
Base.@propagate_inbounds shunt_admittance(params::NamedTuple{(:gl, :bl)}, i) = (params.gl[i], params.bl[i])

"""
    fixed_shunt_models(sys) -> Vector{Symbol}

The models of `sys` that put a fixed admittance at a bus (`AbstractShunt`): the network
builders (Y-bus, fast-decoupled B″, DC injections, and others) read every
shunt through it and `shunt_admittance`.
"""
fixed_shunt_models(sys) = Symbol[name for (name, model) in pairs(sys.models)
                                 if get_base_model_type(model) <: AbstractShunt]

# Connectivity declaration - single source of truth for topology AND ExtParams
const SHUNT_TERMINALS = (Terminal(:Bus, :i),)  # Single-attachment device

get_terminals(::Type{<:AbstractShunt}) = SHUNT_TERMINALS

const SHUNT_VARS = [
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
    ),
]

const SHUNT_VARS_CURRENT = [
    ModelVar(
        :v_r,
        AlgebVar(),
        :Bus,
        :v_r,
        :i;
        description="Bus voltage real component",
        units="pu",
    ),
    ModelVar(
        :v_i,
        AlgebVar(),
        :Bus,
        :v_i,
        :i;
        description="Bus voltage imaginary component",
        units="pu",
    ),
]

const SHUNT_VARS_DC = [
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

const SHUNT_observable_vars = []

const SHUNT_RESIDUALS = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :p,
        :i,
        description="Active power injection from shunt",
    ),
    ModelResidual(
        :q,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :q,
        :i,
        description="Reactive power injection from shunt",
    )
]

const SHUNT_RESIDUALS_CURRENT = [
    ModelResidual(
        :i_r,
        AlgebRes(),
        SharedRes(:add),
        :Bus,
        :i_r,
        :i,
        description="Real current injection from shunt",
    ),
    ModelResidual(
        :i_i,
        AlgebRes(),
        SharedRes(:add),
        :Bus,
        :i_i,
        :i,
        description="Imaginary current injection from shunt",
    )
]

const SHUNT_RESIDUALS_DC = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes(),
        :Bus,
        :p,
        :i,
        description="Active power injection from shunt (DC approximation)",
    )
]

# Shunt parameters for per-unit conversion
# NOTE: gl and bl are in MW/Mvar at 1.0 pu voltage, need PowerNominal conversion
const SHUNT_PARAMS = [
    ModelParam(:gl, PowerNominal, description="Active component", units="MW"),
    ModelParam(:bl, PowerNominal, description="Reactive component", units="MVAR"),
]

# Shunt observations (computed from voltage and admittance)
const SHUNT_OBSERVATIONS = [
    ModelObservable(:p_injection, ActivePower, units="pu", description="Active power injection", default_output=true),
    ModelObservable(:q_injection, ReactivePower, units="pu", description="Reactive power injection", default_output=true),
]

# Shunt has no device ratings
const SHUNT_RATINGS = ModelRating[]

# Shunt references bus for base voltage - uses SHUNT_TERMINALS
const SHUNT_EXT_PARAMS = [
    ExtParam(:basekv, SHUNT_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
    ExtParam(:ide, SHUNT_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type (1=PQ, 2=PV, 3=Slack)"),
    ExtParam(:vm, SHUNT_TERMINALS[1], :vm, Dimensionless(), PostPerUnit(),
             required=true, description="Case bus voltage magnitude (pu), the operating voltage of the DC model")
]

"""
    shunt_metadata(name::Symbol, params) -> ModelMetadata

Metadata of a shunt model (`AbstractShunt`) named `name` with the per-unit parameters
`params`: the variables, residuals, observations and external parameters every shunt
model shares, in each formulation.
"""
function shunt_metadata(name::Symbol, params)
    metadata = ModelMetadata(
        name = name,
        vars = SHUNT_VARS,
        observable_vars = SHUNT_observable_vars,
        residuals = SHUNT_RESIDUALS,
        params = params,
        observations = SHUNT_OBSERVATIONS,
        ratings = SHUNT_RATINGS,
        ext_params = SHUNT_EXT_PARAMS,
        status_field = :status,
    )

    # Register CurrentInjection formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = SHUNT_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = SHUNT_RESIDUALS_CURRENT

    # Register DC Power Flow formulation overrides
    metadata.vars_by_formulation[DCPowerFlow] = SHUNT_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = SHUNT_RESIDUALS_DC

    return metadata
end

model_metadata(::Type{Shunt}) = shunt_metadata(:Shunt, SHUNT_PARAMS)

# Declare formulation support
supports_formulation(::Type{<:AbstractShunt}, ::Type{PowerInjection}) = true
supports_formulation(::Type{<:AbstractShunt}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{<:AbstractShunt}, ::Type{DCPowerFlow}) = true

### === Model Format Adapters ===

supports_format(::Type{Shunt}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Shunt}, ::Type{PSSE}) = (:fixed_shunts,)

function parse_model(::Type{Shunt}, raw, ::FormatSupport{PSSE})
    shunts = raw.fixed_shunts

    # Create StructArray directly from named tuple of vectors
    return StructArray{Shunt{Float64}}((
        i = shunts.i,
        id = shunts.id,
        status = shunts.status,
        gl = shunts.gl,
        bl = shunts.bl,
    ))
end


### === Residual Calculation === ###


# ============================================================================
# Shunt Physics (Reusable for Residuals and Reports)
# ============================================================================

"""
    calc_shunt_power(v, gl, bl) -> (p=..., q=...)

Power of a fixed shunt at voltage magnitude `v`, from its admittance gl + j*bl.

# Returns
Named tuple `(p=..., q=...)` in pu:
- `p = gl*v²`: active power the conductance consumes (typically gl = 0)
- `q = bl*v²`: reactive power the shunt supplies to the bus; positive for a
  capacitor (bl > 0), negative for a reactor (bl < 0)

The residual kernels turn these into injections into the bus (positive =
injection): P injection `-p`, Q injection `q`.

# Arguments
- `v`: Voltage magnitude (pu)
- `gl`: Active conductance (pu)
- `bl`: Reactive susceptance (pu); positive = capacitor, negative = reactor

# Example
Capacitor bank: bl=+0.84 pu at v=1.0 → q=+0.84 (supplies reactive power)
"""
@inline function calc_shunt_power(v, gl, bl)
    v2 = v * v
    p = gl * v2   # Conductance power consumption
    q = bl * v2   # Capacitor (bl>0) injects positive Q
    return (p=p, q=q)
end

# ============================================================================
# Reporting Interface - Vectorized Observation Computation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:AbstractShunt}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all instances of a shunt model (`Shunt`, `SwitchedShunt`)
in the solution (vectorized).

# Supported Variables
- `:p_injection`: Active power injection (pu)
- `:q_injection`: Reactive power injection (pu)

# Implementation Notes
Returns injections into the bus at the solved voltage magnitudes, with the same sign
as the residual kernel (positive = injection), from the admittance `g + j·b` of each
instance (`shunt_admittance`):
- `:p_injection` = -g * |V|²  (the conductance consumes; typically zero)
- `:q_injection` = b * |V|²   (capacitors (b > 0) positive, reactors (b < 0) negative)

Voltages are extracted from solution via AddressManager (same pattern as PQ model).
"""
function compute_observations(sol::PowerFlowSolution, ::Type{M}, var::Symbol) where {M <: AbstractShunt}
    name = nameof(M)
    shunt = sol.prob.sys.models[name]
    var in (:p_injection, :q_injection) || error("Unknown $name output variable: $var")
    params = shunt_admittance_params(shunt)

    # Terminal bus voltages (any formulation); out-of-service shunts inject nothing
    v_vals, _ = terminal_voltages(sol, name, :i)
    n = length(shunt)
    result = Vector{Float64}(undef, n)

    if sol.prob.formulation isa DCPowerFlow
        var == :p_injection || error("$name observation :$var is not defined in the DCPowerFlow formulation")
        @inbounds for i in 1:n
            g, _ = shunt_admittance(params, i)
            result[i] = -shunt.status[i] * dc_shunt_power(g, v_vals[i])
        end
        return result
    end

    @inbounds for i in 1:n
        power = calc_shunt_power(v_vals[i], shunt_admittance(params, i)...)
        # Injection into the bus, as in the residual kernel
        result[i] = shunt.status[i] * (var == :p_injection ? -power.p : power.q)
    end

    return result
end

# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    shunt_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for shunt residuals (every `AbstractShunt` model).

A shunt injects power proportional to voltage squared (positive = injection into
the bus, see `calc_shunt_power`), from its admittance g + j·b (`shunt_admittance`):
- P injection = -g * v²  (the conductance consumes; typically g = 0)
- Q injection = b * v²   (capacitor, b > 0, injects; reactor, b < 0, absorbs)

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:p, :q)}`: Preallocated output buffers
- `params::NamedTuple`: Shunt admittance parameters (pu), `shunt_admittance_params`
  (`(gl, bl)` of a fixed shunt, `(binit,)` of a switched shunt)
- `vars::NamedTuple{(:v,)}`: Voltage magnitude from workspace
- `active::AbstractVector{Bool}`: Mask for in-service shunts
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)

# Example
Bus with BL=+84 Mvar (capacitor) @ 100 MVA base → bl_pu = +0.84 pu.
At V=1.0 pu: Q injection = +0.84 pu.
"""
@inline function shunt_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple,
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        # Branchless masking - only active status (bus type filtering handled by scatter :set mode)
        mask = active[i]

        # Calculate power using extracted physics function
        power = calc_shunt_power(vars.v[i], shunt_admittance(params, i)...)

        # Sign convention: positive = injection into bus
        # P: conductance consumes, so injection is negative
        # Q: capacitor (bl>0) injects positive Q (already correct in calc_shunt_power)
        outputs.p[i] = mask * (-power.p)
        outputs.q[i] = mask * power.q
    end

    nothing
end

"""
    calc_shunt_current(v_r, v_i, gl, bl) -> (i_r=..., i_i=...)

Calculate fixed shunt current injection.

# Current Injection Model
Current injection from admittance: I = Y*V where Y = gl + j*bl

Expanding with V = v_r + j*v_i and Y = gl + j*bl:
- I_r = gl*v_r - bl*v_i  (real component)
- I_i = gl*v_i + bl*v_r  (imaginary component)

# PSS/E Convention
- BL positive = capacitor (supplies reactive power)
- BL negative = reactor (absorbs reactive power)

# ANDES Convention
- Positive current = injection into bus
- Negative current = withdrawal from bus

# Arguments
- `v_r`: Voltage real component (pu)
- `v_i`: Voltage imaginary component (pu)
- `gl`: Active conductance (pu)
- `bl`: Reactive susceptance (pu, PSS/E convention)

# Returns
Named tuple `(i_r=..., i_i=...)` with current in pu (ANDES convention).

# Example
Capacitor bank: bl=+0.84 pu at v=1.0+j0.0 → I_r=0, I_i=-0.84 (capacitive current, ANDES convention)
"""
@inline function calc_shunt_current(v_r, v_i, gl, bl)
    # Complex multiplication: I = Y*V = (gl + j*bl)*(v_r + j*v_i)
    # = gl*v_r + j*gl*v_i + j*bl*v_r + j²*bl*v_i
    # = (gl*v_r - bl*v_i) + j*(gl*v_i + bl*v_r)
    i_r = gl * v_r - bl * v_i
    i_i = gl * v_i + bl * v_r
    return (i_r=i_r, i_i=i_i)
end

"""
    shunt_current_kernel!(outputs, params, external, vars, active, ::Type{AT})

Pure computation kernel for shunt current residuals (every `AbstractShunt` model).

A shunt injects current proportional to voltage: I = Y*V where Y = g + j*b
(`shunt_admittance`)

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:i_r, :i_i)}`: Preallocated output buffers
- `params::NamedTuple`: Shunt admittance parameters (pu), `shunt_admittance_params`
- `external::NamedTuple`: External parameters (ide for bus type)
- `vars::NamedTuple{(:v_r, :v_i)}`: Rectangular voltage from workspace
- `active::AbstractVector{Bool}`: Mask for in-service shunts
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)

# Physics
PSS/E convention for BL (positive = capacitor, negative = reactor)
Current injection: I = Y*V = (gl + j*bl)*(v_r + j*v_i)
- I_r = gl*v_r - bl*v_i
- I_i = gl*v_i + bl*v_r

# Bus Type Handling
Shunt contributes current to ALL non-slack buses:
- PQ buses (ide=1): Receive i_r and i_i contributions
- PV buses (ide=2): Receive i_r and i_i contributions (PV will gather and replace later)
- Slack buses (ide=3): No contributions (voltage set directly)

# Example
Bus with BL=+84 Mvar (capacitor, PSS/E convention) @ 100 MVA base → bl_pu = +0.84 pu
At V=1.0+j0.0 pu: I_r=0, I_i=+0.84 pu (capacitive current injection)
"""
@inline function shunt_current_kernel!(
    outputs::NamedTuple{(:i_r, :i_i)},
    params::NamedTuple,
    external::NamedTuple,
    vars::NamedTuple{(:v_r, :v_i)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.i_r)
        # Single mask: contribute to ALL non-slack buses
        # PV buses will receive contributions and later replace with control equations
        mask = active[i] & (external.ide[i] != BUS_SLACK)

        # Calculate current using extracted physics function
        current = calc_shunt_current(vars.v_r[i], vars.v_i[i], shunt_admittance(params, i)...)

        # Sign convention: Y*V gives current INTO shunt (to ground), negate for current INTO bus
        outputs.i_r[i] = mask * (-current.i_r)
        outputs.i_i[i] = mask * (-current.i_i)
    end

    nothing
end

"""
    shunt_dc_kernel!(outputs, params, external, active, ::Type{AT})

Pure computation kernel for shunt DC power flow residuals (every `AbstractShunt` model).

In DC power flow, a shunt draws constant active power through its conductance g at the
case bus voltage: P = g * vm² (`dc_shunt_power`). Its susceptance (reactive power) is
ignored in the DC formulation, so a switched shunt (a pure susceptance) draws none.

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:p,)}`: Preallocated output buffer
- `params::NamedTuple`: Shunt admittance parameters (pu), `shunt_admittance_params`
- `external::NamedTuple`: External parameters (ide for bus type)
- `active::AbstractVector{Bool}`: Mask for in-service shunts
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)

# Physics
DC approximation for shunt active power:
- P = g * vm² at the case bus voltage (constant)
- Typically g ≈ 0 for most shunts (primarily reactive elements); zero for a switched shunt

# Bus Type Handling
Shunt contributes active power to ALL non-slack buses (same as AC formulation).
"""
@inline function shunt_dc_kernel!(
    outputs::NamedTuple{(:p,)},
    params::NamedTuple,
    external::NamedTuple,
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        # Mask out slack buses (same as AC formulation)
        mask = active[i] & (external.ide[i] != BUS_SLACK)

        # DC model: conductance draws g·vm² at the case voltage; injection is negative
        g, _ = shunt_admittance(params, i)
        outputs.p[i] = mask * (-dc_shunt_power(g, external.vm[i]))
    end

    nothing
end

# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, shunt::StructArray{<:AbstractShunt}, external::NamedTuple) where {T, AT}

Calculate the residual contributions of a shunt model (`Shunt`, `SwitchedShunt`).

A shunt is a voltage-dependent admittance g + j·b (`shunt_admittance`) that injects power
proportional to v² (positive = injection into the bus, see `shunt_residual_kernel!`):
- P injection = -g * v²  (the conductance consumes; typically zero)
- Q injection = b * v²
  - Positive b = capacitive (supplies reactive power)
  - Negative b = inductive (consumes reactive power)

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `shunt`: Shunt model vector (any `AbstractShunt`)
- `external`: External parameters (ide for slack bus masking)

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `ws.active_mask` for out-of-service shunts
- Masks slack buses using `external.ide`
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: `ws.var_vals.v` (voltage magnitude)
- Writes: `ws.res_vals.p`, `ws.res_vals.q` (power residuals)
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    shunt::StructArray{<:AbstractShunt},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    shunt_residual_kernel!(
        # Outputs: residual buffers
        (p = ws.res_vals.p, q = ws.res_vals.q),

        # Parameters: shunt admittance
        shunt_admittance_params(shunt),

        # Variables: voltage from workspace
        (v = ws.var_vals.v,),

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, shunt::StructArray{<:AbstractShunt}, external::NamedTuple, ::CurrentInjection) where {T, AT}

Calculate the current injection residual contributions of a shunt model (`Shunt`, `SwitchedShunt`).

A shunt draws current proportional to voltage, I = Y*V with Y = g + j*b (`shunt_admittance`),
injected into the bus with the opposite sign (see `shunt_current_kernel!`)

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `shunt`: Shunt model vector (any `AbstractShunt`)
- `external`: External parameters (ide for bus type masking)
- `::CurrentInjection`: Formulation tag for dispatch

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `ws.active_mask` for out-of-service shunts
- Masks slack buses using `external.ide`
- Contributes to ALL non-slack buses (including PV buses)
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: `ws.var_vals.v_r`, `ws.var_vals.v_i` (rectangular voltage)
- Writes: `ws.res_vals.i_r`, `ws.res_vals.i_i` (current residuals)
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    shunt::StructArray{<:AbstractShunt},
    external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    shunt_current_kernel!(
        # Outputs: current residual buffers
        (i_r = ws.res_vals.i_r, i_i = ws.res_vals.i_i),

        # Parameters: shunt admittance
        shunt_admittance_params(shunt),

        # External parameters: bus type for slack masking
        external,

        # Variables: rectangular voltage from workspace
        (v_r = ws.var_vals.v_r, v_i = ws.var_vals.v_i),

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, shunt::StructArray{<:AbstractShunt}, external::NamedTuple, ::DCPowerFlow) where {T, AT}

Calculate the DC power flow residual contributions of a shunt model (`Shunt`, `SwitchedShunt`).

In DC power flow, a shunt contributes only the constant active power its conductance draws at
the case bus voltage (P = g * vm², `dc_shunt_power`); reactive power is ignored, so a switched
shunt contributes none.

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `shunt`: Shunt model vector (any `AbstractShunt`)
- `external`: External parameters (ide for bus type masking)
- `::DCPowerFlow`: Formulation tag for dispatch

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `ws.active_mask` for out-of-service shunts
- Masks slack buses using `external.ide`
- Typically g ≈ 0 for most shunts (primarily reactive elements)
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- No variables read (DC uses constant power)
- Writes: `ws.res_vals.p` (active power residual only)
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    shunt::StructArray{<:AbstractShunt},
    external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    shunt_dc_kernel!(
        # Outputs: active power residual buffer only
        (p = ws.res_vals.p,),

        # Parameters: shunt admittance (only its conductance is used in DC)
        shunt_admittance_params(shunt),

        # External parameters: bus type for slack masking
        external,

        # Metadata: active status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

