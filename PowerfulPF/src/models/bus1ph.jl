export Bus, BusVec


"""
    $TYPEDEF

Network bus data records.

# Fields

$TYPEDFIELDS
"""
Base.@kwdef mutable struct Bus{Tv} <: AbstractModel
    "Bus number (1 to 999997)."
    i::Int32

    """
    Alphanumeric identifier assigned to bus 'I'.

    The name may be up to twelve characters and must be enclosed in single
    quotes. NAME may contain any combination of blanks, uppercase letters,
    numbers and special characters, but the first character must not be a minus
    sign.
    """
    name::String15

    "Bus base voltage in kV."
    basekv::Tv

    """
    Bus type code:
    - 1 - load bus or other bus without any generator boundary condition.
    - 2 - generator or plant bus either regulating voltage or with a fixed reactive power (Mvar). A generator that reaches its reactive power limit will no longer control voltage but rather hold reactive power at its limit.
    - 3 - swing bus or slack bus. It has no power or reactive limits and regulates voltage at a fixed reference angle.
    - 4 - disconnected or isolated bus.
    """
    ide::Int8

    "Area number. 1 through the maximum number of areas at the current size level."
    area::Int16

    "Zone number. 1 through the maximum number of zones at the current size level."
    zone::Int16

    "Owner number. 1 through the maximum number of owners at the current size level."
    owner::Int16

    "Bus voltage magnitude in pu."
    vm::Tv

    "Bus voltage angle in degrees."
    va::Tv

    "Normal voltage magnitude high limit in pu."
    nvhi::Tv

    "Normal voltage magnitude low limit in pu."
    nvlo::Tv

    "Emergency voltage magnitude high limit in pu."
    evhi::Tv

    "Emergency voltage magnitude low limit in pu."
    evlo::Tv
end

@register_model Bus

# Type aliases for StructArray-based vectorized types
"""
    BusVec{T}

Alias for `StructArray{Bus{T}}`, the struct-of-arrays container holding all
[`Bus`](@ref) instances of a system.
"""
const BusVec{T} = StructArray{Bus{T}}

get_key_name(::Type{Bus}) = :i

# Base quantity traits - Bus owns base voltage
base_quantity_role(::Type{Bus}, ::Type{BaseVoltage}) = OwnsBase()

# Define Bus model variables (Power Injection formulation - polar coordinates)
const BUS1PH_VARS = [
    ModelVar(
        :theta,
        AlgebVar(),
        description="Bus voltage angle",
        units="rad"
    ),

    ModelVar(
        :v,
        AlgebVar(),
        description="Bus voltage magnitude",
        units="pu"
    )
]

# Current Injection formulation - rectangular coordinates
const BUS1PH_VARS_CURRENT = [
    ModelVar(
        :v_r,
        AlgebVar(),
        description="Bus voltage real part",
        units="pu"
    ),

    ModelVar(
        :v_i,
        AlgebVar(),
        description="Bus voltage imaginary part",
        units="pu"
    )
]

# DC Power Flow formulation - angles only
const BUS1PH_VARS_DC = [
    ModelVar(
        :theta,
        AlgebVar(),
        description="Bus voltage angle",
        units="rad"
    )
]

const BUS1PH_observable_vars = [:theta, :v]
const BUS1PH_observable_vars_CURRENT = [:v_r, :v_i]
const BUS1PH_observable_vars_DC = [:theta]

# Bus observations (output variables from DAE solution)
const BUS_OBSERVATIONS = [
    ModelObservable(:theta, AngleDegrees; description="Voltage angle", units="rad", default_output=true),
    ModelObservable(:v, VoltageNominal; description="Voltage magnitude", units="pu", default_output=true),
]

# Define residual variables (Power Injection formulation)
const BUS1PH_RESIDUALS = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes();
        description="Active power balance equation",
    ),
    ModelResidual(
        :q,
        AlgebRes(),
        SharedRes();
        description="Reactive power balance equation",
    )
]

# Current Injection formulation
const BUS1PH_RESIDUALS_CURRENT = [
    ModelResidual(
        :i_r,
        AlgebRes(),
        SharedRes();
        description="Real current balance equation",
    ),
    ModelResidual(
        :i_i,
        AlgebRes(),
        SharedRes();
        description="Imaginary current balance equation",
    )
]

# DC Power Flow formulation - active power only
const BUS1PH_RESIDUALS_DC = [
    ModelResidual(
        :p,
        AlgebRes(),
        SharedRes();
        description="Active power balance equation (DC approximation)",
    )
]

# Bus parameters
const BUS_PARAMS = [
    ModelParam(:vm, NoPU; description="Voltage magnitude setpoint", units="pu"),
    ModelParam(:va, AngleDegrees; description="Voltage angle setpoint", units="deg")
]

# Bus has no device ratings (it IS the reference)
const BUS_RATINGS = ModelRating[]

# Bus has no external parameter references (it provides them to others)
const BUS_EXT_PARAMS = ExtParam[]

function model_metadata(::Type{Bus})
    metadata = ModelMetadata(
        name = :Bus,
        vars = BUS1PH_VARS,
        observable_vars = BUS1PH_observable_vars,
        residuals = BUS1PH_RESIDUALS,
        params = BUS_PARAMS,
        observations = BUS_OBSERVATIONS,
        ratings = BUS_RATINGS,
        ext_params = BUS_EXT_PARAMS,
    )

    # Register Current Injection formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = BUS1PH_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = BUS1PH_RESIDUALS_CURRENT

    # Register DC Power Flow formulation overrides
    metadata.vars_by_formulation[DCPowerFlow] = BUS1PH_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = BUS1PH_RESIDUALS_DC

    return metadata
end

# Declare formulation support
supports_formulation(::Type{Bus}, ::Type{PowerInjection}) = true
supports_formulation(::Type{Bus}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{Bus}, ::Type{DCPowerFlow}) = true

# ============================================================================
# Group Interface Implementation - :bus_group
# ============================================================================

"""
Bus implements the :bus_group interface, providing voltage variables (v, theta)
that other models (generators, loads) can connect to. This allows generators
to use `Terminal(:bus_group, :i)` instead of `Terminal(:Bus, :i)`, enabling
connection to any bus-like model (Bus, or another model implementing `:bus_group`).
"""
PowerfulCore.implements_group(::Type{<:Bus}) = :bus_group

"""
Map group canonical variable names to Bus's actual variable names.
For Bus, the mapping is identity: :v → :v, :theta → :theta. The port voltage of the voltage form
(`InjectionForms`) is the bus voltage: a bus realizes that form algebraically, I = (E - V)/Z.
"""
function PowerfulCore.group_var_mapping(::Type{<:Bus})
    Dict(
        :v => :v,           # Voltage magnitude
        :theta => :theta,   # Voltage angle
        :v_port => :v,      # Port voltage of the voltage form: the bus voltage
        :theta_port => :theta,
    )
end

# Map group canonical residual names to Bus residuals for bus_group
PowerfulCore.group_res_mapping(::Type{<:Bus}) = Dict{Symbol, Symbol}(
    :p => :p,  # Active power balance
    :q => :q   # Reactive power balance
)

### === Format Support === ###

supports_format(::Type{Bus}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Bus}, ::Type{PSSE}) = (:buses,)

function parse_model(::Type{Bus}, raw, ::FormatSupport{PSSE})
    buses = raw.buses

    # Create StructArray directly from named tuple of vectors
    return StructArray{Bus{Float64}}((
        i = buses.i,
        name = buses.name,
        basekv = buses.basekv,
        ide = buses.ide,
        area = buses.area,
        zone = buses.zone,
        owner = buses.owner,
        vm = buses.vm,
        va = buses.va,
        nvhi = buses.nvhi,
        nvlo = buses.nvlo,
        evhi = buses.evhi,
        evlo = buses.evlo,
    ))
end

# ============================================================================
# Residual Calculation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:Bus}, var::Symbol)

Bus voltage magnitude `:v` (pu) and angle `:theta` (rad) in any formulation, from
[`bus_voltages`](@ref). `DCPowerFlow` does not solve the magnitudes: `:v` is not
defined there.
"""
function compute_observations(sol::PowerFlowSolution, ::Type{<:Bus}, var::Symbol)
    v, theta = bus_voltages(sol)
    var === :theta && return theta
    var === :v || error("Unknown Bus output variable: $var")
    sol.prob.formulation isa DCPowerFlow &&
        error("Bus observation :v is not defined in the DCPowerFlow formulation")
    return v
end

# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    bus_residual_kernel!(outputs, params, vars, ide, ::Type{AT})

Pure computation kernel for Bus residuals.

**Bus model is purely passive** - all residuals are zero for all bus types.
The Bus model only provides the shared residual buffers (p, q) that other models
contribute to:
  - Branch models add power flow contributions
  - PQ load models add consumption
  - PV generator models add generation and handle slack bus constraints

**Note**: Slack bus equations (θ = θ_ref, V = V_ref) are implemented in the
PV model, NOT in the Bus model.

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:p, :q)}`: Preallocated output buffers (all zeros)
- `params::NamedTuple{(:vm, :va)}`: Bus voltage setpoints (unused - for consistency)
- `vars::NamedTuple{(:theta, :v)}`: Current state (unused - for consistency)
- `ide::AbstractVector{Int8}`: Bus type (unused - for consistency)
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuArray for GPU)
"""
@inline function bus_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:vm, :va)},
    vars::NamedTuple{(:theta, :v)},
    ide::AbstractVector{Int8},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        # Bus model is purely passive - all residuals are zero
        # Slack equations have been moved to PV model
        outputs.p[i] = zero(eltype(outputs.p))
        outputs.q[i] = zero(eltype(outputs.q))
    end

    nothing
end

# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, bus::StructArray{<:Bus}, external::NamedTuple, ::PowerInjection) where {T, AT}

Calculate Bus residual contributions for PowerInjection formulation.

**Bus model is purely passive** - contributes zero to all residuals for all bus types.
The Bus model only owns the shared residual buffers (p, q) that other models write to:
  - Branch models add power flow contributions
  - PQ load models add power consumption
  - PV generator models add power generation and implement slack bus constraints

**Important**: Slack bus equations (θ = θ_ref, V = V_ref) are implemented in the
PV model, NOT in the Bus model. The Bus model simply provides zero contribution
regardless of bus type.

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `bus`: Bus model vector with numerical fields
- `external`: External parameters (unused for Bus model)

# Notes
- Zero allocation - operates on pre-allocated buffers
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Writes: ws.res_vals.p = 0, ws.res_vals.q = 0 (passive model)
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    bus::StructArray{<:Bus},
    external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    # Bus doesn't use external (ide is its own parameter)
    # Note: Bus doesn't use active_mask - it initializes residuals for all buses.
    # The active_mask parameter is for API consistency.
    bus_residual_kernel!(
        # Outputs: residual buffers
        (p = ws.res_vals.p, q = ws.res_vals.q),

        # Parameters: voltage reference
        (vm = bus.vm, va = bus.va),

        # Variables: current state from workspace
        (theta = ws.var_vals.theta, v = ws.var_vals.v),

        # Bus type indicator
        bus.ide,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

# CurrentInjection formulation - Bus is still passive (zero residuals)
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    _bus::StructArray{<:Bus},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Bus model is passive in all formulations - zero all residuals
    # Note: active_mask not used - for API consistency
    fill!(ws.res_vals.i_r, zero(T))
    fill!(ws.res_vals.i_i, zero(T))
    nothing
end

# DC Power Flow formulation - Bus is still passive (zero residuals)
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    _bus::StructArray{<:Bus},
    _external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Bus model is passive in all formulations - zero all residuals
    # Note: active_mask not used - for API consistency
    fill!(ws.res_vals.p, zero(T))
    nothing
end

"""
    set_algebraic_guess!(y, am, sys, ::Type{Bus}, formulation)

Populate the algebraic guess for bus voltages based on the chosen formulation.
"""
function set_algebraic_guess!(
    y::AbstractVector,
    am::AddressManager,
    sys::SystemModel,
    ::Type{Bus},
    formulation::FormulationTrait
)
    haskey(sys.models, :Bus) || return
    haskey(am.assignments, :Bus) || return
    vars = am.assignments[:Bus].by_type[AlgebVar()]
    bus = sys.models.Bus

    if formulation isa PowerInjection
        if haskey(vars, :v) && haskey(vars, :theta)
            y[vars[:v]] .= bus.vm
            y[vars[:theta]] .= bus.va  # already radians post per-unit
        end
    elseif formulation isa CurrentInjection
        if haskey(vars, :v_r) && haskey(vars, :v_i)
            y[vars[:v_r]] .= bus.vm .* cos.(bus.va)
            y[vars[:v_i]] .= bus.vm .* sin.(bus.va)
        end
    end

    return nothing
end

