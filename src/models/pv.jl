export PV, PVVec, PVNum, PVNumVec

abstract type AbstractGenerator{Tv} <: AbstractModel end

"""
    $TYPEDEF

PV generator data records.

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct PV{Tv} <: AbstractGenerator{Tv}
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
    Machine status: `1` for in-service, `0` for out-of-service. Default: 1
    """
    stat::Bool

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
const PVVec{T} = StructArray{PV{T}}
const PVNum{T} = PV{T}  # No longer needed, but kept for compatibility
const PVNumVec{T} = StructArray{PV{T}}  # Same as PVVec

get_key_name(::Type{T}) where T<:PV = [:i, :id]

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

const PV_OUTPUT_VARS = []

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
        SharedRes(),
        :Bus,
        :q,  # PV voltage control REPLACES the Q equation at PV buses
        :i,
        description="Voltage magnitude control equation (replaces Q balance)",
    )
]

function model_metadata(::Type{PV}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :PV,
        vars = PV_VARS,
        output_vars = PV_OUTPUT_VARS,
        residuals = PV_RESIDUALS,
        layout = layout
    )
end

### === Model Format Adapters ===

supports_format(::Type{PV}, ::Type{PSSE}) = FormatSupport{PSSE}()

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
        stat = generators.stat,
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

using ..PowerCore: ResidualWorkspace

# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    pv_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for PV generator residuals.

PV generators have two control equations:
1. **Active power injection**: P_inj = Pg
2. **Voltage control**: V - Vs = 0 (replaces Q balance equation)

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:p, :v_ctrl)}`: Preallocated output buffers
- `params::NamedTuple{(:pg, :vs)}`: Generator parameters (active power, voltage setpoint)
- `vars::NamedTuple{(:v,)}`: Voltage magnitude from workspace
- `active::AbstractVector{Bool}`: Mask for in-service generators
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)
"""
@inline function pv_residual_kernel!(
    outputs::NamedTuple{(:p, :v_ctrl)},
    params::NamedTuple{(:pg, :vs)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        # Branchless computation using mask
        mask = active[i]

        # 1. Active power injection (positive into bus)
        outputs.p[i] = mask * params.pg[i]

        # 2. Voltage control equation: V - Vs = 0
        # This REPLACES the Q balance equation at the PV bus
        outputs.v_ctrl[i] = mask * (vars.v[i] - params.vs[i])
    end

    nothing
end

# ============================================================================
# GPU Kernel Implementation
# ============================================================================

"""
CUDA kernel for PV generator residuals - device function executed on GPU.
"""
function pv_residual_cuda_kernel!(
    outputs_p, outputs_v_ctrl,
    params_pg, params_vs,
    vars_v,
    active,
    n
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if i <= n
        # Branchless computation using mask
        mask = active[i] ? one(eltype(outputs_p)) : zero(eltype(outputs_p))

        @inbounds begin
            # 1. Active power injection
            outputs_p[i] = mask * params_pg[i]

            # 2. Voltage control equation
            outputs_v_ctrl[i] = mask * (vars_v[i] - params_vs[i])
        end
    end

    return nothing
end

"""
    pv_residual_kernel!(outputs, params, vars, active, ::Type{<:CuArray})

GPU dispatch for PV residual kernel using CUDA.
"""
@inline function pv_residual_kernel!(
    outputs::NamedTuple{(:p, :v_ctrl)},
    params::NamedTuple{(:pg, :vs)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:CuArray}  # GPU dispatch
)
    n = length(outputs.p)

    # GPU launch configuration
    threads = 256
    blocks = cld(n, threads)

    # Launch CUDA kernel (async on default stream)
    CUDA.@cuda threads=threads blocks=blocks pv_residual_cuda_kernel!(
        outputs.p, outputs.v_ctrl,
        params.pg, params.vs,
        vars.v,
        active,
        n
    )

    # No synchronize - kernels execute in order on default stream
    nothing
end

# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================

"""
    calc_residuals!(ws::ResidualWorkspace{T, AT}, pv::StructArray{<:PV{T}}) where {T, AT}

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
- Respects `ws.active_mask` for out-of-service generators
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: ws.var_vals.v (current voltage)
- Writes: ws.res_vals.p (P injection), ws.res_vals.v_ctrl (V control)
"""
@inline function calc_residuals!(
    ws::ResidualWorkspace{T, AT},
    pv::StructArray{<:PV}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    pv_residual_kernel!(
        # Outputs: residual buffers
        (p = ws.res_vals.p, v_ctrl = ws.res_vals.v_ctrl),

        # Parameters: generator settings
        (pg = pv.pg, vs = pv.vs),

        # Variables: voltage from workspace
        (v = ws.var_vals.v,),

        # Metadata: active status
        ws.active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end


### === Begin Tests === ###

@testitem "PV" begin
    using PowerFlowData
    using Powerful
    using Powerful.Models
    using Powerful.PowerCore

    CASE = PowerFlowData.parse_network(
        joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    )

    pv = parse_model(PV, CASE)
    @test pv isa PVVec

end

@testitem "PV calc_residuals! Unit Test" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    @testset "PV Residual Calculation with Real Data" begin
        # Load real test data
        case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
        case = PowerFlowData.parse_network(case_path)

        # Parse models
        bus = parse_model(Bus, case)
        pv = parse_model(PV, case)

        # Create system with auto_setup (modern pattern)
        am = AddressManager()
        sys = SystemModel(am, (Bus = bus, PV = pv))
        dae = DAE(Float64, sys)

        # Get metadata for workspace creation
        pv_metadata = model_metadata(PV)

        # Create workspace (Vector for CPU)
        ws = ResidualWorkspace(Float64, pv_metadata, am, pv, Vector)

        # Gather variables
        gather_variables!(ws, dae)

        # Zero residuals
        zero_residuals!(ws)

        # Calculate residuals
        calc_residuals!(ws, pv)

        # Verify P residuals equal Pg for active generators
        @testset "Active Power Injection" begin
            for i in 1:ws.n_instances
                if ws.active_mask[i]
                    @test ws.res_vals.p[i] ≈ pv.pg[i]
                end
            end
        end

        # Verify voltage control residuals are V - Vs
        @testset "Voltage Control Equation" begin
            for i in 1:ws.n_instances
                if ws.active_mask[i]
                    v_actual = ws.var_vals.v[i]
                    v_setpoint = pv.vs[i]
                    expected_residual = v_actual - v_setpoint

                    @test ws.res_vals.v_ctrl[i] ≈ expected_residual

                    # Residual should be close to zero if initialized well
                    # (voltages initialized from data)
                    @test abs(ws.res_vals.v_ctrl[i]) < 0.1  # Within 0.1 pu
                end
            end
        end

        # Verify inactive generators have zero residuals
        @testset "Inactive Generators" begin
            for i in 1:ws.n_instances
                if !ws.active_mask[i]
                    @test ws.res_vals.p[i] == 0.0
                    @test ws.res_vals.v_ctrl[i] == 0.0
                end
            end
        end
    end
end