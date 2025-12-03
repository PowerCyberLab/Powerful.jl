export PQ, PQVec, PQNum, PQNumVec

"""
    $TYPEDEF

Static PQ load data records.

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct PQ{Tv} <: AbstractLoad{Tv}
    
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
    Negative for inductive load, positive for capacitive load
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
const PQVec{T} = StructArray{PQ{T}}
const PQNum{T} = PQ{T}  # No longer needed, but kept for compatibility
const PQNumVec{T} = StructArray{PQ{T}}  # Same as PQVec

get_key_name(::Type{T}) where T<:PQ = [:i, :id]

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

const PQ_OUTPUT_VARS = []

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

function model_metadata(::Type{PQ}; layout::LayoutStrategy=ContiguousVariables())
    ModelMetadata(
        name = :PQ,
        vars = PQ_VARS,
        output_vars = PQ_OUTPUT_VARS,
        residuals = PQ_RESIDUALS,
        layout = layout
    )
end

### === Model Format Adapters ===

supports_format(::Type{PQ}, ::Type{PSSE}) = FormatSupport{PSSE}()

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
        yq = loads.yq,
        owner = loads.owner,
        scale = loads.scale,
        intrpt = loads.intrpt,
    ))
end


### === Residual Calculation === ###

using ..PowerCore: ResidualWorkspace

# ============================================================================
# New Design: Pure Computation Kernel (CPU/GPU Ready)
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
- P_load = pl + ip*v + yp*v²
- Q_load = ql + iq*v + yq*v²

# Arguments
- `outputs::NamedTuple{(:p, :q)}`: Preallocated output buffers
- `params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq)}`: Load parameters
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
    (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq),
    (v = voltage_vector,),
    active_mask,
    Vector{Float64}
)

# GPU dispatch (automatically uses CUDA kernel)
pq_residual_kernel!(
    (p = p_res_gpu, q = q_res_gpu),
    (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq),
    (v = voltage_vector_gpu,),
    active_mask_gpu,
    CuVector{Float64}
)
```
"""
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    # Pure vectorized computation - no addressing logic
    @inbounds @simd for i in eachindex(outputs.p)
        # Branchless: multiply by active mask instead of conditional
        # This enables SIMD and GPU vectorization
        mask = active[i]

        vi = vars.v[i]
        v2 = vi * vi

        # Load consumes power, so negative sign for injection into bus
        # ZIP model: constant P + constant I * v + constant Z * v²
        outputs.p[i] = -mask * (params.pl[i] + params.ip[i] * vi + params.yp[i] * v2)
        outputs.q[i] = -mask * (params.ql[i] + params.iq[i] * vi + params.yq[i] * v2)
    end

    nothing
end

# ============================================================================
# GPU Kernel Implementation
# ============================================================================

"""
    pq_residual_cuda_kernel!(outputs_p, outputs_q, params_pl, params_ql, params_ip,
                             params_iq, params_yp, params_yq, vars_v, active, n)

CUDA kernel for PQ load residuals - device function executed on GPU.

Implements identical logic to `pq_residual_kernel!` but optimized for GPU execution.
Each thread processes one load instance.

# Thread/Block Configuration
- Each thread handles one load (index `i`)
- Recommended: 256 threads per block
- Blocks: `cld(n_loads, 256)`

# Arguments
- All arguments are GPU arrays (CuArray or CuDeviceArray)
- Same computation as CPU kernel: ZIP load model with branchless masking
"""
function pq_residual_cuda_kernel!(
    outputs_p, outputs_q,
    params_pl, params_ql, params_ip, params_iq, params_yp, params_yq,
    vars_v,
    active,
    n
)
    # Global thread index
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if i <= n
        # Branchless: multiply by active mask (same as CPU version)
        mask = active[i] ? one(eltype(outputs_p)) : zero(eltype(outputs_p))

        vi = vars_v[i]
        v2 = vi * vi

        # ZIP load model (identical to CPU)
        @inbounds outputs_p[i] = -mask * (params_pl[i] + params_ip[i] * vi + params_yp[i] * v2)
        @inbounds outputs_q[i] = -mask * (params_ql[i] + params_iq[i] * vi + params_yq[i] * v2)
    end

    return nothing
end

"""
    pq_residual_kernel!(outputs, params, vars, active, ::Type{<:CuArray})

GPU dispatch for PQ residual kernel using CUDA.

Launches CUDA kernel with appropriate thread/block configuration.
All arrays (params, vars, outputs, active) should already be on GPU.

This method is selected when the array type parameter is a CuArray subtype.
"""
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:CuArray}  # GPU dispatch
)
    n = length(outputs.p)

    # GPU launch configuration
    threads = 256
    blocks = cld(n, threads)

    # Launch CUDA kernel (async on default stream)
    CUDA.@cuda threads=threads blocks=blocks pq_residual_cuda_kernel!(
        outputs.p, outputs.q,
        params.pl, params.ql, params.ip, params.iq, params.yp, params.yq,
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
    calc_residuals!(ws::ResidualWorkspace{T, AT}, pq::StructArray{<:PQ{T}}) where {T, AT}

Calculate PQ load residual contributions using ZIP load model.

Implements the composite load model with constant MVA (P), constant current (I),
and constant impedance (Z) components:
- P_load = pl + ip*v + yp*v²
- Q_load = ql + iq*v + yq*v²

Residual convention: Positive = injection, so load consumption is negative.
Only active (in-service) loads contribute to residuals.

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `pq`: PQ load model vector with numerical fields

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `ws.active_mask` for out-of-service loads
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: `ws.var_vals.v` (voltage magnitudes)
- Writes: `ws.res_vals.p`, `ws.res_vals.q` (power residuals)
"""
@inline function calc_residuals!(
    ws::ResidualWorkspace{T, AT},
    pq::StructArray{<:PQ}
) where {T, AT}
    # Data adapter - extract arrays and call kernel
    # AT parameter dispatches to CPU (Vector) or GPU (CuArray) implementation
    pq_residual_kernel!(
        (p = ws.res_vals.p, q = ws.res_vals.q),
        (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq),
        (v = ws.var_vals.v,),
        ws.active_mask,
        AT
    )
end


### === Begin Tests === ###

@testitem "PQ Parsing" begin
    using PowerFlowData
    using Powerful
    using Powerful.Models
    using Powerful.PowerCore

    case = PowerFlowData.parse_network(
        joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    )

    pq = parse_model(PQ, case)
    @test pq isa PQVec
end

@testitem "PQ Residual Kernel (Pure Computation)" begin
    using Powerful.Models
    using Powerful.Models: pq_residual_kernel!

    # Test the pure kernel in isolation (no workspace, no DAE, no addressing)
    n = 5

    # Create test data
    p_res = zeros(n)
    q_res = zeros(n)

    # Load parameters (ZIP model coefficients)
    pl = [10.0, 20.0, 30.0, 40.0, 50.0]
    ql = [5.0, 10.0, 15.0, 20.0, 25.0]
    ip = [0.1, 0.2, 0.3, 0.4, 0.5]
    iq = [0.05, 0.1, 0.15, 0.2, 0.25]
    yp = [0.01, 0.02, 0.03, 0.04, 0.05]
    yq = [0.005, 0.01, 0.015, 0.02, 0.025]

    # Voltage values
    v = [1.0, 1.05, 0.98, 1.02, 0.95]

    # Active mask (load 3 is out of service)
    active = [true, true, false, true, true]

    # Call pure kernel
    pq_residual_kernel!(
        (p = p_res, q = q_res),
        (pl = pl, ql = ql, ip = ip, iq = iq, yp = yp, yq = yq),
        (v = v,),
        active,
        Vector{Float64}
    )

    # Verify results
    # Load 1: v=1.0, active=true
    @test p_res[1] ≈ -(pl[1] + ip[1] * v[1] + yp[1] * v[1]^2)
    @test q_res[1] ≈ -(ql[1] + iq[1] * v[1] + yq[1] * v[1]^2)

    # Load 3: v=0.98, active=false (should be zero)
    @test p_res[3] == 0.0
    @test q_res[3] == 0.0

    # Load 5: v=0.95, active=true
    @test p_res[5] ≈ -(pl[5] + ip[5] * v[5] + yp[5] * v[5]^2)
    @test q_res[5] ≈ -(ql[5] + iq[5] * v[5] + yq[5] * v[5]^2)

    # All loads consume power (negative residuals for active loads)
    @test all(p_res[active] .< 0)
    @test all(q_res[active] .< 0)
end

