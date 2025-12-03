export Branch, BranchVec, BranchNum, BranchNumVec

abstract type AbstractBranch{T} <: AbstractModel end

"""
	$TYPEDEF

In PSS/E and similar tools, the basic transmission line model is an Equivalent
Pi model connected between network buses.

Data for shunt equipment units, such as reactors, which are connected to and
switched with the line, are entered in the same data record.

!!! note "Transformers" Branches to be modeled as transformers are not specified
	in this data category; # rather, they are specified in the

# Fields
$TYPEDFIELDS
"""
struct Branch{Tv} <: AbstractBranch{Tv}
	"Branch `from bus` number"
	i::Int32

	"Branch `to bus` number. Negative value designates it as the metered end"
	j::Int32

	"""
	One- or two-character branch circuit identifier. Default is '1'
	"""
	ckt::String3

	"""
	Branch resistance in pu
	"""
	r::Tv

	"""
	Branch reactance in pu
	Must be nonzero
	"""
	x::Tv

	"""
	Total branch charging susceptance in pu
	"""
	b::Tv

	"""
	First loading rating in MVA.
	Zero means branch excluded from loading examination
	"""
	rate_a::Tv

	"""
	Second loading rating in MVA
	"""
	rate_b::Tv

	"""
	Third loading rating in MVA
	"""
	rate_c::Tv

	"""
	Real part of line shunt admittance at bus i end in pu
	"""
	gi::Tv

	"""
	Imaginary part of line shunt admittance at bus i end in pu
	"""
	bi::Tv

	"""
	Real part of line shunt admittance at bus j end in pu
	"""
	gj::Tv

	"""
	Imaginary part of line shunt admittance at bus j end in pu
	"""
	bj::Tv

	"""
	Branch status
	1=in-service, 0=out-of-service
	"""
	st::Bool

	"""
	Metering end indicator
	"""
	met::Int8

	"""
	Line length in user-selected units
	"""
	len::Tv

	"""
	First owner number
	"""
	o1::Int16

	"""
	Fraction of total ownership for first owner
	"""
	f1::Tv

	"""
	Second owner number
	"""
	o2::Union{Missing, Int16}

	"""
	Fraction of total ownership for second owner
	"""
	f2::Union{Missing, Tv}

	"""
	Third owner number
	"""
	o3::Union{Missing, Int16}

	"""
	Fraction of total ownership for third owner
	"""
	f3::Union{Missing, Tv}

	"""
	Fourth owner number
	"""
	o4::Union{Missing, Int16}

	"""
	Fraction of total ownership for fourth owner
	"""
	f4::Union{Missing, Tv}
end

@register_model(Branch)

# Type aliases for StructArray-based vectorized types
const BranchVec{T} = StructArray{Branch{T}}

get_key_name(::Type{Branch}) = (:i, :j, :ckt)

const BRANCH_VARS = [
	ModelVar(
		:theta_i,
		AlgebVar(),
		:Bus,
		:theta,
		:i,
		description = "Branch voltage angle at bus i",
		units = "rad",
	),
	ModelVar(
		:theta_j,
		AlgebVar(),
		:Bus,
		:theta,
		:j,
		description = "Branch voltage angle at bus j",
		units = "rad",
	),
	ModelVar(
		:v_i,
		AlgebVar(),
		:Bus,
		:v,
		:i,
		description = "Branch voltage magnitude at bus i",
		units = "pu",
	),
	ModelVar(
		:v_j,
		AlgebVar(),
		:Bus,
		:v,
		:j,
		description = "Branch voltage magnitude at bus j",
		units = "pu",
	),
]

const BRANCH_OUTPUT_VARS = []

const BRANCH_RESIDUALS = [
	ModelResidual(
		:p_i,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:p,
		:i,
		description = "Active power flow",
	),
	ModelResidual(
		:p_j,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:p,
		:j,
		description = "Active power flow",
	),
	ModelResidual(
		:q_i,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:q,
		:i,
		description = "Reactive power flow",
	),
	ModelResidual(
		:q_j,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:q,
		:j,
		description = "Reactive power flow",
	),
]

"""
Get metadata for branch model
"""
function model_metadata(::Type{Branch}; layout::LayoutStrategy = ContiguousVariables())
	ModelMetadata(
		name = :Branch,
		vars = BRANCH_VARS,
		output_vars = BRANCH_OUTPUT_VARS,
		residuals = BRANCH_RESIDUALS,
		layout = layout,
	)
end

### === Model Format Adapters ===

supports_format(::Type{Branch}, ::Type{PSSE}) = FormatSupport{PSSE}()

function parse_model(::Type{Branch}, raw, ::FormatSupport{PSSE})
	branches = raw.branches

	# Create StructArray directly from named tuple of vectors
	return StructArray{Branch{Float64}}((
		i = branches.i,
		j = branches.j,
		ckt = branches.ckt,
		r = branches.r,
		x = branches.x,
		b = branches.b,
		rate_a = branches.rate_a,
		rate_b = branches.rate_b,
		rate_c = branches.rate_c,
		gi = branches.gi,
		bi = branches.bi,
		gj = branches.gj,
		bj = branches.bj,
		st = branches.st,
		met = branches.met,
		len = branches.len,
		o1 = branches.o1,
		f1 = branches.f1,
		o2 = branches.o2,
		f2 = branches.f2,
		o3 = branches.o3,
		f3 = branches.f3,
		o4 = branches.o4,
		f4 = branches.f4,
	))
end


### === Residual Calculation === ###

using ..PowerCore: ResidualWorkspace

# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    branch_residual_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for branch (transmission line) power flow residuals.

Implements the PI equivalent circuit model for AC transmission lines:
- Series impedance: z = r + jx
- Shunt charging: b (split equally at both ends for PI model)
- Shunt conductance: g (usually zero)

Power flow equations at each end (i → j):
- P_i = g_i * v_i² + v_i * v_j * (g*cos(θ_ij) + b*sin(θ_ij))
- Q_i = -b_i * v_i² - v_i * v_j * (g*sin(θ_ij) - b*cos(θ_ij))
- P_j = g_j * v_j² + v_i * v_j * (g*cos(θ_ji) + b*sin(θ_ji))
- Q_j = -b_j * v_j² - v_i * v_j * (g*sin(θ_ji) - b*cos(θ_ji))

Where θ_ij = θ_i - θ_j, and g + jb = 1/(r + jx) (series admittance).

# Arguments
- `outputs::NamedTuple{(:p_i, :p_j, :q_i, :q_j)}`: Power flow residuals
- `params::NamedTuple{(:r, :x, :b, :gi, :bi, :gj, :bj)}`: Branch parameters
- `vars::NamedTuple{(:theta_i, :theta_j, :v_i, :v_j)}`: Bus voltages
- `active::AbstractVector{Bool}`: Branch in-service status
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)
"""
@inline function branch_residual_kernel!(
    outputs::NamedTuple{(:p_i, :p_j, :q_i, :q_j)},
    params::NamedTuple{(:r, :x, :b, :gi, :bi, :gj, :bj)},
    vars::NamedTuple{(:theta_i, :theta_j, :v_i, :v_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p_i)
        mask = active[i]

        # Extract parameters
        r = params.r[i]
        x = params.x[i]
        b_charging = params.b[i]
        gi = params.gi[i]
        bi = params.bi[i]
        gj = params.gj[i]
        bj = params.bj[i]

        # Extract voltages
        theta_i = vars.theta_i[i]
        theta_j = vars.theta_j[i]
        v_i = vars.v_i[i]
        v_j = vars.v_j[i]

        # Series admittance: y = 1/(r + jx) = (r - jx)/(r² + x²)
        z_sq = r * r + x * x
        g_series = r / z_sq
        b_series = -x / z_sq

        # Angle difference
        theta_ij = theta_i - theta_j
        cos_theta = cos(theta_ij)
        sin_theta = sin(theta_ij)

        # Shunt admittance at each end (half charging + line shunt)
        g_shunt_i = gi
        b_shunt_i = bi + b_charging / 2
        g_shunt_j = gj
        b_shunt_j = bj + b_charging / 2

        # Power flow from i to j (positive = flow out of bus i)
        # P_i = shunt loss + series flow
        p_i_flow = g_shunt_i * v_i * v_i +
                   v_i * v_j * (g_series * cos_theta + b_series * sin_theta)

        # Q_i = shunt generation - series flow
        q_i_flow = -b_shunt_i * v_i * v_i -
                    v_i * v_j * (g_series * sin_theta - b_series * cos_theta)

        # Power flow from j to i (angle reversed)
        cos_theta_ji = cos_theta  # cos(-θ) = cos(θ)
        sin_theta_ji = -sin_theta  # sin(-θ) = -sin(θ)

        p_j_flow = g_shunt_j * v_j * v_j +
                   v_j * v_i * (g_series * cos_theta_ji + b_series * sin_theta_ji)

        q_j_flow = -b_shunt_j * v_j * v_j -
                    v_j * v_i * (g_series * sin_theta_ji - b_series * cos_theta_ji)

        # Branches inject power OUT of buses (negative consumption)
        # Apply mask for out-of-service branches
        outputs.p_i[i] = -mask * p_i_flow
        outputs.q_i[i] = -mask * q_i_flow
        outputs.p_j[i] = -mask * p_j_flow
        outputs.q_j[i] = -mask * q_j_flow
    end

    nothing
end

# ============================================================================
# GPU Kernel Implementation
# ============================================================================

"""
CUDA kernel for branch power flow residuals - device function executed on GPU.
"""
function branch_residual_cuda_kernel!(
    outputs_p_i, outputs_p_j, outputs_q_i, outputs_q_j,
    params_r, params_x, params_b, params_gi, params_bi, params_gj, params_bj,
    vars_theta_i, vars_theta_j, vars_v_i, vars_v_j,
    active,
    n
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if i <= n
        mask = active[i] ? one(eltype(outputs_p_i)) : zero(eltype(outputs_p_i))

        @inbounds begin
            # Extract parameters
            r = params_r[i]
            x = params_x[i]
            b_charging = params_b[i]
            gi = params_gi[i]
            bi = params_bi[i]
            gj = params_gj[i]
            bj = params_bj[i]

            # Extract voltages
            theta_i = vars_theta_i[i]
            theta_j = vars_theta_j[i]
            v_i = vars_v_i[i]
            v_j = vars_v_j[i]

            # Series admittance
            z_sq = r * r + x * x
            g_series = r / z_sq
            b_series = -x / z_sq

            # Angle difference
            theta_ij = theta_i - theta_j
            cos_theta = CUDA.cos(theta_ij)
            sin_theta = CUDA.sin(theta_ij)

            # Shunt admittance
            g_shunt_i = gi
            b_shunt_i = bi + b_charging / 2
            g_shunt_j = gj
            b_shunt_j = bj + b_charging / 2

            # Power flow i → j
            p_i_flow = g_shunt_i * v_i * v_i +
                       v_i * v_j * (g_series * cos_theta + b_series * sin_theta)
            q_i_flow = -b_shunt_i * v_i * v_i -
                        v_i * v_j * (g_series * sin_theta - b_series * cos_theta)

            # Power flow j → i
            cos_theta_ji = cos_theta
            sin_theta_ji = -sin_theta

            p_j_flow = g_shunt_j * v_j * v_j +
                       v_j * v_i * (g_series * cos_theta_ji + b_series * sin_theta_ji)
            q_j_flow = -b_shunt_j * v_j * v_j -
                        v_j * v_i * (g_series * sin_theta_ji - b_series * cos_theta_ji)

            # Apply mask and store (negative for injection)
            outputs_p_i[i] = -mask * p_i_flow
            outputs_q_i[i] = -mask * q_i_flow
            outputs_p_j[i] = -mask * p_j_flow
            outputs_q_j[i] = -mask * q_j_flow
        end
    end

    return nothing
end

"""
    branch_residual_kernel!(outputs, params, vars, active, ::Type{<:CuArray})

GPU dispatch for branch residual kernel using CUDA.
"""
@inline function branch_residual_kernel!(
    outputs::NamedTuple{(:p_i, :p_j, :q_i, :q_j)},
    params::NamedTuple{(:r, :x, :b, :gi, :bi, :gj, :bj)},
    vars::NamedTuple{(:theta_i, :theta_j, :v_i, :v_j)},
    active::AbstractVector{Bool},
    ::Type{<:CuArray}  # GPU dispatch
)
    n = length(outputs.p_i)

    # GPU launch configuration
    threads = 256
    blocks = cld(n, threads)

    # Launch CUDA kernel (async on default stream)
    CUDA.@cuda threads=threads blocks=blocks branch_residual_cuda_kernel!(
        outputs.p_i, outputs.p_j, outputs.q_i, outputs.q_j,
        params.r, params.x, params.b, params.gi, params.bi, params.gj, params.bj,
        vars.theta_i, vars.theta_j, vars.v_i, vars.v_j,
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
    calc_residuals!(ws::ResidualWorkspace{T, AT}, branch::StructArray{<:Branch{T}}) where {T, AT}

Calculate branch (transmission line) power flow residuals.

Implements AC power flow through transmission lines using the PI equivalent circuit.
Each branch contributes power injections at both the "from" bus (i) and "to" bus (j).

**Automatically dispatches to CPU or GPU implementation** based on the workspace
array type parameter `AT`.

# Physics

The PI model represents a transmission line as:
- Series impedance z = r + jx
- Shunt charging susceptance b (split equally at both ends)
- Optional shunt conductance g at each end (line losses)

Power flows are calculated using voltage magnitudes and angle differences,
then injected into the connected buses' power balance equations.

# Arguments
- `ws`: Pre-allocated workspace with variable values and residual buffers
- `branch`: Branch model vector with line parameters

# Notes
- Zero allocation - operates on pre-allocated buffers
- Respects `ws.active_mask` for out-of-service branches
- CPU: Uses `@inbounds` and `@simd` for vectorization
- GPU: Launches CUDA kernel with 256 threads per block
- Reads: `ws.var_vals.{theta_i, theta_j, v_i, v_j}` (bus voltages)
- Writes: `ws.res_vals.{p_i, p_j, q_i, q_j}` (power injections)
"""
@inline function calc_residuals!(
    ws::ResidualWorkspace{T, AT},
    branch::StructArray{<:Branch}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    branch_residual_kernel!(
        # Outputs: power flow residuals at both ends
        (p_i = ws.res_vals.p_i, p_j = ws.res_vals.p_j,
         q_i = ws.res_vals.q_i, q_j = ws.res_vals.q_j),

        # Parameters: line impedance and shunt elements
        (r = branch.r, x = branch.x, b = branch.b,
         gi = branch.gi, bi = branch.bi,
         gj = branch.gj, bj = branch.bj),

        # Variables: bus voltages at both ends
        (theta_i = ws.var_vals.theta_i, theta_j = ws.var_vals.theta_j,
         v_i = ws.var_vals.v_i, v_j = ws.var_vals.v_j),

        # Metadata: branch status
        ws.active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end


### === Begin Tests === ###

@testitem "Branch" begin
	using PowerFlowData
	using Powerful
	using Powerful.Models
	using Powerful.PowerCore

	case = PowerFlowData.parse_network(joinpath(pkgdir(Powerful), "cases", "ieee14.raw"))
	branches = parse_model(Branch, case)
	@test branches isa BranchVec
end
