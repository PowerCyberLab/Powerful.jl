export Branch, BranchVec
export calc_branch_power_flows, calc_branch_currents

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
struct Branch{Tv} <: AbstractModel
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
	Branch status: true = in service
	"""
	status::Bool

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
"""
    BranchVec{T}

Alias for `StructArray{Branch{T}}`, the struct-of-arrays container holding all
[`Branch`](@ref) instances of a system.
"""
const BranchVec{T} = StructArray{Branch{T}}

get_key_name(::Type{Branch}) = [:i, :j, :ckt]

# Base quantity traits - Branch inherits base voltage from both connected buses
base_quantity_role(::Type{Branch}, ::Type{BaseVoltage}) = InheritsBase()

# Connectivity declaration - single source of truth for topology AND ExtParams
const BRANCH_TERMINALS = (
    Terminal(:Bus, :i),  # From bus
    Terminal(:Bus, :j)   # To bus
)

get_terminals(::Type{Branch}) = BRANCH_TERMINALS

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

const BRANCH_observable_vars = Symbol[]  # No DAE variables to report (only computed observations)

# CurrentInjection formulation - rectangular coordinates
const BRANCH_VARS_CURRENT = [
	ModelVar(
		:v_r_i,
		AlgebVar(),
		:Bus,
		:v_r,
		:i,
		description = "Branch voltage real part at bus i",
		units = "pu",
	),
	ModelVar(
		:v_i_i,
		AlgebVar(),
		:Bus,
		:v_i,
		:i,
		description = "Branch voltage imaginary part at bus i",
		units = "pu",
	),
	ModelVar(
		:v_r_j,
		AlgebVar(),
		:Bus,
		:v_r,
		:j,
		description = "Branch voltage real part at bus j",
		units = "pu",
	),
	ModelVar(
		:v_i_j,
		AlgebVar(),
		:Bus,
		:v_i,
		:j,
		description = "Branch voltage imaginary part at bus j",
		units = "pu",
	),
]

# DC Power Flow formulation - angles only
const BRANCH_VARS_DC = [
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
]

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

# CurrentInjection formulation residuals
const BRANCH_RESIDUALS_CURRENT = [
	ModelResidual(
		:i_r_i,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:i_r,
		:i,
		description = "Real current injection at bus i",
	),
	ModelResidual(
		:i_i_i,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:i_i,
		:i,
		description = "Imaginary current injection at bus i",
	),
	ModelResidual(
		:i_r_j,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:i_r,
		:j,
		description = "Real current injection at bus j",
	),
	ModelResidual(
		:i_i_j,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:i_i,
		:j,
		description = "Imaginary current injection at bus j",
	),
]

# DC Power Flow formulation residuals - active power only
const BRANCH_RESIDUALS_DC = [
	ModelResidual(
		:p_i,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:p,
		:i,
		description = "Active power flow (DC approximation)",
	),
	ModelResidual(
		:p_j,
		AlgebRes(),
		SharedRes(),
		:Bus,
		:p,
		:j,
		description = "Active power flow (DC approximation)",
	),
]

# Branch parameters (already in system per-unit from PSS/E)
# All impedance and admittance values (r, x, b, g) are on system MVA base
const BRANCH_PARAMS = ModelParam[]

# Branch observations (computed from state and parameters)
const BRANCH_OBSERVATIONS = [
    ModelObservable(
        :p_from,
        ActivePower,
        units = "pu",
        description = "Active power flow from bus i (computed)",
        default_output = true,
    ),
    ModelObservable(
        :q_from,
        ReactivePower,
        units = "pu",
        description = "Reactive power flow from bus i (computed)",
        default_output = true,
    ),
    ModelObservable(
        :p_to,
        ActivePower,
        units = "pu",
        description = "Active power flow from bus j (computed)",
        default_output = true,
    ),
    ModelObservable(
        :q_to,
        ReactivePower,
        units = "pu",
        description = "Reactive power flow from bus j (computed)",
        default_output = true,
    ),
    ModelObservable(
        :p_loss,
        ActivePower,
        units = "pu",
        description = "Active power loss in branch (computed)",
        default_output = false,
    ),
    ModelObservable(
        :q_loss,
        ReactivePower,
        units = "pu",
        description = "Reactive power loss in branch (computed)",
        default_output = false,
    ),
]

# Branch has no device ratings
const BRANCH_RATINGS = ModelRating[]

# Branch references bus types for slack masking - uses BRANCH_TERMINALS
const BRANCH_EXT_PARAMS = [
    ExtParam(:i_ide, BRANCH_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type at terminal i"),
    ExtParam(:j_ide, BRANCH_TERMINALS[2], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type at terminal j")
]

"""
Get metadata for branch model
"""
function model_metadata(::Type{Branch})
	metadata = ModelMetadata(
		name = :Branch,
		vars = BRANCH_VARS,
		observable_vars = BRANCH_observable_vars,
		residuals = BRANCH_RESIDUALS,
		params = BRANCH_PARAMS,
		observations = BRANCH_OBSERVATIONS,
		ratings = BRANCH_RATINGS,
		ext_params = BRANCH_EXT_PARAMS,
		status_field = :status,
	)

	# Register CurrentInjection formulation overrides
	metadata.vars_by_formulation[CurrentInjection] = BRANCH_VARS_CURRENT
	metadata.residuals_by_formulation[CurrentInjection] = BRANCH_RESIDUALS_CURRENT

	# Register DC Power Flow formulation overrides
	metadata.vars_by_formulation[DCPowerFlow] = BRANCH_VARS_DC
	metadata.residuals_by_formulation[DCPowerFlow] = BRANCH_RESIDUALS_DC

	return metadata
end

# Declare formulation support
supports_formulation(::Type{Branch}, ::Type{PowerInjection}) = true
supports_formulation(::Type{Branch}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{Branch}, ::Type{DCPowerFlow}) = true

### === Model Format Adapters ===

supports_format(::Type{Branch}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Branch}, ::Type{PSSE}) = (:branches,)

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
		status = branches.st,
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


# ============================================================================
# Physics Functions (Reusable, CPU/GPU Compatible)
# ============================================================================

"""
    calc_branch_power_flows(v_i, theta_i, v_j, theta_j, r, x, b, gi, bi, gj, bj) -> NamedTuple

Calculate branch power flows using PI equivalent circuit model.

Pure physics function - works on CPU and GPU (device-compatible).
Called by residual kernels and report generation.

# Arguments
- `v_i, theta_i`: Voltage magnitude and angle at "from" bus (bus i)
- `v_j, theta_j`: Voltage magnitude and angle at "to" bus (bus j)
- `r, x`: Series resistance and reactance (pu)
- `b`: Total line charging susceptance (pu)
- `gi, bi`: Shunt conductance and susceptance at bus i
- `gj, bj`: Shunt conductance and susceptance at bus j

# Returns
NamedTuple with:
- `p_i`: Active power flow from bus i (pu, positive = flowing out)
- `q_i`: Reactive power flow from bus i (pu)
- `p_j`: Active power flow from bus j (pu)
- `q_j`: Reactive power flow from bus j (pu)

# Physics
PI model equations:
- Series admittance: y = 1/(r + jx) = g + jb
- Shunt admittance split between ends: b/2 at each end
- Power flow: P + jQ = V·conj(I) where I = y·(V_i - V_j) + y_shunt·V
"""
@inline function calc_branch_power_flows(v_i, theta_i, v_j, theta_j, r, x, b_charging, gi, bi, gj, bj)
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

    # Power flow from i (positive = flowing out of bus i)
    p_i = g_shunt_i * v_i * v_i +
          g_series * v_i * v_i -
          v_i * v_j * (g_series * cos_theta + b_series * sin_theta)

    q_i = -b_shunt_i * v_i * v_i -
           b_series * v_i * v_i -
           v_i * v_j * (g_series * sin_theta - b_series * cos_theta)

    # Power flow from j (angle reversed: cos(-θ) = cos(θ), sin(-θ) = -sin(θ))
    p_j = g_shunt_j * v_j * v_j +
          g_series * v_j * v_j -
          v_j * v_i * (g_series * cos_theta - b_series * sin_theta)

    q_j = -b_shunt_j * v_j * v_j -
           b_series * v_j * v_j -
           v_j * v_i * (g_series * (-sin_theta) - b_series * cos_theta)

    return (p_i=p_i, q_i=q_i, p_j=p_j, q_j=q_j)
end

"""
    calc_branch_currents(v_r_i, v_i_i, v_r_j, v_i_j, r, x, b_charging, gi, bi, gj, bj) -> NamedTuple

Calculate branch current injections using PI equivalent circuit model with rectangular voltages.

Pure physics function - works on CPU and GPU (device-compatible).
Called by current injection residual kernels and reports.

# Arguments
- `v_r_i, v_i_i`: Real and imaginary voltage at "from" bus (bus i) [pu]
- `v_r_j, v_i_j`: Real and imaginary voltage at "to" bus (bus j) [pu]
- `r, x`: Series resistance and reactance [pu]
- `b_charging`: Total line charging susceptance [pu]
- `gi, bi`: Shunt conductance and susceptance at bus i [pu]
- `gj, bj`: Shunt conductance and susceptance at bus j [pu]

# Returns
NamedTuple with:
- `i_r_i, i_i_i`: Real and imaginary current injection at bus i [pu]
- `i_r_j, i_i_j`: Real and imaginary current injection at bus j [pu]

# Physics
PI model with admittance matrix Y:
```
I = Y * V
where V = [V_i; V_j], I = [I_i; I_j]

Y = [[Y_ii, Y_ij],
     [Y_ji, Y_jj]]

Y_series = 1/(r + jx) = g_series + j*b_series
Y_ii = Y_series + Y_shunt_i
Y_ij = -Y_series
Y_ji = -Y_series
Y_jj = Y_series + Y_shunt_j
```

Current injection I_i = Y_ii * V_i + Y_ij * V_j (positive = flowing into bus from branch)

# Note
Sign convention: Current injection is positive flowing INTO the bus.
This matches the network convention where I = Y*V naturally.
"""
@inline function calc_branch_currents(
    v_r_i, v_i_i, v_r_j, v_i_j,
    r, x, b_charging, gi, bi, gj, bj
)
    # Series admittance: Y_series = 1/(r + jx) = (r - jx)/(r² + x²)
    z_sq = r * r + x * x
    g_series = r / z_sq
    b_series = -x / z_sq

    # Shunt admittance at each end (half charging + line shunt)
    g_shunt_i = gi
    b_shunt_i = bi + b_charging / 2
    g_shunt_j = gj
    b_shunt_j = bj + b_charging / 2

    # Admittance matrix elements
    # Y_ii = Y_series + Y_shunt_i = (g_series + g_shunt_i) + j*(b_series + b_shunt_i)
    y_ii_r = g_series + g_shunt_i
    y_ii_i = b_series + b_shunt_i

    # Y_ij = -Y_series = -g_series - j*b_series
    y_ij_r = -g_series
    y_ij_i = -b_series

    # Y_ji = -Y_series (symmetric network)
    y_ji_r = -g_series
    y_ji_i = -b_series

    # Y_jj = Y_series + Y_shunt_j
    y_jj_r = g_series + g_shunt_j
    y_jj_i = b_series + b_shunt_j

    # Current injection at bus i: I_i = Y_ii * V_i + Y_ij * V_j
    # Complex multiplication: (a + j*b) * (c + j*d) = (ac - bd) + j*(ad + bc)
    # I_i = (y_ii_r + j*y_ii_i) * (v_r_i + j*v_i_i) + (y_ij_r + j*y_ij_i) * (v_r_j + j*v_i_j)
    i_r_i = y_ii_r * v_r_i - y_ii_i * v_i_i + y_ij_r * v_r_j - y_ij_i * v_i_j
    i_i_i = y_ii_r * v_i_i + y_ii_i * v_r_i + y_ij_r * v_i_j + y_ij_i * v_r_j

    # Current injection at bus j: I_j = Y_ji * V_i + Y_jj * V_j
    i_r_j = y_ji_r * v_r_i - y_ji_i * v_i_i + y_jj_r * v_r_j - y_jj_i * v_i_j
    i_i_j = y_ji_r * v_i_i + y_ji_i * v_r_i + y_jj_r * v_i_j + y_jj_i * v_r_j

    return (i_r_i=i_r_i, i_i_i=i_i_i, i_r_j=i_r_j, i_i_j=i_i_j)
end

# ============================================================================
# Reporting Interface - Vectorized Observation Computation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:Branch}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all branches in the solution (vectorized).

# Supported Variables
- `:p_from`: Active power flow from bus i (pu)
- `:q_from`: Reactive power flow from bus i (pu)
- `:p_to`: Active power flow from bus j (pu)
- `:q_to`: Reactive power flow from bus j (pu)
- `:p_loss`: Active power loss = p_from + p_to (pu)
- `:q_loss`: Reactive power loss = q_from + q_to (pu)

# Implementation Notes
Reuses `calc_branch_power_flows` physics function.
"""
function compute_observations(sol::PowerFlowSolution, ::Type{<:Branch}, var::Symbol)
    branches = sol.prob.sys.models.Branch
    var in (:p_from, :q_from, :p_to, :q_to, :p_loss, :q_loss) ||
        error("Unknown Branch output variable: $var")

    # Terminal bus voltages (any formulation); an out-of-service branch carries no flow
    v_i, theta_i = terminal_voltages(sol, :Branch, :i)
    v_j, theta_j = terminal_voltages(sol, :Branch, :j)
    n = length(branches)

    if sol.prob.formulation isa DCPowerFlow
        # Lossless DC flow, as in the DC kernel
        p_from = [branches.status[k] * dc_susceptance(branches.x[k]) * (theta_i[k] - theta_j[k]) for k in 1:n]
        var === :p_from && return p_from
        var === :p_to && return -p_from
        var === :p_loss && return zeros(n)
        error("Branch observation :$var is not defined in the DCPowerFlow formulation")
    end

    # Preallocate result arrays
    p_i = Vector{Float64}(undef, n)
    q_i = Vector{Float64}(undef, n)
    p_j = Vector{Float64}(undef, n)
    q_j = Vector{Float64}(undef, n)

    # Compute power flows using reusable physics kernel
    @inbounds for k in 1:n
        flows = calc_branch_power_flows(
            v_i[k], theta_i[k], v_j[k], theta_j[k],
            branches.r[k], branches.x[k], branches.b[k],
            branches.gi[k], branches.bi[k], branches.gj[k], branches.bj[k]
        )
        on = branches.status[k]
        p_i[k] = on * flows.p_i
        q_i[k] = on * flows.q_i
        p_j[k] = on * flows.p_j
        q_j[k] = on * flows.q_j
    end

    # Return requested variable
    var === :p_from && return p_i
    var === :q_from && return q_i
    var === :p_to && return p_j
    var === :q_to && return q_j
    var === :p_loss && return @. p_i + p_j
    return @. q_i + q_j  # :q_loss
end

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
- P_i = g_i*v_i² + g*v_i² - v_i*v_j*(g*cos(θ_ij) + b*sin(θ_ij))
- Q_i = -b_i*v_i² - b*v_i² - v_i*v_j*(g*sin(θ_ij) - b*cos(θ_ij))
- P_j = g_j*v_j² + g*v_j² - v_i*v_j*(g*cos(θ_ji) + b*sin(θ_ji))
- Q_j = -b_j*v_j² - b*v_j² - v_i*v_j*(g*sin(θ_ji) - b*cos(θ_ji))

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
        # Branchless masking - only active status (bus type filtering handled by scatter :set mode)
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

        # Calculate power flows using extracted physics function
        flows = calc_branch_power_flows(v_i, theta_i, v_j, theta_j, r, x, b_charging, gi, bi, gj, bj)

        # Sign convention: power flowing INTO bus is positive (injection)
        # flows.p_i/q_i = power OUT of bus i, so negate for injection
        outputs.p_i[i] = mask * (-flows.p_i)
        outputs.q_i[i] = mask * (-flows.q_i)
        outputs.p_j[i] = mask * (-flows.p_j)
        outputs.q_j[i] = mask * (-flows.q_j)
    end

    nothing
end

# ============================================================================
# Model-Specific calc_residuals! (Data Adapter)
# ============================================================================

"""
    calc_residuals!(ws::ResidualBuffers{T, AT}, branch::StructArray{<:Branch}) where {T, AT}

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
    ws::ResidualBuffers{T, AT},
    branch::StructArray{<:Branch},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
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
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )

end

# CurrentInjection formulation - uses rectangular voltages and current balance
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    branch::StructArray{<:Branch},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool},
    force_accumulate::Bool=false
) where {T, AT}
    # Data adapter - extract arrays and call current injection kernel
    branch_current_kernel!(
        # Outputs: current injection residuals at both ends
        (i_r_i = ws.res_vals.i_r_i, i_i_i = ws.res_vals.i_i_i,
         i_r_j = ws.res_vals.i_r_j, i_i_j = ws.res_vals.i_i_j),

        # Parameters: line impedance and shunt elements
        (r = branch.r, x = branch.x, b = branch.b,
         gi = branch.gi, bi = branch.bi,
         gj = branch.gj, bj = branch.bj),

        # External parameters: bus types for slack masking
        _external,

        # Variables: rectangular bus voltages at both ends
        (v_r_i = ws.var_vals.v_r_i, v_i_i = ws.var_vals.v_i_i,
         v_r_j = ws.var_vals.v_r_j, v_i_j = ws.var_vals.v_i_j),

        # Metadata: branch status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

# ============================================================================
# Current Injection Kernel (CPU)
# ============================================================================

"""
    branch_current_kernel!(outputs, params, vars, active, ::Type{<:Vector})

Compute branch current injections for CurrentInjection formulation.

Uses rectangular voltage coordinates (V_r, V_i) and computes current balance
at each bus. The current injection from branch into bus i is I_i = Y_ii * V_i + Y_ij * V_j.

Sign convention: Current INJECTED INTO bus is positive (consistent with power injection).
"""
@inline function branch_current_kernel!(
    outputs::NamedTuple{(:i_r_i, :i_i_i, :i_r_j, :i_i_j)},
    params::NamedTuple{(:r, :x, :b, :gi, :bi, :gj, :bj)},
    external::NamedTuple,
    vars::NamedTuple{(:v_r_i, :v_i_i, :v_r_j, :v_i_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.i_r_i)
        # Branchless masking - don't inject currents to slack buses (ide == 3)
        mask_i = active[k] & (external.i_ide[k] != BUS_SLACK)
        mask_j = active[k] & (external.j_ide[k] != BUS_SLACK)

        # Get currents from helper function
        currents = calc_branch_currents(
            vars.v_r_i[k], vars.v_i_i[k], vars.v_r_j[k], vars.v_i_j[k],
            params.r[k], params.x[k], params.b[k],
            params.gi[k], params.bi[k], params.gj[k], params.bj[k]
        )

        # Sign convention: Y*V gives current INTO element, negate for current INTO bus
        outputs.i_r_i[k] = mask_i * (-currents.i_r_i)
        outputs.i_i_i[k] = mask_i * (-currents.i_i_i)
        outputs.i_r_j[k] = mask_j * (-currents.i_r_j)
        outputs.i_i_j[k] = mask_j * (-currents.i_i_j)
    end

    nothing
end

# ============================================================================
# DC Power Flow Kernel (CPU)
# ============================================================================

"""
    branch_dc_kernel!(outputs, params, vars, active, ::Type{<:Vector})

Compute branch active power flows for DC Power Flow formulation.

Uses linearized DC power flow approximation:
- All voltage magnitudes assumed at 1.0 pu
- Resistance neglected (r << x)
- Small angle approximation: sin(theta_i - theta_j) ≈ theta_i - theta_j
- P_ij = (theta_i - theta_j) / x

Sign convention: Power FLOWING OUT of bus is positive (ANDES convention).

# Arguments
- `outputs::NamedTuple{(:p_i, :p_j)}`: Active power flow residuals
- `params::NamedTuple{(:x,)}`: Branch reactance only (resistance ignored in DC)
- `external::NamedTuple`: External parameters (bus types for masking)
- `vars::NamedTuple{(:theta_i, :theta_j)}`: Voltage angles at both ends
- `active::AbstractVector{Bool}`: Branch in-service status
- `::Type{<:Vector}`: CPU dispatch marker
"""
@inline function branch_dc_kernel!(
    outputs::NamedTuple{(:p_i, :p_j)},
    params::NamedTuple{(:x,)},
    vars::NamedTuple{(:theta_i, :theta_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for k in eachindex(outputs.p_i)
        # Branchless masking - only active status
        # Bus type filtering handled by scatter mode (:add vs :set)
        mask = active[k]

        # DC power flow: P = (theta_i - theta_j) / x
        # Note: Resistance neglected, voltage magnitudes = 1.0 pu
        x = params.x[k]
        theta_diff = vars.theta_i[k] - vars.theta_j[k]

        # Power flow from i to j (positive = power leaving bus i)
        p_flow = dc_susceptance(x) * theta_diff

        # Sign convention: power flowing INTO bus is positive (injection)
        # p_flow = power from i to j, so at bus i it leaves (negative), at bus j it enters (positive)
        outputs.p_i[k] = mask * (-p_flow)
        outputs.p_j[k] = mask * p_flow
    end

    nothing
end

# DC Power Flow calc_residuals! dispatcher
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    branch::StructArray{<:Branch},
    _external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and call DC kernel
    branch_dc_kernel!(
        # Outputs: active power flows only (no reactive power in DC)
        (p_i = ws.res_vals.p_i, p_j = ws.res_vals.p_j),

        # Parameters: reactance only (resistance ignored in DC approximation)
        (x = branch.x,),

        # Variables: voltage angles only (magnitudes assumed 1.0 pu)
        (theta_i = ws.var_vals.theta_i, theta_j = ws.var_vals.theta_j),

        # Metadata: branch status
        active_mask,

        # Dispatch parameter: CPU vs GPU
        AT
    )
end

