export Transformer2W, Transformer2WVec
export Transformer3W, Transformer3WVec
export transformer_residual_kernel!, transformer_residual_kernel_current_injection!, transformer_current_kernel!

"""
	$TYPEDEF

This struct represents a two-winding transformer, which is a subset of the data required for
three-winding transformers.

# Fields

$TYPEDFIELDS
"""
Base.@kwdef struct Transformer2W{Tv} <: AbstractModel
	"Bus number of the first winding bus"
	i::Int32

	"Bus number of the second winding bus"
	j::Int32

	"One- or two-character uppercase nonblank alphanumeric transformer circuit identifier"
	ckt::String3

	"Winding data I/O code (1: pu, 2: kV)"
	cw::Int8

	"Impedance data I/O code (1: system base, 2: specified MVA, 3: load loss)"
	cz::Int8

	"Magnetizing admittance I/O code (1: pu, 2: physical units)"
	cm::Int8

	"Magnetizing conductance (cm=1) or no load loss in watts (cm=2)"
	mag1::Tv

	"Magnetizing susceptance (cm=1) or exciting current in pu (cm=2)"
	mag2::Tv

	"Nonmetered end code (1: winding 1, 2: winding 2)"
	nmetr::Int8

	"Transformer name (12 chars max)"
	name::String15

	"Initial transformer status (true = in service)"
	status::Bool

	"First owner number"
	o1::Int16

	"First owner fraction"
	f1::Tv

	"Second owner number"
	o2::Union{Missing, Int16}

	"Second owner fraction"
	f2::Union{Missing, Tv}

	"Third owner number"
	o3::Union{Missing, Int16}

	"Third owner fraction"
	f3::Union{Missing, Tv}

	"Fourth owner number"
	o4::Union{Missing, Int16}

	"Fourth owner fraction"
	f4::Union{Missing, Tv}

	"Vector group identifier"
	vecgrp::Union{Missing, String15}

	"Resistance between buses i and j"
	r1_2::Tv

	"Reactance between buses i and j"
	x1_2::Tv

	"Base MVA for winding 1-2 impedance"
	sbase1_2::Tv

	"Winding 1 off-nominal turns ratio or voltage"
	windv1::Tv

	"Winding 1 nominal voltage in kV"
	nomv1::Tv

	"Winding 1 phase shift angle in degrees"
	ang1::Tv

	"First rating in MVA"
	rata1::Tv

	"Second rating in MVA"
	ratb1::Tv

	"Third rating in MVA"
	ratc1::Tv

	"Control mode for winding 1"
	cod1::Int8

	"Bus number controlled by winding 1"
	cont1::Int32

	"Upper limit on winding 1 adjustment"
	rma1::Tv

	"Lower limit on winding 1 adjustment"
	rmi1::Tv

	"Upper voltage limit for controlled bus"
	vma1::Tv

	"Lower voltage limit for controlled bus"
	vmi1::Tv

	"Number of tap positions"
	ntp1::Int16

	"Impedance correction table number"
	tab1::Int64

	"Load drop compensation resistance"
	cr1::Tv

	"Load drop compensation reactance"
	cx1::Tv

	"Winding connection angle"
	cnxa1::Union{Missing, Tv}

	"Winding 2 off-nominal turns ratio or voltage"
	windv2::Tv

	"Winding 2 nominal voltage in kV"
	nomv2::Tv
end

# ============================================================================
# Three-Winding Transformer
# ============================================================================

"""
	$TYPEDEF

This struct represents a three-winding transformer as parsed from PSS/E data.
It is a transient data structure that gets converted to three 2-winding
transformers with a virtual bus in a star configuration.

Note: This struct does not define residuals - it's only used for parsing
and conversion. Power flow computations use the equivalent 2W transformers.

# Fields

$TYPEDFIELDS
"""
struct Transformer3W{Tv} <: AbstractModel
	"Bus number of the first winding bus"
	i::Int32

	"Bus number of the second winding bus"
	j::Int32

	"Bus number of the third winding bus"
	k::Int32

	"One- or two-character uppercase nonblank alphanumeric transformer circuit identifier"
	ckt::String3

	"Winding data I/O code (1: pu, 2: kV)"
	cw::Int8

	"Impedance data I/O code (1: system base, 2: specified MVA, 3: load loss)"
	cz::Int8

	"Magnetizing admittance I/O code (1: pu, 2: physical units)"
	cm::Int8

	"Magnetizing conductance (cm=1) or no load loss in watts (cm=2)"
	mag1::Tv

	"Magnetizing susceptance (cm=1) or exciting current in pu (cm=2)"
	mag2::Tv

	"Nonmetered end code (1: winding 1, 2: winding 2, 3: winding 3)"
	nmetr::Int8

	"Transformer name (12 chars max)"
	name::String15

	"Initial transformer status (0 = out of service, 1 = in service, 2 = winding 2 out of service, 3 = winding 3 out of service, 4 = only winding 1 in service)"
	stat::Int8

	"First owner number"
	o1::Int16

	"First owner fraction"
	f1::Tv

	"Second owner number"
	o2::Union{Missing, Int16}

	"Second owner fraction"
	f2::Union{Missing, Tv}

	"Third owner number"
	o3::Union{Missing, Int16}

	"Third owner fraction"
	f3::Union{Missing, Tv}

	"Fourth owner number"
	o4::Union{Missing, Int16}

	"Fourth owner fraction"
	f4::Union{Missing, Tv}

	"Vector group identifier"
	vecgrp::Union{Missing, String15}

	# Winding 1-2 (i-j) impedance
	"Resistance between buses i and j"
	r1_2::Tv

	"Reactance between buses i and j"
	x1_2::Tv

	"Base MVA for winding 1-2 impedance"
	sbase1_2::Tv

	# Winding 2-3 (j-k) impedance
	"Resistance between buses j and k"
	r2_3::Tv

	"Reactance between buses j and k"
	x2_3::Tv

	"Base MVA for winding 2-3 impedance"
	sbase2_3::Tv

	# Winding 3-1 (k-i) impedance
	"Resistance between buses k and i"
	r3_1::Tv

	"Reactance between buses k and i"
	x3_1::Tv

	"Base MVA for winding 3-1 impedance"
	sbase3_1::Tv

	"Voltage magnitude at star point"
	vmstar::Tv

	"Angle at star point"
	anstar::Tv

	# Winding 1 parameters
	"Winding 1 off-nominal turns ratio or voltage"
	windv1::Tv

	"Winding 1 nominal voltage in kV"
	nomv1::Tv

	"Winding 1 phase shift angle in degrees"
	ang1::Tv

	"Winding 1 first rating in MVA"
	rata1::Tv

	"Winding 1 second rating in MVA"
	ratb1::Tv

	"Winding 1 third rating in MVA"
	ratc1::Tv

	"Winding 1 control mode"
	cod1::Int8

	"Winding 1 controlled bus number"
	cont1::Int32

	"Winding 1 upper limit on adjustment"
	rma1::Tv

	"Winding 1 lower limit on adjustment"
	rmi1::Tv

	"Winding 1 upper voltage limit for controlled bus"
	vma1::Tv

	"Winding 1 lower voltage limit for controlled bus"
	vmi1::Tv

	"Winding 1 number of tap positions"
	ntp1::Int16

	"Winding 1 impedance correction table number"
	tab1::Int64

	"Winding 1 load drop compensation resistance"
	cr1::Tv

	"Winding 1 load drop compensation reactance"
	cx1::Tv

	"Winding 1 connection angle"
	cnxa1::Union{Missing, Tv}

	# Winding 2 parameters
	"Winding 2 off-nominal turns ratio or voltage"
	windv2::Tv

	"Winding 2 nominal voltage in kV"
	nomv2::Tv

	"Winding 2 phase shift angle in degrees"
	ang2::Tv

	"Winding 2 first rating in MVA"
	rata2::Tv

	"Winding 2 second rating in MVA"
	ratb2::Tv

	"Winding 2 third rating in MVA"
	ratc2::Tv

	"Winding 2 control mode"
	cod2::Int8

	"Winding 2 controlled bus number"
	cont2::Int32

	"Winding 2 upper limit on adjustment"
	rma2::Tv

	"Winding 2 lower limit on adjustment"
	rmi2::Tv

	"Winding 2 upper voltage limit for controlled bus"
	vma2::Tv

	"Winding 2 lower voltage limit for controlled bus"
	vmi2::Tv

	"Winding 2 number of tap positions"
	ntp2::Int16

	"Winding 2 impedance correction table number"
	tab2::Int64

	"Winding 2 load drop compensation resistance"
	cr2::Tv

	"Winding 2 load drop compensation reactance"
	cx2::Tv

	"Winding 2 connection angle"
	cnxa2::Union{Missing, Tv}

	# Winding 3 parameters
	"Winding 3 off-nominal turns ratio or voltage"
	windv3::Tv

	"Winding 3 nominal voltage in kV"
	nomv3::Tv

	"Winding 3 phase shift angle in degrees"
	ang3::Tv

	"Winding 3 first rating in MVA"
	rata3::Tv

	"Winding 3 second rating in MVA"
	ratb3::Tv

	"Winding 3 third rating in MVA"
	ratc3::Tv

	"Winding 3 control mode"
	cod3::Int8

	"Winding 3 controlled bus number"
	cont3::Int32

	"Winding 3 upper limit on adjustment"
	rma3::Tv

	"Winding 3 lower limit on adjustment"
	rmi3::Tv

	"Winding 3 upper voltage limit for controlled bus"
	vma3::Tv

	"Winding 3 lower voltage limit for controlled bus"
	vmi3::Tv

	"Winding 3 number of tap positions"
	ntp3::Int16

	"Winding 3 impedance correction table number"
	tab3::Int64

	"Winding 3 load drop compensation resistance"
	cr3::Tv

	"Winding 3 load drop compensation reactance"
	cx3::Tv

	"Winding 3 connection angle"
	cnxa3::Union{Missing, Tv}
end

# NOTE: Transformer3W is NOT registered as a model because it's only a data structure
# used for conversion to Transformer2W with virtual buses. It's handled in load_models.

# Type aliases
"""
    Transformer3WVec{T}

Alias for `StructArray{Transformer3W{T}}`, the struct-of-arrays container holding all
[`Transformer3W`](@ref) instances of a system.
"""
const Transformer3WVec{T} = StructArray{Transformer3W{T}}

get_key_name(::Type{T}) where T <: Transformer3W = [:i, :j, :k, :ckt]

# Note: No ModelMetadata for Transformer3W - it's converted before residual calculation

### === Format Support for 3W Transformers ===

supports_format(::Type{Transformer3W}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Transformer3W}, ::Type{PSSE}) = (:transformers,)

function parse_model(::Type{Transformer3W}, raw, ::FormatSupport{PSSE})
	transformers = raw.transformers

	# Find indices of three-winding transformers (k != 0)
	idx_3w = findall(x -> x != 0, transformers.k)

	if isempty(idx_3w)
		# Return empty StructArray - create from empty vector
		return StructArray(Transformer3W{Float64}[])
	end

	# Create StructArray from PSS/E data
	return StructArray{Transformer3W{Float64}}((
		i = transformers.i[idx_3w],
		j = transformers.j[idx_3w],
		k = transformers.k[idx_3w],
		ckt = transformers.ckt[idx_3w],
		cw = transformers.cw[idx_3w],
		cz = transformers.cz[idx_3w],
		cm = transformers.cm[idx_3w],
		mag1 = transformers.mag1[idx_3w],
		mag2 = transformers.mag2[idx_3w],
		nmetr = transformers.nmetr[idx_3w],
		name = transformers.name[idx_3w],
		stat = transformers.stat[idx_3w],
		o1 = transformers.o1[idx_3w],
		f1 = transformers.f1[idx_3w],
		o2 = transformers.o2[idx_3w],
		f2 = transformers.f2[idx_3w],
		o3 = transformers.o3[idx_3w],
		f3 = transformers.f3[idx_3w],
		o4 = transformers.o4[idx_3w],
		f4 = transformers.f4[idx_3w],
		vecgrp = transformers.vecgrp[idx_3w],
		r1_2 = transformers.r1_2[idx_3w],
		x1_2 = transformers.x1_2[idx_3w],
		sbase1_2 = transformers.sbase1_2[idx_3w],
		r2_3 = transformers.r2_3[idx_3w],
		x2_3 = transformers.x2_3[idx_3w],
		sbase2_3 = transformers.sbase2_3[idx_3w],
		r3_1 = transformers.r3_1[idx_3w],
		x3_1 = transformers.x3_1[idx_3w],
		sbase3_1 = transformers.sbase3_1[idx_3w],
		vmstar = transformers.vmstar[idx_3w],
		anstar = transformers.anstar[idx_3w],
		windv1 = transformers.windv1[idx_3w],
		nomv1 = transformers.nomv1[idx_3w],
		ang1 = transformers.ang1[idx_3w],
		rata1 = transformers.rata1[idx_3w],
		ratb1 = transformers.ratb1[idx_3w],
		ratc1 = transformers.ratc1[idx_3w],
		cod1 = transformers.cod1[idx_3w],
		cont1 = transformers.cont1[idx_3w],
		rma1 = transformers.rma1[idx_3w],
		rmi1 = transformers.rmi1[idx_3w],
		vma1 = transformers.vma1[idx_3w],
		vmi1 = transformers.vmi1[idx_3w],
		ntp1 = transformers.ntp1[idx_3w],
		tab1 = transformers.tab1[idx_3w],
		cr1 = transformers.cr1[idx_3w],
		cx1 = transformers.cx1[idx_3w],
		cnxa1 = transformers.cnxa1[idx_3w],
		windv2 = transformers.windv2[idx_3w],
		nomv2 = transformers.nomv2[idx_3w],
		ang2 = transformers.ang2[idx_3w],
		rata2 = transformers.rata2[idx_3w],
		ratb2 = transformers.ratb2[idx_3w],
		ratc2 = transformers.ratc2[idx_3w],
		cod2 = transformers.cod2[idx_3w],
		cont2 = transformers.cont2[idx_3w],
		rma2 = transformers.rma2[idx_3w],
		rmi2 = transformers.rmi2[idx_3w],
		vma2 = transformers.vma2[idx_3w],
		vmi2 = transformers.vmi2[idx_3w],
		ntp2 = transformers.ntp2[idx_3w],
		tab2 = transformers.tab2[idx_3w],
		cr2 = transformers.cr2[idx_3w],
		cx2 = transformers.cx2[idx_3w],
		cnxa2 = transformers.cnxa2[idx_3w],
		windv3 = transformers.windv3[idx_3w],
		nomv3 = transformers.nomv3[idx_3w],
		ang3 = transformers.ang3[idx_3w],
		rata3 = transformers.rata3[idx_3w],
		ratb3 = transformers.ratb3[idx_3w],
		ratc3 = transformers.ratc3[idx_3w],
		cod3 = transformers.cod3[idx_3w],
		cont3 = transformers.cont3[idx_3w],
		rma3 = transformers.rma3[idx_3w],
		rmi3 = transformers.rmi3[idx_3w],
		vma3 = transformers.vma3[idx_3w],
		vmi3 = transformers.vmi3[idx_3w],
		ntp3 = transformers.ntp3[idx_3w],
		tab3 = transformers.tab3[idx_3w],
		cr3 = transformers.cr3[idx_3w],
		cx3 = transformers.cx3[idx_3w],
		cnxa3 = transformers.cnxa3[idx_3w],
	))
end

# ============================================================================
# End of Three-Winding Transformer
# ============================================================================

@register_model(Transformer2W)

# Type aliases for StructArray-based vectorized types
"""
    Transformer2WVec{T}

Alias for `StructArray{Transformer2W{T}}`, the struct-of-arrays container holding all
[`Transformer2W`](@ref) instances of a system.
"""
const Transformer2WVec{T} = StructArray{Transformer2W{T}}

get_key_name(::Type{T}) where T <: Transformer2W = [:i, :j, :ckt]

# Base quantity traits - Transformer inherits base voltage from both connected buses
base_quantity_role(::Type{T}, ::Type{BaseVoltage}) where T <: Transformer2W = InheritsBase()

# Connectivity declaration - single source of truth for topology AND ExtParams
const TRANSFORMER2W_TERMINALS = (
    Terminal(:Bus, :i),  # Winding 1 bus
    Terminal(:Bus, :j)   # Winding 2 bus
)

get_terminals(::Type{T}) where T <: Transformer2W = TRANSFORMER2W_TERMINALS

const TRANSFORMER2W_VARS = [
	ModelVar(:theta_i, AlgebVar(), :Bus, :theta, :i,
		description = "Transformer voltage angle at bus i",
		units = "rad",
	),
	ModelVar(:theta_j, AlgebVar(), :Bus, :theta, :j,
		description = "Transformer voltage angle at bus j",
		units = "rad",
	),
	ModelVar(:v_i, AlgebVar(), :Bus, :v, :i,
		description = "Transformer voltage magnitude at bus i",
		units = "pu",
	),
	ModelVar(:v_j, AlgebVar(), :Bus, :v, :j,
		description = "Transformer voltage magnitude at bus j",
		units = "pu",
	),
]

const TRANSFORMER2W_observable_vars = []

# CurrentInjection formulation - rectangular coordinates
const TRANSFORMER2W_VARS_CURRENT = [
	ModelVar(:v_r_i, AlgebVar(), :Bus, :v_r, :i,
		description = "Transformer voltage real part at bus i",
		units = "pu",
	),
	ModelVar(:v_i_i, AlgebVar(), :Bus, :v_i, :i,
		description = "Transformer voltage imaginary part at bus i",
		units = "pu",
	),
	ModelVar(:v_r_j, AlgebVar(), :Bus, :v_r, :j,
		description = "Transformer voltage real part at bus j",
		units = "pu",
	),
	ModelVar(:v_i_j, AlgebVar(), :Bus, :v_i, :j,
		description = "Transformer voltage imaginary part at bus j",
		units = "pu",
	),
]

const TRANSFORMER2W_VARS_DC = [
	ModelVar(:theta_i, AlgebVar(), :Bus, :theta, :i,
		description = "Transformer voltage angle at bus i",
		units = "rad",
	),
	ModelVar(:theta_j, AlgebVar(), :Bus, :theta, :j,
		description = "Transformer voltage angle at bus j",
		units = "rad",
	),
]

const TRANSFORMER2W_RESIDUALS = [
	ModelResidual(:p_i, AlgebRes(), SharedRes(), :Bus, :p, :i,
		description = "Active power flow to bus i",
	),
	ModelResidual(:p_j, AlgebRes(), SharedRes(), :Bus, :p, :j,
		description = "Active power flow to bus j",
	),
	ModelResidual(:q_i, AlgebRes(), SharedRes(), :Bus, :q, :i,
		description = "Reactive power flow to bus i",
	),
	ModelResidual(:q_j, AlgebRes(), SharedRes(), :Bus, :q, :j,
		description = "Reactive power flow to bus j",
	),
]

# CurrentInjection formulation residuals
const TRANSFORMER2W_RESIDUALS_CURRENT = [
	ModelResidual(:i_r_i, AlgebRes(), SharedRes(), :Bus, :i_r, :i,
		description = "Real current flow to bus i",
	),
	ModelResidual(:i_i_i, AlgebRes(), SharedRes(), :Bus, :i_i, :i,
		description = "Imaginary current flow to bus i",
	),
	ModelResidual(:i_r_j, AlgebRes(), SharedRes(), :Bus, :i_r, :j,
		description = "Real current flow to bus j",
	),
	ModelResidual(:i_i_j, AlgebRes(), SharedRes(), :Bus, :i_i, :j,
		description = "Imaginary current flow to bus j",
	),
]

# DC Power Flow formulation residuals
const TRANSFORMER2W_RESIDUALS_DC = [
	ModelResidual(:p_i, AlgebRes(), SharedRes(), :Bus, :p, :i,
		description = "Active power flow to bus i (DC approximation)",
	),
	ModelResidual(:p_j, AlgebRes(), SharedRes(), :Bus, :p, :j,
		description = "Active power flow to bus j (DC approximation)",
	),
]

# Transformer parameters for per-unit conversion
# NOTE: Transformer per-unit conversion is complex due to cz/cw/cm flags:
# - cz: impedance data code (1=system base, 2=specified MVA, 3=load loss)
# - cw: winding data code (1=pu, 2=kV)
# - cm: magnetizing admittance code (1=pu, 2=physical units)
# The generic conversion system may need special handling for transformers.
const TRANSFORMER2W_PARAMS = [
	ModelParam(:r1_2, ImpedancePU, description="Winding 1-2 resistance", units="pu"),
	ModelParam(:x1_2, ImpedancePU, description="Winding 1-2 reactance", units="pu"),
]

# Transformer observations (computed from bus voltages and parameters)
const TRANSFORMER2W_OBSERVATIONS = [
	ModelObservable(:p_from, ActivePower, units = "pu",
		description = "Active power flow from bus i (computed)", default_output = true),
	ModelObservable(:q_from, ReactivePower, units = "pu",
		description = "Reactive power flow from bus i (computed)", default_output = true),
	ModelObservable(:p_to, ActivePower, units = "pu",
		description = "Active power flow from bus j (computed)", default_output = true),
	ModelObservable(:q_to, ReactivePower, units = "pu",
		description = "Reactive power flow from bus j (computed)", default_output = true),
	ModelObservable(:p_loss, ActivePower, units = "pu",
		description = "Active power loss in transformer (computed)", default_output = false),
	ModelObservable(:q_loss, ReactivePower, units = "pu",
		description = "Reactive power loss in transformer (computed)", default_output = false),
]

# Transformer device rating (MVA base for winding impedance when cz=2)
const TRANSFORMER2W_RATINGS = [
	ModelRating(:sbase1_2, NominalPower, description="Winding 1-2 MVA base")
]

# Transformer references both bus base voltages - uses TRANSFORMER2W_TERMINALS
# Following ANDES convention: primary terminal :i used for per-unit conversion
const TRANSFORMER2W_EXT_PARAMS = [
	ExtParam(:i_basekv, TRANSFORMER2W_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
		required=true, description="Bus i base voltage"),
	ExtParam(:j_basekv, TRANSFORMER2W_TERMINALS[2], :basekv, BaseVoltage(), PrePerUnit(),
		required=true, description="Bus j base voltage"),
	ExtParam(:i_ide, TRANSFORMER2W_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
		required=true, description="Bus type at terminal i"),
	ExtParam(:j_ide, TRANSFORMER2W_TERMINALS[2], :ide, Dimensionless(), PrePerUnit(),
		required=true, description="Bus type at terminal j")
]

function model_metadata(::Type{Transformer2W})
	metadata = ModelMetadata(
		name = :Transformer2W,
		vars = TRANSFORMER2W_VARS,
		observable_vars = TRANSFORMER2W_observable_vars,
		residuals = TRANSFORMER2W_RESIDUALS,
		params = TRANSFORMER2W_PARAMS,
		observations = TRANSFORMER2W_OBSERVATIONS,
		ratings = TRANSFORMER2W_RATINGS,
		ext_params = TRANSFORMER2W_EXT_PARAMS,
		status_field = :status,
	)

	# Register CurrentInjection formulation overrides
	metadata.vars_by_formulation[CurrentInjection] = TRANSFORMER2W_VARS_CURRENT
	metadata.residuals_by_formulation[CurrentInjection] = TRANSFORMER2W_RESIDUALS_CURRENT

	# Register DC Power Flow formulation overrides
	metadata.vars_by_formulation[DCPowerFlow] = TRANSFORMER2W_VARS_DC
	metadata.residuals_by_formulation[DCPowerFlow] = TRANSFORMER2W_RESIDUALS_DC

	return metadata
end

# Declare formulation support
supports_formulation(::Type{Transformer2W}, ::Type{PowerInjection}) = true
supports_formulation(::Type{Transformer2W}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{Transformer2W}, ::Type{DCPowerFlow}) = true

### === Model Format Adapters ===

supports_format(::Type{Transformer2W}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{Transformer2W}, ::Type{PSSE}) = (:transformers,)

function parse_model(::Type{Transformer2W}, raw, ::FormatSupport{PSSE})
	transformers = raw.transformers

	# Find indices of two-winding transformers (k == 0)
	# Three-winding transformers (k != 0) are parsed separately by Transformer3W
	idx_2w = findall(x -> x == 0, transformers.k)

	# Extract raw r/x values and cz codes
	cz_vec = transformers.cz[idx_2w]
	r_raw = transformers.r1_2[idx_2w]
	x_raw = transformers.x1_2[idx_2w]
	sbase_winding = transformers.sbase1_2[idx_2w]

	# Convert CZ=3 (load loss format) to CZ=2 equivalent (pu on winding base)
	# CZ=3: r = load loss in watts, x = impedance magnitude |Z| in pu on winding base
	# Need to convert to: r, x in pu on winding base (same as CZ=2)
	r_converted = similar(r_raw)
	x_converted = similar(x_raw)
	cz_converted = copy(cz_vec)  # Copy so we can update converted entries
	for i in eachindex(cz_vec)
		if cz_vec[i] == 3
			r_converted[i], x_converted[i] = cz3_to_winding_pu(r_raw[i], x_raw[i], sbase_winding[i])
			cz_converted[i] = Int8(2)  # Now in pu on winding base (CZ=2 format)
		else
			# CZ=1 or CZ=2: no pre-conversion needed (handled by per-unit system)
			r_converted[i] = r_raw[i]
			x_converted[i] = x_raw[i]
		end
	end

	# Convert windv to CW=1: tap ratio in pu of the winding bus base voltage
	# (winding 1 at bus i, winding 2 at bus j; see tap_ratio_pu)
	cw_vec = transformers.cw[idx_2w]
	windv1_raw = transformers.windv1[idx_2w]
	windv2_raw = transformers.windv2[idx_2w]
	nomv1_raw = transformers.nomv1[idx_2w]
	nomv2_raw = transformers.nomv2[idx_2w]
	basekv = Dict(zip(raw.buses.i, raw.buses.basekv))
	bus_i = transformers.i[idx_2w]
	bus_j = transformers.j[idx_2w]

	windv1_converted = similar(windv1_raw)
	windv2_converted = similar(windv2_raw)
	cw_converted = fill(Int8(1), length(cw_vec))  # Now in per-unit of the bus base voltage
	for i in eachindex(cw_vec)
		windv1_converted[i] = tap_ratio_pu(windv1_raw[i], nomv1_raw[i], cw_vec[i], basekv[bus_i[i]])
		windv2_converted[i] = tap_ratio_pu(windv2_raw[i], nomv2_raw[i], cw_vec[i], basekv[bus_j[i]])
	end

	# Convert magnetizing data to CM=1: pu on the system MVA base and winding 1 bus voltage base
	cm_vec = transformers.cm[idx_2w]
	mag1_raw = transformers.mag1[idx_2w]
	mag2_raw = transformers.mag2[idx_2w]
	mag1_converted = similar(mag1_raw)
	mag2_converted = similar(mag2_raw)
	for i in eachindex(cm_vec)
		mag1_converted[i], mag2_converted[i] = convert_magnetizing_to_pu(
			mag1_raw[i], mag2_raw[i], cm_vec[i], sbase_winding[i], raw.caseid.sbase,
			basekv[bus_i[i]], nomv1_raw[i])
	end

	# Create StructArray directly from named tuple of vectors
	return StructArray{Transformer2W{Float64}}((
		i = transformers.i[idx_2w],
		j = transformers.j[idx_2w],
		ckt = transformers.ckt[idx_2w],
		cw = cw_converted,  # CW=2/3 converted to CW=1
		cz = cz_converted,
		cm = fill(Int8(1), length(cm_vec)),  # CM=2 converted to CM=1
		mag1 = mag1_converted,
		mag2 = mag2_converted,
		nmetr = transformers.nmetr[idx_2w],
		name = transformers.name[idx_2w],
		status = transformers.stat[idx_2w] .!= 0,  # PSS/E STAT code to Bool
		o1 = transformers.o1[idx_2w],
		f1 = transformers.f1[idx_2w],
		o2 = transformers.o2[idx_2w],
		f2 = transformers.f2[idx_2w],
		o3 = transformers.o3[idx_2w],
		f3 = transformers.f3[idx_2w],
		o4 = transformers.o4[idx_2w],
		f4 = transformers.f4[idx_2w],
		vecgrp = transformers.vecgrp[idx_2w],
		r1_2 = r_converted,
		x1_2 = x_converted,
		sbase1_2 = sbase_winding,
		windv1 = windv1_converted,  # pu of bus base voltage
		nomv1 = nomv1_raw,
		ang1 = transformers.ang1[idx_2w],
		rata1 = transformers.rata1[idx_2w],
		ratb1 = transformers.ratb1[idx_2w],
		ratc1 = transformers.ratc1[idx_2w],
		cod1 = transformers.cod1[idx_2w],
		cont1 = transformers.cont1[idx_2w],
		rma1 = transformers.rma1[idx_2w],
		rmi1 = transformers.rmi1[idx_2w],
		vma1 = transformers.vma1[idx_2w],
		vmi1 = transformers.vmi1[idx_2w],
		ntp1 = transformers.ntp1[idx_2w],
		tab1 = transformers.tab1[idx_2w],
		cr1 = transformers.cr1[idx_2w],
		cx1 = transformers.cx1[idx_2w],
		cnxa1 = transformers.cnxa1[idx_2w],
		windv2 = windv2_converted,  # pu of bus base voltage
		nomv2 = nomv2_raw,
	))
end


### === Residual Calculation === ###


# ============================================================================
# Pure Computation Kernel (CPU/GPU Ready)
# ============================================================================

"""
    transformer_residual_kernel_current_injection!(outputs, params, vars, active, ::Type{AT})

CPU kernel for PowerInjection formulation using current-based calculation.

NOTE: Despite the name, this kernel is used for PowerInjection formulation.
It computes power residuals (P, Q) by first calculating current injections,
then converting them to power via S = V * conj(I).

This approach is more robust for transformers with taps and phase shifts than
the direct power formulation (used in `transformer_residual_kernel!`).

Use `transformer_current_kernel!` for the CurrentInjection formulation.

# Arguments
- `outputs::NamedTuple{(:p_i, :p_j, :q_i, :q_j)}`: Power flow residuals
- `params::NamedTuple{(:r, :x, :mag1, :mag2, :windv1, :windv2, :ang1)}`: Raw transformer parameters
- `vars::NamedTuple{(:theta_i, :theta_j, :v_i, :v_j)}`: Bus voltages in polar coordinates
- `active::AbstractVector{Bool}`: Transformer in-service status
- `::Type{AT}`: Array type for dispatch (Vector for CPU)
"""
@inline function transformer_residual_kernel_current_injection!(
    outputs::NamedTuple{(:p_i, :p_j, :q_i, :q_j)},
    params::NamedTuple{(:r, :x, :mag1, :mag2, :windv1, :windv2, :ang1)},
    vars::NamedTuple{(:theta_i, :theta_j, :v_i, :v_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p_i)
        # Branchless masking - only active status (bus type filtering handled by scatter :set mode)
        mask = active[i]

        # Calculate current injections using extracted physics function
        currents = calc_transformer_currents(
            vars.v_i[i], vars.theta_i[i], vars.v_j[i], vars.theta_j[i],
            params.r[i], params.x[i], params.mag1[i], params.mag2[i],
            params.windv1[i], params.windv2[i], params.ang1[i]
        )

        # Extract values for readability
        v_i = vars.v_i[i]
        v_j = vars.v_j[i]
        theta_ij = vars.theta_i[i] - vars.theta_j[i]
        cos_theta_ij = cos(theta_ij)
        sin_theta_ij = sin(theta_ij)

        # Convert currents to power: S = V * conj(I)
        # Bus i (reference, θ_i=0): P = v*Re(I), Q = -v*Im(I)
        p_i_inj = v_i * currents.i_i_real
        q_i_inj = -v_i * currents.i_i_imag

        # Bus j (non-reference, θ_j=-θ_ij): Must account for voltage phase!
        # S_j = V_j * conj(I_j) = v_j*e^(jθ_j) * (i_j_real - j*i_j_imag)
        # With θ_j = -θ_ij:
        # P_j = v_j * (i_j_real*cos(θ_ij) - i_j_imag*sin(θ_ij))
        # Q_j = v_j * (-i_j_real*sin(θ_ij) - i_j_imag*cos(θ_ij))
        p_j_inj = v_j * (currents.i_j_real * cos_theta_ij - currents.i_j_imag * sin_theta_ij)
        q_j_inj = v_j * (-currents.i_j_real * sin_theta_ij - currents.i_j_imag * cos_theta_ij)

        # Sign convention: power flowing INTO bus is positive (injection)
        # p_i_inj etc. compute power OUT, so negate for injection
        outputs.p_i[i] = mask * (-p_i_inj)
        outputs.q_i[i] = mask * (-q_i_inj)
        outputs.p_j[i] = mask * (-p_j_inj)
        outputs.q_j[i] = mask * (-q_j_inj)
    end

    nothing
end

# ============================================================================
# Data Adapter Layer
# ============================================================================

# ============================================================================
# NEW IMPLEMENTATION: Direct Power Formulation (following ANDES)
# ============================================================================

# ============================================================================
# Transformer Physics (Reusable for CPU/GPU Kernels and Reports)
# ============================================================================

"""
    calc_transformer_currents(v_i, theta_i, v_j, theta_j, r, x, g_mag, b_mag,
                               windv1, windv2, ang1) -> NamedTuple

Calculate transformer current injections using π-model with tap and phase shift.

# Transformer Model
- Series admittance: y_series = 1/(r + jx) = (r - jx)/(r² + x²)
- Magnetizing admittance: y_mag = g_mag + j*b_mag, a shunt at the winding 1 bus (bus i)
  on the bus side of the ideal transformer, so it is not scaled by the tap (RAW model)
- Tap ratio: tap = windv1/windv2 (applied to winding 1 side)
- Phase shift: φ = ang1 (degrees → radians, applied to winding 1)

# Arguments
- `v_i, theta_i`: Voltage magnitude and angle at bus i (winding 1)
- `v_j, theta_j`: Voltage magnitude and angle at bus j (winding 2)
- `r, x`: Series resistance and reactance (pu on system base)
- `g_mag, b_mag`: Magnetizing conductance and susceptance (pu on system base)
- `windv1, windv2`: Winding voltage ratios
- `ang1`: Phase shift angle at winding 1 (degrees)

# Returns
Named tuple `(i_i_real=..., i_i_imag=..., i_j_real=..., i_j_imag=...)` with current injections in pu.
- Current components are in rectangular form (real and imaginary parts)
- Positive current: flowing OUT of the bus (into the network)

# Note
This is a pure physics function compatible with CPU, GPU, and reporting.
Reusable for both PowerInjection (converts I→P) and CurrentInjection formulations.
"""
@inline function calc_transformer_currents(
    v_i, theta_i, v_j, theta_j,
    r, x, g_mag, b_mag,
    windv1, windv2, ang1
)
    # Tap ratio: t = windv1/windv2
    t = windv1 / windv2

    # Phase shift (degrees to radians)
    T = typeof(ang1)
    phi = ang1 * (T(π) / T(180))
    cos_phi = cos(phi)
    sin_phi = sin(phi)

    # Series admittance: y = 1/(r + jx) = (r - jx)/(r² + x²)
    z_sq = r * r + x * x
    g_series = r / z_sq
    b_series = -x / z_sq

    # Angle difference
    theta_ij = theta_i - theta_j
    cos_theta_ij = cos(theta_ij)
    sin_theta_ij = sin(theta_ij)

    # Admittance matrix elements (π-model with tap ratio t on bus i side)
    t_sq = t * t

    # Y_ii = (g + jb) / t² + (g_mag + jb_mag)   (magnetizing shunt on the bus side of the tap)
    g_ii = g_series / t_sq + g_mag
    b_ii = b_series / t_sq + b_mag

    # Y_ij = -(g + jb)/t * e^(jφ)
    #      = -[(g*cos(φ) - b*sin(φ)) + j(b*cos(φ) + g*sin(φ))] / t
    g_ij = -(g_series * cos_phi - b_series * sin_phi) / t
    b_ij = -(b_series * cos_phi + g_series * sin_phi) / t

    # Y_ji = -(g + jb)/t * e^(-jφ)
    #      = -[(g*cos(φ) + b*sin(φ)) + j(b*cos(φ) - g*sin(φ))] / t
    g_ji = -(g_series * cos_phi + b_series * sin_phi) / t
    b_ji = -(b_series * cos_phi - g_series * sin_phi) / t

    # Y_jj = (g + jb)   (no magnetizing admittance on winding 2)
    g_jj = g_series
    b_jj = b_series

    # Current injections: I = Y * V
    # I_i = Y_ii*V_i + Y_ij*V_j where V_i = v_i∠θ_i, V_j = v_j∠θ_j
    i_i_real = g_ii * v_i + v_j * (g_ij * cos_theta_ij + b_ij * sin_theta_ij)
    i_i_imag = b_ii * v_i + v_j * (b_ij * cos_theta_ij - g_ij * sin_theta_ij)

    # I_j = Y_ji*V_i + Y_jj*V_j
    i_j_real = v_i * g_ji + v_j * (g_jj * cos_theta_ij + b_jj * sin_theta_ij)
    i_j_imag = v_i * b_ji + v_j * (b_jj * cos_theta_ij - g_jj * sin_theta_ij)

    return (i_i_real=i_i_real, i_i_imag=i_i_imag, i_j_real=i_j_real, i_j_imag=i_j_imag)
end

"""
    calc_transformer_power_flows(v_i, theta_i, v_j, theta_j, r, x, g_mag, b_mag,
                                  windv1, windv2, ang1) -> NamedTuple

Calculate transformer power flow using ANDES formulation with split π-model.

# Transformer Model
- Series admittance: y_series = 1/(r + jx) = (r - jx)/(r² + x²)
- Magnetizing admittance: y_mag = g_mag + j*b_mag, a shunt at the winding 1 bus (bus i)
  on the bus side of the ideal transformer, so it is not scaled by the tap (RAW model)
- Tap ratio: tap = windv1/windv2 (applied to winding 1 side)
- Phase shift: φ = ang1 (degrees → radians, applied to winding 1)

# Arguments
- `v_i, theta_i`: Voltage magnitude and angle at bus i (winding 1)
- `v_j, theta_j`: Voltage magnitude and angle at bus j (winding 2)
- `r, x`: Series resistance and reactance (pu on system base)
- `g_mag, b_mag`: Magnetizing conductance and susceptance (pu on system base)
- `windv1, windv2`: Winding voltage ratios
- `ang1`: Phase shift angle at winding 1 (degrees)

# Returns
Named tuple `(p_i=..., q_i=..., p_j=..., q_j=...)` with power flows in pu.
- Positive values: power flowing OUT of the bus (ANDES convention)

# Note
This is a pure physics function compatible with CPU, GPU, and reporting.
The magnetizing admittance gh/bh is a shunt at bus i only (see the model description above).
"""
@inline function calc_transformer_power_flows(
    v_i, theta_i, v_j, theta_j,
    r, x, g_mag, b_mag,
    windv1, windv2, ang1
)
    # Tap ratio (windv1/windv2)
    tap = windv1 / windv2
    itap = one(tap) / tap
    itap2 = itap * itap

    # Phase shift (degrees to radians)
    T = typeof(ang1)
    phi = ang1 * (T(π) / T(180))

    # Series admittance: y = 1/(r + jx) = (r - jx)/(r² + x²)
    z_sq = r * r + x * x
    ghk = r / z_sq
    bhk = -x / z_sq

    # Magnetizing admittance: shunt at bus i on the bus side of the tap (not scaled)
    gh = g_mag
    bh = b_mag

    # Angle difference
    delta = theta_i - theta_j - phi
    cos_delta = cos(delta)
    sin_delta = sin(delta)

    # P_i = v1² * (gh + ghk*itap2) - v1*v2*(ghk*cos(delta) + bhk*sin(delta)) * itap
    p_i = v_i * v_i * (gh + ghk * itap2) - v_i * v_j * (ghk * cos_delta + bhk * sin_delta) * itap

    # Q_i = -v1² * (bh + bhk*itap2) - v1*v2*(ghk*sin(delta) - bhk*cos(delta)) * itap
    q_i = -v_i * v_i * (bh + bhk * itap2) - v_i * v_j * (ghk * sin_delta - bhk * cos_delta) * itap

    # P_j = v2² * ghk - v1*v2*(ghk*cos(delta) - bhk*sin(delta)) * itap
    p_j = v_j * v_j * ghk - v_i * v_j * (ghk * cos_delta - bhk * sin_delta) * itap

    # Q_j = -v2² * bhk + v1*v2*(ghk*sin(delta) + bhk*cos(delta)) * itap
    q_j = -v_j * v_j * bhk + v_i * v_j * (ghk * sin_delta + bhk * cos_delta) * itap

    return (p_i=p_i, q_i=q_i, p_j=p_j, q_j=q_j)
end

"""
CPU kernel for transformer using ANDES' direct power formulation.
This avoids the current computation and directly computes power flows.
"""
@inline function transformer_residual_kernel!(
    outputs::NamedTuple{(:p_i, :p_j, :q_i, :q_j)},
    params::NamedTuple{(:r, :x, :mag1, :mag2, :windv1, :windv2, :ang1)},
    vars::NamedTuple{(:theta_i, :theta_j, :v_i, :v_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p_i)
        # Branchless masking - only active status
        # Bus type filtering handled by scatter mode (:add vs :set)
        mask = active[i]

        # Calculate power flows using extracted physics function
        flows = calc_transformer_power_flows(
            vars.v_i[i], vars.theta_i[i], vars.v_j[i], vars.theta_j[i],
            params.r[i], params.x[i], params.mag1[i], params.mag2[i],
            params.windv1[i], params.windv2[i], params.ang1[i]
        )

        # Sign convention: power flowing INTO bus is positive (injection)
        # flows.p_i etc. compute power OUT, so negate for injection
        outputs.p_i[i] = mask * (-flows.p_i)
        outputs.q_i[i] = mask * (-flows.q_i)
        outputs.p_j[i] = mask * (-flows.p_j)
        outputs.q_j[i] = mask * (-flows.q_j)
    end

    nothing
end

# ============================================================================
# Data Adapter
# ============================================================================

"""
    calc_residuals!(ws::ResidualBuffers, model::StructArray{<:Transformer2W}, external::NamedTuple)

Compute transformer power flow residuals and scatter to system residual vector.
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    model::StructArray{<:Transformer2W},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and dispatch to kernel
    # Do NOT compute derived parameters here - breaks ForwardDiff!
    transformer_residual_kernel!(
        # Outputs: power flow residuals
        (p_i = ws.res_vals.p_i, p_j = ws.res_vals.p_j,
         q_i = ws.res_vals.q_i, q_j = ws.res_vals.q_j),

        # Parameters: pass raw model fields, compute tap/mag inside kernel
        (r = model.r1_2, x = model.x1_2,
         mag1 = model.mag1, mag2 = model.mag2,
         windv1 = model.windv1, windv2 = model.windv2,
         ang1 = model.ang1),

        # Variables: bus voltages
        (theta_i = ws.var_vals.theta_i, theta_j = ws.var_vals.theta_j,
         v_i = ws.var_vals.v_i, v_j = ws.var_vals.v_j),

        active_mask,
        AT  # Dispatches to CPU (Vector) or GPU (CuArray)
    )
end

# ============================================================================
# CurrentInjection Formulation Kernel
# ============================================================================

"""
CPU kernel for CurrentInjection formulation using rectangular voltages.

Computes current injection residuals directly from rectangular voltage coordinates.
Uses the same physics as PowerInjection but outputs currents instead of power.
"""
@inline function transformer_current_kernel!(
    outputs::NamedTuple{(:i_r_i, :i_i_i, :i_r_j, :i_i_j)},
    params::NamedTuple{(:r, :x, :mag1, :mag2, :windv1, :windv2, :ang1)},
    vars::NamedTuple{(:v_r_i, :v_i_i, :v_r_j, :v_i_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.i_r_i)
        # Branchless masking - only active status
        # Bus type filtering handled by scatter mode (:add vs :set)
        mask = active[i]

        # Convert rectangular voltages to polar
        v_r_i = vars.v_r_i[i]
        v_i_i = vars.v_i_i[i]
        v_r_j = vars.v_r_j[i]
        v_i_j = vars.v_i_j[i]

        v_i = sqrt(v_r_i * v_r_i + v_i_i * v_i_i)
        theta_i = atan(v_i_i, v_r_i)
        v_j = sqrt(v_r_j * v_r_j + v_i_j * v_i_j)
        theta_j = atan(v_i_j, v_r_j)

        # Calculate current injections using physics function
        # Note: calc_transformer_currents returns currents in LOCAL frame (bus i at angle 0)
        currents = calc_transformer_currents(
            v_i, theta_i, v_j, theta_j,
            params.r[i], params.x[i], params.mag1[i], params.mag2[i],
            params.windv1[i], params.windv2[i], params.ang1[i]
        )

        # Rotate currents from local frame to GLOBAL rectangular frame
        # Local frame: bus i at angle 0, bus j at angle -(theta_i - theta_j)
        # Global frame: bus i at theta_i, bus j at theta_j
        # Rotation: Rotate entire local frame by theta_i
        # I_global = I_local * e^(j*theta_i)
        cos_theta_i = cos(theta_i)
        sin_theta_i = sin(theta_i)

        # Both buses rotate by theta_i (the reference rotation of the local frame)
        i_r_i_global = currents.i_i_real * cos_theta_i - currents.i_i_imag * sin_theta_i
        i_i_i_global = currents.i_i_real * sin_theta_i + currents.i_i_imag * cos_theta_i
        i_r_j_global = currents.i_j_real * cos_theta_i - currents.i_j_imag * sin_theta_i
        i_i_j_global = currents.i_j_real * sin_theta_i + currents.i_j_imag * cos_theta_i

        # Sign convention: Y*V gives current INTO element, negate for current INTO bus
        outputs.i_r_i[i] = mask * (-i_r_i_global)
        outputs.i_i_i[i] = mask * (-i_i_i_global)
        outputs.i_r_j[i] = mask * (-i_r_j_global)
        outputs.i_i_j[i] = mask * (-i_i_j_global)
    end

    nothing
end

# ============================================================================
# DC Power Flow Formulation Kernel
# ============================================================================

"""
    transformer_dc_kernel!(outputs, params, vars, active, ::Type{AT})

Pure computation kernel for transformer DC power flow residuals.

In DC power flow, transformers are modeled as reactances with linear power flow
(see `dc_susceptance`):
- P = (theta_i - theta_j) / (x · tap), tap = windv1/windv2
- Voltage magnitudes assumed at 1.0 pu, phase shifts ignored
- Resistance negligible (r << x)

Dispatches to CPU or GPU implementation based on array type parameter.

# Arguments
- `outputs::NamedTuple{(:p_i, :p_j)}`: Preallocated output buffers
- `params::NamedTuple{(:x, :windv1, :windv2)}`: Transformer reactance (pu) and winding ratios
- `vars::NamedTuple{(:theta_i, :theta_j)}`: Voltage angles from workspace
- `active::AbstractVector{Bool}`: Mask for in-service transformers
- `::Type{AT}`: Array type for dispatch (Vector for CPU, CuVector for GPU)

# Physics
DC approximation for transformer active power flow:
- P_i = (theta_i - theta_j) / (x · tap)  (power flow from i to j)
- P_j = -P_i (power conservation)
- All voltage magnitudes = 1.0 pu (phase shifts ignored)
"""
@inline function transformer_dc_kernel!(
    outputs::NamedTuple{(:p_i, :p_j)},
    params::NamedTuple{(:x, :windv1, :windv2)},
    vars::NamedTuple{(:theta_i, :theta_j)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p_i)
        # Branchless masking - only active status
        # Bus type filtering handled by scatter mode (:add vs :set)
        mask = active[i]

        # DC approximation: P = (theta_i - theta_j) / (x · tap)
        b = dc_susceptance(params.x[i], params.windv1[i] / params.windv2[i])
        p_flow = b * (vars.theta_i[i] - vars.theta_j[i])

        # Sign convention: power flowing INTO bus is positive (injection)
        # p_flow = power from i to j, so at bus i it leaves (negative), at bus j it enters (positive)
        outputs.p_i[i] = mask * (-p_flow)
        outputs.p_j[i] = mask * p_flow
    end

    nothing
end

"""
    calc_residuals!(ws::ResidualBuffers, model::StructArray{<:Transformer2W}, external::NamedTuple, ::CurrentInjection)

Compute transformer current injection residuals for CurrentInjection formulation.
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    model::StructArray{<:Transformer2W},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and call current injection kernel
    transformer_current_kernel!(
        # Outputs: current injection residuals
        (i_r_i = ws.res_vals.i_r_i, i_i_i = ws.res_vals.i_i_i,
         i_r_j = ws.res_vals.i_r_j, i_i_j = ws.res_vals.i_i_j),

        # Parameters: pass raw model fields
        (r = model.r1_2, x = model.x1_2,
         mag1 = model.mag1, mag2 = model.mag2,
         windv1 = model.windv1, windv2 = model.windv2,
         ang1 = model.ang1),

        # Variables: rectangular voltages
        (v_r_i = ws.var_vals.v_r_i, v_i_i = ws.var_vals.v_i_i,
         v_r_j = ws.var_vals.v_r_j, v_i_j = ws.var_vals.v_i_j),

        # Active mask
        active_mask,

        # Dispatch type
        AT
    )
end

"""
    calc_residuals!(ws::ResidualBuffers, model::StructArray{<:Transformer2W}, external::NamedTuple, ::DCPowerFlow)

Compute transformer DC power flow residuals.

In DC power flow, transformers are modeled as pure reactances with linear power flow.
Off-nominal tap ratios and phase shifts are ignored.
"""
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    model::StructArray{<:Transformer2W},
    _external::NamedTuple,
    ::DCPowerFlow,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    # Data adapter - extract arrays and call DC kernel
    transformer_dc_kernel!(
        # Outputs: active power residuals only
        (p_i = ws.res_vals.p_i, p_j = ws.res_vals.p_j),

        # Parameters: reactance and turns ratio (resistance ignored in DC)
        (x = model.x1_2, windv1 = model.windv1, windv2 = model.windv2),

        # Variables: voltage angles only
        (theta_i = ws.var_vals.theta_i, theta_j = ws.var_vals.theta_j),

        # Active mask
        active_mask,

        # Dispatch type
        AT
    )
end


# ============================================================================
# Reporting Interface - Vectorized Observation Computation
# ============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{<:Transformer2W}, var::Symbol) -> Vector{Float64}

Compute observation `var` for all two-winding transformers in the solution (vectorized).

# Supported Variables
- `:p_from`, `:q_from`: Power flow out of bus i into the transformer (pu)
- `:p_to`, `:q_to`: Power flow out of bus j into the transformer (pu)
- `:p_loss`, `:q_loss`: Losses = from + to (pu)

# Implementation Notes
Reuses `calc_transformer_power_flows`, the physics of the residual kernel; in
`DCPowerFlow` the flow is that of the DC kernel.
"""
function compute_observations(sol::PowerFlowSolution, ::Type{<:Transformer2W}, var::Symbol)
    xf = sol.prob.sys.models.Transformer2W
    var in (:p_from, :q_from, :p_to, :q_to, :p_loss, :q_loss) ||
        error("Unknown Transformer2W output variable: $var")

    # Terminal bus voltages (any formulation); an out-of-service transformer carries no flow
    v_i, theta_i = terminal_voltages(sol, :Transformer2W, :i)
    v_j, theta_j = terminal_voltages(sol, :Transformer2W, :j)
    n = length(xf)

    if sol.prob.formulation isa DCPowerFlow
        # Lossless DC flow, as in the DC kernel
        p_from = [xf.status[k] * dc_susceptance(xf.x1_2[k], xf.windv1[k] / xf.windv2[k]) *
                  (theta_i[k] - theta_j[k]) for k in 1:n]
        var === :p_from && return p_from
        var === :p_to && return -p_from
        var === :p_loss && return zeros(n)
        error("Transformer2W observation :$var is not defined in the DCPowerFlow formulation")
    end

    p_i = Vector{Float64}(undef, n)
    q_i = Vector{Float64}(undef, n)
    p_j = Vector{Float64}(undef, n)
    q_j = Vector{Float64}(undef, n)

    @inbounds for k in 1:n
        flows = calc_transformer_power_flows(
            v_i[k], theta_i[k], v_j[k], theta_j[k],
            xf.r1_2[k], xf.x1_2[k], xf.mag1[k], xf.mag2[k],
            xf.windv1[k], xf.windv2[k], xf.ang1[k]
        )
        on = xf.status[k]
        p_i[k] = on * flows.p_i
        q_i[k] = on * flows.q_i
        p_j[k] = on * flows.p_j
        q_j[k] = on * flows.q_j
    end

    var === :p_from && return p_i
    var === :q_from && return q_i
    var === :p_to && return p_j
    var === :q_to && return q_j
    var === :p_loss && return @. p_i + p_j
    return @. q_i + q_j  # :q_loss
end
