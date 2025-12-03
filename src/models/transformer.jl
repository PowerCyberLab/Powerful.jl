export Transformer2W, Transformer2WVec
export Transformer2WNum, Transformer2WNumVec

abstract type AbstractTransformer{T} <: AbstractModel end

"""
	$TYPEDEF

This struct represents a two-winding transformer, which is a subset of the data required for 
three-winding transformers.

# Fields

$TYPEDFIELDS
"""
struct Transformer2W{Tv} <: AbstractTransformer{Tv}
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

	"Initial transformer status (true=in-service)"
	stat::Bool

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

@register_model(Transformer2W)

# Type aliases for StructArray-based vectorized types
const Transformer2WVec{T} = StructArray{Transformer2W{T}}

get_key_name(::Type{T}) where T <: Transformer2W = [:i, :j, :ckt]

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

const TRANSFORMER2W_OUTPUT_VARS = []

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


function model_metadata(::Type{Transformer2W}; layout::LayoutStrategy = ContiguousVariables())
	ModelMetadata(
		name = :Transformer2W,
		vars = TRANSFORMER2W_VARS,
		output_vars = TRANSFORMER2W_OUTPUT_VARS,
		residuals = TRANSFORMER2W_RESIDUALS,
		layout = layout,
	)
end

### === Model Format Adapters ===

supports_format(::Type{Transformer2W}, ::Type{PSSE}) = FormatSupport{PSSE}()

function parse_model(::Type{Transformer2W}, raw, ::FormatSupport{PSSE})
	transformers = raw.transformers

	# Find indices of two-winding transformers
	idx_2w = findall(x -> x == 0, transformers.k)

	# Warn if there are three-winding transformers
	n_3w = length(transformers.k) - length(idx_2w)
	if n_3w > 0
		@warn "Found $n_3w three-winding transformer(s) which will be ignored"
	end

	# Create StructArray directly from named tuple of vectors
	return StructArray{Transformer2W{Float64}}((
		i = transformers.i[idx_2w],
		j = transformers.j[idx_2w],
		ckt = transformers.ckt[idx_2w],
		cw = transformers.cw[idx_2w],
		cz = transformers.cz[idx_2w],
		cm = transformers.cm[idx_2w],
		mag1 = transformers.mag1[idx_2w],
		mag2 = transformers.mag2[idx_2w],
		nmetr = transformers.nmetr[idx_2w],
		name = transformers.name[idx_2w],
		stat = transformers.stat[idx_2w],
		o1 = transformers.o1[idx_2w],
		f1 = transformers.f1[idx_2w],
		o2 = transformers.o2[idx_2w],
		f2 = transformers.f2[idx_2w],
		o3 = transformers.o3[idx_2w],
		f3 = transformers.f3[idx_2w],
		o4 = transformers.o4[idx_2w],
		f4 = transformers.f4[idx_2w],
		vecgrp = transformers.vecgrp[idx_2w],
		r1_2 = transformers.r1_2[idx_2w],
		x1_2 = transformers.x1_2[idx_2w],
		sbase1_2 = transformers.sbase1_2[idx_2w],
		windv1 = transformers.windv1[idx_2w],
		nomv1 = transformers.nomv1[idx_2w],
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
		windv2 = transformers.windv2[idx_2w],
		nomv2 = transformers.nomv2[idx_2w],
	))
end


### === Begin Tests === ###

@testitem "Transformer2W" begin
	using PowerFlowData
	using Powerful
	using Powerful.Models
	using Powerful.PowerCore

	case = PowerFlowData.parse_network(
		joinpath(pkgdir(Powerful), "cases", "ieee14.raw"),
	)

	transformers = parse_model(Transformer2W, case)
	@test transformers isa Transformer2WVec
	@test length(transformers) == 4
end
