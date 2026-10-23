"""GENCLS: Classical Generator Model (Swing Equation)

Hand-maintained: edit this file directly. It is not regenerated; it started as output of
a DSL compiler from `models/generator/gencls.toml`.
"""


using InlineStrings
using OrderedCollections: OrderedDict

export GENCLS

struct GENCLS{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number (external ID from DYR)"
    BUS::Int32
    "Generator ID (1-2 char alphanumeric)"
    ID::String3

    # === Parameters ===
    H::T  # Inertia constant
    D::T  # Damping coefficient
    xd_p::T  # d-axis transient reactance
    r_a::T  # Armature resistance

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function GENCLS(; BUS::Int32, ID, H, D, xd_p, r_a, status::Bool=true)
    T = typeof(H)
    GENCLS{T}(BUS, ID, T(H), T(D), T(xd_p), T(r_a), status)
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:GENCLS}) = true
get_dyr_model_name(::Type{<:GENCLS}) = "GENCLS"


# ============================================================================
# Group Interface Implementation - :generators
# ============================================================================

# Declare that GENCLS implements the :generators group interface
implements_group(::Type{<:GENCLS}) = :generators

# Map group canonical variable names to GENCLS's actual variable names
function group_var_mapping(::Type{<:GENCLS})
    Dict(
        :E_fd => :E_q, # GENCLS uses E_q
        :E_q => :E_q, # Same name
        :P_m => :P_m, # Same name
        :delta => :delta, # Same name
        :omega => :omega, # Same name
    )
end

# Map group canonical residual names to GENCLS's actual residual names
function group_res_mapping(::Type{<:GENCLS})
    Dict(
        :g_E_fd => :g_E_q, # GENCLS uses g_E_q
        :g_E_q => :g_E_q, # Same name
        :g_P_m => :g_P_m, # Same name
    )
end

# Link specification - which models this links to and how
function get_replacement_map(::Type{<:GENCLS})
    OrderedDict(
        :Slack => (
            source_fields = (:BUS, :ID,),
            target_fields = (:i, :id,),
        ),
        :PV => (
            source_fields = (:BUS, :ID,),
            target_fields = (:i, :id,),
        ),
    )
end

# Data inheritance specification - what to copy from linked models
function get_data_inheritance(::Type{<:GENCLS})
    OrderedDict(
        :Slack => OrderedDict(:xd_p => :zx, :r_a => :zr),
        :PV => OrderedDict(:xd_p => :zx, :r_a => :zr),
    )
end

# ============================================================================
# Model Registration
# ============================================================================

supports_format(::Type{<:GENCLS}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model GENCLS

const GENCLS_TERMINALS = (
    Terminal(:bus_group, :BUS),
)

get_terminals(::Type{<:GENCLS}) = GENCLS_TERMINALS

# External parameter declarations
const GENCLS_EXT_PARAMS = [
    ExtParam(:bus_basekv, GENCLS_TERMINALS[1], :basekv,
        BaseVoltage(), PrePerUnit(),
        required=true, description="Base voltage at generator bus"),
]

get_key_name(::Type{<:GENCLS}) = [:BUS, :ID]

function model_metadata(::Type{GENCLS})
    ModelMetadata(
        name = :GENCLS,
        vars = [
            ModelVar(:delta, StateVar();
                description="Rotor angle (rad)"),
            ModelVar(:omega, StateVar();
                description="Rotor speed (pu)", mass_coefficient=MassCoefficient(m -> 2.0 .* m.H)),
            ModelVar(:E_q, AlgebVar();
                description="Internal EMF (control input - receives contribution from exciter)"),
            ModelVar(:P_m, AlgebVar();
                description="Mechanical power input (control input - receives contribution from governor)"),
            ModelVar(:v, AlgebVar(), :bus_group, :v_port, :BUS;
                description="Port voltage magnitude (external): the bus voltage, or behind a bridge the voltage behind which the delivered current flows", units="pu"),
            ModelVar(:theta, AlgebVar(), :bus_group, :theta_port, :BUS;
                description="Port voltage angle (external)", units="rad")
        ],
        residuals = [
        ModelResidual(:ddelta, StateRes(), PrivateRes();
            description="Rotor angle (rad)"),
        ModelResidual(:domega, StateRes(), PrivateRes();
            description="Rotor speed (pu)"),
        ModelResidual(:g_E_q, AlgebRes(), SharedRes(:add);
            description="Internal EMF (control input - receives contribution from exciter)"),
        ModelResidual(:g_P_m, AlgebRes(), SharedRes(:add);
            description="Mechanical power input (control input - receives contribution from governor)"),
        ModelResidual(:p_inj, AlgebRes(), SharedRes(:add), :bus_group, :p, :BUS;
            description="Active power injection into bus_group (negative of P_e)"),
        ModelResidual(:q_inj, AlgebRes(), SharedRes(:add), :bus_group, :q, :BUS;
            description="Reactive power injection into bus_group (negative of Q_e)"),
        ModelResidual(:ip_inj, AlgebRes(), SharedRes(:add), :bus_group, :ip, :BUS;
            description="Active current injection into bus_group, in phase with the bus voltage (system base)"),
        ModelResidual(:iq_inj, AlgebRes(), SharedRes(:add), :bus_group, :iq, :BUS;
            description="Reactive current injection into bus_group, in quadrature to the bus voltage (system base)"),
        ModelResidual(:jr_inj, AlgebRes(), SharedRes(:add), :bus_group, :jr, :BUS;
            description="Short-circuit current E/Z into bus_group, real part in the network frame (system base): the voltage form"),
        ModelResidual(:ji_inj, AlgebRes(), SharedRes(:add), :bus_group, :ji, :BUS;
            description="Short-circuit current E/Z into bus_group, imaginary part in the network frame (system base)")
    ],
        ratings = [ModelRating(:mbase, NominalPower; from = LinkedDevice(),
            description="MVA base of the generator this machine replaces; its parameters, states, currents and mechanical power are per unit on it")],
        derived_params = [RatingRatio(:Mratio, NominalPower;
            description="Machine base over the system base (Sn/Sb), where the machine meets the network")],
        ext_params = GENCLS_EXT_PARAMS,
        # An EMF behind its subtransient (transient) impedance: the voltage form, J = E/Z; a bus that does not
        # realize it takes the current, I = (E - V)/Z, or the power p = v ip, q = v iq
        injection = InjectionForms(:voltage, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj], voltage = [:jr_inj, :ji_inj]),
        init_inputs = [
        InitInput(:V_bus, :from_solution, :Bus, :v;
            description="Bus voltage magnitude"),
        InitInput(:theta_bus, :from_solution, :Bus, :theta;
            description="Bus voltage angle"),
        InitInput(:p_gen_generator, :from_observation, :Generator, :p_gen;
            description="Active power from linked PV/Slack generator"),
        InitInput(:q_gen_generator, :from_observation, :Generator, :q_gen;
            description="Reactive power from linked PV/Slack generator")
    ],
        observable_vars = [:P_e, :Q_e],
        initialized_params = [
        InitParam(:E_q0, :E_q; description="Internal EMF (control input - receives contribution from exciter)"),
        InitParam(:P_m0, :P_m; description="Mechanical power input (control input - receives contribution from governor)")
    ],
        status_field = :status,
    )
end

@inline function gencls_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, ::Type{AT}
) where AT
    omega_base = params.omega_base
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        zd_p_sq = params.r_a[i] * params.r_a[i] + params.xd_p[i] * params.xd_p[i]

        vd = vars.v[i] * sin(vars.delta[i] - vars.theta[i])
        vq = vars.v[i] * cos(vars.delta[i] - vars.theta[i])
        Id = (params.xd_p[i] * (vars.E_q[i] - vq) - params.r_a[i] * vd) / zd_p_sq
        Iq = (params.r_a[i] * (vars.E_q[i] - vq) + params.xd_p[i] * vd) / zd_p_sq
        P_e = vd * Id + vq * Iq
        pe = P_e
        # === Residual Assignments ===
        outputs.ddelta[i] = mask * (omega_base * (vars.omega[i] - 1.0))
        outputs.domega[i] = mask * (vars.P_m[i] - P_e - params.D[i] * (vars.omega[i] - 1.0))
        outputs.g_E_q[i] = mask * (params.E_q0[i] - vars.E_q[i])
        outputs.g_P_m[i] = mask * (params.P_m0[i] - vars.P_m[i])
        outputs.P_e[i] = mask * (params.Mratio[i] * pe)
        outputs.Q_e[i] = mask * (params.Mratio[i] * (vq * Id - vd * Iq))
        # Sign convention: positive = injection into bus
        outputs.p_inj[i] = mask * (params.Mratio[i] * pe)
        outputs.q_inj[i] = mask * (params.Mratio[i] * (vq * Id - vd * Iq))
        # The current port: I = (Id + j Iq) e^{j(delta - pi/2)} = (ip - j iq) e^{j theta}, so
        # v ip = P_e, v iq = Q_e (derived with SymPy)
        sa, ca = sincos(vars.delta[i] - vars.theta[i])
        outputs.ip_inj[i] = mask * (params.Mratio[i] * (Id * sa + Iq * ca))
        outputs.iq_inj[i] = mask * (params.Mratio[i] * (Id * ca - Iq * sa))
        # The voltage form: the EMF E' = E_q e^{j delta} behind
        # Z = r + jX, as its short-circuit current J = Mratio E/Z on the system base (derived with SymPy)
        sd, cd = sincos(vars.delta[i])
        Er, Ei = vars.E_q[i] * cd, vars.E_q[i] * sd
        zz = params.r_a[i] * params.r_a[i] + params.xd_p[i] * params.xd_p[i]
        outputs.jr_inj[i] = mask * (params.Mratio[i] * (Er * params.r_a[i] + Ei * params.xd_p[i]) / zz)
        outputs.ji_inj[i] = mask * (params.Mratio[i] * (Ei * params.r_a[i] - Er * params.xd_p[i]) / zz)
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:GENCLS},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    omega_base = ws.system.omega_base
    gencls_residual_kernel!(
        (ddelta = ws.res_vals.ddelta, domega = ws.res_vals.domega, g_E_q = ws.res_vals.g_E_q, g_P_m = ws.res_vals.g_P_m, P_e = ws.res_vals.P_e, Q_e = ws.res_vals.Q_e, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj, ip_inj = ws.res_vals.ip_inj, iq_inj = ws.res_vals.iq_inj, jr_inj = ws.res_vals.jr_inj, ji_inj = ws.res_vals.ji_inj),
        (H = model.H, D = model.D, xd_p = model.xd_p, r_a = model.r_a, E_q0 = ws.setpoints.E_q0, P_m0 = ws.setpoints.P_m0, omega_base = omega_base, Mratio = external.Mratio),
        (delta = ws.var_vals.delta, omega = ws.var_vals.omega, E_q = ws.var_vals.E_q, P_m = ws.var_vals.P_m, v = ws.var_vals.v, theta = ws.var_vals.theta),
        active_mask,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{GENCLS{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, GENCLS)

    v = data.V_bus
    theta = data.theta_bus
    # The network powers are per unit on the system base; the machine is on its own base
    P_gen = data.p_gen_generator ./ data.Mratio
    Q_gen = data.q_gen_generator ./ data.Mratio

    # === Pre-allocate Result Arrays ===
    delta = Vector{T}(undef, n)
    omega = fill(one(T), n)
    E_q = Vector{T}(undef, n)
    P_m = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # Compute init-only services
        Ir = (P_gen[i] * cos(theta[i]) + Q_gen[i] * sin(theta[i])) / v[i]
        Ii = (P_gen[i] * sin(theta[i]) - Q_gen[i] * cos(theta[i])) / v[i]
        Er = v[i] * cos(theta[i]) + model.r_a[i] * Ir - model.xd_p[i] * Ii
        Ei = v[i] * sin(theta[i]) + model.xd_p[i] * Ir + model.r_a[i] * Ii
        delta0 = atan(Ei, Er)
        vd0 = v[i] * sin(delta0 - theta[i])
        vq0 = v[i] * cos(delta0 - theta[i])
        Id0 = (P_gen[i] * sin(delta0 - theta[i]) + Q_gen[i] * cos(delta0 - theta[i])) / v[i]
        Iq0 = (P_gen[i] * cos(delta0 - theta[i]) - Q_gen[i] * sin(delta0 - theta[i])) / v[i]
        E_q_init = sqrt(Er * Er + Ei * Ei)

        # State initial values
        delta[i] = delta0

        # Algebraic initial values
        E_q[i] = E_q_init
        P_m[i] = P_gen[i]
    end

    # === Return Initial Values ===
    return (delta=delta, omega=omega, E_q=E_q, P_m=P_m)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{GENCLS}, dyr_data, sys_pf, links) -> StructArray{GENCLS}

Parse GENCLS model from DYR data with linking to power flow models.
"""
function parse_dynamic_model(
    ::Type{GENCLS},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    links::ModelLinks
)
    n = length(links)  # One instance per linked DYR row
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{GENCLS{T}}(undef, n)

    # Get data inheritance mapping
    data_inheritance = get_data_inheritance(GENCLS)

    for idx in 1:n
        dyr_row = dyr_data[links.source_idx[idx]]
        pf_model = sys_pf.models[links.target_type[idx]][links.target_idx[idx]]

        # Get inheritance mapping for this target type
        inherit_map = data_inheritance[links.target_type[idx]]

        # Build inherited fields
        inherited_fields = Dict{Symbol, Any}()
        for (dst_field, src_field) in inherit_map
            inherited_fields[dst_field] = T(getproperty(pf_model, src_field))
        end

        instances[idx] = GENCLS(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            H = T(dyr_row.H),
            D = T(dyr_row.D),
            xd_p = haskey(inherited_fields, :xd_p) ? inherited_fields[:xd_p] : T(dyr_row.xd_p),
            r_a = haskey(inherited_fields, :r_a) ? inherited_fields[:r_a] : T(dyr_row.r_a)
        )
    end

    return StructArray(instances)
end