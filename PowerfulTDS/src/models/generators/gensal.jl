"""GENSAL: Salient Pole Generator with Quadratic Saturation

Hand-maintained: edit this file directly. It is not regenerated; it started as output of
a DSL compiler from `models/generator/gensal.toml`.
"""


using InlineStrings
using OrderedCollections: OrderedDict

export GENSAL

struct GENSAL{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number where generator is connected"
    BUS::Int32
    "Generator identifier (1-2 characters)"
    ID::String3

    # === Parameters ===
    Td10::T  # d-axis transient open-circuit time constant
    Td20::T  # d-axis subtransient open-circuit time constant
    Tq20::T  # q-axis subtransient open-circuit time constant
    H::T  # Inertia constant
    D::T  # Damping coefficient
    Xd::T  # d-axis synchronous reactance
    Xq::T  # q-axis synchronous reactance
    Xd1::T  # d-axis transient reactance
    Xd2::T  # d-axis subtransient reactance
    Xl::T  # Leakage reactance
    S10::T  # Saturation factor at 1.0 pu flux
    S12::T  # Saturation factor at 1.2 pu flux
    ra::T  # Armature resistance

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function GENSAL(; BUS::Int32, ID, Td10, Td20, Tq20, H, D, Xd, Xq, Xd1, Xd2, Xl, S10, S12, ra, status::Bool=true)
    T = typeof(Td10)
    GENSAL{T}(BUS, ID, T(Td10), T(Td20), T(Tq20), T(H), T(D), T(Xd), T(Xq), T(Xd1), T(Xd2), T(Xl), T(S10), T(S12), T(ra), status)
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:GENSAL}) = true
get_dyr_model_name(::Type{<:GENSAL}) = "GENSAL"


# ============================================================================
# Group Interface Implementation - :generators
# ============================================================================

# Declare that GENSAL implements the :generators group interface
implements_group(::Type{<:GENSAL}) = :generators

# Map group canonical variable names to GENSAL's actual variable names
function group_var_mapping(::Type{<:GENSAL})
    Dict(
        :E_fd => :E_fd, # Same name
        :P_m => :P_m, # Same name
        :delta => :delta, # Same name
        :e1q => :e1q, # Same name
        :omega => :omega, # Same name
        :psi1d => :psi1d, # Same name
        :psi2q => :psi2q, # Same name
    )
end

# Map group canonical residual names to GENSAL's actual residual names
function group_res_mapping(::Type{<:GENSAL})
    Dict(
        :g_E_fd => :g_E_fd, # Same name
        :g_P_m => :g_P_m, # Same name
    )
end

# Link specification - which models this links to and how
function get_replacement_map(::Type{<:GENSAL})
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
function get_data_inheritance(::Type{<:GENSAL})
    OrderedDict(
        :Slack => OrderedDict(:ra => :zr),
        :PV => OrderedDict(:ra => :zr),
    )
end

# ============================================================================
# Model Registration
# ============================================================================


supports_format(::Type{<:GENSAL}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model GENSAL

# Terminal declarations - where GENSAL connects
# Use :bus_group (not :Bus) so GENSAL connects to any bus model implementing it (e.g. Bus)
const GENSAL_TERMINALS = (
    Terminal(:bus_group, :BUS),
)

get_terminals(::Type{<:GENSAL}) = GENSAL_TERMINALS


get_key_name(::Type{<:GENSAL}) = [:BUS, :ID]

function model_metadata(::Type{GENSAL})
    ModelMetadata(
        name = :GENSAL,
        vars = [
            ModelVar(:delta, StateVar();
                description="Rotor angle (rad)"),
            ModelVar(:omega, StateVar();
                description="Rotor speed (pu)", mass_coefficient=MassCoefficient(m -> 2.0 .* m.H)),
            ModelVar(:e1q, StateVar();
                description="q-axis transient voltage"),
            ModelVar(:psi1d, StateVar();
                description="d-axis transient flux"),
            ModelVar(:psi2q, StateVar();
                description="q-axis subtransient flux"),
            ModelVar(:Id, AlgebVar();
                description="d-axis current"),
            ModelVar(:Iq, AlgebVar();
                description="q-axis current"),
            ModelVar(:E_fd, AlgebVar();
                description="Field voltage"),
            ModelVar(:P_m, AlgebVar();
                description="Mechanical power"),
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
        ModelResidual(:de1q, StateRes(), PrivateRes();
            description="q-axis transient voltage"),
        ModelResidual(:dpsi1d, StateRes(), PrivateRes();
            description="d-axis transient flux"),
        ModelResidual(:dpsi2q, StateRes(), PrivateRes();
            description="q-axis subtransient flux"),
        ModelResidual(:g_Id, AlgebRes(), PrivateRes();
            description="d-axis current"),
        ModelResidual(:g_Iq, AlgebRes(), PrivateRes();
            description="q-axis current"),
        ModelResidual(:g_E_fd, AlgebRes(), SharedRes(:add);
            description="Field voltage"),
        ModelResidual(:g_P_m, AlgebRes(), SharedRes(:add);
            description="Mechanical power"),
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
        # An EMF behind its subtransient (transient) impedance: the voltage form, J = E/Z; a bus that does not
        # realize it takes the current, I = (E - V)/Z, or the power p = v ip, q = v iq
        injection = InjectionForms(:voltage, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj], voltage = [:jr_inj, :ji_inj]),
        ext_params = ExtParam[],
        init_inputs = [
        InitInput(:V_bus, :from_solution, :Bus, :v;
            description="Bus voltage magnitude"),
        InitInput(:theta_bus, :from_solution, :Bus, :theta;
            description="Bus voltage angle"),
        InitInput(:p_gen_generator, :from_observation, :Generator, :p_gen;
            description="Active power from power flow"),
        InitInput(:q_gen_generator, :from_observation, :Generator, :q_gen;
            description="Reactive power from power flow")
    ],
        observable_vars = [:Pe, :Qe],
        initialized_params = [
        InitParam(:E_fd0, :E_fd; description="Field voltage"),
        InitParam(:P_m0, :P_m; description="Mechanical power")
    ],
        status_field = :status,
        evaluation_pass = 1,
    )
end

@inline function gensal_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, ::Type{AT}
) where AT
    omega_base = params.omega_base
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        # The salient-pole machine has no q-axis subtransient reactance of its own: X''q = X''d (PSS/E GENSAL)
        Xq2 = params.Xd2[i]
        K1 = (params.Xd2[i] - params.Xl[i]) / (params.Xd1[i] - params.Xl[i])
        K2 = (params.Xd1[i] - params.Xd2[i]) / (params.Xd1[i] - params.Xl[i])
        SAT0_S12_safe = ifelse(params.S12[i] < 1e-10, 1.0, params.S12[i])
        SAT0_ratio = sqrt(params.S10[i] * 1.2 / SAT0_S12_safe)
        SAT0_sqrt12 = sqrt(1.2)
        SAT0_denom = ifelse(abs(SAT0_ratio - SAT0_sqrt12) < 1e-10, 1e-10, SAT0_ratio - SAT0_sqrt12)
        SAT0_A = SAT0_sqrt12 * (SAT0_sqrt12 * SAT0_ratio - 1.0) / SAT0_denom
        SAT0_B = params.S10[i] / max((1.0 - SAT0_A) * (1.0 - SAT0_A), 1e-10)

        vd = vars.v[i] * sin(vars.delta[i] - vars.theta[i])
        vq = vars.v[i] * cos(vars.delta[i] - vars.theta[i])
        psi2d = K1 * vars.e1q[i] + K2 * vars.psi1d[i]
        # --- SAT: QuadraticSaturation ---
        SAT_S12_safe = ifelse((params.S12[i]) < 1e-10, 1.0, (params.S12[i]))
        SAT_ratio = sqrt((params.S10[i]) * 1.2 / SAT_S12_safe)
        SAT_sqrt12 = sqrt(1.2)
        SAT_denom = ifelse(abs(SAT_ratio - SAT_sqrt12) < 1e-10, 1e-10, SAT_ratio - SAT_sqrt12)
        SAT_SAT_A = SAT_sqrt12 * (SAT_sqrt12 * SAT_ratio - 1.0) / SAT_denom
        SAT_SAT_B = (params.S10[i]) / max((1.0 - SAT_SAT_A) * (1.0 - SAT_SAT_A), 1e-10)
        SAT_y = ifelse((vars.e1q[i]) >= SAT_SAT_A, ((vars.e1q[i]) - SAT_SAT_A) * ((vars.e1q[i]) - SAT_SAT_A) * SAT_SAT_B / (vars.e1q[i]), 0.0)
        SAT_SAT_A = SAT_SAT_A
        SAT_SAT_B = SAT_SAT_B

        Se = SAT_y
        det = params.ra[i] * params.ra[i] + params.Xd2[i] * Xq2
        Id_calc = (params.ra[i] * (-vars.psi2q[i] - vd) + Xq2 * (psi2d - vq)) / det
        Iq_calc = (params.ra[i] * (psi2d - vq) - params.Xd2[i] * (-vars.psi2q[i] - vd)) / det
        te = psi2d * vars.Iq[i] - vars.psi2q[i] * vars.Id[i] + vars.Id[i] * vars.Iq[i] * (Xq2 - params.Xd2[i])
        Pe = vd * vars.Id[i] + vq * vars.Iq[i]
        Qe = vq * vars.Id[i] - vd * vars.Iq[i]
        # === Residual Assignments ===
        outputs.ddelta[i] = mask * (omega_base * (vars.omega[i] - 1.0))
        outputs.domega[i] = mask * (vars.P_m[i] - te - params.D[i] * (vars.omega[i] - 1.0))
        outputs.de1q[i] = mask * ((vars.E_fd[i] - vars.e1q[i] * (1.0 + Se) - (params.Xd[i] - params.Xd1[i]) * vars.Id[i]) / params.Td10[i])
        outputs.dpsi1d[i] = mask * ((vars.e1q[i] - vars.psi1d[i] - (params.Xd1[i] - params.Xl[i]) * vars.Id[i]) / params.Td20[i])
        outputs.dpsi2q[i] = mask * ((-vars.psi2q[i] - (params.Xq[i] - Xq2) * vars.Iq[i]) / params.Tq20[i])
        outputs.g_Id[i] = mask * (Id_calc - vars.Id[i])
        outputs.g_Iq[i] = mask * (Iq_calc - vars.Iq[i])
        outputs.g_E_fd[i] = mask * (params.E_fd0[i] - vars.E_fd[i])
        outputs.g_P_m[i] = mask * (params.P_m0[i] - vars.P_m[i])
        outputs.Pe[i] = mask * (params.Mratio[i] * Pe)
        outputs.Qe[i] = mask * (params.Mratio[i] * Qe)
        # Sign convention: positive = injection into bus
        outputs.p_inj[i] = mask * (params.Mratio[i] * Pe)
        outputs.q_inj[i] = mask * (params.Mratio[i] * Qe)
        # The current port: I = (Id + j Iq) e^{j(delta - pi/2)} = (ip - j iq) e^{j theta}, so
        # v ip = Pe, v iq = Qe (derived with SymPy)
        sa, ca = sincos(vars.delta[i] - vars.theta[i])
        outputs.ip_inj[i] = mask * (params.Mratio[i] * (vars.Id[i] * sa + vars.Iq[i] * ca))
        outputs.iq_inj[i] = mask * (params.Mratio[i] * (vars.Id[i] * ca - vars.Iq[i] * sa))
        # The voltage form: the EMF E'' = (-psi2q + j psi2d) e^{j(delta - pi/2)} behind
        # Z = r + jX, as its short-circuit current J = Mratio E/Z on the system base (derived with SymPy)
        sd, cd = sincos(vars.delta[i])
        Er, Ei = psi2d * cd - vars.psi2q[i] * sd, psi2d * sd + vars.psi2q[i] * cd
        zz = params.ra[i] * params.ra[i] + params.Xd2[i] * params.Xd2[i]
        outputs.jr_inj[i] = mask * (params.Mratio[i] * (Er * params.ra[i] + Ei * params.Xd2[i]) / zz)
        outputs.ji_inj[i] = mask * (params.Mratio[i] * (Ei * params.ra[i] - Er * params.Xd2[i]) / zz)
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:GENSAL},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    omega_base = ws.system.omega_base
    gensal_residual_kernel!(
        (ddelta = ws.res_vals.ddelta, domega = ws.res_vals.domega, de1q = ws.res_vals.de1q, dpsi1d = ws.res_vals.dpsi1d, dpsi2q = ws.res_vals.dpsi2q, g_Id = ws.res_vals.g_Id, g_Iq = ws.res_vals.g_Iq, g_E_fd = ws.res_vals.g_E_fd, g_P_m = ws.res_vals.g_P_m, Pe = ws.res_vals.Pe, Qe = ws.res_vals.Qe, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj, ip_inj = ws.res_vals.ip_inj, iq_inj = ws.res_vals.iq_inj, jr_inj = ws.res_vals.jr_inj, ji_inj = ws.res_vals.ji_inj),
        (Td10 = model.Td10, Td20 = model.Td20, Tq20 = model.Tq20, H = model.H, D = model.D, Xd = model.Xd, Xq = model.Xq, Xd1 = model.Xd1, Xd2 = model.Xd2, Xl = model.Xl, S10 = model.S10, S12 = model.S12, ra = model.ra, E_fd0 = ws.setpoints.E_fd0, P_m0 = ws.setpoints.P_m0, omega_base = omega_base, Mratio = external.Mratio),
        (delta = ws.var_vals.delta, omega = ws.var_vals.omega, e1q = ws.var_vals.e1q, psi1d = ws.var_vals.psi1d, psi2q = ws.var_vals.psi2q, Id = ws.var_vals.Id, Iq = ws.var_vals.Iq, E_fd = ws.var_vals.E_fd, P_m = ws.var_vals.P_m, v = ws.var_vals.v, theta = ws.var_vals.theta),
        active_mask,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{GENSAL{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, GENSAL)

    v = data.V_bus
    theta = data.theta_bus
    # The network powers are per unit on the system base; the machine is on its own base
    P_gen = data.p_gen_generator ./ data.Mratio
    Q_gen = data.q_gen_generator ./ data.Mratio

    # === Pre-allocate Result Arrays ===
    delta = Vector{T}(undef, n)
    omega = fill(one(T), n)
    e1q = Vector{T}(undef, n)
    psi1d = Vector{T}(undef, n)
    psi2q = Vector{T}(undef, n)
    Id = Vector{T}(undef, n)
    Iq = Vector{T}(undef, n)
    E_fd = Vector{T}(undef, n)
    P_m = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # Compute init-only services
        Xq2 = model.Xd2[i]
        SAT0_S12_safe = ifelse(model.S12[i] < 1e-10, 1.0, model.S12[i])
        SAT0_ratio = sqrt(model.S10[i] * 1.2 / SAT0_S12_safe)
        SAT0_sqrt12 = sqrt(1.2)
        SAT0_denom = ifelse(abs(SAT0_ratio - SAT0_sqrt12) < 1e-10, 1e-10, SAT0_ratio - SAT0_sqrt12)
        SAT0_A = SAT0_sqrt12 * (SAT0_sqrt12 * SAT0_ratio - 1.0) / SAT0_denom
        SAT0_B = model.S10[i] / max((1.0 - SAT0_A) * (1.0 - SAT0_A), 1e-10)
        It_re = (P_gen[i] * cos(theta[i]) + Q_gen[i] * sin(theta[i])) / v[i]
        It_im = (P_gen[i] * sin(theta[i]) - Q_gen[i] * cos(theta[i])) / v[i]
        It_abs = sqrt(It_re * It_re + It_im * It_im)
        Eq_re = v[i] * cos(theta[i]) + model.ra[i] * It_re - model.Xq[i] * It_im
        Eq_im = v[i] * sin(theta[i]) + model.ra[i] * It_im + model.Xq[i] * It_re
        delta0 = atan(Eq_im, Eq_re)
        Id0 = It_re * sin(delta0) - It_im * cos(delta0)
        Iq0 = It_re * cos(delta0) + It_im * sin(delta0)
        Vd0 = v[i] * sin(delta0 - theta[i])
        Vq0 = v[i] * cos(delta0 - theta[i])
        psi2q0 = -Vd0 - model.ra[i] * Id0 + Xq2 * Iq0
        psi2d0 = Vq0 + model.ra[i] * Iq0 + model.Xd2[i] * Id0
        e1q0 = psi2d0 + (model.Xd1[i] - model.Xd2[i]) * Id0
        psi1d0 = e1q0 - (model.Xd1[i] - model.Xl[i]) * Id0
        Se0 = ifelse(abs(e1q0) >= SAT0_A, (abs(e1q0) - SAT0_A) * (abs(e1q0) - SAT0_A) * SAT0_B / abs(e1q0), 0.0)
        E_fd0 = e1q0 * (1.0 + Se0) + (model.Xd[i] - model.Xd1[i]) * Id0
        tm0 = psi2d0 * Iq0 - psi2q0 * Id0

        # State initial values
        delta[i] = delta0
        e1q[i] = e1q0
        psi1d[i] = psi1d0
        psi2q[i] = psi2q0

        # Algebraic initial values
        Id[i] = Id0
        Iq[i] = Iq0
        E_fd[i] = E_fd0
        P_m[i] = tm0
    end

    # === Return Initial Values ===
    return (delta=delta, omega=omega, e1q=e1q, psi1d=psi1d, psi2q=psi2q, Id=Id, Iq=Iq, E_fd=E_fd, P_m=P_m)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{GENSAL}, dyr_data, sys_pf, links) -> StructArray{GENSAL}

Parse GENSAL model from DYR data with linking to power flow models.
"""
function parse_dynamic_model(
    ::Type{GENSAL},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    links::ModelLinks
)
    n = length(links)  # One instance per linked DYR row
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{GENSAL{T}}(undef, n)

    # Get data inheritance mapping
    data_inheritance = get_data_inheritance(GENSAL)

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

        instances[idx] = GENSAL(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            Td10 = T(dyr_row.Td10),
            Td20 = T(dyr_row.Td20),
            Tq20 = T(dyr_row.Tq20),
            H = T(dyr_row.H),
            D = T(dyr_row.D),
            Xd = T(dyr_row.Xd),
            Xq = T(dyr_row.Xq),
            Xd1 = T(dyr_row.Xd1),
            Xd2 = T(dyr_row.Xd2),
            Xl = T(dyr_row.Xl),
            S10 = T(dyr_row.S10),
            S12 = T(dyr_row.S12),
            ra = haskey(inherited_fields, :ra) ? inherited_fields[:ra] : T(dyr_row.ra)
        )
    end

    return StructArray(instances)
end