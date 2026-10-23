"""GENTPJ1: Round Rotor Generator with Quadratic Saturation and Kis term

Hand-maintained: edit this file directly. It is not regenerated; it started as output of
a DSL compiler from `models/generator/gentpj1.toml`.
"""


using InlineStrings
using OrderedCollections: OrderedDict

export GENTPJ1

struct GENTPJ1{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number where generator is connected"
    BUS::Int32
    "Generator identifier (1-2 characters)"
    ID::String3

    # === Parameters ===
    Td10::T  # d-axis transient open-circuit time constant
    Td20::T  # d-axis subtransient open-circuit time constant
    Tq10::T  # q-axis transient open-circuit time constant
    Tq20::T  # q-axis subtransient open-circuit time constant
    H::T  # Inertia constant
    D::T  # Damping coefficient
    Xd::T  # d-axis synchronous reactance
    Xq::T  # q-axis synchronous reactance
    Xd1::T  # d-axis transient reactance
    Xq1::T  # q-axis transient reactance
    Xd2::T  # d-axis subtransient reactance
    Xq2::T  # q-axis subtransient reactance
    Xl::T  # Leakage reactance
    S10::T  # Saturation factor at 1.0 pu flux
    S12::T  # Saturation factor at 1.2 pu flux
    Kis::T  # Saturation input constant
    ra::T  # Armature resistance

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function GENTPJ1(; BUS::Int32, ID, Td10, Td20, Tq10, Tq20, H, D, Xd, Xq, Xd1, Xq1, Xd2, Xq2, Xl, S10, S12, Kis, ra, status::Bool=true)
    T = typeof(Td10)
    GENTPJ1{T}(BUS, ID, T(Td10), T(Td20), T(Tq10), T(Tq20), T(H), T(D), T(Xd), T(Xq), T(Xd1), T(Xq1), T(Xd2), T(Xq2), T(Xl), T(S10), T(S12), T(Kis), T(ra), status)
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:GENTPJ1}) = true
get_dyr_model_name(::Type{<:GENTPJ1}) = "GENTPJ1"


# ============================================================================
# Group Interface Implementation - :generators
# ============================================================================

# Declare that GENTPJ1 implements the :generators group interface
implements_group(::Type{<:GENTPJ1}) = :generators

# Map group canonical variable names to GENTPJ1's actual variable names
function group_var_mapping(::Type{<:GENTPJ1})
    Dict(
        :E_fd => :E_fd, # Same name
        :P_m => :P_m, # Same name
        :delta => :delta, # Same name
        :e1d => :e1d, # Same name
        :e1q => :e1q, # Same name
        :e2d => :e2d, # Same name
        :e2q => :e2q, # Same name
        :omega => :omega, # Same name
    )
end

# Map group canonical residual names to GENTPJ1's actual residual names
function group_res_mapping(::Type{<:GENTPJ1})
    Dict(
        :g_E_fd => :g_E_fd, # Same name
        :g_P_m => :g_P_m, # Same name
    )
end

# Link specification - which models this links to and how
function get_replacement_map(::Type{<:GENTPJ1})
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
function get_data_inheritance(::Type{<:GENTPJ1})
    OrderedDict(
        :Slack => OrderedDict(:ra => :zr),
        :PV => OrderedDict(:ra => :zr),
    )
end

# ============================================================================
# Model Registration
# ============================================================================


supports_format(::Type{<:GENTPJ1}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model GENTPJ1

# Terminal declarations - where GENTPJ1 connects
# Use :bus_group (not :Bus) so GENTPJ1 connects to any bus model implementing it (e.g. Bus)
const GENTPJ1_TERMINALS = (
    Terminal(:bus_group, :BUS),
)

get_terminals(::Type{<:GENTPJ1}) = GENTPJ1_TERMINALS


get_key_name(::Type{<:GENTPJ1}) = [:BUS, :ID]

function model_metadata(::Type{GENTPJ1})
    ModelMetadata(
        name = :GENTPJ1,
        vars = [
            ModelVar(:delta, StateVar();
                description="Rotor angle (rad)"),
            ModelVar(:omega, StateVar();
                description="Rotor speed (pu)", mass_coefficient=MassCoefficient(m -> 2.0 .* m.H)),
            ModelVar(:e1q, StateVar();
                description="q-axis transient voltage"),
            ModelVar(:e1d, StateVar();
                description="d-axis transient voltage"),
            ModelVar(:e2d, StateVar();
                description="d-axis subtransient voltage"),
            ModelVar(:e2q, StateVar();
                description="q-axis subtransient voltage"),
            ModelVar(:Id, AlgebVar();
                description="d-axis current"),
            ModelVar(:Iq, AlgebVar();
                description="q-axis current"),
            ModelVar(:E_fd, AlgebVar();
                description="Field voltage (exciter output)"),
            ModelVar(:P_m, AlgebVar();
                description="Mechanical power (governor output)"),
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
        ModelResidual(:de1d, StateRes(), PrivateRes();
            description="d-axis transient voltage"),
        ModelResidual(:de2d, StateRes(), PrivateRes();
            description="d-axis subtransient voltage"),
        ModelResidual(:de2q, StateRes(), PrivateRes();
            description="q-axis subtransient voltage"),
        ModelResidual(:g_Id, AlgebRes(), PrivateRes();
            description="d-axis current"),
        ModelResidual(:g_Iq, AlgebRes(), PrivateRes();
            description="q-axis current"),
        ModelResidual(:g_E_fd, AlgebRes(), SharedRes(:add);
            description="Field voltage (exciter output)"),
        ModelResidual(:g_P_m, AlgebRes(), SharedRes(:add);
            description="Mechanical power (governor output)"),
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
        InitParam(:E_fd0, :E_fd; description="Field voltage (exciter output)"),
        InitParam(:P_m0, :P_m; description="Mechanical power (governor output)")
    ],
        status_field = :status,
        evaluation_pass = 1,
    )
end

@inline function gentpj1_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, ::Type{AT}
) where AT
    omega_base = params.omega_base
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        gd1 = (params.Xd2[i] - params.Xl[i]) / (params.Xd1[i] - params.Xl[i])
        gq1 = (params.Xq2[i] - params.Xl[i]) / (params.Xq1[i] - params.Xl[i])
        gd2 = (params.Xd1[i] - params.Xd2[i]) / ((params.Xd1[i] - params.Xl[i]) * (params.Xd1[i] - params.Xl[i]))
        gq2 = (params.Xq1[i] - params.Xq2[i]) / ((params.Xq1[i] - params.Xl[i]) * (params.Xq1[i] - params.Xl[i]))
        gqd = (params.Xq[i] - params.Xl[i]) / (params.Xd[i] - params.Xl[i])
        Se_q_ratio = params.Xq[i] / params.Xd[i]
        SAT0_S12_safe = ifelse(params.S12[i] < 1e-10, 1.0, params.S12[i])
        SAT0_ratio = sqrt(params.S10[i] * 1.2 / SAT0_S12_safe)
        SAT0_sqrt12 = sqrt(1.2)
        SAT0_denom = ifelse(abs(SAT0_ratio - SAT0_sqrt12) < 1e-10, 1e-10, SAT0_ratio - SAT0_sqrt12)
        SAT0_A = SAT0_sqrt12 * (SAT0_sqrt12 * SAT0_ratio - 1.0) / SAT0_denom
        SAT0_B = params.S10[i] / max((1.0 - SAT0_A) * (1.0 - SAT0_A), 1e-10)

        vd = vars.v[i] * sin(vars.delta[i] - vars.theta[i])
        vq = vars.v[i] * cos(vars.delta[i] - vars.theta[i])
        psi2d = gd1 * vars.e1q[i] + gd2 * (params.Xd1[i] - params.Xl[i]) * vars.e2d[i]
        psi2q = gq1 * vars.e1d[i] + (1 - gq1) * vars.e2q[i]
        psi2 = sqrt(psi2d * psi2d + psi2q * psi2q)
        It_mag = sqrt(vars.Id[i] * vars.Id[i] + vars.Iq[i] * vars.Iq[i])
        psi_sat_input = psi2 + params.Kis[i] * It_mag
        # --- SAT: QuadraticSaturation ---
        SAT_S12_safe = ifelse((params.S12[i]) < 1e-10, 1.0, (params.S12[i]))
        SAT_ratio = sqrt((params.S10[i]) * 1.2 / SAT_S12_safe)
        SAT_sqrt12 = sqrt(1.2)
        SAT_denom = ifelse(abs(SAT_ratio - SAT_sqrt12) < 1e-10, 1e-10, SAT_ratio - SAT_sqrt12)
        SAT_SAT_A = SAT_sqrt12 * (SAT_sqrt12 * SAT_ratio - 1.0) / SAT_denom
        SAT_SAT_B = (params.S10[i]) / max((1.0 - SAT_SAT_A) * (1.0 - SAT_SAT_A), 1e-10)
        SAT_y = ifelse((psi_sat_input) >= SAT_SAT_A, ((psi_sat_input) - SAT_SAT_A) * ((psi_sat_input) - SAT_SAT_A) * SAT_SAT_B / (psi_sat_input), 0.0)
        SAT_SAT_A = SAT_SAT_A
        SAT_SAT_B = SAT_SAT_B

        Se = SAT_y
        Se_q = Se * Se_q_ratio
        det = params.ra[i] * params.ra[i] + params.Xd2[i] * params.Xq2[i]
        Id_calc = (params.ra[i] * (psi2q - vd) + params.Xq2[i] * (psi2d - vq)) / det
        Iq_calc = (params.ra[i] * (psi2d - vq) - params.Xd2[i] * (psi2q - vd)) / det
        XadIfd = vars.e1q[i] + (params.Xd[i] - params.Xd1[i]) * (gd1 * vars.Id[i] - gd2 * vars.e2d[i] + gd2 * vars.e1q[i]) + Se * psi2d
        XaqI1q = vars.e1d[i] + (params.Xq[i] - params.Xq1[i]) * (gq2 * vars.e1d[i] - gq2 * vars.e2q[i] - gq1 * vars.Iq[i]) + Se_q * psi2q
        te = psi2d * vars.Iq[i] + psi2q * vars.Id[i]
        Pe = vd * vars.Id[i] + vq * vars.Iq[i]
        Qe = vq * vars.Id[i] - vd * vars.Iq[i]
        # === Residual Assignments ===
        outputs.ddelta[i] = mask * (omega_base * (vars.omega[i] - 1.0))
        outputs.domega[i] = mask * (vars.P_m[i] - te - params.D[i] * (vars.omega[i] - 1.0))
        outputs.de1q[i] = mask * ((-XadIfd + vars.E_fd[i]) / params.Td10[i])
        outputs.de1d[i] = mask * (-XaqI1q / params.Tq10[i])
        outputs.de2d[i] = mask * ((-vars.e2d[i] + vars.e1q[i] - (params.Xd1[i] - params.Xl[i]) * vars.Id[i]) / params.Td20[i])
        outputs.de2q[i] = mask * ((-vars.e2q[i] + vars.e1d[i] + (params.Xq1[i] - params.Xl[i]) * vars.Iq[i]) / params.Tq20[i])
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
        # The voltage form: the EMF E'' = (psi2q + j psi2d) e^{j(delta - pi/2)} (Xq2 = Xd2: one impedance) behind
        # Z = r + jX, as its short-circuit current J = Mratio E/Z on the system base (derived with SymPy)
        sd, cd = sincos(vars.delta[i])
        Er, Ei = psi2d * cd + psi2q * sd, psi2d * sd - psi2q * cd
        zz = params.ra[i] * params.ra[i] + params.Xd2[i] * params.Xd2[i]
        outputs.jr_inj[i] = mask * (params.Mratio[i] * (Er * params.ra[i] + Ei * params.Xd2[i]) / zz)
        outputs.ji_inj[i] = mask * (params.Mratio[i] * (Ei * params.ra[i] - Er * params.Xd2[i]) / zz)
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:GENTPJ1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    omega_base = ws.system.omega_base
    gentpj1_residual_kernel!(
        (ddelta = ws.res_vals.ddelta, domega = ws.res_vals.domega, de1q = ws.res_vals.de1q, de1d = ws.res_vals.de1d, de2d = ws.res_vals.de2d, de2q = ws.res_vals.de2q, g_Id = ws.res_vals.g_Id, g_Iq = ws.res_vals.g_Iq, g_E_fd = ws.res_vals.g_E_fd, g_P_m = ws.res_vals.g_P_m, Pe = ws.res_vals.Pe, Qe = ws.res_vals.Qe, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj, ip_inj = ws.res_vals.ip_inj, iq_inj = ws.res_vals.iq_inj, jr_inj = ws.res_vals.jr_inj, ji_inj = ws.res_vals.ji_inj),
        (Td10 = model.Td10, Td20 = model.Td20, Tq10 = model.Tq10, Tq20 = model.Tq20, H = model.H, D = model.D, Xd = model.Xd, Xq = model.Xq, Xd1 = model.Xd1, Xq1 = model.Xq1, Xd2 = model.Xd2, Xq2 = model.Xq2, Xl = model.Xl, S10 = model.S10, S12 = model.S12, Kis = model.Kis, ra = model.ra, E_fd0 = ws.setpoints.E_fd0, P_m0 = ws.setpoints.P_m0, omega_base = omega_base, Mratio = external.Mratio),
        (delta = ws.var_vals.delta, omega = ws.var_vals.omega, e1q = ws.var_vals.e1q, e1d = ws.var_vals.e1d, e2d = ws.var_vals.e2d, e2q = ws.var_vals.e2q, Id = ws.var_vals.Id, Iq = ws.var_vals.Iq, E_fd = ws.var_vals.E_fd, P_m = ws.var_vals.P_m, v = ws.var_vals.v, theta = ws.var_vals.theta),
        active_mask,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{GENTPJ1{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, GENTPJ1)

    v = data.V_bus
    theta = data.theta_bus
    # The network powers are per unit on the system base; the machine is on its own base
    P_gen = data.p_gen_generator ./ data.Mratio
    Q_gen = data.q_gen_generator ./ data.Mratio

    # === Pre-allocate Result Arrays ===
    delta = Vector{T}(undef, n)
    omega = fill(one(T), n)
    e1q = Vector{T}(undef, n)
    e1d = Vector{T}(undef, n)
    e2d = Vector{T}(undef, n)
    e2q = Vector{T}(undef, n)
    Id = Vector{T}(undef, n)
    Iq = Vector{T}(undef, n)
    E_fd = Vector{T}(undef, n)
    P_m = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # Compute init-only services
        Se_q_ratio = model.Xq[i] / model.Xd[i]
        SAT0_S12_safe = ifelse(model.S12[i] < 1e-10, 1.0, model.S12[i])
        SAT0_ratio = sqrt(model.S10[i] * 1.2 / SAT0_S12_safe)
        SAT0_sqrt12 = sqrt(1.2)
        SAT0_denom = ifelse(abs(SAT0_ratio - SAT0_sqrt12) < 1e-10, 1e-10, SAT0_ratio - SAT0_sqrt12)
        SAT0_A = SAT0_sqrt12 * (SAT0_sqrt12 * SAT0_ratio - 1.0) / SAT0_denom
        SAT0_B = model.S10[i] / max((1.0 - SAT0_A) * (1.0 - SAT0_A), 1e-10)
        It_re = (P_gen[i] * cos(theta[i]) + Q_gen[i] * sin(theta[i])) / v[i]
        It_im = (P_gen[i] * sin(theta[i]) - Q_gen[i] * cos(theta[i])) / v[i]
        It = complex(It_re, It_im)
        Vt = v[i] * cis(theta[i])
        Kis_It = model.Kis[i] * abs(It)
        # The start from the residual kernel (derived with SymPy): the stator rows give the
        # subtransient fluxes psi2q = vd + ra Id - Xq2 Iq, psi2d = vq + ra Iq + Xd2 Id, the de2d,
        # de2q rows and the fluxes the states, and the de1d row is Im[Et e^(-j delta)] = 0 with
        # Et = (1 + Se_q) (V + (ra + j Xq2) I) + j (Xq - Xq2) I. The saturation Se = SAT(|psi2| +
        # Kis |I|) depends on delta through |psi2| unless Xq2 = Xd2: delta and Se are the fixed
        # point of the two
        Se0 = zero(T)
        delta0 = Id0 = Iq0 = psi2d0 = psi2q0 = zero(T)
        converged = false
        for _ in 1:100
            Et = (1 + Se0 * Se_q_ratio) * (Vt + complex(model.ra[i], model.Xq2[i]) * It) +
                 im * (model.Xq[i] - model.Xq2[i]) * It
            delta0 = angle(Et)
            Idq = It * im * cis(-delta0)
            Id0, Iq0 = real(Idq), imag(Idq)
            psi2d0 = v[i] * cos(delta0 - theta[i]) + model.ra[i] * Iq0 + model.Xd2[i] * Id0
            psi2q0 = v[i] * sin(delta0 - theta[i]) + model.ra[i] * Id0 - model.Xq2[i] * Iq0
            u0 = sqrt(psi2d0 * psi2d0 + psi2q0 * psi2q0) + Kis_It
            Se_next = ifelse(u0 >= SAT0_A, (u0 - SAT0_A) * (u0 - SAT0_A) * SAT0_B / u0, zero(T))
            converged = abs(Se_next - Se0) <= 1e-14
            Se0 = Se_next
            converged && break
        end
        converged || error("GENTPJ1 $(model.BUS[i]) $(model.ID[i]): the start of the saturation does not " *
                           "converge in 100 iterations (Se = $Se0)")
        tm0 = psi2d0 * Iq0 + psi2q0 * Id0
        E_fd0 = (1 + Se0) * psi2d0 + (model.Xd[i] - model.Xd2[i]) * Id0
        e1q0 = psi2d0 + (model.Xd1[i] - model.Xd2[i]) * Id0
        e1d0 = psi2q0 - (model.Xq1[i] - model.Xq2[i]) * Iq0
        e2d0 = e1q0 - (model.Xd1[i] - model.Xl[i]) * Id0
        e2q0 = e1d0 + (model.Xq1[i] - model.Xl[i]) * Iq0

        # State initial values
        delta[i] = delta0
        e1q[i] = e1q0
        e1d[i] = e1d0
        e2d[i] = e2d0
        e2q[i] = e2q0

        # Algebraic initial values
        Id[i] = Id0
        Iq[i] = Iq0
        E_fd[i] = E_fd0
        P_m[i] = tm0
    end

    # === Return Initial Values ===
    return (delta=delta, omega=omega, e1q=e1q, e1d=e1d, e2d=e2d, e2q=e2q, Id=Id, Iq=Iq, E_fd=E_fd, P_m=P_m)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{GENTPJ1}, dyr_data, sys_pf, links) -> StructArray{GENTPJ1}

Parse GENTPJ1 model from DYR data with linking to power flow models.
"""
function parse_dynamic_model(
    ::Type{GENTPJ1},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    links::ModelLinks
)
    n = length(links)  # One instance per linked DYR row
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{GENTPJ1{T}}(undef, n)

    # Get data inheritance mapping
    data_inheritance = get_data_inheritance(GENTPJ1)

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

        instances[idx] = GENTPJ1(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            Td10 = T(dyr_row.Td10),
            Td20 = T(dyr_row.Td20),
            Tq10 = T(dyr_row.Tq10),
            Tq20 = T(dyr_row.Tq20),
            H = T(dyr_row.H),
            D = T(dyr_row.D),
            Xd = T(dyr_row.Xd),
            Xq = T(dyr_row.Xq),
            Xd1 = T(dyr_row.Xd1),
            Xq1 = T(dyr_row.Xq1),
            Xd2 = T(dyr_row.Xd2),
            Xq2 = T(dyr_row.Xq2),
            Xl = T(dyr_row.Xl),
            S10 = T(dyr_row.S10),
            S12 = T(dyr_row.S12),
            Kis = T(dyr_row.Kis),
            ra = haskey(inherited_fields, :ra) ? inherited_fields[:ra] : T(dyr_row.ra)
        )
    end

    return StructArray(instances)
end