"""REGCA1: Renewable Energy Generator/Converter Model A (REGC_A)

DO NOT EDIT - changes will be overwritten.
"""


using InlineStrings
using PowerfulCore: Automaton, Mode, Transition, InitialMode, Condition, Atom, AllOf, AnyOf, Negation

export REGCA1
export AbstractREGCA1

abstract type AbstractREGCA1 <: AbstractModel end

struct REGCA1{T} <: AbstractREGCA1
    # ===== Indexing Fields =====
    "Bus number where the renewable generator is connected"
    BUS::Int32
    "Generator identifier (1-2 characters)"
    ID::String3

    # === Parameters ===
    Lvplsw::T  # Low voltage power logic switch: 1 LVPL acts, 0 it does not
    Tg::T  # Converter time constant of both current lags (s0, s1)
    Rrpwr::T  # Up-ramp rate limit of the active current s1 (pu/s), on every rise of s1 while LVPL does not hold it
    Brkpt::T  # LVPL breakpoint: above it LVPL does not limit s1
    Zerox::T  # LVPL zero crossing: below it LVPL holds s1 at zero
    Lvpl1::T  # LVPL at Brkpt
    Volim::T  # Vtmax: the voltage the high-voltage reactive current management keeps Vt at or below
    Lvpnt1::T  # Low-voltage active current management: voltage above which all of s1 passes
    Lvpnt0::T  # Low-voltage active current management: voltage below which none of s1 passes
    Iolim::T  # qmin: the floor of the reactive power q of the high-voltage reactive current management
    Tfltr::T  # Time constant of the voltage filter s2 that LVPL reads
    Khv::T  # accel: acceleration factor of the iterative network solution of the high-voltage reactive current management; kept, it has no effect on the simultaneous solution
    Iqrmax::T  # Up rate limit of the reactive current s0 (pu/s), acting where Qgen0 > 0
    Iqrmin::T  # Down rate limit of the reactive current s0 (pu/s), acting where Qgen0 < 0
    Accel::T  # Acceleration factor of the record, which the specification does not describe; kept, it has no effect on the simultaneous solution

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function REGCA1(; BUS::Int32, ID, Lvplsw, Tg, Rrpwr, Brkpt, Zerox, Lvpl1, Volim, Lvpnt1, Lvpnt0, Iolim, Tfltr, Khv, Iqrmax, Iqrmin, Accel, status::Bool=true)
    T = typeof(Lvplsw)
    check_parameter_domains(REGCA1{T}(BUS, ID, T(Lvplsw), T(Tg), T(Rrpwr), T(Brkpt), T(Zerox), T(Lvpl1), T(Volim), T(Lvpnt1), T(Lvpnt0), T(Iolim), T(Tfltr), T(Khv), T(Iqrmax), T(Iqrmin), T(Accel), status),
                            "keyword construction")
end

# Parameter domains, checked where data enters
function check_parameter_domains(x::REGCA1, from::AbstractString)
    (x.Lvplsw in (0.0, 1.0,)) || throw(ArgumentError("REGCA1 ($from): Lvplsw = $(x.Lvplsw) is not {0, 1}"))
    x.Zerox < x.Brkpt || throw(ArgumentError("REGCA1 ($from) does not satisfy its constraint lvpl_breakpoints: LVPL: Zerox is below Brkpt ([E18] Appendix F, Brkpt > Zerox)"))
    x.Lvpnt0 < x.Lvpnt1 || throw(ArgumentError("REGCA1 ($from) does not satisfy its constraint lvacm_breakpoints: low-voltage active current management: Lvpnt0 is below Lvpnt1 ([E18] Appendix F, Lvpnt1 > Lvpnt0)"))
    x.Rrpwr > 0 || throw(ArgumentError("REGCA1 ($from) does not satisfy its constraint rrpwr_positive: the up-ramp rate limit Rrpwr is positive ([E18] Appendix F, rrpwr > 0)"))
    x.Iqrmax > 0 || throw(ArgumentError("REGCA1 ($from) does not satisfy its constraint iqrmax_positive: the up rate limit Iqrmax is positive ([E18] Appendix F, iqrmax > 0)"))
    x.Iqrmin < 0 || throw(ArgumentError("REGCA1 ($from) does not satisfy its constraint iqrmin_negative: the down rate limit Iqrmin is negative ([E18] Appendix F, iqrmin < 0)"))
    x.Iolim < 0 || throw(ArgumentError("REGCA1 ($from) does not satisfy its constraint qmin_negative: qmin, Iolim, is negative ([E18] Appendix F, qmin < 0)"))
    (x.Brkpt - x.Zerox != 0) || throw(ArgumentError("REGCA1 ($from): Brkpt - Zerox is zero, and the model divides by it (services.kLVPL)"))
    (x.Tfltr != 0) || throw(ArgumentError("REGCA1 ($from): Tfltr is zero, and the model divides by it (graph[0].Vflt_rate)"))
    (x.Tg != 0) || throw(ArgumentError("REGCA1 ($from): Tg is zero, and the model divides by it (graph[1].Ip_raw)"))
    (x.Lvpnt1 - x.Lvpnt0 != 0) || throw(ArgumentError("REGCA1 ($from): Lvpnt1 - Lvpnt0 is zero, and the model divides by it (graph[3].multp)"))
    return x
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:REGCA1}) = true
get_dyr_model_name(::Type{<:REGCA1}) = "REGCA1"


# ============================================================================
# Group Interface Implementation - :ren_gen
# ============================================================================

# Declare that REGCA1 implements the :ren_gen group interface
implements_group(::Type{<:REGCA1}) = :ren_gen

# Map group canonical variable names to REGCA1's actual variable names
function group_var_mapping(::Type{<:REGCA1})
    Dict(
        :Ip => :Ip, # Same name
        :Ipcmd => :Ipcmd, # Same name
        :Iq_lag => :Iq_lag, # Same name
        :Iqcmd => :Iqcmd, # Same name
        :Pe => :Pe, # Same name
        :Qe => :Qe, # Same name
        :Vflt => :Vflt, # Same name
    )
end

# Map group canonical residual names to REGCA1's actual residual names
function group_res_mapping(::Type{<:REGCA1})
    Dict(
        :g_Ipcmd => :g_Ipcmd, # Same name
        :g_Iqcmd => :g_Iqcmd, # Same name
    )
end

# Link specification - which models this links to and how
function get_replacement_map(::Type{<:REGCA1})
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

# ============================================================================
# Model Registration
# ============================================================================

supports_format(::Type{<:REGCA1}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model REGCA1

# Terminal declarations - where REGCA1 connects
const REGCA1_TERMINALS = (
    Terminal(:bus_group, :BUS),
)

get_terminals(::Type{<:REGCA1}) = REGCA1_TERMINALS

get_key_name(::Type{<:REGCA1}) = [:BUS, :ID]

function model_metadata(::Type{REGCA1})
    ModelMetadata(
        name = :REGCA1,
        vars = [
            ModelVar(:Vflt, StateVar();
                description="Filtered terminal voltage (s2), which LVPL reads"),
            ModelVar(:Ip, StateVar();
                description="Active current (s1): free, its own rate, its rise limited at Rrpwr; held within LVPL (automaton Ip)"),
            ModelVar(:Iq_lag, StateVar();
                description="Reactive current (s0), positive injecting; its recovery rate limited by Iqrmax or Iqrmin"),
            ModelVar(:Ipcmd, AlgebVar();
                description="Active current command (receives the electrical controller's contribution)"),
            ModelVar(:Iqcmd, AlgebVar();
                description="Reactive current command, positive absorbing (receives the electrical controller's contribution)"),
            ModelVar(:dq, AlgebVar();
                description="Reduction of the reactive power by the high-voltage reactive current management"),
            ModelVar(:Pe, AlgebVar();
                description="Active power injected into the bus, Mbase"),
            ModelVar(:Qe, AlgebVar();
                description="Reactive power injected into the bus, q of the high-voltage reactive current management, Mbase"),
            ModelVar(:v, AlgebVar(), :bus_group, :v, :BUS;
                description="Terminal voltage magnitude Vt (external)", units="pu")
        ],
        residuals = [
        ModelResidual(:dVflt, StateRes(), PrivateRes();
            description="Filtered terminal voltage (s2), which LVPL reads"),
        ModelResidual(:dIp, StateRes(), PrivateRes();
            description="Active current (s1): free, its own rate, its rise limited at Rrpwr; held within LVPL (automaton Ip)"),
        ModelResidual(:dIq_lag, StateRes(), PrivateRes();
            description="Reactive current (s0), positive injecting; its recovery rate limited by Iqrmax or Iqrmin"),
        ModelResidual(:g_Ipcmd, AlgebRes(), SharedRes(:add);
            description="Active current command (receives the electrical controller's contribution)"),
        ModelResidual(:g_Iqcmd, AlgebRes(), SharedRes(:add);
            description="Reactive current command, positive absorbing (receives the electrical controller's contribution)"),
        ModelResidual(:g_dq, AlgebRes(), PrivateRes();
            description="Reduction of the reactive power by the high-voltage reactive current management"),
        ModelResidual(:g_Pe, AlgebRes(), PrivateRes();
            description="Active power injected into the bus, Mbase"),
        ModelResidual(:g_Qe, AlgebRes(), PrivateRes();
            description="Reactive power injected into the bus, q of the high-voltage reactive current management, Mbase"),
        ModelResidual(:ip_inj, AlgebRes(), SharedRes(:add), :bus_group, :ip, :BUS;
            description="Active current injection into the bus, system base"),
        ModelResidual(:iq_inj, AlgebRes(), SharedRes(:add), :bus_group, :iq, :BUS;
            description="Reactive current injection into the bus, system base: Qe / v, defined at v = 0"),
        ModelResidual(:p_inj, AlgebRes(), SharedRes(:add), :bus_group, :p, :BUS;
            description="p injection into the bus, derived from the current form the model declares"),
        ModelResidual(:q_inj, AlgebRes(), SharedRes(:add), :bus_group, :q, :BUS;
            description="q injection into the bus, derived from the current form the model declares")
    ],
        ratings = [ModelRating(:mbase, NominalPower; from = LinkedDevice(),
            description="Rating of the static device this model replaces; the model's per-unit data are on it")],
        derived_params = [RatingRatio(:Mratio, NominalPower;
            description="Converter rating Mbase over the system base Sbase")],
        ext_params = ExtParam[],
        init_inputs = [
        InitInput(:v, :from_solution, :Bus, :v;
            description="Terminal voltage from the power flow solution"),
        InitInput(:p_gen, :from_observation, :Generator, :p_gen;
            description="Active power from the power flow"),
        InitInput(:q_gen, :from_observation, :Generator, :q_gen;
            description="Reactive power from the power flow")
    ],
        observable_vars = Symbol[],
        initialized_params = [
        InitParam(:Qgen0, :Qe; description="Initial reactive power, whose sign chooses the rate limit of s0"),
        InitParam(:Ipcmd0, :Ipcmd; description="Active current command (receives the electrical controller's contribution)"),
        InitParam(:Iqcmd0, :Iqcmd; description="Reactive current command, positive absorbing (receives the electrical controller's contribution)")
    ],
        injection = InjectionForms(:current, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj]),
        automata = (
            Automaton(:iq_recovery;
                modes = [Mode(:upward; holds = [:up_limited => 1.0, :down_limited => 0.0]),
                    Mode(:downward; holds = [:up_limited => 0.0, :down_limited => 1.0]),
                    Mode(:neither; holds = [:up_limited => 0.0, :down_limited => 0.0])],
                transitions = Transition[],
                initial = [InitialMode(:upward; when = Atom(:Qgen0, :>)),
                    InitialMode(:downward; when = Atom(:Qgen0, :<)),
                    InitialMode(:neither)]),
            Automaton(:hv;
                modes = [Mode(:clear; holds = [:hv_clear => 1.0, :hv_clamp => 0.0, :hv_floor => 0.0]),
                    Mode(:clamp; holds = [:hv_clear => 0.0, :hv_clamp => 1.0, :hv_floor => 0.0]),
                    Mode(:floor; holds = [:hv_clear => 0.0, :hv_clamp => 0.0, :hv_floor => 1.0])],
                transitions = [Transition(:clear, :clamp; when = AllOf(Condition[Atom(:v_minus_Volim, :>), Atom(:hv_value, :>)]), continuity = [:v_minus_Volim => [:dq]]),
                    Transition(:clamp, :clear; when = Atom(:dq, :<), crossing_continuity = [:dq => [:v_minus_Volim]]),
                    Transition(:clamp, :floor; when = AllOf(Condition[Atom(:dq, :>=), Atom(:hv_value_2, :<)]), continuity = [:hv_value_2 => [:v_minus_Volim]]),
                    Transition(:floor, :clamp; when = AllOf(Condition[Atom(:dq, :>=), Atom(:v_minus_Volim, :<)]), continuity = [:v_minus_Volim => [:hv_value_2]]),
                    Transition(:floor, :clear; when = Atom(:dq, :<), continuity = [:dq => [:hv_value, :hv_value_2]])],
                initial = [InitialMode(:clear)]),
            Automaton(:Ip;
                modes = [Mode(:free; holds = [:Ip_gate => 1.0, :Ip_at_upper_1 => 0.0]),
                    Mode(:upper_0; holds = [:Ip_gate => 0.0, :Ip_at_upper_1 => 0.0], states = [:Ip => 0.0]),
                    Mode(:upper_1; holds = [:Ip_gate => 0.0, :Ip_at_upper_1 => 1.0], states = [:Ip => :Ip_upper_1])],
                transitions = [Transition(:free, :upper_0; when = AllOf(Condition[Atom(:Vflt_minus_Zerox, :<=), Atom(:Lvplsw, :>), AnyOf(Condition[Atom(:Ip, :>), AllOf(Condition[Atom(:Ip, :>=), Atom(:Ip_rate, :>)])])])),
                    Transition(:upper_0, :free; when = AnyOf(Condition[Atom(:Ip, :<), Atom(:Ip_rate, :<), Negation(AllOf(Condition[Atom(:Vflt_minus_Zerox, :<=), Atom(:Lvplsw, :>)]))]), continuity = [:Vflt_minus_Zerox => [:Ip, :Ip_value]]),
                    Transition(:free, :upper_1; when = AllOf(Condition[Atom(:Lvplsw, :>), Atom(:Vflt_minus_Zerox, :>), Atom(:Vflt_minus_Brkpt, :<), AnyOf(Condition[Atom(:Ip_value, :>), AllOf(Condition[Atom(:Ip_value_2, :>), Atom(:Ip_value, :>=)])])])),
                    Transition(:upper_1, :free; when = AnyOf(Condition[Atom(:Ip_value_2, :<), Atom(:Ip_value, :<), Negation(AllOf(Condition[Atom(:Lvplsw, :>), Atom(:Vflt_minus_Zerox, :>), Atom(:Vflt_minus_Brkpt, :<)]))]), continuity = [:Ip_value_2 => [:Ip_value], :Vflt_minus_Zerox => [:Ip, :Ip_value], :Vflt_minus_Brkpt => [:Ip_value]])],
                initial = [InitialMode(:free)],
                at_rest = [:Ip_rate, :Ip_value_2]),
            Automaton(:multp_clamp;
                modes = [Mode(:high; holds = [:multp_clamp_high => 1.0, :multp_clamp_low => 0.0, :multp_clamp_mid => 0.0]),
                    Mode(:low; holds = [:multp_clamp_high => 0.0, :multp_clamp_low => 1.0, :multp_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:multp_clamp_high => 0.0, :multp_clamp_low => 0.0, :multp_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:multp_clamp_value, :<), Atom(:multp_clamp_value_2, :<=)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:multp_clamp_value, :>=), Atom(:multp_clamp_value_2, :<=)])),
                    Transition(:low, :high; when = Atom(:multp_clamp_value_2, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:multp_clamp_value, :>=), Atom(:multp_clamp_value_2, :<=)])),
                    Transition(:mid, :high; when = Atom(:multp_clamp_value_2, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:multp_clamp_value, :<), Atom(:multp_clamp_value_2, :<=)]))],
                initial = [InitialMode(:high; when = Atom(:multp_clamp_value_2, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:multp_clamp_value, :<), Atom(:multp_clamp_value_2, :<=)])),
                    InitialMode(:mid)]),
            Automaton(:Ip_rate_min;
                modes = [Mode(:second; holds = [:Ip_rate_min_second => 1.0, :Ip_rate_min_first => 0.0]),
                    Mode(:first; holds = [:Ip_rate_min_second => 0.0, :Ip_rate_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Ip_rate_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:Ip_rate_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:Ip_rate_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:Qe_max;
                modes = [Mode(:second; holds = [:Qe_max_second => 1.0, :Qe_max_first => 0.0]),
                    Mode(:first; holds = [:Qe_max_second => 0.0, :Qe_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:hv_value, :<=)),
                    Transition(:first, :second; when = Atom(:hv_value, :>))],
                initial = [InitialMode(:second; when = Atom(:hv_value, :>)),
                    InitialMode(:first)]),
            Automaton(:Iq_lag_rate_min;
                modes = [Mode(:second; holds = [:Iq_lag_rate_min_second => 1.0, :Iq_lag_rate_min_first => 0.0]),
                    Mode(:first; holds = [:Iq_lag_rate_min_second => 0.0, :Iq_lag_rate_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Iq_lag_rate_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:Iq_lag_rate_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:Iq_lag_rate_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:Iq_lag_rate_max;
                modes = [Mode(:second; holds = [:Iq_lag_rate_max_second => 1.0, :Iq_lag_rate_max_first => 0.0]),
                    Mode(:first; holds = [:Iq_lag_rate_max_second => 0.0, :Iq_lag_rate_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Iq_lag_rate_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:Iq_lag_rate_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:Iq_lag_rate_max_value, :>)),
                    InitialMode(:first)]),
        ),
        status_field = :status,
        evaluation_pass = 1,
    )
end

@inline function regca1_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t, ::Type{AT}
) where AT
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        kLVPL = params.Lvpl1[i] / (params.Brkpt[i] - params.Zerox[i])

        Vflt_rate = (vars.v[i] - vars.Vflt[i]) / params.Tfltr[i]
        Ip_raw = (vars.Ipcmd[i] - vars.Ip[i]) / params.Tg[i]
        Iq_raw = (-vars.Iq_lag[i] - vars.Iqcmd[i]) / params.Tg[i]
        multp = params.multp_clamp_high[i] + params.multp_clamp_mid[i] * (vars.v[i] - params.Lvpnt0[i]) / (params.Lvpnt1[i] - params.Lvpnt0[i])
        Ip_rate = params.Ip_rate_min_first[i] * Ip_raw + params.Ip_rate_min_second[i] * params.Rrpwr[i]
        # === Residual Assignments ===
        outputs.dVflt[i] = mask * (Vflt_rate)
        outputs.dIp[i] = mask * (params.Ip_gate[i] * Ip_rate + params.Ip_at_upper_1[i] * Vflt_rate * kLVPL)
        outputs.dIq_lag[i] = mask * (Iq_raw - params.Iq_lag_rate_max_second[i] * params.up_limited[i] * (Iq_raw - params.Iqrmax[i]) - params.Iq_lag_rate_min_second[i] * params.down_limited[i] * (Iq_raw - params.Iqrmin[i]))
        outputs.g_Ipcmd[i] = mask * (params.Ipcmd0[i] - vars.Ipcmd[i])
        outputs.g_Iqcmd[i] = mask * (params.Iqcmd0[i] - vars.Iqcmd[i])
        outputs.g_dq[i] = mask * (params.hv_clamp[i] * (params.Volim[i] - vars.v[i]) + params.hv_floor[i] * (vars.Iq_lag[i] * vars.v[i] - params.Iolim[i] - vars.dq[i]) - vars.dq[i] * params.hv_clear[i])
        outputs.g_Pe[i] = mask * (vars.Ip[i] * multp * vars.v[i] - vars.Pe[i])
        outputs.g_Qe[i] = mask * (params.Iolim[i] * params.hv_floor[i] + params.hv_clamp[i] * (vars.Iq_lag[i] * vars.v[i] - vars.dq[i]) + params.hv_clear[i] * (params.Iolim[i] * params.Qe_max_first[i] + vars.Iq_lag[i] * params.Qe_max_second[i] * vars.v[i]) - vars.Qe[i])
        outputs.ip_inj[i] = mask * (vars.Ip[i] * params.Mratio[i] * multp)
        outputs.iq_inj[i] = mask * (params.Mratio[i] * ifelse(vars.v[i] > 0, vars.Qe[i] / vars.v[i], vars.Iq_lag[i]))
        outputs.p_inj[i] = mask * (vars.Ip[i] * params.Mratio[i] * multp * vars.v[i])
        outputs.q_inj[i] = mask * (params.Mratio[i] * vars.v[i] * ifelse(vars.v[i] > 0, vars.Qe[i] / vars.v[i], vars.Iq_lag[i]))
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REGCA1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    regca1_residual_kernel!(
        (dVflt = ws.res_vals.dVflt, dIp = ws.res_vals.dIp, dIq_lag = ws.res_vals.dIq_lag, g_Ipcmd = ws.res_vals.g_Ipcmd, g_Iqcmd = ws.res_vals.g_Iqcmd, g_dq = ws.res_vals.g_dq, g_Pe = ws.res_vals.g_Pe, g_Qe = ws.res_vals.g_Qe, ip_inj = ws.res_vals.ip_inj, iq_inj = ws.res_vals.iq_inj, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj),
        (Lvplsw = model.Lvplsw, Tg = model.Tg, Rrpwr = model.Rrpwr, Brkpt = model.Brkpt, Zerox = model.Zerox, Lvpl1 = model.Lvpl1, Volim = model.Volim, Lvpnt1 = model.Lvpnt1, Lvpnt0 = model.Lvpnt0, Iolim = model.Iolim, Tfltr = model.Tfltr, Khv = model.Khv, Iqrmax = model.Iqrmax, Iqrmin = model.Iqrmin, Accel = model.Accel, Qgen0 = ws.setpoints.Qgen0, Ipcmd0 = ws.setpoints.Ipcmd0, Iqcmd0 = ws.setpoints.Iqcmd0, up_limited = ws.automata.discrete.up_limited, down_limited = ws.automata.discrete.down_limited, hv_clear = ws.automata.discrete.hv_clear, hv_clamp = ws.automata.discrete.hv_clamp, hv_floor = ws.automata.discrete.hv_floor, Ip_gate = ws.automata.discrete.Ip_gate, Ip_at_upper_1 = ws.automata.discrete.Ip_at_upper_1, multp_clamp_high = ws.automata.discrete.multp_clamp_high, multp_clamp_low = ws.automata.discrete.multp_clamp_low, multp_clamp_mid = ws.automata.discrete.multp_clamp_mid, Ip_rate_min_second = ws.automata.discrete.Ip_rate_min_second, Ip_rate_min_first = ws.automata.discrete.Ip_rate_min_first, Qe_max_second = ws.automata.discrete.Qe_max_second, Qe_max_first = ws.automata.discrete.Qe_max_first, Iq_lag_rate_min_second = ws.automata.discrete.Iq_lag_rate_min_second, Iq_lag_rate_min_first = ws.automata.discrete.Iq_lag_rate_min_first, Iq_lag_rate_max_second = ws.automata.discrete.Iq_lag_rate_max_second, Iq_lag_rate_max_first = ws.automata.discrete.Iq_lag_rate_max_first, Mratio = external.Mratio),
        (Vflt = ws.var_vals.Vflt, Ip = ws.var_vals.Ip, Iq_lag = ws.var_vals.Iq_lag, Ipcmd = ws.var_vals.Ipcmd, Iqcmd = ws.var_vals.Iqcmd, dq = ws.var_vals.dq, Pe = ws.var_vals.Pe, Qe = ws.var_vals.Qe, v = ws.var_vals.v),
        active_mask,
        t,
        AT
    )
end

@inline function regca1_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.Qgen0)
        kLVPL = params.Lvpl1[i] / (params.Brkpt[i] - params.Zerox[i])
        Vflt_rate = (vars.v[i] - vars.Vflt[i]) / params.Tfltr[i]
        Ip_raw = (vars.Ipcmd[i] - vars.Ip[i]) / params.Tg[i]
        Iq_raw = (-vars.Iq_lag[i] - vars.Iqcmd[i]) / params.Tg[i]
        Ip_rate = params.Ip_rate_min_first[i] * Ip_raw + params.Ip_rate_min_second[i] * params.Rrpwr[i]
        # === Automaton Values ===
        values.Qgen0[i] = params.Qgen0[i]
        values.v_minus_Volim[i] = vars.v[i] - params.Volim[i]
        values.hv_value[i] = vars.Iq_lag[i] * vars.v[i] - params.Iolim[i]
        values.dq[i] = vars.dq[i]
        values.hv_value_2[i] = vars.Iq_lag[i] * vars.v[i] - params.Iolim[i] - vars.dq[i]
        values.Ip_upper_1[i] = kLVPL * (vars.Vflt[i] - params.Zerox[i])
        values.Vflt_minus_Zerox[i] = vars.Vflt[i] - params.Zerox[i]
        values.Lvplsw[i] = params.Lvplsw[i]
        values.Ip[i] = vars.Ip[i]
        values.Ip_rate[i] = Ip_rate
        values.Vflt_minus_Brkpt[i] = vars.Vflt[i] - params.Brkpt[i]
        values.Ip_value_2[i] = Ip_rate - Vflt_rate * kLVPL
        values.multp_clamp_value[i] = (vars.v[i] - params.Lvpnt0[i]) / (params.Lvpnt1[i] - params.Lvpnt0[i])
        values.multp_clamp_value_2[i] = (vars.v[i] - params.Lvpnt0[i]) / (params.Lvpnt1[i] - params.Lvpnt0[i]) - 1
        values.Ip_rate_min_value[i] = params.Rrpwr[i] - Ip_raw
        values.Iq_lag_rate_min_value[i] = Iq_raw - params.Iqrmin[i]
        values.Iq_lag_rate_max_value[i] = Iq_raw - params.Iqrmax[i]
        values.Ip_value[i] = vars.Ip[i] - values.Ip_upper_1[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REGCA1},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    regca1_automaton_values_kernel!(
        (Qgen0 = ws.automata.values.Qgen0, v_minus_Volim = ws.automata.values.v_minus_Volim, hv_value = ws.automata.values.hv_value, dq = ws.automata.values.dq, hv_value_2 = ws.automata.values.hv_value_2, Ip_upper_1 = ws.automata.values.Ip_upper_1, Vflt_minus_Zerox = ws.automata.values.Vflt_minus_Zerox, Lvplsw = ws.automata.values.Lvplsw, Ip = ws.automata.values.Ip, Ip_rate = ws.automata.values.Ip_rate, Vflt_minus_Brkpt = ws.automata.values.Vflt_minus_Brkpt, Ip_value = ws.automata.values.Ip_value, Ip_value_2 = ws.automata.values.Ip_value_2, multp_clamp_value = ws.automata.values.multp_clamp_value, multp_clamp_value_2 = ws.automata.values.multp_clamp_value_2, Ip_rate_min_value = ws.automata.values.Ip_rate_min_value, Iq_lag_rate_min_value = ws.automata.values.Iq_lag_rate_min_value, Iq_lag_rate_max_value = ws.automata.values.Iq_lag_rate_max_value),
        (Lvplsw = model.Lvplsw, Tg = model.Tg, Rrpwr = model.Rrpwr, Brkpt = model.Brkpt, Zerox = model.Zerox, Lvpl1 = model.Lvpl1, Volim = model.Volim, Lvpnt1 = model.Lvpnt1, Lvpnt0 = model.Lvpnt0, Iolim = model.Iolim, Tfltr = model.Tfltr, Khv = model.Khv, Iqrmax = model.Iqrmax, Iqrmin = model.Iqrmin, Accel = model.Accel, Qgen0 = ws.setpoints.Qgen0, Ipcmd0 = ws.setpoints.Ipcmd0, Iqcmd0 = ws.setpoints.Iqcmd0, up_limited = ws.automata.discrete.up_limited, down_limited = ws.automata.discrete.down_limited, hv_clear = ws.automata.discrete.hv_clear, hv_clamp = ws.automata.discrete.hv_clamp, hv_floor = ws.automata.discrete.hv_floor, Ip_gate = ws.automata.discrete.Ip_gate, Ip_at_upper_1 = ws.automata.discrete.Ip_at_upper_1, multp_clamp_high = ws.automata.discrete.multp_clamp_high, multp_clamp_low = ws.automata.discrete.multp_clamp_low, multp_clamp_mid = ws.automata.discrete.multp_clamp_mid, Ip_rate_min_second = ws.automata.discrete.Ip_rate_min_second, Ip_rate_min_first = ws.automata.discrete.Ip_rate_min_first, Qe_max_second = ws.automata.discrete.Qe_max_second, Qe_max_first = ws.automata.discrete.Qe_max_first, Iq_lag_rate_min_second = ws.automata.discrete.Iq_lag_rate_min_second, Iq_lag_rate_min_first = ws.automata.discrete.Iq_lag_rate_min_first, Iq_lag_rate_max_second = ws.automata.discrete.Iq_lag_rate_max_second, Iq_lag_rate_max_first = ws.automata.discrete.Iq_lag_rate_max_first, Mratio = external.Mratio),
        (Vflt = ws.var_vals.Vflt, Ip = ws.var_vals.Ip, Iq_lag = ws.var_vals.Iq_lag, Ipcmd = ws.var_vals.Ipcmd, Iqcmd = ws.var_vals.Iqcmd, dq = ws.var_vals.dq, Pe = ws.var_vals.Pe, Qe = ws.var_vals.Qe, v = ws.var_vals.v),
        t,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{REGCA1{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, REGCA1)

    v = T.(data.v)
    p_gen = T.(data.p_gen)
    q_gen = T.(data.q_gen)
    Mratio = data.Mratio

    # === Pre-allocate Result Arrays ===
    Vflt = Vector{T}(undef, n)
    Ip = Vector{T}(undef, n)
    Iq_lag = Vector{T}(undef, n)
    Ipcmd = Vector{T}(undef, n)
    Iqcmd = Vector{T}(undef, n)
    dq = fill(zero(T), n)
    Pe = Vector{T}(undef, n)
    Qe = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        kLVPL = model.Lvpl1[i] / (model.Brkpt[i] - model.Zerox[i])

        ifelse(p_gen[i] == 0, true, v[i] > model.Lvpnt0[i]) || throw(ArgumentError("REGCA1 instance $i cannot start: it requires that the power-flow voltage is above Lvpnt0, or the power flow's active power is zero: at or below Lvpnt0 the low-voltage active current management passes no active current (Appendix A, A-2), and none carries a nonzero active power (check active_current_passes)"))
        q_gen[i] / Mratio[i] >= model.Iolim[i] || throw(ArgumentError("REGCA1 instance $i cannot start: it requires that the power flow's reactive power, on the converter's base, is at least qmin (Iolim), the floor of the high-voltage reactive current management (Appendix A, A-1) (check reactive_power_above_floor)"))
        ifelse(v[i] <= model.Volim[i], true, q_gen[i] / Mratio[i] <= model.Iolim[i]) || throw(ArgumentError("REGCA1 instance $i cannot start: it requires that the power-flow voltage is at most Volim, unless the reactive power is at its floor: above it the high-voltage reactive current management holds the voltage at Volim (Appendix A, A-1), and the power flow does not (check voltage_within_clamp)"))

        # State initial values
        Vflt[i] = v[i]
        Ip[i] = ifelse(p_gen[i] == 0, 0, p_gen[i] / Mratio[i] / v[i] / clamp((v[i] - model.Lvpnt0[i]) / (model.Lvpnt1[i] - model.Lvpnt0[i]), 0, 1))
        ifelse(ifelse(Ip[i] <= 0, true, !ifelse(Vflt[i] <= model.Zerox[i], model.Lvplsw[i] > 0, false)), ifelse(Ip[i] <= kLVPL * (Vflt[i] - model.Zerox[i]), true, !ifelse(model.Lvplsw[i] > 0, ifelse(Vflt[i] > model.Zerox[i], Vflt[i] < model.Brkpt[i], false), false)), false) || throw(ArgumentError("REGCA1 instance $i cannot start: it requires that Ip starts within its upper bound (check Ip_within_upper)"))
        Iq_lag[i] = q_gen[i] / Mratio[i] / v[i]

        # Algebraic initial values
        Ipcmd[i] = Ip[i]
        Iqcmd[i] = -Iq_lag[i]
        Pe[i] = p_gen[i] / Mratio[i]
        Qe[i] = q_gen[i] / Mratio[i]
    end

    # === Return Initial Values ===
    return (Vflt=Vflt, Ip=Ip, Iq_lag=Iq_lag, Ipcmd=Ipcmd, Iqcmd=Iqcmd, dq=dq, Pe=Pe, Qe=Qe)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{REGCA1}, dyr_data, sys_pf, links) -> StructArray{REGCA1}

Parse REGCA1 model from DYR data with linking to power flow models.
"""
function parse_dynamic_model(
    ::Type{REGCA1},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    links::ModelLinks
)
    n = length(links)
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{REGCA1{T}}(undef, n)

    for idx in 1:n
        dyr_row = dyr_data[links.source_idx[idx]]
        instances[idx] = check_parameter_domains(REGCA1(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            Lvplsw = T(dyr_row.Lvplsw),
            Tg = T(dyr_row.Tg),
            Rrpwr = T(dyr_row.Rrpwr),
            Brkpt = T(dyr_row.Brkpt),
            Zerox = T(dyr_row.Zerox),
            Lvpl1 = T(dyr_row.Lvpl1),
            Volim = T(dyr_row.Volim),
            Lvpnt1 = T(dyr_row.Lvpnt1),
            Lvpnt0 = T(dyr_row.Lvpnt0),
            Iolim = T(dyr_row.Iolim),
            Tfltr = T(dyr_row.Tfltr),
            Khv = T(dyr_row.Khv),
            Iqrmax = T(dyr_row.Iqrmax),
            Iqrmin = T(dyr_row.Iqrmin),
            Accel = T(dyr_row.Accel)
        ), "dyr record $(idx)")
    end

    return StructArray(instances)
end
