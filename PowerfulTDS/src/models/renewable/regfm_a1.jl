"""REGFM_A1: Droop-controlled grid-forming converter (REGFM_A1)

DO NOT EDIT - changes will be overwritten.
"""


using InlineStrings
using PowerfulCore: Automaton, Mode, Transition, InitialMode, Condition, Atom, AllOf

export REGFM_A1

struct REGFM_A1{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number where the converter is connected"
    BUS::Int32
    "Converter identifier (1-2 characters)"
    ID::String3

    # === Parameters ===
    XL::T  # Coupling reactance, inverter base
    mq::T  # Q-V droop gain
    kpv::T  # Proportional gain of the voltage controller (VFlag = 1)
    kiv::T  # Integral gain of the voltage controller (VFlag = 1)
    Emax::T  # Upper limit of the output of the voltage loop
    Emin::T  # Lower limit of the output of the voltage loop
    mp::T  # P-f droop gain
    Pmax::T  # Upper limit of the inverter active power output
    Pmin::T  # Lower limit of the inverter active power output (negative for energy storage)
    kppmax::T  # Proportional gain of the overload mitigation controller
    kipmax::T  # Integral gain of the overload mitigation controller
    w0::T  # Rated angular frequency
    TPf::T  # Time constant of the low-pass filter for P measurement
    TQf::T  # Time constant of the low-pass filter for Q measurement
    TVf::T  # Time constant of the low-pass filter for V measurement
    Qmax::T  # Upper limit of the inverter reactive power output
    Qmin::T  # Lower limit of the inverter reactive power output
    VFlag::T  # Voltage control mode: 1 regulates the terminal voltage V through the PI, 0 the internal voltage E directly
    QVFlag::T  # 1: the plant controller changes Vref and the start sets Qref = 0; 0: it changes Qref and the start sets Qref = Qinv
    kpqmax::T  # Proportional gain of the Qmax and Qmin controller (3 for VFlag = 1, 0.1 for VFlag = 0)
    kiqmax::T  # Integral gain of the Qmax and Qmin controller (20 for VFlag = 1, 10 for VFlag = 0)
    ImaxF::T  # Inverter maximum transient output current

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function REGFM_A1(; BUS::Int32, ID, XL, mq, kpv, kiv, Emax, Emin, mp, Pmax, Pmin, kppmax, kipmax, w0, TPf, TQf, TVf, Qmax, Qmin, VFlag, QVFlag, kpqmax, kiqmax, ImaxF, status::Bool=true)
    T = typeof(XL)
    check_parameter_domains(REGFM_A1{T}(BUS, ID, T(XL), T(mq), T(kpv), T(kiv), T(Emax), T(Emin), T(mp), T(Pmax), T(Pmin), T(kppmax), T(kipmax), T(w0), T(TPf), T(TQf), T(TVf), T(Qmax), T(Qmin), T(VFlag), T(QVFlag), T(kpqmax), T(kiqmax), T(ImaxF), status),
                            "keyword construction")
end

# Parameter domains, checked where data enters
function check_parameter_domains(x::REGFM_A1, from::AbstractString)
    (x.VFlag in (0.0, 1.0,)) || throw(ArgumentError("REGFM_A1 ($from): VFlag = $(x.VFlag) is not {0, 1}"))
    (x.QVFlag in (0.0, 1.0,)) || throw(ArgumentError("REGFM_A1 ($from): QVFlag = $(x.QVFlag) is not {0, 1}"))
    x.TPf >= 0.0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.P_f.T_min: Lag admits T >= 0"))
    x.TQf >= 0.0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.Q_f.T_min: Lag admits T >= 0"))
    x.TVf >= 0.0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.V_f.T_min: Lag admits T >= 0"))
    x.Emin <= x.Emax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.XV.y_limits_ordered: control.XV.y: the lower limit of control.XV.y is at most its upper limit"))
    x.XL > 0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.coupling_reactance: the coupling reactance XL is positive ([P23] Table 1: 0.05 pu to 0.25 pu)"))
    x.ImaxF > 0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.current_limit: the maximum transient current ImaxF is positive ([P23] Table 1: 1.5 pu to 3 pu)"))
    ifelse(x.mp >= 0, x.mq >= 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.droop_gains: the droop gains are not negative: a negative gain is positive feedback in the P-f or the Q-V droop ([P23] Table 1)"))
    x.w0 > 0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.rated_frequency: the rated angular frequency w0 is positive"))
    x.Emin < x.Emax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.voltage_limits: Emin is below Emax: the voltage loop's limits are ordered ([P23] Table 1)"))
    x.Pmin < x.Pmax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.active_power_limits: Pmin is below Pmax: the active power limits are ordered ([P23] Table 1)"))
    x.Qmin < x.Qmax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.reactive_power_limits: Qmin is below Qmax: the reactive power limits are ordered ([P23] Table 1)"))
    ifelse(x.kipmax > 0, x.kppmax >= 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.overload_gains: the overload mitigation controller's integral gain is positive and its proportional gain not negative: the limits it holds its integrators at act on the sign of the integral's rate"))
    ifelse(x.kiqmax > 0, x.kpqmax >= 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.reactive_limit_gains: the Qmax/Qmin controller's integral gain is positive and its proportional gain not negative: the limits it holds its integrators at act on the sign of the integral's rate"))
    ifelse(x.kiv >= 0, x.kpv >= 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint control.voltage_gains: the voltage controller's gains are not negative ([P23] Table 1)"))
    (x.XL != 0) || throw(ArgumentError("REGFM_A1 ($from): XL is zero, and the model divides by it (parts.control.intermediates.Id_p.expr)"))
    (x.TPf != 0) || throw(ArgumentError("REGFM_A1 ($from): TPf is zero, and the model divides by it (parts.control.parts.P_f.states.y.derivative)"))
    (x.TQf != 0) || throw(ArgumentError("REGFM_A1 ($from): TQf is zero, and the model divides by it (parts.control.parts.Q_f.states.y.derivative)"))
    (x.TVf != 0) || throw(ArgumentError("REGFM_A1 ($from): TVf is zero, and the model divides by it (parts.control.parts.V_f.states.y.derivative)"))
    return x
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:REGFM_A1}) = true
get_dyr_model_name(::Type{<:REGFM_A1}) = "REGFM_A1"


# ============================================================================
# Group Interface Implementation - :ren_gfm
# ============================================================================

# Declare that REGFM_A1 implements the :ren_gfm group interface
implements_group(::Type{<:REGFM_A1}) = :ren_gfm

# Map group canonical variable names to REGFM_A1's actual variable names
function group_var_mapping(::Type{<:REGFM_A1})
    Dict(
        :Pe => :Pe, # Same name
        :Pref => :control_Pref, # REGFM_A1 uses control_Pref
        :Qe => :Qe, # Same name
        :Qref => :control_Qref, # REGFM_A1 uses control_Qref
        :Vref => :control_Vref, # REGFM_A1 uses control_Vref
        :control_DELTA_y => :control_DELTA_y, # Same name
        :control_P_f_y => :control_P_f_y, # Same name
        :control_Pref => :control_Pref, # Same name
        :control_Q_f_y => :control_Q_f_y, # Same name
        :control_Qref => :control_Qref, # Same name
        :control_V_f_y => :control_V_f_y, # Same name
        :control_Vref => :control_Vref, # Same name
        :control_XPMAX_y => :control_XPMAX_y, # Same name
        :control_XPMIN_y => :control_XPMIN_y, # Same name
        :control_XQMAX_y => :control_XQMAX_y, # Same name
        :control_XQMIN_y => :control_XQMIN_y, # Same name
        :control_XV_y => :control_XV_y, # Same name
    )
end

# Map group canonical residual names to REGFM_A1's actual residual names
function group_res_mapping(::Type{<:REGFM_A1})
    Dict(
        :g_Pref => :g_control_Pref, # REGFM_A1 uses g_control_Pref
        :g_Qref => :g_control_Qref, # REGFM_A1 uses g_control_Qref
        :g_Vref => :g_control_Vref, # REGFM_A1 uses g_control_Vref
        :g_control_Pref => :g_control_Pref, # Same name
        :g_control_Qref => :g_control_Qref, # Same name
        :g_control_Vref => :g_control_Vref, # Same name
    )
end

# Link specification - which models this links to and how
function get_replacement_map(::Type{<:REGFM_A1})
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

supports_format(::Type{<:REGFM_A1}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model REGFM_A1

# Terminal declarations - where REGFM_A1 connects
const REGFM_A1_TERMINALS = (
    Terminal(:bus_group, :BUS),
)

get_terminals(::Type{<:REGFM_A1}) = REGFM_A1_TERMINALS

get_key_name(::Type{<:REGFM_A1}) = [:BUS, :ID]

function model_metadata(::Type{REGFM_A1})
    ModelMetadata(
        name = :REGFM_A1,
        vars = [
            ModelVar(:control_P_f_y, StateVar();
                description="Filtered active power, inverter base (eq. 4)"),
            ModelVar(:control_Q_f_y, StateVar();
                description="Filtered reactive power, inverter base (eq. 5)"),
            ModelVar(:control_V_f_y, StateVar();
                description="Filtered terminal voltage (eq. 6)"),
            ModelVar(:control_XPMAX_y, StateVar();
                description="Integrator of the Pmax branch of the overload mitigation control, not above 0"),
            ModelVar(:control_XPMIN_y, StateVar();
                description="Integrator of the Pmin branch of the overload mitigation control, not below 0"),
            ModelVar(:control_XQMAX_y, StateVar();
                description="Integrator of the Qmax controller, not above 0"),
            ModelVar(:control_XQMIN_y, StateVar();
                description="Integrator of the Qmin controller, not below 0"),
            ModelVar(:control_DELTA_y, StateVar();
                description="Angle of the internal voltage from the phasor reference, the integral of the frequency deviation"),
            ModelVar(:control_XV_y, StateVar();
                description="Integrator of the voltage controller, within [Emin, Emax]"),
            ModelVar(:control_Idroop, AlgebVar();
                description="Magnitude of the droop current of eq. 7, inverter base"),
            ModelVar(:control_Pref, AlgebVar();
                description="Active power reference, inverter base (the plant controller changes it)"),
            ModelVar(:control_Qref, AlgebVar();
                description="Reactive power reference, inverter base (the plant controller changes it where QVFlag = 0)"),
            ModelVar(:control_Vref, AlgebVar();
                description="Voltage set point (the plant controller changes it where QVFlag = 1); the start makes the droop output the internal voltage the power flow needs"),
            ModelVar(:Pe, AlgebVar();
                description="Active power the converter delivers at the terminal, inverter base"),
            ModelVar(:Qe, AlgebVar();
                description="Reactive power the converter delivers at the terminal, inverter base"),
            ModelVar(:v, AlgebVar(), :bus_group, :v, :BUS;
                description="Terminal voltage magnitude V (external)", units="pu"),
            ModelVar(:theta, AlgebVar(), :bus_group, :theta, :BUS;
                description="Terminal voltage angle (external)", units="rad")
        ],
        residuals = [
        ModelResidual(:dcontrol_P_f_y, StateRes(), PrivateRes();
            description="Filtered active power, inverter base (eq. 4)"),
        ModelResidual(:dcontrol_Q_f_y, StateRes(), PrivateRes();
            description="Filtered reactive power, inverter base (eq. 5)"),
        ModelResidual(:dcontrol_V_f_y, StateRes(), PrivateRes();
            description="Filtered terminal voltage (eq. 6)"),
        ModelResidual(:dcontrol_XPMAX_y, StateRes(), PrivateRes();
            description="Integrator of the Pmax branch of the overload mitigation control, not above 0"),
        ModelResidual(:dcontrol_XPMIN_y, StateRes(), PrivateRes();
            description="Integrator of the Pmin branch of the overload mitigation control, not below 0"),
        ModelResidual(:dcontrol_XQMAX_y, StateRes(), PrivateRes();
            description="Integrator of the Qmax controller, not above 0"),
        ModelResidual(:dcontrol_XQMIN_y, StateRes(), PrivateRes();
            description="Integrator of the Qmin controller, not below 0"),
        ModelResidual(:dcontrol_DELTA_y, StateRes(), PrivateRes();
            description="Angle of the internal voltage from the phasor reference, the integral of the frequency deviation"),
        ModelResidual(:dcontrol_XV_y, StateRes(), PrivateRes();
            description="Integrator of the voltage controller, within [Emin, Emax]"),
        ModelResidual(:g_control_Idroop, AlgebRes(), PrivateRes();
            description="Magnitude of the droop current of eq. 7, inverter base"),
        ModelResidual(:g_control_Pref, AlgebRes(), SharedRes(:add);
            description="Active power reference, inverter base (the plant controller changes it)"),
        ModelResidual(:g_control_Qref, AlgebRes(), SharedRes(:add);
            description="Reactive power reference, inverter base (the plant controller changes it where QVFlag = 0)"),
        ModelResidual(:g_control_Vref, AlgebRes(), SharedRes(:add);
            description="Voltage set point (the plant controller changes it where QVFlag = 1); the start makes the droop output the internal voltage the power flow needs"),
        ModelResidual(:g_Pe, AlgebRes(), PrivateRes();
            description="Active power the converter delivers at the terminal, inverter base"),
        ModelResidual(:g_Qe, AlgebRes(), PrivateRes();
            description="Reactive power the converter delivers at the terminal, inverter base"),
        ModelResidual(:ip_inj, AlgebRes(), SharedRes(:add), :bus_group, :ip, :BUS;
            description="Active current injection into the bus"),
        ModelResidual(:iq_inj, AlgebRes(), SharedRes(:add), :bus_group, :iq, :BUS;
            description="Reactive current injection into the bus"),
        ModelResidual(:p_inj, AlgebRes(), SharedRes(:add), :bus_group, :p, :BUS;
            description="p injection into the bus, derived from the current form the model declares"),
        ModelResidual(:q_inj, AlgebRes(), SharedRes(:add), :bus_group, :q, :BUS;
            description="q injection into the bus, derived from the current form the model declares")
    ],
        ratings = [ModelRating(:mbase, NominalPower; from = LinkedDevice(),
            description="Rating of the static device this model replaces; the model's per-unit data are on it")],
        derived_params = [RatingRatio(:Mratio, NominalPower;
            description="Inverter rating Mbase over the system base Sbase")],
        ext_params = ExtParam[],
        init_inputs = [
        InitInput(:v, :from_solution, :Bus, :v;
            description="Terminal voltage magnitude from the power flow solution"),
        InitInput(:theta, :from_solution, :Bus, :theta;
            description="Terminal voltage angle from the power flow solution"),
        InitInput(:p_gen, :from_observation, :Generator, :p_gen;
            description="Active power from the power flow"),
        InitInput(:q_gen, :from_observation, :Generator, :q_gen;
            description="Reactive power from the power flow")
    ],
        observable_vars = Symbol[],
        initialized_params = [
        InitParam(:control_Pref0, :control_Pref; description="Active power reference, inverter base (the plant controller changes it)"),
        InitParam(:control_Qref0, :control_Qref; description="Reactive power reference, inverter base (the plant controller changes it where QVFlag = 0)"),
        InitParam(:control_Vref0, :control_Vref; description="Voltage set point (the plant controller changes it where QVFlag = 1); the start makes the droop output the internal voltage the power flow needs")
    ],
        injection = InjectionForms(:current, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj]),
        variable_limits = [VariableLimit(:control_XPMAX_y, nothing, m -> 0), VariableLimit(:control_XPMIN_y, m -> 0, nothing), VariableLimit(:control_XQMAX_y, nothing, m -> 0), VariableLimit(:control_XQMIN_y, m -> 0, nothing), VariableLimit(:control_XV_y, :Emin, :Emax), VariableLimit(:control_P_f_y, :Pmin, :Pmax), VariableLimit(:control_Q_f_y, :Qmin, :Qmax), VariableLimit(:control_Idroop, nothing, :ImaxF)],
        automata = (
            Automaton(:control_XPMAX_y;
                modes = [Mode(:free; holds = [:control_XPMAX_y_gate => 1.0]),
                    Mode(:upper; holds = [:control_XPMAX_y_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:control_XPMAX_y, :>=), Atom(:control_XPMAX_y_rate, :>)])),
                    Transition(:upper, :free; when = Atom(:control_XPMAX_y_rate, :<))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:control_XPMAX_y, :>=), Atom(:control_XPMAX_y_rate, :>)])),
                    InitialMode(:free)]),
            Automaton(:control_XPMIN_y;
                modes = [Mode(:free; holds = [:control_XPMIN_y_gate => 1.0]),
                    Mode(:lower; holds = [:control_XPMIN_y_gate => 0.0])],
                transitions = [Transition(:free, :lower; when = AllOf(Condition[Atom(:control_XPMIN_y, :<=), Atom(:control_XPMIN_y_rate, :<)])),
                    Transition(:lower, :free; when = Atom(:control_XPMIN_y_rate, :>))],
                initial = [InitialMode(:lower; when = AllOf(Condition[Atom(:control_XPMIN_y, :<=), Atom(:control_XPMIN_y_rate, :<)])),
                    InitialMode(:free)]),
            Automaton(:control_XQMAX_y;
                modes = [Mode(:free; holds = [:control_XQMAX_y_gate => 1.0]),
                    Mode(:upper; holds = [:control_XQMAX_y_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:control_XQMAX_y, :>=), Atom(:control_XQMAX_y_rate, :>)])),
                    Transition(:upper, :free; when = Atom(:control_XQMAX_y_rate, :<))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:control_XQMAX_y, :>=), Atom(:control_XQMAX_y_rate, :>)])),
                    InitialMode(:free)]),
            Automaton(:control_XQMIN_y;
                modes = [Mode(:free; holds = [:control_XQMIN_y_gate => 1.0]),
                    Mode(:lower; holds = [:control_XQMIN_y_gate => 0.0])],
                transitions = [Transition(:free, :lower; when = AllOf(Condition[Atom(:control_XQMIN_y, :<=), Atom(:control_XQMIN_y_rate, :<)])),
                    Transition(:lower, :free; when = Atom(:control_XQMIN_y_rate, :>))],
                initial = [InitialMode(:lower; when = AllOf(Condition[Atom(:control_XQMIN_y, :<=), Atom(:control_XQMIN_y_rate, :<)])),
                    InitialMode(:free)]),
            Automaton(:control_XV_y;
                modes = [Mode(:free; holds = [:control_XV_y_gate => 1.0]),
                    Mode(:upper; holds = [:control_XV_y_gate => 0.0]),
                    Mode(:lower; holds = [:control_XV_y_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:control_XV_y_minus_Emax, :>=), Atom(:control_XV_y_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:control_XV_y_minus_Emin, :<=), Atom(:control_XV_y_rate, :<)])),
                    Transition(:upper, :free; when = Atom(:control_XV_y_rate, :<)),
                    Transition(:lower, :free; when = Atom(:control_XV_y_rate, :>))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:control_XV_y_minus_Emax, :>=), Atom(:control_XV_y_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:control_XV_y_minus_Emin, :<=), Atom(:control_XV_y_rate, :<)])),
                    InitialMode(:free)]),
            Automaton(:control_current_limit;
                modes = [Mode(:droop; holds = [:control_current_limited => 0.0]),
                    Mode(:limited; holds = [:control_current_limited => 1.0])],
                transitions = [Transition(:droop, :limited; when = Atom(:control_Idroop_minus_ImaxF, :>=)),
                    Transition(:limited, :droop; when = Atom(:control_Idroop_minus_ImaxF, :<))],
                initial = [InitialMode(:droop)]),
            Automaton(:control_Id_floor_max;
                modes = [Mode(:second; holds = [:control_Id_floor_max_second => 1.0, :control_Id_floor_max_first => 0.0]),
                    Mode(:first; holds = [:control_Id_floor_max_second => 0.0, :control_Id_floor_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:control_Id_floor_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:control_Id_floor_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:control_Id_floor_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:control_pmax_out_min;
                modes = [Mode(:second; holds = [:control_pmax_out_min_second => 1.0, :control_pmax_out_min_first => 0.0]),
                    Mode(:first; holds = [:control_pmax_out_min_second => 0.0, :control_pmax_out_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:control_pmax_out_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:control_pmax_out_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:control_pmax_out_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:control_pmin_out_max;
                modes = [Mode(:second; holds = [:control_pmin_out_max_second => 1.0, :control_pmin_out_max_first => 0.0]),
                    Mode(:first; holds = [:control_pmin_out_max_second => 0.0, :control_pmin_out_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:control_pmin_out_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:control_pmin_out_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:control_pmin_out_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:control_qmax_out_min;
                modes = [Mode(:second; holds = [:control_qmax_out_min_second => 1.0, :control_qmax_out_min_first => 0.0]),
                    Mode(:first; holds = [:control_qmax_out_min_second => 0.0, :control_qmax_out_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:control_qmax_out_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:control_qmax_out_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:control_qmax_out_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:control_qmin_out_max;
                modes = [Mode(:second; holds = [:control_qmin_out_max_second => 1.0, :control_qmin_out_max_first => 0.0]),
                    Mode(:first; holds = [:control_qmin_out_max_second => 0.0, :control_qmin_out_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:control_qmin_out_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:control_qmin_out_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:control_qmin_out_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:control_e_cmd_limited_clamp;
                modes = [Mode(:high; holds = [:control_e_cmd_limited_clamp_high => 1.0, :control_e_cmd_limited_clamp_low => 0.0, :control_e_cmd_limited_clamp_mid => 0.0]),
                    Mode(:low; holds = [:control_e_cmd_limited_clamp_high => 0.0, :control_e_cmd_limited_clamp_low => 1.0, :control_e_cmd_limited_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:control_e_cmd_limited_clamp_high => 0.0, :control_e_cmd_limited_clamp_low => 0.0, :control_e_cmd_limited_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:control_e_cmd_limited_clamp_value, :<=), Atom(:control_e_cmd_limited_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:control_e_cmd_limited_clamp_value, :<=), Atom(:control_e_cmd_limited_clamp_value_2, :>=)])),
                    Transition(:low, :high; when = Atom(:control_e_cmd_limited_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:control_e_cmd_limited_clamp_value, :<=), Atom(:control_e_cmd_limited_clamp_value_2, :>=)])),
                    Transition(:mid, :high; when = Atom(:control_e_cmd_limited_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:control_e_cmd_limited_clamp_value, :<=), Atom(:control_e_cmd_limited_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:control_e_cmd_limited_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:control_e_cmd_limited_clamp_value, :<=), Atom(:control_e_cmd_limited_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:control_e_pi_limited_clamp;
                modes = [Mode(:high; holds = [:control_e_pi_limited_clamp_high => 1.0, :control_e_pi_limited_clamp_low => 0.0, :control_e_pi_limited_clamp_mid => 0.0]),
                    Mode(:low; holds = [:control_e_pi_limited_clamp_high => 0.0, :control_e_pi_limited_clamp_low => 1.0, :control_e_pi_limited_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:control_e_pi_limited_clamp_high => 0.0, :control_e_pi_limited_clamp_low => 0.0, :control_e_pi_limited_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:control_e_pi_limited_clamp_value, :<=), Atom(:control_e_pi_limited_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:control_e_pi_limited_clamp_value, :<=), Atom(:control_e_pi_limited_clamp_value_2, :>=)])),
                    Transition(:low, :high; when = Atom(:control_e_pi_limited_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:control_e_pi_limited_clamp_value, :<=), Atom(:control_e_pi_limited_clamp_value_2, :>=)])),
                    Transition(:mid, :high; when = Atom(:control_e_pi_limited_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:control_e_pi_limited_clamp_value, :<=), Atom(:control_e_pi_limited_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:control_e_pi_limited_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:control_e_pi_limited_clamp_value, :<=), Atom(:control_e_pi_limited_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
        ),
        status_field = :status,
        evaluation_pass = 1,
    )
end

@inline function regfm_a1_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t, ::Type{AT}
) where AT
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        control_p_ref = vars.control_Pref[i] + 0
        control_q_ref = vars.control_Qref[i] + 0
        control_v_ref = vars.control_Vref[i] + 0
        control_cos_theta = cos(vars.theta[i])
        control_sin_theta = sin(vars.theta[i])
        control_Id_floor = params.ImaxF[i] * params.control_Id_floor_max_second[i] + vars.control_Idroop[i] * params.control_Id_floor_max_first[i]
        control_pmax_err = params.Pmax[i] - vars.control_P_f_y[i]
        control_pmin_err = params.Pmin[i] - vars.control_P_f_y[i]
        control_qmax_err = params.Qmax[i] - vars.control_Q_f_y[i]
        control_qmin_err = params.Qmin[i] - vars.control_Q_f_y[i]
        control_p_dev = control_p_ref - vars.control_P_f_y[i]
        control_v_x = vars.v[i] * control_cos_theta
        control_v_y = vars.v[i] * control_sin_theta
        control_Id_ratio = params.ImaxF[i] / control_Id_floor
        control_pmax_rate = params.kipmax[i] * control_pmax_err
        control_pmin_rate = params.kipmax[i] * control_pmin_err
        control_qmax_rate = params.kiqmax[i] * control_qmax_err
        control_qmin_rate = params.kiqmax[i] * control_qmin_err
        control_limit_excess = 1 - control_Id_ratio
        control_limit_cut = params.control_current_limited[i] * control_limit_excess
        control_XPMAX_y_rate = 1.0 * control_pmax_rate / 1.0
        control_XPMIN_y_rate = 1.0 * control_pmin_rate / 1.0
        control_XQMAX_y_rate = 1.0 * control_qmax_rate / 1.0
        control_XQMIN_y_rate = 1.0 * control_qmin_rate / 1.0
        control_pmax_sum = params.kppmax[i] * control_pmax_err + vars.control_XPMAX_y[i]
        control_pmin_sum = params.kppmax[i] * control_pmin_err + vars.control_XPMIN_y[i]
        control_qmax_sum = params.kpqmax[i] * control_qmax_err + vars.control_XQMAX_y[i]
        control_qmin_sum = params.kpqmax[i] * control_qmin_err + vars.control_XQMIN_y[i]
        control_limit_scale = 1 - control_limit_cut
        control_pmax_out = 0 * params.control_pmax_out_min_second[i] + control_pmax_sum * params.control_pmax_out_min_first[i]
        control_pmin_out = 0 * params.control_pmin_out_max_second[i] + control_pmin_sum * params.control_pmin_out_max_first[i]
        control_qmax_out = 0 * params.control_qmax_out_min_second[i] + control_qmax_sum * params.control_qmax_out_min_first[i]
        control_qmin_out = 0 * params.control_qmin_out_max_second[i] + control_qmin_sum * params.control_qmin_out_max_first[i]
        control_delta_rate = params.w0[i] * (params.mp[i] * control_p_dev + control_pmax_out + control_pmin_out)
        control_Vcmd = control_v_ref + params.mq[i] * (control_q_ref - vars.control_Q_f_y[i]) + control_qmax_out + control_qmin_out
        control_Verr = control_Vcmd - vars.control_V_f_y[i]
        control_e_cmd_limited = params.Emax[i] * params.control_e_cmd_limited_clamp_high[i] + params.Emin[i] * params.control_e_cmd_limited_clamp_low[i] + control_Vcmd * params.control_e_cmd_limited_clamp_mid[i]
        control_xv_rate = params.VFlag[i] * params.kiv[i] * control_Verr
        control_a = vars.control_DELTA_y[i] - vars.theta[i]
        control_cos_delta = cos(vars.control_DELTA_y[i])
        control_sin_delta = sin(vars.control_DELTA_y[i])
        control_sin_a = sin(control_a)
        control_cos_a = cos(control_a)
        control_XV_y_rate = 1.0 * control_xv_rate / 1.0
        control_e_pi = params.kpv[i] * control_Verr + vars.control_XV_y[i]
        control_e_pi_limited = params.Emax[i] * params.control_e_pi_limited_clamp_high[i] + params.Emin[i] * params.control_e_pi_limited_clamp_low[i] + control_e_pi * params.control_e_pi_limited_clamp_mid[i]
        control_Edroop = ifelse(params.VFlag[i] > 0.5, control_e_pi_limited, control_e_cmd_limited)
        control_Ed_sin = control_Edroop * control_sin_a
        control_Ed_cos = control_Edroop * control_cos_a
        control_Ed_x = control_Edroop * control_cos_delta
        control_Ed_y = control_Edroop * control_sin_delta
        control_dE_p = control_Ed_sin
        control_dE_q = control_Ed_cos - vars.v[i]
        control_dE_x = control_Ed_x - control_v_x
        control_dE_y = control_Ed_y - control_v_y
        control_Id_p = control_dE_p / params.XL[i]
        control_Id_q = control_dE_q / params.XL[i]
        control_dE_x_sq = control_dE_x ^ 2
        control_dE_y_sq = control_dE_y ^ 2
        control_dE2 = control_dE_x_sq + control_dE_y_sq
        control_ip_t = control_limit_scale * control_Id_p
        control_iq_t = control_limit_scale * control_Id_q
        control_dE_root = sqrt(control_dE2)
        control_dE_mag = ifelse(control_dE2 > 0, control_dE_root, 0)
        control_Id_mag = control_dE_mag / params.XL[i]
        # === Residual Assignments ===
        outputs.dcontrol_P_f_y[i] = mask * ((1.0 * vars.Pe[i] - vars.control_P_f_y[i]) / params.TPf[i])
        outputs.dcontrol_Q_f_y[i] = mask * ((1.0 * vars.Qe[i] - vars.control_Q_f_y[i]) / params.TQf[i])
        outputs.dcontrol_V_f_y[i] = mask * ((1.0 * vars.v[i] - vars.control_V_f_y[i]) / params.TVf[i])
        outputs.dcontrol_XPMAX_y[i] = mask * (params.control_XPMAX_y_gate[i] * control_XPMAX_y_rate)
        outputs.dcontrol_XPMIN_y[i] = mask * (params.control_XPMIN_y_gate[i] * control_XPMIN_y_rate)
        outputs.dcontrol_XQMAX_y[i] = mask * (params.control_XQMAX_y_gate[i] * control_XQMAX_y_rate)
        outputs.dcontrol_XQMIN_y[i] = mask * (params.control_XQMIN_y_gate[i] * control_XQMIN_y_rate)
        outputs.dcontrol_DELTA_y[i] = mask * (1.0 * control_delta_rate / 1.0)
        outputs.dcontrol_XV_y[i] = mask * (params.control_XV_y_gate[i] * control_XV_y_rate)
        outputs.g_control_Idroop[i] = mask * (control_Id_mag - vars.control_Idroop[i])
        outputs.g_control_Pref[i] = mask * (params.control_Pref0[i] - vars.control_Pref[i])
        outputs.g_control_Qref[i] = mask * (params.control_Qref0[i] - vars.control_Qref[i])
        outputs.g_control_Vref[i] = mask * (params.control_Vref0[i] - vars.control_Vref[i])
        outputs.g_Pe[i] = mask * (vars.v[i] * control_ip_t - vars.Pe[i])
        outputs.g_Qe[i] = mask * (vars.v[i] * control_iq_t - vars.Qe[i])
        outputs.ip_inj[i] = mask * (params.Mratio[i] * control_ip_t)
        outputs.iq_inj[i] = mask * (params.Mratio[i] * control_iq_t)
        outputs.p_inj[i] = mask * (vars.v[i] * (params.Mratio[i] * control_ip_t))
        outputs.q_inj[i] = mask * (vars.v[i] * (params.Mratio[i] * control_iq_t))
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REGFM_A1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    regfm_a1_residual_kernel!(
        (dcontrol_P_f_y = ws.res_vals.dcontrol_P_f_y, dcontrol_Q_f_y = ws.res_vals.dcontrol_Q_f_y, dcontrol_V_f_y = ws.res_vals.dcontrol_V_f_y, dcontrol_XPMAX_y = ws.res_vals.dcontrol_XPMAX_y, dcontrol_XPMIN_y = ws.res_vals.dcontrol_XPMIN_y, dcontrol_XQMAX_y = ws.res_vals.dcontrol_XQMAX_y, dcontrol_XQMIN_y = ws.res_vals.dcontrol_XQMIN_y, dcontrol_DELTA_y = ws.res_vals.dcontrol_DELTA_y, dcontrol_XV_y = ws.res_vals.dcontrol_XV_y, g_control_Idroop = ws.res_vals.g_control_Idroop, g_control_Pref = ws.res_vals.g_control_Pref, g_control_Qref = ws.res_vals.g_control_Qref, g_control_Vref = ws.res_vals.g_control_Vref, g_Pe = ws.res_vals.g_Pe, g_Qe = ws.res_vals.g_Qe, ip_inj = ws.res_vals.ip_inj, iq_inj = ws.res_vals.iq_inj, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj),
        (XL = model.XL, mq = model.mq, kpv = model.kpv, kiv = model.kiv, Emax = model.Emax, Emin = model.Emin, mp = model.mp, Pmax = model.Pmax, Pmin = model.Pmin, kppmax = model.kppmax, kipmax = model.kipmax, w0 = model.w0, TPf = model.TPf, TQf = model.TQf, TVf = model.TVf, Qmax = model.Qmax, Qmin = model.Qmin, VFlag = model.VFlag, QVFlag = model.QVFlag, kpqmax = model.kpqmax, kiqmax = model.kiqmax, ImaxF = model.ImaxF, control_Pref0 = ws.setpoints.control_Pref0, control_Qref0 = ws.setpoints.control_Qref0, control_Vref0 = ws.setpoints.control_Vref0, control_current_limited = ws.automata.discrete.control_current_limited, control_XPMAX_y_gate = ws.automata.discrete.control_XPMAX_y_gate, control_XPMIN_y_gate = ws.automata.discrete.control_XPMIN_y_gate, control_XQMAX_y_gate = ws.automata.discrete.control_XQMAX_y_gate, control_XQMIN_y_gate = ws.automata.discrete.control_XQMIN_y_gate, control_XV_y_gate = ws.automata.discrete.control_XV_y_gate, control_Id_floor_max_second = ws.automata.discrete.control_Id_floor_max_second, control_Id_floor_max_first = ws.automata.discrete.control_Id_floor_max_first, control_pmax_out_min_second = ws.automata.discrete.control_pmax_out_min_second, control_pmax_out_min_first = ws.automata.discrete.control_pmax_out_min_first, control_pmin_out_max_second = ws.automata.discrete.control_pmin_out_max_second, control_pmin_out_max_first = ws.automata.discrete.control_pmin_out_max_first, control_qmax_out_min_second = ws.automata.discrete.control_qmax_out_min_second, control_qmax_out_min_first = ws.automata.discrete.control_qmax_out_min_first, control_qmin_out_max_second = ws.automata.discrete.control_qmin_out_max_second, control_qmin_out_max_first = ws.automata.discrete.control_qmin_out_max_first, control_e_cmd_limited_clamp_high = ws.automata.discrete.control_e_cmd_limited_clamp_high, control_e_cmd_limited_clamp_low = ws.automata.discrete.control_e_cmd_limited_clamp_low, control_e_cmd_limited_clamp_mid = ws.automata.discrete.control_e_cmd_limited_clamp_mid, control_e_pi_limited_clamp_high = ws.automata.discrete.control_e_pi_limited_clamp_high, control_e_pi_limited_clamp_low = ws.automata.discrete.control_e_pi_limited_clamp_low, control_e_pi_limited_clamp_mid = ws.automata.discrete.control_e_pi_limited_clamp_mid, Mratio = external.Mratio),
        (control_P_f_y = ws.var_vals.control_P_f_y, control_Q_f_y = ws.var_vals.control_Q_f_y, control_V_f_y = ws.var_vals.control_V_f_y, control_XPMAX_y = ws.var_vals.control_XPMAX_y, control_XPMIN_y = ws.var_vals.control_XPMIN_y, control_XQMAX_y = ws.var_vals.control_XQMAX_y, control_XQMIN_y = ws.var_vals.control_XQMIN_y, control_DELTA_y = ws.var_vals.control_DELTA_y, control_XV_y = ws.var_vals.control_XV_y, control_Idroop = ws.var_vals.control_Idroop, control_Pref = ws.var_vals.control_Pref, control_Qref = ws.var_vals.control_Qref, control_Vref = ws.var_vals.control_Vref, Pe = ws.var_vals.Pe, Qe = ws.var_vals.Qe, v = ws.var_vals.v, theta = ws.var_vals.theta),
        active_mask,
        t,
        AT
    )
end

@inline function regfm_a1_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.control_XPMAX_y)
        control_q_ref = vars.control_Qref[i] + 0
        control_v_ref = vars.control_Vref[i] + 0
        control_pmax_err = params.Pmax[i] - vars.control_P_f_y[i]
        control_pmin_err = params.Pmin[i] - vars.control_P_f_y[i]
        control_qmax_err = params.Qmax[i] - vars.control_Q_f_y[i]
        control_qmin_err = params.Qmin[i] - vars.control_Q_f_y[i]
        control_pmax_rate = params.kipmax[i] * control_pmax_err
        control_pmin_rate = params.kipmax[i] * control_pmin_err
        control_qmax_rate = params.kiqmax[i] * control_qmax_err
        control_qmin_rate = params.kiqmax[i] * control_qmin_err
        control_XPMAX_y_rate = 1.0 * control_pmax_rate / 1.0
        control_XPMIN_y_rate = 1.0 * control_pmin_rate / 1.0
        control_XQMAX_y_rate = 1.0 * control_qmax_rate / 1.0
        control_XQMIN_y_rate = 1.0 * control_qmin_rate / 1.0
        control_pmax_sum = params.kppmax[i] * control_pmax_err + vars.control_XPMAX_y[i]
        control_pmin_sum = params.kppmax[i] * control_pmin_err + vars.control_XPMIN_y[i]
        control_qmax_sum = params.kpqmax[i] * control_qmax_err + vars.control_XQMAX_y[i]
        control_qmin_sum = params.kpqmax[i] * control_qmin_err + vars.control_XQMIN_y[i]
        control_qmax_out = 0 * params.control_qmax_out_min_second[i] + control_qmax_sum * params.control_qmax_out_min_first[i]
        control_qmin_out = 0 * params.control_qmin_out_max_second[i] + control_qmin_sum * params.control_qmin_out_max_first[i]
        control_Vcmd = control_v_ref + params.mq[i] * (control_q_ref - vars.control_Q_f_y[i]) + control_qmax_out + control_qmin_out
        control_Verr = control_Vcmd - vars.control_V_f_y[i]
        control_xv_rate = params.VFlag[i] * params.kiv[i] * control_Verr
        control_XV_y_rate = 1.0 * control_xv_rate / 1.0
        control_e_pi = params.kpv[i] * control_Verr + vars.control_XV_y[i]
        # === Automaton Values ===
        values.control_XPMAX_y[i] = vars.control_XPMAX_y[i]
        values.control_XPMAX_y_rate[i] = control_XPMAX_y_rate
        values.control_XPMIN_y[i] = vars.control_XPMIN_y[i]
        values.control_XPMIN_y_rate[i] = control_XPMIN_y_rate
        values.control_XQMAX_y[i] = vars.control_XQMAX_y[i]
        values.control_XQMAX_y_rate[i] = control_XQMAX_y_rate
        values.control_XQMIN_y[i] = vars.control_XQMIN_y[i]
        values.control_XQMIN_y_rate[i] = control_XQMIN_y_rate
        values.control_XV_y_minus_Emax[i] = vars.control_XV_y[i] - params.Emax[i]
        values.control_XV_y_rate[i] = control_XV_y_rate
        values.control_XV_y_minus_Emin[i] = vars.control_XV_y[i] - params.Emin[i]
        values.control_Idroop_minus_ImaxF[i] = vars.control_Idroop[i] - params.ImaxF[i]
        values.control_Id_floor_max_value[i] = params.ImaxF[i] - vars.control_Idroop[i]
        values.control_pmax_out_min_value[i] = 0 - control_pmax_sum
        values.control_pmin_out_max_value[i] = 0 - control_pmin_sum
        values.control_qmax_out_min_value[i] = 0 - control_qmax_sum
        values.control_qmin_out_max_value[i] = 0 - control_qmin_sum
        values.control_e_cmd_limited_clamp_value[i] = control_Vcmd - params.Emax[i]
        values.control_e_cmd_limited_clamp_value_2[i] = control_Vcmd - params.Emin[i]
        values.control_e_pi_limited_clamp_value[i] = control_e_pi - params.Emax[i]
        values.control_e_pi_limited_clamp_value_2[i] = control_e_pi - params.Emin[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REGFM_A1},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    regfm_a1_automaton_values_kernel!(
        (control_XPMAX_y = ws.automata.values.control_XPMAX_y, control_XPMAX_y_rate = ws.automata.values.control_XPMAX_y_rate, control_XPMIN_y = ws.automata.values.control_XPMIN_y, control_XPMIN_y_rate = ws.automata.values.control_XPMIN_y_rate, control_XQMAX_y = ws.automata.values.control_XQMAX_y, control_XQMAX_y_rate = ws.automata.values.control_XQMAX_y_rate, control_XQMIN_y = ws.automata.values.control_XQMIN_y, control_XQMIN_y_rate = ws.automata.values.control_XQMIN_y_rate, control_XV_y_minus_Emax = ws.automata.values.control_XV_y_minus_Emax, control_XV_y_rate = ws.automata.values.control_XV_y_rate, control_XV_y_minus_Emin = ws.automata.values.control_XV_y_minus_Emin, control_Idroop_minus_ImaxF = ws.automata.values.control_Idroop_minus_ImaxF, control_Id_floor_max_value = ws.automata.values.control_Id_floor_max_value, control_pmax_out_min_value = ws.automata.values.control_pmax_out_min_value, control_pmin_out_max_value = ws.automata.values.control_pmin_out_max_value, control_qmax_out_min_value = ws.automata.values.control_qmax_out_min_value, control_qmin_out_max_value = ws.automata.values.control_qmin_out_max_value, control_e_cmd_limited_clamp_value = ws.automata.values.control_e_cmd_limited_clamp_value, control_e_cmd_limited_clamp_value_2 = ws.automata.values.control_e_cmd_limited_clamp_value_2, control_e_pi_limited_clamp_value = ws.automata.values.control_e_pi_limited_clamp_value, control_e_pi_limited_clamp_value_2 = ws.automata.values.control_e_pi_limited_clamp_value_2),
        (XL = model.XL, mq = model.mq, kpv = model.kpv, kiv = model.kiv, Emax = model.Emax, Emin = model.Emin, mp = model.mp, Pmax = model.Pmax, Pmin = model.Pmin, kppmax = model.kppmax, kipmax = model.kipmax, w0 = model.w0, TPf = model.TPf, TQf = model.TQf, TVf = model.TVf, Qmax = model.Qmax, Qmin = model.Qmin, VFlag = model.VFlag, QVFlag = model.QVFlag, kpqmax = model.kpqmax, kiqmax = model.kiqmax, ImaxF = model.ImaxF, control_Pref0 = ws.setpoints.control_Pref0, control_Qref0 = ws.setpoints.control_Qref0, control_Vref0 = ws.setpoints.control_Vref0, control_current_limited = ws.automata.discrete.control_current_limited, control_XPMAX_y_gate = ws.automata.discrete.control_XPMAX_y_gate, control_XPMIN_y_gate = ws.automata.discrete.control_XPMIN_y_gate, control_XQMAX_y_gate = ws.automata.discrete.control_XQMAX_y_gate, control_XQMIN_y_gate = ws.automata.discrete.control_XQMIN_y_gate, control_XV_y_gate = ws.automata.discrete.control_XV_y_gate, control_Id_floor_max_second = ws.automata.discrete.control_Id_floor_max_second, control_Id_floor_max_first = ws.automata.discrete.control_Id_floor_max_first, control_pmax_out_min_second = ws.automata.discrete.control_pmax_out_min_second, control_pmax_out_min_first = ws.automata.discrete.control_pmax_out_min_first, control_pmin_out_max_second = ws.automata.discrete.control_pmin_out_max_second, control_pmin_out_max_first = ws.automata.discrete.control_pmin_out_max_first, control_qmax_out_min_second = ws.automata.discrete.control_qmax_out_min_second, control_qmax_out_min_first = ws.automata.discrete.control_qmax_out_min_first, control_qmin_out_max_second = ws.automata.discrete.control_qmin_out_max_second, control_qmin_out_max_first = ws.automata.discrete.control_qmin_out_max_first, control_e_cmd_limited_clamp_high = ws.automata.discrete.control_e_cmd_limited_clamp_high, control_e_cmd_limited_clamp_low = ws.automata.discrete.control_e_cmd_limited_clamp_low, control_e_cmd_limited_clamp_mid = ws.automata.discrete.control_e_cmd_limited_clamp_mid, control_e_pi_limited_clamp_high = ws.automata.discrete.control_e_pi_limited_clamp_high, control_e_pi_limited_clamp_low = ws.automata.discrete.control_e_pi_limited_clamp_low, control_e_pi_limited_clamp_mid = ws.automata.discrete.control_e_pi_limited_clamp_mid, Mratio = external.Mratio),
        (control_P_f_y = ws.var_vals.control_P_f_y, control_Q_f_y = ws.var_vals.control_Q_f_y, control_V_f_y = ws.var_vals.control_V_f_y, control_XPMAX_y = ws.var_vals.control_XPMAX_y, control_XPMIN_y = ws.var_vals.control_XPMIN_y, control_XQMAX_y = ws.var_vals.control_XQMAX_y, control_XQMIN_y = ws.var_vals.control_XQMIN_y, control_DELTA_y = ws.var_vals.control_DELTA_y, control_XV_y = ws.var_vals.control_XV_y, control_Idroop = ws.var_vals.control_Idroop, control_Pref = ws.var_vals.control_Pref, control_Qref = ws.var_vals.control_Qref, control_Vref = ws.var_vals.control_Vref, Pe = ws.var_vals.Pe, Qe = ws.var_vals.Qe, v = ws.var_vals.v, theta = ws.var_vals.theta),
        t,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{REGFM_A1{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, REGFM_A1)

    v = T.(data.v)
    theta = T.(data.theta)
    p_gen = T.(data.p_gen)
    q_gen = T.(data.q_gen)
    Mratio = data.Mratio

    # === Pre-allocate Result Arrays ===
    control_P_f_y = Vector{T}(undef, n)
    control_Q_f_y = Vector{T}(undef, n)
    control_V_f_y = Vector{T}(undef, n)
    control_XPMAX_y = fill(zero(T), n)
    control_XPMIN_y = fill(zero(T), n)
    control_XQMAX_y = fill(zero(T), n)
    control_XQMIN_y = fill(zero(T), n)
    control_DELTA_y = Vector{T}(undef, n)
    control_XV_y = Vector{T}(undef, n)
    control_Idroop = Vector{T}(undef, n)
    control_Pref = Vector{T}(undef, n)
    control_Qref = Vector{T}(undef, n)
    control_Vref = Vector{T}(undef, n)
    Pe = Vector{T}(undef, n)
    Qe = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        v[i] > 0 || throw(ArgumentError("REGFM_A1 instance $i cannot start: it requires that the power-flow terminal voltage is positive: the start divides the powers by it (check control.terminal_voltage_positive)"))

        # State initial values
        control_P_f_y[i] = p_gen[i] / Mratio[i]
        control_Q_f_y[i] = q_gen[i] / Mratio[i]
        control_V_f_y[i] = v[i]
        control_DELTA_y[i] = theta[i] + atan((model.XL[i] * 1 * control_P_f_y[i] - 0 * control_Q_f_y[i]) / v[i], v[i] + (0 * control_P_f_y[i] + model.XL[i] * 1 * control_Q_f_y[i]) / v[i])
        control_XV_y[i] = sqrt((v[i] + (0 * control_P_f_y[i] + model.XL[i] * 1 * control_Q_f_y[i]) / v[i]) ^ 2 + ((model.XL[i] * 1 * control_P_f_y[i] - 0 * control_Q_f_y[i]) / v[i]) ^ 2)

        # Algebraic initial values
        control_Pref[i] = p_gen[i] / Mratio[i] + 0
        control_Qref[i] = ifelse(model.QVFlag[i] == 0, control_Q_f_y[i], 0)
        control_Vref[i] = ifelse(model.VFlag[i] == 1, v[i], control_XV_y[i]) - model.mq[i] * (ifelse(model.QVFlag[i] == 0, control_Q_f_y[i], 0) - control_Q_f_y[i])
        control_Idroop[i] = sqrt(0.0 ^ 2 + (model.XL[i] * 1) ^ 2) * sqrt(control_P_f_y[i] ^ 2 + control_Q_f_y[i] ^ 2) / v[i] / model.XL[i]
        Pe[i] = p_gen[i] / Mratio[i]
        Qe[i] = q_gen[i] / Mratio[i]
    end

    # === Return Initial Values ===
    return (control_P_f_y=control_P_f_y, control_Q_f_y=control_Q_f_y, control_V_f_y=control_V_f_y, control_XPMAX_y=control_XPMAX_y, control_XPMIN_y=control_XPMIN_y, control_XQMAX_y=control_XQMAX_y, control_XQMIN_y=control_XQMIN_y, control_DELTA_y=control_DELTA_y, control_XV_y=control_XV_y, control_Idroop=control_Idroop, control_Pref=control_Pref, control_Qref=control_Qref, control_Vref=control_Vref, Pe=Pe, Qe=Qe)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{REGFM_A1}, dyr_data, sys_pf, links) -> StructArray{REGFM_A1}

Parse REGFM_A1 model from DYR data with linking to power flow models.
"""
function parse_dynamic_model(
    ::Type{REGFM_A1},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    links::ModelLinks
)
    n = length(links)
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{REGFM_A1{T}}(undef, n)

    for idx in 1:n
        dyr_row = dyr_data[links.source_idx[idx]]
        instances[idx] = check_parameter_domains(REGFM_A1(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            XL = T(dyr_row.XL),
            mq = T(dyr_row.mq),
            kpv = T(dyr_row.kpv),
            kiv = T(dyr_row.kiv),
            Emax = T(dyr_row.Emax),
            Emin = T(dyr_row.Emin),
            mp = T(dyr_row.mp),
            Pmax = T(dyr_row.Pmax),
            Pmin = T(dyr_row.Pmin),
            kppmax = T(dyr_row.kppmax),
            kipmax = T(dyr_row.kipmax),
            w0 = T(dyr_row.w0),
            TPf = T(dyr_row.TPf),
            TQf = T(dyr_row.TQf),
            TVf = T(dyr_row.TVf),
            Qmax = T(dyr_row.Qmax),
            Qmin = T(dyr_row.Qmin),
            VFlag = T(dyr_row.VFlag),
            QVFlag = T(dyr_row.QVFlag),
            kpqmax = T(dyr_row.kpqmax),
            kiqmax = T(dyr_row.kiqmax),
            ImaxF = T(dyr_row.ImaxF)
        ), "dyr record $(idx)")
    end

    return StructArray(instances)
end
