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
    x.Emin <= x.Emax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint xv_limits_ordered: xv: the lower limit of xv is at most its upper limit"))
    x.XL > 0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint coupling_reactance: the coupling reactance XL is positive ([P23] Table 1: 0.05 pu to 0.25 pu)"))
    x.ImaxF > 0 || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint current_limit: the maximum transient current ImaxF is positive ([P23] Table 1: 1.5 pu to 3 pu)"))
    x.Emin < x.Emax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint voltage_limits: Emin is below Emax: the voltage loop's limits are ordered ([P23] Table 1)"))
    x.Pmin < x.Pmax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint active_power_limits: Pmin is below Pmax: the active power limits are ordered ([P23] Table 1)"))
    x.Qmin < x.Qmax || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint reactive_power_limits: Qmin is below Qmax: the reactive power limits are ordered ([P23] Table 1)"))
    ifelse(x.kppmax >= 0, x.kipmax > 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint overload_gains: the overload mitigation controller's integral gain is positive and its proportional gain not negative: the limits it holds its integrators at act on the sign of the integral's rate"))
    ifelse(x.kpqmax >= 0, x.kiqmax > 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint reactive_limit_gains: the Qmax/Qmin controller's integral gain is positive and its proportional gain not negative: the limits it holds its integrators at act on the sign of the integral's rate"))
    ifelse(x.kiv >= 0, x.kpv >= 0, false) || throw(ArgumentError("REGFM_A1 ($from) does not satisfy its constraint voltage_gains: the voltage controller's gains are not negative ([P23] Table 1)"))
    (x.XL != 0) || throw(ArgumentError("REGFM_A1 ($from): XL is zero, and the model divides by it (graph[1].Id_p)"))
    (x.TPf != 0) || throw(ArgumentError("REGFM_A1 ($from): TPf is zero, and the model divides by it (states.Pinv.derivative)"))
    (x.TQf != 0) || throw(ArgumentError("REGFM_A1 ($from): TQf is zero, and the model divides by it (states.Qinv.derivative)"))
    (x.TVf != 0) || throw(ArgumentError("REGFM_A1 ($from): TVf is zero, and the model divides by it (states.Vinv.derivative)"))
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
        :Pinv => :Pinv, # Same name
        :Pref => :Pref, # Same name
        :Qe => :Qe, # Same name
        :Qinv => :Qinv, # Same name
        :Qref => :Qref, # Same name
        :Vinv => :Vinv, # Same name
        :Vref => :Vref, # Same name
        :delta => :delta, # Same name
        :xpmax => :xpmax, # Same name
        :xpmin => :xpmin, # Same name
        :xqmax => :xqmax, # Same name
        :xqmin => :xqmin, # Same name
        :xv => :xv, # Same name
    )
end

# Map group canonical residual names to REGFM_A1's actual residual names
function group_res_mapping(::Type{<:REGFM_A1})
    Dict(
        :g_Pref => :g_Pref, # Same name
        :g_Qref => :g_Qref, # Same name
        :g_Vref => :g_Vref, # Same name
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
            ModelVar(:Pinv, StateVar();
                description="Filtered active power, inverter base (eq. 4)"),
            ModelVar(:Qinv, StateVar();
                description="Filtered reactive power, inverter base (eq. 5)"),
            ModelVar(:Vinv, StateVar();
                description="Filtered terminal voltage (eq. 6)"),
            ModelVar(:delta, StateVar();
                description="Angle of the internal voltage, the integral of the frequency deviation (Fig. 3a)"),
            ModelVar(:xpmax, StateVar();
                description="Integrator of the Pmax branch of the overload mitigation control, not above 0"),
            ModelVar(:xpmin, StateVar();
                description="Integrator of the Pmin branch of the overload mitigation control, not below 0"),
            ModelVar(:xqmax, StateVar();
                description="Integrator of the Qmax controller, not above 0"),
            ModelVar(:xqmin, StateVar();
                description="Integrator of the Qmin controller, not below 0"),
            ModelVar(:xv, StateVar();
                description="Integrator of the voltage controller, within [Emin, Emax]; idle where VFlag = 0"),
            ModelVar(:Pref, AlgebVar();
                description="Active power reference, inverter base (the plant controller changes it)"),
            ModelVar(:Qref, AlgebVar();
                description="Reactive power reference, inverter base (the plant controller changes it where QVFlag = 0)"),
            ModelVar(:Vref, AlgebVar();
                description="Voltage set point (the plant controller changes it where QVFlag = 1); the start makes the droop output the internal voltage the power flow needs"),
            ModelVar(:Vcmd, AlgebVar();
                description="Voltage command of the Q-V droop with the Qmax/Qmin controller (Fig. 3b)"),
            ModelVar(:Edroop, AlgebVar();
                description="Output of the droop control, within [Emin, Emax] (Fig. 3b)"),
            ModelVar(:Idroop, AlgebVar();
                description="Magnitude of the droop current of eq. 7, inverter base"),
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
        ModelResidual(:dPinv, StateRes(), PrivateRes();
            description="Filtered active power, inverter base (eq. 4)"),
        ModelResidual(:dQinv, StateRes(), PrivateRes();
            description="Filtered reactive power, inverter base (eq. 5)"),
        ModelResidual(:dVinv, StateRes(), PrivateRes();
            description="Filtered terminal voltage (eq. 6)"),
        ModelResidual(:ddelta, StateRes(), PrivateRes();
            description="Angle of the internal voltage, the integral of the frequency deviation (Fig. 3a)"),
        ModelResidual(:dxpmax, StateRes(), PrivateRes();
            description="Integrator of the Pmax branch of the overload mitigation control, not above 0"),
        ModelResidual(:dxpmin, StateRes(), PrivateRes();
            description="Integrator of the Pmin branch of the overload mitigation control, not below 0"),
        ModelResidual(:dxqmax, StateRes(), PrivateRes();
            description="Integrator of the Qmax controller, not above 0"),
        ModelResidual(:dxqmin, StateRes(), PrivateRes();
            description="Integrator of the Qmin controller, not below 0"),
        ModelResidual(:dxv, StateRes(), PrivateRes();
            description="Integrator of the voltage controller, within [Emin, Emax]; idle where VFlag = 0"),
        ModelResidual(:g_Pref, AlgebRes(), SharedRes(:add);
            description="Active power reference, inverter base (the plant controller changes it)"),
        ModelResidual(:g_Qref, AlgebRes(), SharedRes(:add);
            description="Reactive power reference, inverter base (the plant controller changes it where QVFlag = 0)"),
        ModelResidual(:g_Vref, AlgebRes(), SharedRes(:add);
            description="Voltage set point (the plant controller changes it where QVFlag = 1); the start makes the droop output the internal voltage the power flow needs"),
        ModelResidual(:g_Vcmd, AlgebRes(), PrivateRes();
            description="Voltage command of the Q-V droop with the Qmax/Qmin controller (Fig. 3b)"),
        ModelResidual(:g_Edroop, AlgebRes(), PrivateRes();
            description="Output of the droop control, within [Emin, Emax] (Fig. 3b)"),
        ModelResidual(:g_Idroop, AlgebRes(), PrivateRes();
            description="Magnitude of the droop current of eq. 7, inverter base"),
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
        InitParam(:Pref0, :Pref; description="Active power reference, inverter base (the plant controller changes it)"),
        InitParam(:Qref0, :Qref; description="Reactive power reference, inverter base (the plant controller changes it where QVFlag = 0)"),
        InitParam(:Vref0, :Vref; description="Voltage set point (the plant controller changes it where QVFlag = 1); the start makes the droop output the internal voltage the power flow needs")
    ],
        injection = InjectionForms(:current, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj]),
        variable_limits = [VariableLimit(:xpmax, nothing, m -> 0), VariableLimit(:xpmin, m -> 0, nothing), VariableLimit(:xqmax, nothing, m -> 0), VariableLimit(:xqmin, m -> 0, nothing), VariableLimit(:xv, :Emin, :Emax), VariableLimit(:Pinv, :Pmin, :Pmax), VariableLimit(:Qinv, :Qmin, :Qmax), VariableLimit(:Edroop, :Emin, :Emax), VariableLimit(:Idroop, nothing, :ImaxF)],
        automata = (
            Automaton(:xpmax;
                modes = [Mode(:free; holds = [:xpmax_gate => 1.0]),
                    Mode(:upper; holds = [:xpmax_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:xpmax, :>=), Atom(:xpmax_rate, :>)])),
                    Transition(:upper, :free; when = Atom(:xpmax_rate, :<))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:xpmax, :>=), Atom(:xpmax_rate, :>)])),
                    InitialMode(:free)]),
            Automaton(:xpmin;
                modes = [Mode(:free; holds = [:xpmin_gate => 1.0]),
                    Mode(:lower; holds = [:xpmin_gate => 0.0])],
                transitions = [Transition(:free, :lower; when = AllOf(Condition[Atom(:xpmin, :<=), Atom(:xpmin_rate, :<)])),
                    Transition(:lower, :free; when = Atom(:xpmin_rate, :>))],
                initial = [InitialMode(:lower; when = AllOf(Condition[Atom(:xpmin, :<=), Atom(:xpmin_rate, :<)])),
                    InitialMode(:free)]),
            Automaton(:xqmax;
                modes = [Mode(:free; holds = [:xqmax_gate => 1.0]),
                    Mode(:upper; holds = [:xqmax_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:xqmax, :>=), Atom(:xqmax_rate, :>)])),
                    Transition(:upper, :free; when = Atom(:xqmax_rate, :<))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:xqmax, :>=), Atom(:xqmax_rate, :>)])),
                    InitialMode(:free)]),
            Automaton(:xqmin;
                modes = [Mode(:free; holds = [:xqmin_gate => 1.0]),
                    Mode(:lower; holds = [:xqmin_gate => 0.0])],
                transitions = [Transition(:free, :lower; when = AllOf(Condition[Atom(:xqmin, :<=), Atom(:xqmin_rate, :<)])),
                    Transition(:lower, :free; when = Atom(:xqmin_rate, :>))],
                initial = [InitialMode(:lower; when = AllOf(Condition[Atom(:xqmin, :<=), Atom(:xqmin_rate, :<)])),
                    InitialMode(:free)]),
            Automaton(:xv;
                modes = [Mode(:free; holds = [:xv_gate => 1.0]),
                    Mode(:upper; holds = [:xv_gate => 0.0]),
                    Mode(:lower; holds = [:xv_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:xv_minus_Emax, :>=), Atom(:xv_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:xv_minus_Emin, :<=), Atom(:xv_rate, :<)])),
                    Transition(:upper, :free; when = Atom(:xv_rate, :<)),
                    Transition(:lower, :free; when = Atom(:xv_rate, :>))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:xv_minus_Emax, :>=), Atom(:xv_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:xv_minus_Emin, :<=), Atom(:xv_rate, :<)])),
                    InitialMode(:free)],
                at_rest = [:xv_rate]),
            Automaton(:current_limit;
                modes = [Mode(:droop; holds = [:current_limited => 0.0]),
                    Mode(:limited; holds = [:current_limited => 1.0])],
                transitions = [Transition(:droop, :limited; when = Atom(:Idroop_minus_ImaxF, :>=)),
                    Transition(:limited, :droop; when = Atom(:Idroop_minus_ImaxF, :<))],
                initial = [InitialMode(:droop)]),
            Automaton(:limit_scale_max;
                modes = [Mode(:second; holds = [:limit_scale_max_second => 1.0, :limit_scale_max_first => 0.0]),
                    Mode(:first; holds = [:limit_scale_max_second => 0.0, :limit_scale_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:limit_scale_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:limit_scale_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:limit_scale_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:pmax_out_min;
                modes = [Mode(:second; holds = [:pmax_out_min_second => 1.0, :pmax_out_min_first => 0.0]),
                    Mode(:first; holds = [:pmax_out_min_second => 0.0, :pmax_out_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:pmax_out_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:pmax_out_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:pmax_out_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:pmin_out_max;
                modes = [Mode(:second; holds = [:pmin_out_max_second => 1.0, :pmin_out_max_first => 0.0]),
                    Mode(:first; holds = [:pmin_out_max_second => 0.0, :pmin_out_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:pmin_out_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:pmin_out_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:pmin_out_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:qmax_out_min;
                modes = [Mode(:second; holds = [:qmax_out_min_second => 1.0, :qmax_out_min_first => 0.0]),
                    Mode(:first; holds = [:qmax_out_min_second => 0.0, :qmax_out_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:qmax_out_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:qmax_out_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:qmax_out_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:qmin_out_max;
                modes = [Mode(:second; holds = [:qmin_out_max_second => 1.0, :qmin_out_max_first => 0.0]),
                    Mode(:first; holds = [:qmin_out_max_second => 0.0, :qmin_out_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:qmin_out_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:qmin_out_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:qmin_out_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:Edroop_clamp;
                modes = [Mode(:high; holds = [:Edroop_clamp_high => 1.0, :Edroop_clamp_low => 0.0, :Edroop_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Edroop_clamp_high => 0.0, :Edroop_clamp_low => 1.0, :Edroop_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Edroop_clamp_high => 0.0, :Edroop_clamp_low => 0.0, :Edroop_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Edroop_clamp_value, :<=), Atom(:Edroop_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Edroop_clamp_value_2, :>=), Atom(:Edroop_clamp_value, :<=)])),
                    Transition(:low, :high; when = Atom(:Edroop_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Edroop_clamp_value_2, :>=), Atom(:Edroop_clamp_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:Edroop_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Edroop_clamp_value, :<=), Atom(:Edroop_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:Edroop_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Edroop_clamp_value, :<=), Atom(:Edroop_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:Edroop_clamp_2;
                modes = [Mode(:high; holds = [:Edroop_clamp_2_high => 1.0, :Edroop_clamp_2_low => 0.0, :Edroop_clamp_2_mid => 0.0]),
                    Mode(:low; holds = [:Edroop_clamp_2_high => 0.0, :Edroop_clamp_2_low => 1.0, :Edroop_clamp_2_mid => 0.0]),
                    Mode(:mid; holds = [:Edroop_clamp_2_high => 0.0, :Edroop_clamp_2_low => 0.0, :Edroop_clamp_2_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Edroop_clamp_2_value, :<=), Atom(:Edroop_clamp_2_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Edroop_clamp_2_value_2, :>=), Atom(:Edroop_clamp_2_value, :<=)])),
                    Transition(:low, :high; when = Atom(:Edroop_clamp_2_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Edroop_clamp_2_value_2, :>=), Atom(:Edroop_clamp_2_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:Edroop_clamp_2_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Edroop_clamp_2_value, :<=), Atom(:Edroop_clamp_2_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:Edroop_clamp_2_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Edroop_clamp_2_value, :<=), Atom(:Edroop_clamp_2_value_2, :<)])),
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

        a = vars.delta[i] - vars.theta[i]
        Id_p = vars.Edroop[i] * sin(a) / params.XL[i]
        Id_q = (vars.Edroop[i] * cos(a) - vars.v[i]) / params.XL[i]
        Id2 = Id_p ^ 2 + Id_q ^ 2
        limit_scale = 1 - params.current_limited[i] * (1 - params.ImaxF[i] / (vars.Idroop[i] * params.limit_scale_max_first[i] + params.ImaxF[i] * params.limit_scale_max_second[i]))
        ip_t = Id_p * limit_scale
        iq_t = Id_q * limit_scale
        pmax_out = params.pmax_out_min_second[i] * (vars.xpmax[i] + params.kppmax[i] * (params.Pmax[i] - vars.Pinv[i]))
        pmin_out = params.pmin_out_max_second[i] * (vars.xpmin[i] + params.kppmax[i] * (params.Pmin[i] - vars.Pinv[i]))
        qmax_out = params.qmax_out_min_second[i] * (vars.xqmax[i] + params.kpqmax[i] * (params.Qmax[i] - vars.Qinv[i]))
        qmin_out = params.qmin_out_max_second[i] * (vars.xqmin[i] + params.kpqmax[i] * (params.Qmin[i] - vars.Qinv[i]))
        Verr = vars.Vcmd[i] - vars.Vinv[i]
        xpmax_rate = params.kipmax[i] * (params.Pmax[i] - vars.Pinv[i])
        xpmin_rate = params.kipmax[i] * (params.Pmin[i] - vars.Pinv[i])
        xqmax_rate = params.kiqmax[i] * (params.Qmax[i] - vars.Qinv[i])
        xqmin_rate = params.kiqmax[i] * (params.Qmin[i] - vars.Qinv[i])
        xv_rate = params.VFlag[i] * Verr * params.kiv[i]
        # === Residual Assignments ===
        outputs.dPinv[i] = mask * ((vars.Pe[i] - vars.Pinv[i]) / params.TPf[i])
        outputs.dQinv[i] = mask * ((vars.Qe[i] - vars.Qinv[i]) / params.TQf[i])
        outputs.dVinv[i] = mask * ((vars.v[i] - vars.Vinv[i]) / params.TVf[i])
        outputs.ddelta[i] = mask * (params.w0[i] * (pmax_out + pmin_out + params.mp[i] * (vars.Pref[i] - vars.Pinv[i])))
        outputs.dxpmax[i] = mask * (params.xpmax_gate[i] * xpmax_rate)
        outputs.dxpmin[i] = mask * (params.xpmin_gate[i] * xpmin_rate)
        outputs.dxqmax[i] = mask * (params.xqmax_gate[i] * xqmax_rate)
        outputs.dxqmin[i] = mask * (params.xqmin_gate[i] * xqmin_rate)
        outputs.dxv[i] = mask * (params.xv_gate[i] * xv_rate)
        outputs.g_Pref[i] = mask * (params.Pref0[i] - vars.Pref[i])
        outputs.g_Qref[i] = mask * (params.Qref0[i] - vars.Qref[i])
        outputs.g_Vref[i] = mask * (params.Vref0[i] - vars.Vref[i])
        outputs.g_Vcmd[i] = mask * (vars.Vref[i] + qmax_out + qmin_out + params.mq[i] * (vars.Qref[i] - vars.Qinv[i]) - vars.Vcmd[i])
        outputs.g_Edroop[i] = mask * (params.VFlag[i] * (params.Edroop_clamp_high[i] * params.Emax[i] + params.Edroop_clamp_low[i] * params.Emin[i] + params.Edroop_clamp_mid[i] * (vars.xv[i] + Verr * params.kpv[i])) + (1 - params.VFlag[i]) * (params.Edroop_clamp_2_high[i] * params.Emax[i] + params.Edroop_clamp_2_low[i] * params.Emin[i] + params.Edroop_clamp_2_mid[i] * vars.Vcmd[i]) - vars.Edroop[i])
        outputs.g_Idroop[i] = mask * (ifelse(Id2 > 0, sqrt(Id2), 0) - vars.Idroop[i])
        outputs.g_Pe[i] = mask * (ip_t * vars.v[i] - vars.Pe[i])
        outputs.g_Qe[i] = mask * (iq_t * vars.v[i] - vars.Qe[i])
        outputs.ip_inj[i] = mask * (params.Mratio[i] * ip_t)
        outputs.iq_inj[i] = mask * (params.Mratio[i] * iq_t)
        outputs.p_inj[i] = mask * (params.Mratio[i] * ip_t * vars.v[i])
        outputs.q_inj[i] = mask * (params.Mratio[i] * iq_t * vars.v[i])
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REGFM_A1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    regfm_a1_residual_kernel!(
        (dPinv = ws.res_vals.dPinv, dQinv = ws.res_vals.dQinv, dVinv = ws.res_vals.dVinv, ddelta = ws.res_vals.ddelta, dxpmax = ws.res_vals.dxpmax, dxpmin = ws.res_vals.dxpmin, dxqmax = ws.res_vals.dxqmax, dxqmin = ws.res_vals.dxqmin, dxv = ws.res_vals.dxv, g_Pref = ws.res_vals.g_Pref, g_Qref = ws.res_vals.g_Qref, g_Vref = ws.res_vals.g_Vref, g_Vcmd = ws.res_vals.g_Vcmd, g_Edroop = ws.res_vals.g_Edroop, g_Idroop = ws.res_vals.g_Idroop, g_Pe = ws.res_vals.g_Pe, g_Qe = ws.res_vals.g_Qe, ip_inj = ws.res_vals.ip_inj, iq_inj = ws.res_vals.iq_inj, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj),
        (XL = model.XL, mq = model.mq, kpv = model.kpv, kiv = model.kiv, Emax = model.Emax, Emin = model.Emin, mp = model.mp, Pmax = model.Pmax, Pmin = model.Pmin, kppmax = model.kppmax, kipmax = model.kipmax, w0 = model.w0, TPf = model.TPf, TQf = model.TQf, TVf = model.TVf, Qmax = model.Qmax, Qmin = model.Qmin, VFlag = model.VFlag, QVFlag = model.QVFlag, kpqmax = model.kpqmax, kiqmax = model.kiqmax, ImaxF = model.ImaxF, Pref0 = ws.setpoints.Pref0, Qref0 = ws.setpoints.Qref0, Vref0 = ws.setpoints.Vref0, current_limited = ws.automata.discrete.current_limited, xpmax_gate = ws.automata.discrete.xpmax_gate, xpmin_gate = ws.automata.discrete.xpmin_gate, xqmax_gate = ws.automata.discrete.xqmax_gate, xqmin_gate = ws.automata.discrete.xqmin_gate, xv_gate = ws.automata.discrete.xv_gate, limit_scale_max_second = ws.automata.discrete.limit_scale_max_second, limit_scale_max_first = ws.automata.discrete.limit_scale_max_first, pmax_out_min_second = ws.automata.discrete.pmax_out_min_second, pmax_out_min_first = ws.automata.discrete.pmax_out_min_first, pmin_out_max_second = ws.automata.discrete.pmin_out_max_second, pmin_out_max_first = ws.automata.discrete.pmin_out_max_first, qmax_out_min_second = ws.automata.discrete.qmax_out_min_second, qmax_out_min_first = ws.automata.discrete.qmax_out_min_first, qmin_out_max_second = ws.automata.discrete.qmin_out_max_second, qmin_out_max_first = ws.automata.discrete.qmin_out_max_first, Edroop_clamp_high = ws.automata.discrete.Edroop_clamp_high, Edroop_clamp_low = ws.automata.discrete.Edroop_clamp_low, Edroop_clamp_mid = ws.automata.discrete.Edroop_clamp_mid, Edroop_clamp_2_high = ws.automata.discrete.Edroop_clamp_2_high, Edroop_clamp_2_low = ws.automata.discrete.Edroop_clamp_2_low, Edroop_clamp_2_mid = ws.automata.discrete.Edroop_clamp_2_mid, Mratio = external.Mratio),
        (Pinv = ws.var_vals.Pinv, Qinv = ws.var_vals.Qinv, Vinv = ws.var_vals.Vinv, delta = ws.var_vals.delta, xpmax = ws.var_vals.xpmax, xpmin = ws.var_vals.xpmin, xqmax = ws.var_vals.xqmax, xqmin = ws.var_vals.xqmin, xv = ws.var_vals.xv, Pref = ws.var_vals.Pref, Qref = ws.var_vals.Qref, Vref = ws.var_vals.Vref, Vcmd = ws.var_vals.Vcmd, Edroop = ws.var_vals.Edroop, Idroop = ws.var_vals.Idroop, Pe = ws.var_vals.Pe, Qe = ws.var_vals.Qe, v = ws.var_vals.v, theta = ws.var_vals.theta),
        active_mask,
        t,
        AT
    )
end

@inline function regfm_a1_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.xpmax)
        Verr = vars.Vcmd[i] - vars.Vinv[i]
        xpmax_rate = params.kipmax[i] * (params.Pmax[i] - vars.Pinv[i])
        xpmin_rate = params.kipmax[i] * (params.Pmin[i] - vars.Pinv[i])
        xqmax_rate = params.kiqmax[i] * (params.Qmax[i] - vars.Qinv[i])
        xqmin_rate = params.kiqmax[i] * (params.Qmin[i] - vars.Qinv[i])
        xv_rate = params.VFlag[i] * Verr * params.kiv[i]
        # === Automaton Values ===
        values.xpmax[i] = vars.xpmax[i]
        values.xpmax_rate[i] = xpmax_rate
        values.xpmin[i] = vars.xpmin[i]
        values.xpmin_rate[i] = xpmin_rate
        values.xqmax[i] = vars.xqmax[i]
        values.xqmax_rate[i] = xqmax_rate
        values.xqmin[i] = vars.xqmin[i]
        values.xqmin_rate[i] = xqmin_rate
        values.xv_minus_Emax[i] = vars.xv[i] - params.Emax[i]
        values.xv_rate[i] = xv_rate
        values.xv_minus_Emin[i] = vars.xv[i] - params.Emin[i]
        values.Idroop_minus_ImaxF[i] = vars.Idroop[i] - params.ImaxF[i]
        values.limit_scale_max_value[i] = params.ImaxF[i] - vars.Idroop[i]
        values.pmax_out_min_value[i] = vars.xpmax[i] + params.kppmax[i] * (params.Pmax[i] - vars.Pinv[i])
        values.pmin_out_max_value[i] = vars.xpmin[i] + params.kppmax[i] * (params.Pmin[i] - vars.Pinv[i])
        values.qmax_out_min_value[i] = vars.xqmax[i] + params.kpqmax[i] * (params.Qmax[i] - vars.Qinv[i])
        values.qmin_out_max_value[i] = vars.xqmin[i] + params.kpqmax[i] * (params.Qmin[i] - vars.Qinv[i])
        values.Edroop_clamp_value[i] = vars.xv[i] + Verr * params.kpv[i] - params.Emax[i]
        values.Edroop_clamp_value_2[i] = vars.xv[i] + Verr * params.kpv[i] - params.Emin[i]
        values.Edroop_clamp_2_value[i] = vars.Vcmd[i] - params.Emax[i]
        values.Edroop_clamp_2_value_2[i] = vars.Vcmd[i] - params.Emin[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REGFM_A1},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    regfm_a1_automaton_values_kernel!(
        (xpmax = ws.automata.values.xpmax, xpmax_rate = ws.automata.values.xpmax_rate, xpmin = ws.automata.values.xpmin, xpmin_rate = ws.automata.values.xpmin_rate, xqmax = ws.automata.values.xqmax, xqmax_rate = ws.automata.values.xqmax_rate, xqmin = ws.automata.values.xqmin, xqmin_rate = ws.automata.values.xqmin_rate, xv_minus_Emax = ws.automata.values.xv_minus_Emax, xv_rate = ws.automata.values.xv_rate, xv_minus_Emin = ws.automata.values.xv_minus_Emin, Idroop_minus_ImaxF = ws.automata.values.Idroop_minus_ImaxF, limit_scale_max_value = ws.automata.values.limit_scale_max_value, pmax_out_min_value = ws.automata.values.pmax_out_min_value, pmin_out_max_value = ws.automata.values.pmin_out_max_value, qmax_out_min_value = ws.automata.values.qmax_out_min_value, qmin_out_max_value = ws.automata.values.qmin_out_max_value, Edroop_clamp_value = ws.automata.values.Edroop_clamp_value, Edroop_clamp_value_2 = ws.automata.values.Edroop_clamp_value_2, Edroop_clamp_2_value = ws.automata.values.Edroop_clamp_2_value, Edroop_clamp_2_value_2 = ws.automata.values.Edroop_clamp_2_value_2),
        (XL = model.XL, mq = model.mq, kpv = model.kpv, kiv = model.kiv, Emax = model.Emax, Emin = model.Emin, mp = model.mp, Pmax = model.Pmax, Pmin = model.Pmin, kppmax = model.kppmax, kipmax = model.kipmax, w0 = model.w0, TPf = model.TPf, TQf = model.TQf, TVf = model.TVf, Qmax = model.Qmax, Qmin = model.Qmin, VFlag = model.VFlag, QVFlag = model.QVFlag, kpqmax = model.kpqmax, kiqmax = model.kiqmax, ImaxF = model.ImaxF, Pref0 = ws.setpoints.Pref0, Qref0 = ws.setpoints.Qref0, Vref0 = ws.setpoints.Vref0, current_limited = ws.automata.discrete.current_limited, xpmax_gate = ws.automata.discrete.xpmax_gate, xpmin_gate = ws.automata.discrete.xpmin_gate, xqmax_gate = ws.automata.discrete.xqmax_gate, xqmin_gate = ws.automata.discrete.xqmin_gate, xv_gate = ws.automata.discrete.xv_gate, limit_scale_max_second = ws.automata.discrete.limit_scale_max_second, limit_scale_max_first = ws.automata.discrete.limit_scale_max_first, pmax_out_min_second = ws.automata.discrete.pmax_out_min_second, pmax_out_min_first = ws.automata.discrete.pmax_out_min_first, pmin_out_max_second = ws.automata.discrete.pmin_out_max_second, pmin_out_max_first = ws.automata.discrete.pmin_out_max_first, qmax_out_min_second = ws.automata.discrete.qmax_out_min_second, qmax_out_min_first = ws.automata.discrete.qmax_out_min_first, qmin_out_max_second = ws.automata.discrete.qmin_out_max_second, qmin_out_max_first = ws.automata.discrete.qmin_out_max_first, Edroop_clamp_high = ws.automata.discrete.Edroop_clamp_high, Edroop_clamp_low = ws.automata.discrete.Edroop_clamp_low, Edroop_clamp_mid = ws.automata.discrete.Edroop_clamp_mid, Edroop_clamp_2_high = ws.automata.discrete.Edroop_clamp_2_high, Edroop_clamp_2_low = ws.automata.discrete.Edroop_clamp_2_low, Edroop_clamp_2_mid = ws.automata.discrete.Edroop_clamp_2_mid, Mratio = external.Mratio),
        (Pinv = ws.var_vals.Pinv, Qinv = ws.var_vals.Qinv, Vinv = ws.var_vals.Vinv, delta = ws.var_vals.delta, xpmax = ws.var_vals.xpmax, xpmin = ws.var_vals.xpmin, xqmax = ws.var_vals.xqmax, xqmin = ws.var_vals.xqmin, xv = ws.var_vals.xv, Pref = ws.var_vals.Pref, Qref = ws.var_vals.Qref, Vref = ws.var_vals.Vref, Vcmd = ws.var_vals.Vcmd, Edroop = ws.var_vals.Edroop, Idroop = ws.var_vals.Idroop, Pe = ws.var_vals.Pe, Qe = ws.var_vals.Qe, v = ws.var_vals.v, theta = ws.var_vals.theta),
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
    Pinv = Vector{T}(undef, n)
    Qinv = Vector{T}(undef, n)
    Vinv = Vector{T}(undef, n)
    delta = Vector{T}(undef, n)
    xpmax = fill(zero(T), n)
    xpmin = fill(zero(T), n)
    xqmax = fill(zero(T), n)
    xqmin = fill(zero(T), n)
    xv = Vector{T}(undef, n)
    Pref = Vector{T}(undef, n)
    Qref = Vector{T}(undef, n)
    Vref = Vector{T}(undef, n)
    Vcmd = Vector{T}(undef, n)
    Edroop = Vector{T}(undef, n)
    Idroop = Vector{T}(undef, n)
    Pe = Vector{T}(undef, n)
    Qe = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # State initial values
        Pinv[i] = p_gen[i] / Mratio[i]
        Qinv[i] = q_gen[i] / Mratio[i]
        Vinv[i] = v[i]
        delta[i] = theta[i] + atan(Pinv[i] * model.XL[i] / v[i], v[i] + Qinv[i] * model.XL[i] / v[i])

        # Algebraic initial values
        Edroop[i] = sqrt((v[i] + Qinv[i] * model.XL[i] / v[i]) ^ 2 + Pinv[i] ^ 2 * model.XL[i] ^ 2 / v[i] ^ 2)

        # State initial values
        xv[i] = Edroop[i]

        # Algebraic initial values
        Pref[i] = Pinv[i]
        Qref[i] = ifelse(model.QVFlag[i] == 0, Qinv[i], 0)
        Vref[i] = ifelse(model.VFlag[i] == 1, v[i], Edroop[i]) - model.mq[i] * (Qref[i] - Qinv[i])
        Vcmd[i] = ifelse(model.VFlag[i] == 1, v[i], Edroop[i])
        Idroop[i] = sqrt(Pinv[i] ^ 2 + Qinv[i] ^ 2) / v[i]
        Pe[i] = p_gen[i] / Mratio[i]
        Qe[i] = q_gen[i] / Mratio[i]
    end

    # === Return Initial Values ===
    return (Pinv=Pinv, Qinv=Qinv, Vinv=Vinv, delta=delta, xpmax=xpmax, xpmin=xpmin, xqmax=xqmax, xqmin=xqmin, xv=xv, Pref=Pref, Qref=Qref, Vref=Vref, Vcmd=Vcmd, Edroop=Edroop, Idroop=Idroop, Pe=Pe, Qe=Qe)
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
