"""REECA1: Renewable Energy Electrical Control Model A (REEC_A)

The current limit logic of Appendix B has a square root that is not Lipschitz at the current limit: with Q priority (PQFLAG = 0) Ipmax = min(VDL2, sqrt(Imax^2 - Iqcmd^2)), whose derivative with respect to Iqcmd is unbounded as |Iqcmd| reaches Imax, which it can where VDL1 >= Imax at the filtered voltage; with P priority (PQFLAG = 1) the same holds of Iqmax where VDL2 >= Imax and Ipcmd reaches Imax. It is a property of the specification, kept as written: variable-step integration takes ever smaller steps near that point.

DO NOT EDIT - changes will be overwritten.
"""


using InlineStrings
using PowerfulCore: Automaton, Mode, Transition, InitialMode, Condition, Atom, AllOf, AnyOf, MemoryQuantity

export REECA1

struct REECA1{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number where the renewable generator is connected"
    BUS::Int32
    "Generator identifier"
    ID::String3
    "Remote bus for voltage control; stored, the terminal bus voltage is used"
    VREG::Int32

    # === Parameters ===
    PFFLAG::T  # Power factor control flag: 1 power factor control, 0 Q control
    VFLAG::T  # Voltage control flag: 1 Q control, 0 voltage control
    QFLAG::T  # Reactive power control flag: 1 voltage/Q control, 0 constant pf or Q control
    PFLAG::T  # Power reference flag: 1 multiplied by the speed, 0 not
    PQFLAG::T  # Current limit priority: 0 Q priority, 1 P priority
    Vdip::T  # Voltage below which the voltage dip logic is activated
    Vup::T  # Voltage above which the voltage dip logic is activated
    Trv::T  # Filter time constant for voltage measurement
    dbd1::T  # Deadband of the voltage error, overvoltage side (<= 0)
    dbd2::T  # Deadband of the voltage error, undervoltage side (>= 0)
    Kqv::T  # Gain of the reactive current injection
    Iqh1::T  # Maximum reactive current injection
    Iql1::T  # Minimum reactive current injection
    Vref0::T  # Reference of the voltage error; 0 uses the initial terminal voltage
    Iqfrz::T  # Value Iqinj is held at for Thld seconds after a dip, if Thld > 0
    Thld::T  # Hold of the injection after a dip: > 0 at Iqfrz, < 0 in State 1, 0 none
    Thld2::T  # Time Ipmax is held at its dip value after a dip
    Tp::T  # Filter time constant for electrical power measurement
    QMax::T  # Reactive power limit maximum
    QMin::T  # Reactive power limit minimum
    VMAX::T  # Voltage control maximum
    VMIN::T  # Voltage control minimum
    Kqp::T  # Proportional gain of the Q PI controller (s2)
    Kqi::T  # Integral gain of the Q PI controller (s2)
    Kvp::T  # Proportional gain of the V PI controller (s3)
    Kvi::T  # Integral gain of the V PI controller (s3)
    Vref1::T  # Reference/bias of the inner voltage control (VFLAG = 0)
    Tiq::T  # Time constant of the reactive current lag (s4)
    dPmax::T  # Maximum ramp rate of the power reference (pu/s)
    dPmin::T  # Minimum ramp rate of the power reference (pu/s)
    PMAX::T  # Maximum power order
    PMIN::T  # Minimum power order
    Imax::T  # Maximum total converter current
    Tpord::T  # Time constant of the power order (s5)
    Vq1::T  # VDL1 point 1, voltage
    Iq1::T  # VDL1 point 1, reactive current limit
    Vq2::T  # VDL1 point 2, voltage
    Iq2::T  # VDL1 point 2, reactive current limit
    Vq3::T  # VDL1 point 3, voltage
    Iq3::T  # VDL1 point 3, reactive current limit
    Vq4::T  # VDL1 point 4, voltage
    Iq4::T  # VDL1 point 4, reactive current limit
    Vp1::T  # VDL2 point 1, voltage
    Ip1::T  # VDL2 point 1, active current limit
    Vp2::T  # VDL2 point 2, voltage
    Ip2::T  # VDL2 point 2, active current limit
    Vp3::T  # VDL2 point 3, voltage
    Ip3::T  # VDL2 point 3, active current limit
    Vp4::T  # VDL2 point 4, voltage
    Ip4::T  # VDL2 point 4, active current limit

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function REECA1(; BUS::Int32, ID, VREG::Int32, PFFLAG, VFLAG, QFLAG, PFLAG, PQFLAG, Vdip, Vup, Trv, dbd1, dbd2, Kqv, Iqh1, Iql1, Vref0, Iqfrz, Thld, Thld2, Tp, QMax, QMin, VMAX, VMIN, Kqp, Kqi, Kvp, Kvi, Vref1, Tiq, dPmax, dPmin, PMAX, PMIN, Imax, Tpord, Vq1, Iq1, Vq2, Iq2, Vq3, Iq3, Vq4, Iq4, Vp1, Ip1, Vp2, Ip2, Vp3, Ip3, Vp4, Ip4, status::Bool=true)
    T = typeof(PFFLAG)
    check_parameter_domains(REECA1{T}(BUS, ID, VREG, T(PFFLAG), T(VFLAG), T(QFLAG), T(PFLAG), T(PQFLAG), T(Vdip), T(Vup), T(Trv), T(dbd1), T(dbd2), T(Kqv), T(Iqh1), T(Iql1), T(Vref0), T(Iqfrz), T(Thld), T(Thld2), T(Tp), T(QMax), T(QMin), T(VMAX), T(VMIN), T(Kqp), T(Kqi), T(Kvp), T(Kvi), T(Vref1), T(Tiq), T(dPmax), T(dPmin), T(PMAX), T(PMIN), T(Imax), T(Tpord), T(Vq1), T(Iq1), T(Vq2), T(Iq2), T(Vq3), T(Iq3), T(Vq4), T(Iq4), T(Vp1), T(Ip1), T(Vp2), T(Ip2), T(Vp3), T(Ip3), T(Vp4), T(Ip4), status),
                            "keyword construction")
end

# Parameter domains, checked where data enters
function check_parameter_domains(x::REECA1, from::AbstractString)
    (x.PFFLAG in (0.0, 1.0,)) || throw(ArgumentError("REECA1 ($from): PFFLAG = $(x.PFFLAG) is not {0, 1}"))
    (x.VFLAG in (0.0, 1.0,)) || throw(ArgumentError("REECA1 ($from): VFLAG = $(x.VFLAG) is not {0, 1}"))
    (x.QFLAG in (0.0, 1.0,)) || throw(ArgumentError("REECA1 ($from): QFLAG = $(x.QFLAG) is not {0, 1}"))
    (x.PFLAG in (0.0, 1.0,)) || throw(ArgumentError("REECA1 ($from): PFLAG = $(x.PFLAG) is not {0, 1}"))
    (x.PQFLAG in (0.0, 1.0,)) || throw(ArgumentError("REECA1 ($from): PQFLAG = $(x.PQFLAG) is not {0, 1}"))
    (x.Thld2 >= 0.0) || throw(ArgumentError("REECA1 ($from): Thld2 = $(x.Thld2) is not >= 0"))
    x.VMIN <= x.VMAX || throw(ArgumentError("REECA1 ($from) does not satisfy its constraint PIQ_limits_ordered: PIQ: the lower limit of its output is at most its upper limit"))
    x.PMIN <= x.PMAX || throw(ArgumentError("REECA1 ($from) does not satisfy its constraint s5_y_limits_ordered: s5_y: the lower limit of s5_y is at most its upper limit"))
    ifelse(x.Vq1 < x.Vq2, ifelse(x.Vq2 < x.Vq3, x.Vq3 < x.Vq4, false), false) || throw(ArgumentError("REECA1 ($from) does not satisfy its constraint vdl1_increasing: VDL1: the voltages Vq1 < Vq2 < Vq3 < Vq4 increase"))
    ifelse(x.Vp1 < x.Vp2, ifelse(x.Vp2 < x.Vp3, x.Vp3 < x.Vp4, false), false) || throw(ArgumentError("REECA1 ($from) does not satisfy its constraint vdl2_increasing: VDL2: the voltages Vp1 < Vp2 < Vp3 < Vp4 increase"))
    x.Vdip < x.Vup || @warn "REECA1 ($from) does not satisfy its constraint dip_band: the dip band: Vdip is below Vup ([E18] Appendix F, vup > vdip)"
    x.dbd1 <= x.dbd2 || @warn "REECA1 ($from) does not satisfy its constraint deadband: the deadband: dbd1 is at most dbd2 ([E18] Appendix F, dbd2 >= dbd1)"
    x.Iql1 < x.Iqh1 || @warn "REECA1 ($from) does not satisfy its constraint injection_limits: the injection limits: Iql1 is below Iqh1 ([E18] Appendix F, iqh1 > iql1)"
    x.Iqfrz >= 0 || @warn "REECA1 ($from) does not satisfy its constraint injection_held_value: the held injection Iqfrz is not negative ([E18] Appendix F, iqfrz >= 0)"
    x.QMin < x.QMax || @warn "REECA1 ($from) does not satisfy its constraint reactive_power_limits: the reactive power limits: QMin is below QMax ([E18] Appendix F, qmax > qmin)"
    x.dPmin < x.dPmax || @warn "REECA1 ($from) does not satisfy its constraint ramp_limits: the ramp rates: dPmin is below dPmax ([E18] Appendix F, dpmax > dpmin)"
    x.Imax > 0 || throw(ArgumentError("REECA1 ($from) does not satisfy its constraint current_limit: the converter current limit Imax is positive ([E18] Appendix F, imax > 0)"))
    ifelse(x.Kqi >= 0, ifelse(x.Kqp >= 0, ifelse(x.Kqv >= 0, ifelse(x.Kvi >= 0, x.Kvp >= 0, false), false), false), false) || @warn "REECA1 ($from) does not satisfy its constraint gains: the gains Kqv, Kqp, Kqi, Kvp and Kvi are not negative ([E18] Appendix F)"
    (x.Tpord != 0) || throw(ArgumentError("REECA1 ($from): Tpord is zero, and the model divides by it (states.s5_y.derivative)"))
    (x.Tiq != 0) || throw(ArgumentError("REECA1 ($from): Tiq is zero, and the model divides by it (states.s4_y.derivative)"))
    (x.Trv != 0) || throw(ArgumentError("REECA1 ($from): Trv is zero, and the model divides by it (graph[0].s0)"))
    (x.Tp != 0) || throw(ArgumentError("REECA1 ($from): Tp is zero, and the model divides by it (graph[1].S1)"))
    return x
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:REECA1}) = true
get_dyr_model_name(::Type{<:REECA1}) = "REECA1"
get_controlled_group(::Type{<:REECA1}) = :ren_gen


# ============================================================================
# Group Interface Implementation - :ren_exciter
# ============================================================================

# Declare that REECA1 implements the :ren_exciter group interface
implements_group(::Type{<:REECA1}) = :ren_exciter

# Map group canonical variable names to REECA1's actual variable names
function group_var_mapping(::Type{<:REECA1})
    Dict(
        :Ipcmd_out => :Ipcmd_out, # Same name
        :Iqcmd_out => :Iqcmd_out, # Same name
        :Pref => :Pref, # Same name
        :Qref => :Qref, # Same name
        :pfaref => :pfaref, # Same name
        :pfilt_y => :pfilt_y, # Same name
        :s4_y => :s4_y, # Same name
        :s5_y => :s5_y, # Same name
    )
end

# Map group canonical residual names to REECA1's actual residual names
function group_res_mapping(::Type{<:REECA1})
    Dict(
        :g_Pref => :g_Pref, # Same name
        :g_Qref => :g_Qref, # Same name
        :g_pfaref => :g_pfaref, # Same name
    )
end

# Controllers don't replace anything - they control generators
get_replacement_map(::Type{<:REECA1}) = OrderedDict{Symbol, NamedTuple}()

# No data inheritance for controllers (all parameters from DYR)
function get_data_inheritance(::Type{<:REECA1})
    OrderedDict{Symbol, OrderedDict{Symbol, Symbol}}()
end

# ============================================================================
# Model Registration
# ============================================================================

supports_format(::Type{<:REECA1}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model REECA1

# Terminal declarations - where REECA1 connects
const REECA1_TERMINALS = (
    Terminal(:bus_group, :BUS),
    Terminal(:ren_gen, [:BUS, :ID]),
)

get_terminals(::Type{<:REECA1}) = REECA1_TERMINALS

get_key_name(::Type{<:REECA1}) = [:BUS, :ID]

function model_metadata(::Type{REECA1})
    ModelMetadata(
        name = :REECA1,
        vars = [
            ModelVar(:s4_y, StateVar();
                description="Reactive current of the Q or power factor reference (s4), frozen in a dip"),
            ModelVar(:pfilt_y, StateVar();
                description="Power reference, rate limited by dPmin and dPmax (s6)"),
            ModelVar(:s5_y, StateVar();
                description="Power order Pord (s5), within [Pmin, Pmax], frozen in a dip"),
            ModelVar(:s0_y, StateVar();
                description="Voltage filter (s0): Vt_filt"),
            ModelVar(:S1_y, StateVar();
                description="Electrical power filter (s1)"),
            ModelVar(:PIV_xi, StateVar();
                description="Voltage PI controller (s3), non-windup within [-Iqmax, Iqmax], frozen in a dip"),
            ModelVar(:PIQ_xi, StateVar();
                description="Reactive power PI controller (s2), non-windup within [Vmin, Vmax], frozen in a dip"),
            ModelVar(:Pref, AlgebVar();
                description="Active power reference (receives a plant controller's contribution)"),
            ModelVar(:pfaref, AlgebVar();
                description="Power factor angle reference, at which Pe tan(pfaref) is the reactive reference the start needs"),
            ModelVar(:Qref, AlgebVar();
                description="Reactive power reference Qext (receives a plant controller's contribution)"),
            ModelVar(:Ipcmd_out, AlgebVar();
                description="Active current command, Pord / Vt_filt within [0, Ipmax]"),
            ModelVar(:Iqcmd_out, AlgebVar();
                description="Reactive current command, positive injecting, within [-Iqmax, Iqmax]"),
            ModelVar(:Qerr, AlgebVar();
                description="Reactive power error, the reference within [Qmin, Qmax] less Qgen"),
            ModelVar(:Vsel, AlgebVar();
                description="VFlag switch, within [Vmin, Vmax]: the voltage reference of s3"),
            ModelVar(:Iqinj, AlgebVar();
                description="Reactive current injection of the state machine"),
            ModelVar(:Iqmax, AlgebVar();
                description="Reactive current limit (Appendix B); Iqmin = -Iqmax"),
            ModelVar(:Ipmax, AlgebVar();
                description="Active current limit (Appendix B), held at its dip value for Thld2 s; Ipmin = 0"),
            ModelVar(:v, AlgebVar(), :bus_group, :v, :BUS;
                description="Terminal voltage magnitude (Vt) (external)", units="pu"),
            ModelVar(:Pe, AlgebVar(), :ren_gen, :Pe, [:BUS, :ID];
                description="Electrical power Pe, measured at the converter's terminal (external)", units="pu"),
            ModelVar(:Qgen, AlgebVar(), :ren_gen, :Qe, [:BUS, :ID];
                description="Reactive power generated Qgen, measured at the converter's terminal (external)", units="pu")
        ],
        residuals = [
        ModelResidual(:ds4_y, StateRes(), PrivateRes();
            description="Reactive current of the Q or power factor reference (s4), frozen in a dip"),
        ModelResidual(:dpfilt_y, StateRes(), PrivateRes();
            description="Power reference, rate limited by dPmin and dPmax (s6)"),
        ModelResidual(:ds5_y, StateRes(), PrivateRes();
            description="Power order Pord (s5), within [Pmin, Pmax], frozen in a dip"),
        ModelResidual(:ds0_y, StateRes(), PrivateRes();
            description="Voltage filter (s0): Vt_filt"),
        ModelResidual(:dS1_y, StateRes(), PrivateRes();
            description="Electrical power filter (s1)"),
        ModelResidual(:dPIV_xi, StateRes(), PrivateRes();
            description="Voltage PI controller (s3), non-windup within [-Iqmax, Iqmax], frozen in a dip"),
        ModelResidual(:dPIQ_xi, StateRes(), PrivateRes();
            description="Reactive power PI controller (s2), non-windup within [Vmin, Vmax], frozen in a dip"),
        ModelResidual(:g_Pref, AlgebRes(), SharedRes(:add);
            description="Active power reference (receives a plant controller's contribution)"),
        ModelResidual(:g_pfaref, AlgebRes(), SharedRes(:add);
            description="Power factor angle reference, at which Pe tan(pfaref) is the reactive reference the start needs"),
        ModelResidual(:g_Qref, AlgebRes(), SharedRes(:add);
            description="Reactive power reference Qext (receives a plant controller's contribution)"),
        ModelResidual(:g_Ipcmd_out, AlgebRes(), PrivateRes();
            description="Active current command, Pord / Vt_filt within [0, Ipmax]"),
        ModelResidual(:g_Iqcmd_out, AlgebRes(), PrivateRes();
            description="Reactive current command, positive injecting, within [-Iqmax, Iqmax]"),
        ModelResidual(:g_Qerr, AlgebRes(), PrivateRes();
            description="Reactive power error, the reference within [Qmin, Qmax] less Qgen"),
        ModelResidual(:g_Vsel, AlgebRes(), PrivateRes();
            description="VFlag switch, within [Vmin, Vmax]: the voltage reference of s3"),
        ModelResidual(:g_Iqinj, AlgebRes(), PrivateRes();
            description="Reactive current injection of the state machine"),
        ModelResidual(:g_Iqmax, AlgebRes(), PrivateRes();
            description="Reactive current limit (Appendix B); Iqmin = -Iqmax"),
        ModelResidual(:g_Ipmax, AlgebRes(), PrivateRes();
            description="Active current limit (Appendix B), held at its dip value for Thld2 s; Ipmin = 0"),
        ModelResidual(:Ipcmd_contrib, AlgebRes(), SharedRes(:add), :ren_gen, :g_Ipcmd, [:BUS, :ID];
            description="Active current command contribution (delta from initial)"),
        ModelResidual(:Iqcmd_contrib, AlgebRes(), SharedRes(:add), :ren_gen, :g_Iqcmd, [:BUS, :ID];
            description="Reactive current command contribution, negated into the converter's sign (delta from initial)")
    ],
        ext_params = ExtParam[],
        init_inputs = [
        InitInput(:v, :from_solution, :Bus, :v;
            description="Terminal voltage from the power flow solution"),
        InitInput(:Ipcmd_init, :from_group_init, :ren_gen, :Ipcmd;
            description="Initial active current command of the converter model"),
        InitInput(:Iqcmd_init, :from_group_init, :ren_gen, :Iqcmd;
            description="Initial reactive current command of the converter model (positive absorbs)"),
        InitInput(:Pe_init, :from_group_init, :ren_gen, :Pe;
            description="Initial electrical power at the converter's terminal"),
        InitInput(:Qgen_init, :from_group_init, :ren_gen, :Qe;
            description="Initial reactive power generated at the converter's terminal")
    ],
        observable_vars = Symbol[],
        initialized_params = [
        InitParam(:Ipcmd0, :Ipcmd_out; description="Initial active current command"),
        InitParam(:Iqcmd0, :Iqcmd_out; description="Initial reactive current command (positive injecting)"),
        InitParam(:v0, :s0_y; description="Initial filtered terminal voltage"),
        InitParam(:Pref0, :Pref; description="Active power reference (receives a plant controller's contribution)"),
        InitParam(:pfaref0, :pfaref; description="Power factor angle reference, at which Pe tan(pfaref) is the reactive reference the start needs"),
        InitParam(:Qref0, :Qref; description="Reactive power reference Qext (receives a plant controller's contribution)")
    ],
        variable_limits = [VariableLimit(:PIV_xi, (m, x) -> -x.Iqmax, (m, x) -> x.Iqmax), VariableLimit(:PIQ_xi, :VMIN, :VMAX), VariableLimit(:s5_y, :PMIN, :PMAX), VariableLimit(:Qref, :QMin, :QMax; when = (m, x) -> ((m.VFLAG .== 1) .& (x.s0_y .>= m.Vdip) .& (x.s0_y .<= m.Vup)))],
        automata = (
            Automaton(:PIV;
                modes = [Mode(:free; holds = [:PIV_gate => 1.0, :PIV_side => 0.0]),
                    Mode(:upper; holds = [:PIV_gate => 0.0, :PIV_side => 1.0]),
                    Mode(:lower; holds = [:PIV_gate => 0.0, :PIV_side => -1.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:PIV_unlimited_minus_Iqmax, :>=), Atom(:PIV_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:PIV_rate, :<), Atom(:PIV_value, :<=)])),
                    Transition(:upper, :free; when = AnyOf(Condition[Atom(:PIV_rate, :<=), Atom(:PIV_unlimited_minus_Iqmax, :<)])),
                    Transition(:lower, :free; when = AnyOf(Condition[Atom(:PIV_rate, :>=), Atom(:PIV_value, :>)]))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:PIV_unlimited_minus_Iqmax, :>=), Atom(:PIV_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:PIV_rate, :<), Atom(:PIV_value, :<=)])),
                    InitialMode(:free)],
                at_rest = [:PIV_rate]),
            Automaton(:PIQ;
                modes = [Mode(:free; holds = [:PIQ_gate => 1.0, :PIQ_side => 0.0]),
                    Mode(:upper; holds = [:PIQ_gate => 0.0, :PIQ_side => 1.0]),
                    Mode(:lower; holds = [:PIQ_gate => 0.0, :PIQ_side => -1.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:PIQ_rate, :>), AnyOf(Condition[Atom(:PIQ_unlimited_minus_VMAX, :>=), Atom(:PIV_side, :>)])])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:PIQ_rate, :<), AnyOf(Condition[Atom(:PIQ_unlimited_minus_VMIN, :<=), Atom(:PIV_side, :<)])])),
                    Transition(:upper, :free; when = AnyOf(Condition[Atom(:PIQ_rate, :<=), AllOf(Condition[Atom(:PIV_side, :<=), Atom(:PIQ_unlimited_minus_VMAX, :<)])])),
                    Transition(:lower, :free; when = AnyOf(Condition[Atom(:PIQ_rate, :>=), AllOf(Condition[Atom(:PIV_side, :>=), Atom(:PIQ_unlimited_minus_VMIN, :>)])]))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:PIQ_rate, :>), AnyOf(Condition[Atom(:PIQ_unlimited_minus_VMAX, :>=), Atom(:PIV_side, :>)])])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:PIQ_rate, :<), AnyOf(Condition[Atom(:PIQ_unlimited_minus_VMIN, :<=), Atom(:PIV_side, :<)])])),
                    InitialMode(:free)],
                at_rest = [:PIQ_rate]),
            Automaton(:s5_y;
                modes = [Mode(:free; holds = [:s5_y_gate => 1.0]),
                    Mode(:upper; holds = [:s5_y_gate => 0.0]),
                    Mode(:lower; holds = [:s5_y_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:s5_y_minus_PMAX, :>=), Atom(:s5_y_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:s5_y_minus_PMIN, :<=), Atom(:s5_y_rate, :<)])),
                    Transition(:upper, :free; when = Atom(:s5_y_rate, :<)),
                    Transition(:lower, :free; when = Atom(:s5_y_rate, :>))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:s5_y_minus_PMAX, :>=), Atom(:s5_y_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:s5_y_minus_PMIN, :<=), Atom(:s5_y_rate, :<)])),
                    InitialMode(:free)],
                at_rest = [:s5_y_rate]),
            Automaton(:injection;
                modes = [Mode(:state0; holds = [:Voltage_dip => 0.0, :in_state1 => 0.0, :in_state2 => 0.0]),
                    Mode(:state1; holds = [:Voltage_dip => 1.0, :in_state1 => 1.0, :in_state2 => 0.0]),
                    Mode(:state1_extended; holds = [:Voltage_dip => 0.0, :in_state1 => 1.0, :in_state2 => 0.0]),
                    Mode(:state2; holds = [:Voltage_dip => 0.0, :in_state1 => 0.0, :in_state2 => 1.0])],
                transitions = [Transition(:state0, :state1; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)])),
                    Transition(:state1, :state0; when = AllOf(Condition[Atom(:Thld, :>=), Atom(:v_minus_Vdip, :>=), Atom(:Thld, :<=), Atom(:v_minus_Vup, :<=)])),
                    Transition(:state1, :state2; when = AllOf(Condition[Atom(:v_minus_Vdip, :>=), Atom(:v_minus_Vup, :<=), Atom(:Thld, :>)]), sets = [:t_release => :t]),
                    Transition(:state1, :state1_extended; when = AllOf(Condition[Atom(:v_minus_Vdip, :>=), Atom(:v_minus_Vup, :<=), Atom(:Thld, :<)]), sets = [:t_release => :t]),
                    Transition(:state2, :state0; when = AllOf(Condition[Atom(:v_minus_Vdip, :>=), Atom(:v_minus_Vup, :<=), Atom(:injection_value, :>=)])),
                    Transition(:state1_extended, :state0; when = AllOf(Condition[Atom(:v_minus_Vdip, :>=), Atom(:v_minus_Vup, :<=), Atom(:injection_value_2, :>=)])),
                    Transition(:state2, :state1; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)])),
                    Transition(:state1_extended, :state1; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)]))],
                initial = [InitialMode(:state1; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)])),
                    InitialMode(:state0)]),
            Automaton(:ipmax_hold;
                modes = [Mode(:normal; holds = [:Ipmax_hold => 0.0]),
                    Mode(:dip; holds = [:Ipmax_hold => 0.0]),
                    Mode(:held; holds = [:Ipmax_hold => 1.0])],
                transitions = [Transition(:normal, :dip; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)])),
                    Transition(:dip, :held; when = AllOf(Condition[Atom(:v_minus_Vdip, :>=), Atom(:v_minus_Vup, :<=)]), sets = [:t_release2 => :t, :Ipmax_held => :Ipmax]),
                    Transition(:held, :normal; when = AllOf(Condition[Atom(:v_minus_Vdip, :>=), Atom(:v_minus_Vup, :<=), Atom(:ipmax_hold_value, :>=)])),
                    Transition(:held, :dip; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)]))],
                initial = [InitialMode(:dip; when = AnyOf(Condition[Atom(:v_minus_Vup, :>), Atom(:v_minus_Vdip, :<)])),
                    InitialMode(:normal)]),
            Automaton(:vfloor_max;
                modes = [Mode(:second; holds = [:vfloor_max_second => 1.0, :vfloor_max_first => 0.0]),
                    Mode(:first; holds = [:vfloor_max_second => 0.0, :vfloor_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:vfloor_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:vfloor_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:vfloor_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:PIV_y_clamp;
                modes = [Mode(:high; holds = [:PIV_y_clamp_high => 1.0, :PIV_y_clamp_low => 0.0, :PIV_y_clamp_mid => 0.0]),
                    Mode(:low; holds = [:PIV_y_clamp_high => 0.0, :PIV_y_clamp_low => 1.0, :PIV_y_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:PIV_y_clamp_high => 0.0, :PIV_y_clamp_low => 0.0, :PIV_y_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:PIV_y_clamp_value, :<), Atom(:PIV_y_clamp_value_2, :<=)]), continuity = [:PIV_y_clamp_value_2 => [:PIV_unlimited_minus_Iqmax], :PIV_y_clamp_value => [:PIV_value]], continuity_writers = [:PIV_unlimited_minus_Iqmax => [:injection, :ipmax_hold], :PIV_value => [:injection, :ipmax_hold]]),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:PIV_y_clamp_value, :>=), Atom(:PIV_y_clamp_value_2, :<=)]), continuity = [:PIV_y_clamp_value_2 => [:PIV_unlimited_minus_Iqmax], :PIV_y_clamp_value => [:PIV_value]], continuity_writers = [:PIV_unlimited_minus_Iqmax => [:injection, :ipmax_hold], :PIV_value => [:injection, :ipmax_hold]]),
                    Transition(:low, :high; when = Atom(:PIV_y_clamp_value_2, :>), continuity = [:PIV_y_clamp_value_2 => [:PIV_unlimited_minus_Iqmax]], continuity_writers = [:PIV_unlimited_minus_Iqmax => [:injection, :ipmax_hold]]),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:PIV_y_clamp_value, :>=), Atom(:PIV_y_clamp_value_2, :<=)]), continuity = [:PIV_y_clamp_value_2 => [:PIV_unlimited_minus_Iqmax], :PIV_y_clamp_value => [:PIV_value]], continuity_writers = [:PIV_unlimited_minus_Iqmax => [:injection, :ipmax_hold], :PIV_value => [:injection, :ipmax_hold]]),
                    Transition(:mid, :high; when = Atom(:PIV_y_clamp_value_2, :>), continuity = [:PIV_y_clamp_value_2 => [:PIV_unlimited_minus_Iqmax]], continuity_writers = [:PIV_unlimited_minus_Iqmax => [:injection, :ipmax_hold]]),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:PIV_y_clamp_value, :<), Atom(:PIV_y_clamp_value_2, :<=)]), continuity = [:PIV_y_clamp_value_2 => [:PIV_unlimited_minus_Iqmax], :PIV_y_clamp_value => [:PIV_value]], continuity_writers = [:PIV_unlimited_minus_Iqmax => [:injection, :ipmax_hold], :PIV_value => [:injection, :ipmax_hold]])],
                initial = [InitialMode(:high; when = Atom(:PIV_y_clamp_value_2, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:PIV_y_clamp_value, :<), Atom(:PIV_y_clamp_value_2, :<=)])),
                    InitialMode(:mid)]),
            Automaton(:PIQ_y_clamp;
                modes = [Mode(:high; holds = [:PIQ_y_clamp_high => 1.0, :PIQ_y_clamp_low => 0.0, :PIQ_y_clamp_mid => 0.0]),
                    Mode(:low; holds = [:PIQ_y_clamp_high => 0.0, :PIQ_y_clamp_low => 1.0, :PIQ_y_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:PIQ_y_clamp_high => 0.0, :PIQ_y_clamp_low => 0.0, :PIQ_y_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:PIQ_y_clamp_value, :<=), Atom(:PIQ_y_clamp_value_2, :<)]), continuity = [:PIQ_y_clamp_value => [:PIQ_unlimited_minus_VMAX], :PIQ_y_clamp_value_2 => [:PIQ_unlimited_minus_VMIN]]),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:PIQ_y_clamp_value_2, :>=), Atom(:PIQ_y_clamp_value, :<=)]), continuity = [:PIQ_y_clamp_value => [:PIQ_unlimited_minus_VMAX], :PIQ_y_clamp_value_2 => [:PIQ_unlimited_minus_VMIN]]),
                    Transition(:low, :high; when = Atom(:PIQ_y_clamp_value, :>), continuity = [:PIQ_y_clamp_value => [:PIQ_unlimited_minus_VMAX]]),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:PIQ_y_clamp_value_2, :>=), Atom(:PIQ_y_clamp_value, :<=)]), continuity = [:PIQ_y_clamp_value => [:PIQ_unlimited_minus_VMAX], :PIQ_y_clamp_value_2 => [:PIQ_unlimited_minus_VMIN]]),
                    Transition(:mid, :high; when = Atom(:PIQ_y_clamp_value, :>), continuity = [:PIQ_y_clamp_value => [:PIQ_unlimited_minus_VMAX]]),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:PIQ_y_clamp_value, :<=), Atom(:PIQ_y_clamp_value_2, :<)]), continuity = [:PIQ_y_clamp_value => [:PIQ_unlimited_minus_VMAX], :PIQ_y_clamp_value_2 => [:PIQ_unlimited_minus_VMIN]])],
                initial = [InitialMode(:high; when = Atom(:PIQ_y_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:PIQ_y_clamp_value, :<=), Atom(:PIQ_y_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:Iqv_clamp;
                modes = [Mode(:high; holds = [:Iqv_clamp_high => 1.0, :Iqv_clamp_low => 0.0, :Iqv_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Iqv_clamp_high => 0.0, :Iqv_clamp_low => 1.0, :Iqv_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Iqv_clamp_high => 0.0, :Iqv_clamp_low => 0.0, :Iqv_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Iqv_clamp_value, :<=), Atom(:Iqv_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Iqv_clamp_value_2, :>=), Atom(:Iqv_clamp_value, :<=)])),
                    Transition(:low, :high; when = Atom(:Iqv_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Iqv_clamp_value_2, :>=), Atom(:Iqv_clamp_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:Iqv_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Iqv_clamp_value, :<=), Atom(:Iqv_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:Iqv_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Iqv_clamp_value, :<=), Atom(:Iqv_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:VDL1_vdl1;
                modes = [Mode(:below; holds = [:VDL1_vdl1_below => 1.0, :VDL1_vdl1_segment1 => 0.0, :VDL1_vdl1_segment2 => 0.0, :VDL1_vdl1_segment3 => 0.0, :VDL1_vdl1_above => 0.0]),
                    Mode(:segment1; holds = [:VDL1_vdl1_below => 0.0, :VDL1_vdl1_segment1 => 1.0, :VDL1_vdl1_segment2 => 0.0, :VDL1_vdl1_segment3 => 0.0, :VDL1_vdl1_above => 0.0]),
                    Mode(:segment2; holds = [:VDL1_vdl1_below => 0.0, :VDL1_vdl1_segment1 => 0.0, :VDL1_vdl1_segment2 => 1.0, :VDL1_vdl1_segment3 => 0.0, :VDL1_vdl1_above => 0.0]),
                    Mode(:segment3; holds = [:VDL1_vdl1_below => 0.0, :VDL1_vdl1_segment1 => 0.0, :VDL1_vdl1_segment2 => 0.0, :VDL1_vdl1_segment3 => 1.0, :VDL1_vdl1_above => 0.0]),
                    Mode(:above; holds = [:VDL1_vdl1_below => 0.0, :VDL1_vdl1_segment1 => 0.0, :VDL1_vdl1_segment2 => 0.0, :VDL1_vdl1_segment3 => 0.0, :VDL1_vdl1_above => 1.0])],
                transitions = [Transition(:below, :segment1; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :<)])),
                    Transition(:below, :segment2; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :<)])),
                    Transition(:below, :segment3; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :<=)])),
                    Transition(:below, :above; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :>)])),
                    Transition(:segment1, :below; when = Atom(:VDL1_vdl1_value, :<)),
                    Transition(:segment1, :segment2; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :<)])),
                    Transition(:segment1, :segment3; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :<=)])),
                    Transition(:segment1, :above; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :>)])),
                    Transition(:segment2, :below; when = Atom(:VDL1_vdl1_value, :<)),
                    Transition(:segment2, :segment1; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :<)])),
                    Transition(:segment2, :segment3; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :<=)])),
                    Transition(:segment2, :above; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :>)])),
                    Transition(:segment3, :below; when = Atom(:VDL1_vdl1_value, :<)),
                    Transition(:segment3, :segment1; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :<)])),
                    Transition(:segment3, :segment2; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :<)])),
                    Transition(:segment3, :above; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :>)])),
                    Transition(:above, :below; when = Atom(:VDL1_vdl1_value, :<)),
                    Transition(:above, :segment1; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :<)])),
                    Transition(:above, :segment2; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :<)])),
                    Transition(:above, :segment3; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :<=)]))],
                initial = [InitialMode(:below; when = Atom(:VDL1_vdl1_value, :<)),
                    InitialMode(:segment1; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :<)])),
                    InitialMode(:segment2; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :<)])),
                    InitialMode(:segment3; when = AllOf(Condition[Atom(:VDL1_vdl1_value, :>=), Atom(:VDL1_vdl1_value_2, :>=), Atom(:VDL1_vdl1_value_3, :>=), Atom(:VDL1_vdl1_value_4, :<=)])),
                    InitialMode(:above)]),
            Automaton(:VDL2_vdl2;
                modes = [Mode(:below; holds = [:VDL2_vdl2_below => 1.0, :VDL2_vdl2_segment1 => 0.0, :VDL2_vdl2_segment2 => 0.0, :VDL2_vdl2_segment3 => 0.0, :VDL2_vdl2_above => 0.0]),
                    Mode(:segment1; holds = [:VDL2_vdl2_below => 0.0, :VDL2_vdl2_segment1 => 1.0, :VDL2_vdl2_segment2 => 0.0, :VDL2_vdl2_segment3 => 0.0, :VDL2_vdl2_above => 0.0]),
                    Mode(:segment2; holds = [:VDL2_vdl2_below => 0.0, :VDL2_vdl2_segment1 => 0.0, :VDL2_vdl2_segment2 => 1.0, :VDL2_vdl2_segment3 => 0.0, :VDL2_vdl2_above => 0.0]),
                    Mode(:segment3; holds = [:VDL2_vdl2_below => 0.0, :VDL2_vdl2_segment1 => 0.0, :VDL2_vdl2_segment2 => 0.0, :VDL2_vdl2_segment3 => 1.0, :VDL2_vdl2_above => 0.0]),
                    Mode(:above; holds = [:VDL2_vdl2_below => 0.0, :VDL2_vdl2_segment1 => 0.0, :VDL2_vdl2_segment2 => 0.0, :VDL2_vdl2_segment3 => 0.0, :VDL2_vdl2_above => 1.0])],
                transitions = [Transition(:below, :segment1; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :<)])),
                    Transition(:below, :segment2; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :<)])),
                    Transition(:below, :segment3; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :<=)])),
                    Transition(:below, :above; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :>)])),
                    Transition(:segment1, :below; when = Atom(:VDL2_vdl2_value, :<)),
                    Transition(:segment1, :segment2; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :<)])),
                    Transition(:segment1, :segment3; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :<=)])),
                    Transition(:segment1, :above; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :>)])),
                    Transition(:segment2, :below; when = Atom(:VDL2_vdl2_value, :<)),
                    Transition(:segment2, :segment1; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :<)])),
                    Transition(:segment2, :segment3; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :<=)])),
                    Transition(:segment2, :above; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :>)])),
                    Transition(:segment3, :below; when = Atom(:VDL2_vdl2_value, :<)),
                    Transition(:segment3, :segment1; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :<)])),
                    Transition(:segment3, :segment2; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :<)])),
                    Transition(:segment3, :above; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :>)])),
                    Transition(:above, :below; when = Atom(:VDL2_vdl2_value, :<)),
                    Transition(:above, :segment1; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :<)])),
                    Transition(:above, :segment2; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :<)])),
                    Transition(:above, :segment3; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :<=)]))],
                initial = [InitialMode(:below; when = Atom(:VDL2_vdl2_value, :<)),
                    InitialMode(:segment1; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :<)])),
                    InitialMode(:segment2; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :<)])),
                    InitialMode(:segment3; when = AllOf(Condition[Atom(:VDL2_vdl2_value, :>=), Atom(:VDL2_vdl2_value_2, :>=), Atom(:VDL2_vdl2_value_3, :>=), Atom(:VDL2_vdl2_value_4, :<=)])),
                    InitialMode(:above)]),
            Automaton(:Ipcmd_out_clamp;
                modes = [Mode(:high; holds = [:Ipcmd_out_clamp_high => 1.0, :Ipcmd_out_clamp_low => 0.0, :Ipcmd_out_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Ipcmd_out_clamp_high => 0.0, :Ipcmd_out_clamp_low => 1.0, :Ipcmd_out_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Ipcmd_out_clamp_high => 0.0, :Ipcmd_out_clamp_low => 0.0, :Ipcmd_out_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Ipcmd_out_clamp_value, :<), Atom(:Ipcmd_out_clamp_value_2, :<=)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Ipcmd_out_clamp_value, :>=), Atom(:Ipcmd_out_clamp_value_2, :<=)])),
                    Transition(:low, :high; when = Atom(:Ipcmd_out_clamp_value_2, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Ipcmd_out_clamp_value, :>=), Atom(:Ipcmd_out_clamp_value_2, :<=)])),
                    Transition(:mid, :high; when = Atom(:Ipcmd_out_clamp_value_2, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Ipcmd_out_clamp_value, :<), Atom(:Ipcmd_out_clamp_value_2, :<=)]))],
                initial = [InitialMode(:high; when = Atom(:Ipcmd_out_clamp_value_2, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Ipcmd_out_clamp_value, :<), Atom(:Ipcmd_out_clamp_value_2, :<=)])),
                    InitialMode(:mid)]),
            Automaton(:Iqcmd_out_clamp;
                modes = [Mode(:high; holds = [:Iqcmd_out_clamp_high => 1.0, :Iqcmd_out_clamp_low => 0.0, :Iqcmd_out_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Iqcmd_out_clamp_high => 0.0, :Iqcmd_out_clamp_low => 1.0, :Iqcmd_out_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Iqcmd_out_clamp_high => 0.0, :Iqcmd_out_clamp_low => 0.0, :Iqcmd_out_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Iqcmd_out_clamp_value, :<), Atom(:Iqcmd_out_clamp_value_2, :<=)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Iqcmd_out_clamp_value, :>=), Atom(:Iqcmd_out_clamp_value_2, :<=)])),
                    Transition(:low, :high; when = Atom(:Iqcmd_out_clamp_value_2, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Iqcmd_out_clamp_value, :>=), Atom(:Iqcmd_out_clamp_value_2, :<=)])),
                    Transition(:mid, :high; when = Atom(:Iqcmd_out_clamp_value_2, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Iqcmd_out_clamp_value, :<), Atom(:Iqcmd_out_clamp_value_2, :<=)]))],
                initial = [InitialMode(:high; when = Atom(:Iqcmd_out_clamp_value_2, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Iqcmd_out_clamp_value, :<), Atom(:Iqcmd_out_clamp_value_2, :<=)])),
                    InitialMode(:mid)]),
            Automaton(:Qerr_clamp;
                modes = [Mode(:high; holds = [:Qerr_clamp_high => 1.0, :Qerr_clamp_low => 0.0, :Qerr_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Qerr_clamp_high => 0.0, :Qerr_clamp_low => 1.0, :Qerr_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Qerr_clamp_high => 0.0, :Qerr_clamp_low => 0.0, :Qerr_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Qerr_clamp_value, :<=), Atom(:Qerr_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Qerr_clamp_value_2, :>=), Atom(:Qerr_clamp_value, :<=)])),
                    Transition(:low, :high; when = Atom(:Qerr_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Qerr_clamp_value_2, :>=), Atom(:Qerr_clamp_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:Qerr_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Qerr_clamp_value, :<=), Atom(:Qerr_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:Qerr_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Qerr_clamp_value, :<=), Atom(:Qerr_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:Vsel_clamp;
                modes = [Mode(:high; holds = [:Vsel_clamp_high => 1.0, :Vsel_clamp_low => 0.0, :Vsel_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Vsel_clamp_high => 0.0, :Vsel_clamp_low => 1.0, :Vsel_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Vsel_clamp_high => 0.0, :Vsel_clamp_low => 0.0, :Vsel_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Vsel_clamp_value, :<=), Atom(:Vsel_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Vsel_clamp_value_2, :>=), Atom(:Vsel_clamp_value, :<=)])),
                    Transition(:low, :high; when = Atom(:Vsel_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Vsel_clamp_value_2, :>=), Atom(:Vsel_clamp_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:Vsel_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Vsel_clamp_value, :<=), Atom(:Vsel_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:Vsel_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Vsel_clamp_value, :<=), Atom(:Vsel_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:Iqinj_clamp;
                modes = [Mode(:high; holds = [:Iqinj_clamp_high => 1.0, :Iqinj_clamp_low => 0.0, :Iqinj_clamp_mid => 0.0]),
                    Mode(:low; holds = [:Iqinj_clamp_high => 0.0, :Iqinj_clamp_low => 1.0, :Iqinj_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:Iqinj_clamp_high => 0.0, :Iqinj_clamp_low => 0.0, :Iqinj_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:Iqinj_clamp_value, :<=), Atom(:Iqinj_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:Iqinj_clamp_value_2, :>=), Atom(:Iqinj_clamp_value, :<=)])),
                    Transition(:low, :high; when = Atom(:Iqinj_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:Iqinj_clamp_value_2, :>=), Atom(:Iqinj_clamp_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:Iqinj_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:Iqinj_clamp_value, :<=), Atom(:Iqinj_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:Iqinj_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:Iqinj_clamp_value, :<=), Atom(:Iqinj_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
            Automaton(:Iqmax_max;
                modes = [Mode(:second; holds = [:Iqmax_max_second => 1.0, :Iqmax_max_first => 0.0]),
                    Mode(:first; holds = [:Iqmax_max_second => 0.0, :Iqmax_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Iqmax_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:Iqmax_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:Iqmax_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:Iqmax_min;
                modes = [Mode(:second; holds = [:Iqmax_min_second => 1.0, :Iqmax_min_first => 0.0]),
                    Mode(:first; holds = [:Iqmax_min_second => 0.0, :Iqmax_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Iqmax_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:Iqmax_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:Iqmax_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:Iqmax_min_2;
                modes = [Mode(:second; holds = [:Iqmax_min_2_second => 1.0, :Iqmax_min_2_first => 0.0]),
                    Mode(:first; holds = [:Iqmax_min_2_second => 0.0, :Iqmax_min_2_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Iqmax_min_2_value, :>=)),
                    Transition(:first, :second; when = Atom(:Iqmax_min_2_value, :<))],
                initial = [InitialMode(:second; when = Atom(:Iqmax_min_2_value, :<)),
                    InitialMode(:first)]),
            Automaton(:Ipmax_min;
                modes = [Mode(:second; holds = [:Ipmax_min_second => 1.0, :Ipmax_min_first => 0.0]),
                    Mode(:first; holds = [:Ipmax_min_second => 0.0, :Ipmax_min_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Ipmax_min_value, :>=)),
                    Transition(:first, :second; when = Atom(:Ipmax_min_value, :<))],
                initial = [InitialMode(:second; when = Atom(:Ipmax_min_value, :<)),
                    InitialMode(:first)]),
            Automaton(:Ipmax_max;
                modes = [Mode(:second; holds = [:Ipmax_max_second => 1.0, :Ipmax_max_first => 0.0]),
                    Mode(:first; holds = [:Ipmax_max_second => 0.0, :Ipmax_max_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Ipmax_max_value, :<=)),
                    Transition(:first, :second; when = Atom(:Ipmax_max_value, :>))],
                initial = [InitialMode(:second; when = Atom(:Ipmax_max_value, :>)),
                    InitialMode(:first)]),
            Automaton(:Ipmax_min_2;
                modes = [Mode(:second; holds = [:Ipmax_min_2_second => 1.0, :Ipmax_min_2_first => 0.0]),
                    Mode(:first; holds = [:Ipmax_min_2_second => 0.0, :Ipmax_min_2_first => 1.0])],
                transitions = [Transition(:second, :first; when = Atom(:Ipmax_min_2_value, :>=)),
                    Transition(:first, :second; when = Atom(:Ipmax_min_2_value, :<))],
                initial = [InitialMode(:second; when = Atom(:Ipmax_min_2_value, :<)),
                    InitialMode(:first)]),
            Automaton(:pfilt_y_rate_clamp;
                modes = [Mode(:high; holds = [:pfilt_y_rate_clamp_high => 1.0, :pfilt_y_rate_clamp_low => 0.0, :pfilt_y_rate_clamp_mid => 0.0]),
                    Mode(:low; holds = [:pfilt_y_rate_clamp_high => 0.0, :pfilt_y_rate_clamp_low => 1.0, :pfilt_y_rate_clamp_mid => 0.0]),
                    Mode(:mid; holds = [:pfilt_y_rate_clamp_high => 0.0, :pfilt_y_rate_clamp_low => 0.0, :pfilt_y_rate_clamp_mid => 1.0])],
                transitions = [Transition(:high, :low; when = AllOf(Condition[Atom(:pfilt_y_rate_clamp_value, :<=), Atom(:pfilt_y_rate_clamp_value_2, :<)])),
                    Transition(:high, :mid; when = AllOf(Condition[Atom(:pfilt_y_rate_clamp_value_2, :>=), Atom(:pfilt_y_rate_clamp_value, :<=)])),
                    Transition(:low, :high; when = Atom(:pfilt_y_rate_clamp_value, :>)),
                    Transition(:low, :mid; when = AllOf(Condition[Atom(:pfilt_y_rate_clamp_value_2, :>=), Atom(:pfilt_y_rate_clamp_value, :<=)])),
                    Transition(:mid, :high; when = Atom(:pfilt_y_rate_clamp_value, :>)),
                    Transition(:mid, :low; when = AllOf(Condition[Atom(:pfilt_y_rate_clamp_value, :<=), Atom(:pfilt_y_rate_clamp_value_2, :<)]))],
                initial = [InitialMode(:high; when = Atom(:pfilt_y_rate_clamp_value, :>)),
                    InitialMode(:low; when = AllOf(Condition[Atom(:pfilt_y_rate_clamp_value, :<=), Atom(:pfilt_y_rate_clamp_value_2, :<)])),
                    InitialMode(:mid)]),
        ),
        memory = (
            MemoryQuantity(:t_release; initial = 0.0),
            MemoryQuantity(:t_release2; initial = 0.0),
            MemoryQuantity(:Ipmax_held; initial = 0.0),
        ),
        status_field = :status,
        evaluation_pass = 2,
    )
end

@inline function reeca1_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t, ::Type{AT}
) where AT
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        tpfilt = 0.02
        wg = 1.0

        # --- s0: Lag ---
        s0_dy = ((1.0) * (vars.v[i]) - vars.s0_y[i]) / (params.Trv[i])

        # --- S1: Lag ---
        S1_dy = ((1.0) * (vars.Pe[i]) - vars.S1_y[i]) / (params.Tp[i])

        vfloor = 0.01 * params.vfloor_max_first[i] + vars.s0_y[i] * params.vfloor_max_second[i]
        Qcpf = vars.S1_y[i] * tan(vars.pfaref[i])
        PFsel = params.PFFLAG[i] * Qcpf + vars.Qref[i] * (1 - params.PFFLAG[i])
        moving = 1 - params.Voltage_dip[i]
        PIV_rate = params.Kvi[i] * params.QFLAG[i] * moving * (vars.Vsel[i] - vars.s0_y[i])
        PIV_unlimited = vars.PIV_xi[i] + params.Kvp[i] * params.QFLAG[i] * (vars.Vsel[i] - vars.s0_y[i])
        PIV_y = vars.Iqmax[i] * params.PIV_y_clamp_high[i] + params.PIV_y_clamp_mid[i] * (vars.PIV_xi[i] + params.Kvp[i] * params.QFLAG[i] * (vars.Vsel[i] - vars.s0_y[i])) - vars.Iqmax[i] * params.PIV_y_clamp_low[i]
        PIQ_rate = params.Kqi[i] * vars.Qerr[i] * params.VFLAG[i] * moving
        PIQ_unlimited = vars.PIQ_xi[i] + params.Kqp[i] * vars.Qerr[i] * params.VFLAG[i]
        PIQ_y = params.PIQ_y_clamp_high[i] * params.VMAX[i] + params.PIQ_y_clamp_low[i] * params.VMIN[i] + params.PIQ_y_clamp_mid[i] * (vars.PIQ_xi[i] + params.Kqp[i] * vars.Qerr[i] * params.VFLAG[i])
        Vref0_eff = ifelse(params.Vref0[i] == 0, params.v0[i], params.Vref0[i])
        Verr = Vref0_eff - vars.s0_y[i]
        Iqv = params.Kqv[i] * (Verr - params.Iqv_clamp_high[i] * params.dbd2[i] - params.Iqv_clamp_low[i] * params.dbd1[i] - params.Iqv_clamp_mid[i] * Verr)
        VDL1 = params.Iq1[i] * params.VDL1_vdl1_below[i] + params.Iq4[i] * params.VDL1_vdl1_above[i] + params.VDL1_vdl1_segment1[i] * (params.Iq1[i] + (params.Iq2[i] - params.Iq1[i]) * (vars.s0_y[i] - params.Vq1[i]) / (params.Vq2[i] - params.Vq1[i])) + params.VDL1_vdl1_segment2[i] * (params.Iq2[i] + (params.Iq3[i] - params.Iq2[i]) * (vars.s0_y[i] - params.Vq2[i]) / (params.Vq3[i] - params.Vq2[i])) + params.VDL1_vdl1_segment3[i] * (params.Iq3[i] + (params.Iq4[i] - params.Iq3[i]) * (vars.s0_y[i] - params.Vq3[i]) / (params.Vq4[i] - params.Vq3[i]))
        VDL2 = params.Ip1[i] * params.VDL2_vdl2_below[i] + params.Ip4[i] * params.VDL2_vdl2_above[i] + params.VDL2_vdl2_segment1[i] * (params.Ip1[i] + (params.Ip2[i] - params.Ip1[i]) * (vars.s0_y[i] - params.Vp1[i]) / (params.Vp2[i] - params.Vp1[i])) + params.VDL2_vdl2_segment2[i] * (params.Ip2[i] + (params.Ip3[i] - params.Ip2[i]) * (vars.s0_y[i] - params.Vp2[i]) / (params.Vp3[i] - params.Vp2[i])) + params.VDL2_vdl2_segment3[i] * (params.Ip3[i] + (params.Ip4[i] - params.Ip3[i]) * (vars.s0_y[i] - params.Vp3[i]) / (params.Vp4[i] - params.Vp3[i]))
        Psel = vars.pfilt_y[i] * (1 - params.PFLAG[i]) + params.PFLAG[i] * vars.pfilt_y[i] * wg
        s5_y_rate = moving * (Psel - vars.s5_y[i]) / params.Tpord[i]
        # === Residual Assignments ===
        outputs.ds4_y[i] = mask * (moving * (PFsel / vfloor - vars.s4_y[i]) / params.Tiq[i])
        outputs.dpfilt_y[i] = mask * (params.dPmax[i] * params.pfilt_y_rate_clamp_high[i] + params.dPmin[i] * params.pfilt_y_rate_clamp_low[i] + params.pfilt_y_rate_clamp_mid[i] * (vars.Pref[i] - vars.pfilt_y[i]) / tpfilt)
        outputs.ds5_y[i] = mask * (params.s5_y_gate[i] * s5_y_rate)
        outputs.ds0_y[i] = mask * (s0_dy)
        outputs.dS1_y[i] = mask * (S1_dy)
        outputs.dPIV_xi[i] = mask * (params.PIV_gate[i] * PIV_rate)
        outputs.dPIQ_xi[i] = mask * (params.PIQ_gate[i] * PIQ_rate)
        outputs.g_Pref[i] = mask * (params.Pref0[i] - vars.Pref[i])
        outputs.g_pfaref[i] = mask * (params.pfaref0[i] - vars.pfaref[i])
        outputs.g_Qref[i] = mask * (params.Qref0[i] - vars.Qref[i])
        outputs.g_Ipcmd_out[i] = mask * (params.Ipcmd_out_clamp_high[i] * vars.Ipmax[i] + params.Ipcmd_out_clamp_mid[i] * vars.s5_y[i] / vfloor - vars.Ipcmd_out[i])
        outputs.g_Iqcmd_out[i] = mask * (params.Iqcmd_out_clamp_high[i] * vars.Iqmax[i] + params.Iqcmd_out_clamp_mid[i] * (vars.Iqinj[i] + PIV_y * params.QFLAG[i] + vars.s4_y[i] * (1 - params.QFLAG[i])) - vars.Iqcmd_out[i] - params.Iqcmd_out_clamp_low[i] * vars.Iqmax[i])
        outputs.g_Qerr[i] = mask * (PFsel * params.Qerr_clamp_mid[i] + params.QMax[i] * params.Qerr_clamp_high[i] + params.QMin[i] * params.Qerr_clamp_low[i] - vars.Qerr[i] - vars.Qgen[i])
        outputs.g_Vsel[i] = mask * (params.VMAX[i] * params.Vsel_clamp_high[i] + params.VMIN[i] * params.Vsel_clamp_low[i] + params.Vsel_clamp_mid[i] * (PIQ_y * params.VFLAG[i] + (1 - params.VFLAG[i]) * (PFsel + params.Vref1[i])) - vars.Vsel[i])
        outputs.g_Iqinj[i] = mask * (params.Iqfrz[i] * params.in_state2[i] + params.in_state1[i] * (params.Iqh1[i] * params.Iqinj_clamp_high[i] + params.Iqinj_clamp_low[i] * params.Iql1[i] + params.Iqinj_clamp_mid[i] * Iqv) - vars.Iqinj[i])
        outputs.g_Iqmax[i] = mask * (params.PQFLAG[i] * (params.Iqmax_min_first[i] * VDL1 + params.Iqmax_min_second[i] * sqrt(params.Iqmax_max_second[i]) * sqrt(max(0, params.Imax[i] ^ 2 - vars.Ipcmd_out[i] ^ 2))) + (1 - params.PQFLAG[i]) * (params.Imax[i] * params.Iqmax_min_2_first[i] + params.Iqmax_min_2_second[i] * VDL1) - vars.Iqmax[i])
        outputs.g_Ipmax[i] = mask * (params.Ipmax_held[i] * params.Ipmax_hold[i] + (1 - params.Ipmax_hold[i]) * (params.PQFLAG[i] * (params.Imax[i] * params.Ipmax_min_first[i] + params.Ipmax_min_second[i] * VDL2) + (1 - params.PQFLAG[i]) * (params.Ipmax_min_2_first[i] * VDL2 + params.Ipmax_min_2_second[i] * sqrt(params.Ipmax_max_second[i]) * sqrt(max(0, params.Imax[i] ^ 2 - vars.Iqcmd_out[i] ^ 2)))) - vars.Ipmax[i])
        outputs.Ipcmd_contrib[i] = mask * (vars.Ipcmd_out[i] - params.Ipcmd0[i])
        outputs.Iqcmd_contrib[i] = mask * (params.Iqcmd0[i] - vars.Iqcmd_out[i])
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REECA1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    reeca1_residual_kernel!(
        (ds4_y = ws.res_vals.ds4_y, dpfilt_y = ws.res_vals.dpfilt_y, ds5_y = ws.res_vals.ds5_y, ds0_y = ws.res_vals.ds0_y, dS1_y = ws.res_vals.dS1_y, dPIV_xi = ws.res_vals.dPIV_xi, dPIQ_xi = ws.res_vals.dPIQ_xi, g_Pref = ws.res_vals.g_Pref, g_pfaref = ws.res_vals.g_pfaref, g_Qref = ws.res_vals.g_Qref, g_Ipcmd_out = ws.res_vals.g_Ipcmd_out, g_Iqcmd_out = ws.res_vals.g_Iqcmd_out, g_Qerr = ws.res_vals.g_Qerr, g_Vsel = ws.res_vals.g_Vsel, g_Iqinj = ws.res_vals.g_Iqinj, g_Iqmax = ws.res_vals.g_Iqmax, g_Ipmax = ws.res_vals.g_Ipmax, Ipcmd_contrib = ws.res_vals.Ipcmd_contrib, Iqcmd_contrib = ws.res_vals.Iqcmd_contrib),
        (PFFLAG = model.PFFLAG, VFLAG = model.VFLAG, QFLAG = model.QFLAG, PFLAG = model.PFLAG, PQFLAG = model.PQFLAG, Vdip = model.Vdip, Vup = model.Vup, Trv = model.Trv, dbd1 = model.dbd1, dbd2 = model.dbd2, Kqv = model.Kqv, Iqh1 = model.Iqh1, Iql1 = model.Iql1, Vref0 = model.Vref0, Iqfrz = model.Iqfrz, Thld = model.Thld, Thld2 = model.Thld2, Tp = model.Tp, QMax = model.QMax, QMin = model.QMin, VMAX = model.VMAX, VMIN = model.VMIN, Kqp = model.Kqp, Kqi = model.Kqi, Kvp = model.Kvp, Kvi = model.Kvi, Vref1 = model.Vref1, Tiq = model.Tiq, dPmax = model.dPmax, dPmin = model.dPmin, PMAX = model.PMAX, PMIN = model.PMIN, Imax = model.Imax, Tpord = model.Tpord, Vq1 = model.Vq1, Iq1 = model.Iq1, Vq2 = model.Vq2, Iq2 = model.Iq2, Vq3 = model.Vq3, Iq3 = model.Iq3, Vq4 = model.Vq4, Iq4 = model.Iq4, Vp1 = model.Vp1, Ip1 = model.Ip1, Vp2 = model.Vp2, Ip2 = model.Ip2, Vp3 = model.Vp3, Ip3 = model.Ip3, Vp4 = model.Vp4, Ip4 = model.Ip4, Ipcmd0 = ws.setpoints.Ipcmd0, Iqcmd0 = ws.setpoints.Iqcmd0, v0 = ws.setpoints.v0, Pref0 = ws.setpoints.Pref0, pfaref0 = ws.setpoints.pfaref0, Qref0 = ws.setpoints.Qref0, Voltage_dip = ws.automata.discrete.Voltage_dip, in_state1 = ws.automata.discrete.in_state1, in_state2 = ws.automata.discrete.in_state2, Ipmax_hold = ws.automata.discrete.Ipmax_hold, PIV_gate = ws.automata.discrete.PIV_gate, PIV_side = ws.automata.discrete.PIV_side, PIQ_gate = ws.automata.discrete.PIQ_gate, PIQ_side = ws.automata.discrete.PIQ_side, s5_y_gate = ws.automata.discrete.s5_y_gate, vfloor_max_second = ws.automata.discrete.vfloor_max_second, vfloor_max_first = ws.automata.discrete.vfloor_max_first, PIV_y_clamp_high = ws.automata.discrete.PIV_y_clamp_high, PIV_y_clamp_low = ws.automata.discrete.PIV_y_clamp_low, PIV_y_clamp_mid = ws.automata.discrete.PIV_y_clamp_mid, PIQ_y_clamp_high = ws.automata.discrete.PIQ_y_clamp_high, PIQ_y_clamp_low = ws.automata.discrete.PIQ_y_clamp_low, PIQ_y_clamp_mid = ws.automata.discrete.PIQ_y_clamp_mid, Iqv_clamp_high = ws.automata.discrete.Iqv_clamp_high, Iqv_clamp_low = ws.automata.discrete.Iqv_clamp_low, Iqv_clamp_mid = ws.automata.discrete.Iqv_clamp_mid, VDL1_vdl1_below = ws.automata.discrete.VDL1_vdl1_below, VDL1_vdl1_segment1 = ws.automata.discrete.VDL1_vdl1_segment1, VDL1_vdl1_segment2 = ws.automata.discrete.VDL1_vdl1_segment2, VDL1_vdl1_segment3 = ws.automata.discrete.VDL1_vdl1_segment3, VDL1_vdl1_above = ws.automata.discrete.VDL1_vdl1_above, VDL2_vdl2_below = ws.automata.discrete.VDL2_vdl2_below, VDL2_vdl2_segment1 = ws.automata.discrete.VDL2_vdl2_segment1, VDL2_vdl2_segment2 = ws.automata.discrete.VDL2_vdl2_segment2, VDL2_vdl2_segment3 = ws.automata.discrete.VDL2_vdl2_segment3, VDL2_vdl2_above = ws.automata.discrete.VDL2_vdl2_above, Ipcmd_out_clamp_high = ws.automata.discrete.Ipcmd_out_clamp_high, Ipcmd_out_clamp_low = ws.automata.discrete.Ipcmd_out_clamp_low, Ipcmd_out_clamp_mid = ws.automata.discrete.Ipcmd_out_clamp_mid, Iqcmd_out_clamp_high = ws.automata.discrete.Iqcmd_out_clamp_high, Iqcmd_out_clamp_low = ws.automata.discrete.Iqcmd_out_clamp_low, Iqcmd_out_clamp_mid = ws.automata.discrete.Iqcmd_out_clamp_mid, Qerr_clamp_high = ws.automata.discrete.Qerr_clamp_high, Qerr_clamp_low = ws.automata.discrete.Qerr_clamp_low, Qerr_clamp_mid = ws.automata.discrete.Qerr_clamp_mid, Vsel_clamp_high = ws.automata.discrete.Vsel_clamp_high, Vsel_clamp_low = ws.automata.discrete.Vsel_clamp_low, Vsel_clamp_mid = ws.automata.discrete.Vsel_clamp_mid, Iqinj_clamp_high = ws.automata.discrete.Iqinj_clamp_high, Iqinj_clamp_low = ws.automata.discrete.Iqinj_clamp_low, Iqinj_clamp_mid = ws.automata.discrete.Iqinj_clamp_mid, Iqmax_max_second = ws.automata.discrete.Iqmax_max_second, Iqmax_max_first = ws.automata.discrete.Iqmax_max_first, Iqmax_min_second = ws.automata.discrete.Iqmax_min_second, Iqmax_min_first = ws.automata.discrete.Iqmax_min_first, Iqmax_min_2_second = ws.automata.discrete.Iqmax_min_2_second, Iqmax_min_2_first = ws.automata.discrete.Iqmax_min_2_first, Ipmax_min_second = ws.automata.discrete.Ipmax_min_second, Ipmax_min_first = ws.automata.discrete.Ipmax_min_first, Ipmax_max_second = ws.automata.discrete.Ipmax_max_second, Ipmax_max_first = ws.automata.discrete.Ipmax_max_first, Ipmax_min_2_second = ws.automata.discrete.Ipmax_min_2_second, Ipmax_min_2_first = ws.automata.discrete.Ipmax_min_2_first, pfilt_y_rate_clamp_high = ws.automata.discrete.pfilt_y_rate_clamp_high, pfilt_y_rate_clamp_low = ws.automata.discrete.pfilt_y_rate_clamp_low, pfilt_y_rate_clamp_mid = ws.automata.discrete.pfilt_y_rate_clamp_mid, t_release = ws.automata.memory.t_release, t_release2 = ws.automata.memory.t_release2, Ipmax_held = ws.automata.memory.Ipmax_held),
        (s4_y = ws.var_vals.s4_y, pfilt_y = ws.var_vals.pfilt_y, s5_y = ws.var_vals.s5_y, s0_y = ws.var_vals.s0_y, S1_y = ws.var_vals.S1_y, PIV_xi = ws.var_vals.PIV_xi, PIQ_xi = ws.var_vals.PIQ_xi, Pref = ws.var_vals.Pref, pfaref = ws.var_vals.pfaref, Qref = ws.var_vals.Qref, Ipcmd_out = ws.var_vals.Ipcmd_out, Iqcmd_out = ws.var_vals.Iqcmd_out, Qerr = ws.var_vals.Qerr, Vsel = ws.var_vals.Vsel, Iqinj = ws.var_vals.Iqinj, Iqmax = ws.var_vals.Iqmax, Ipmax = ws.var_vals.Ipmax, v = ws.var_vals.v, Pe = ws.var_vals.Pe, Qgen = ws.var_vals.Qgen),
        active_mask,
        t,
        AT
    )
end

@inline function reeca1_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.PIV_unlimited_minus_Iqmax)
        tpfilt = 0.02
        wg = 1.0
        vfloor = 0.01 * params.vfloor_max_first[i] + vars.s0_y[i] * params.vfloor_max_second[i]
        Qcpf = vars.S1_y[i] * tan(vars.pfaref[i])
        PFsel = params.PFFLAG[i] * Qcpf + vars.Qref[i] * (1 - params.PFFLAG[i])
        moving = 1 - params.Voltage_dip[i]
        PIV_rate = params.Kvi[i] * params.QFLAG[i] * moving * (vars.Vsel[i] - vars.s0_y[i])
        PIV_unlimited = vars.PIV_xi[i] + params.Kvp[i] * params.QFLAG[i] * (vars.Vsel[i] - vars.s0_y[i])
        PIV_y = vars.Iqmax[i] * params.PIV_y_clamp_high[i] + params.PIV_y_clamp_mid[i] * (vars.PIV_xi[i] + params.Kvp[i] * params.QFLAG[i] * (vars.Vsel[i] - vars.s0_y[i])) - vars.Iqmax[i] * params.PIV_y_clamp_low[i]
        PIQ_rate = params.Kqi[i] * vars.Qerr[i] * params.VFLAG[i] * moving
        PIQ_unlimited = vars.PIQ_xi[i] + params.Kqp[i] * vars.Qerr[i] * params.VFLAG[i]
        PIQ_y = params.PIQ_y_clamp_high[i] * params.VMAX[i] + params.PIQ_y_clamp_low[i] * params.VMIN[i] + params.PIQ_y_clamp_mid[i] * (vars.PIQ_xi[i] + params.Kqp[i] * vars.Qerr[i] * params.VFLAG[i])
        Vref0_eff = ifelse(params.Vref0[i] == 0, params.v0[i], params.Vref0[i])
        Verr = Vref0_eff - vars.s0_y[i]
        Iqv = params.Kqv[i] * (Verr - params.Iqv_clamp_high[i] * params.dbd2[i] - params.Iqv_clamp_low[i] * params.dbd1[i] - params.Iqv_clamp_mid[i] * Verr)
        VDL1 = params.Iq1[i] * params.VDL1_vdl1_below[i] + params.Iq4[i] * params.VDL1_vdl1_above[i] + params.VDL1_vdl1_segment1[i] * (params.Iq1[i] + (params.Iq2[i] - params.Iq1[i]) * (vars.s0_y[i] - params.Vq1[i]) / (params.Vq2[i] - params.Vq1[i])) + params.VDL1_vdl1_segment2[i] * (params.Iq2[i] + (params.Iq3[i] - params.Iq2[i]) * (vars.s0_y[i] - params.Vq2[i]) / (params.Vq3[i] - params.Vq2[i])) + params.VDL1_vdl1_segment3[i] * (params.Iq3[i] + (params.Iq4[i] - params.Iq3[i]) * (vars.s0_y[i] - params.Vq3[i]) / (params.Vq4[i] - params.Vq3[i]))
        VDL2 = params.Ip1[i] * params.VDL2_vdl2_below[i] + params.Ip4[i] * params.VDL2_vdl2_above[i] + params.VDL2_vdl2_segment1[i] * (params.Ip1[i] + (params.Ip2[i] - params.Ip1[i]) * (vars.s0_y[i] - params.Vp1[i]) / (params.Vp2[i] - params.Vp1[i])) + params.VDL2_vdl2_segment2[i] * (params.Ip2[i] + (params.Ip3[i] - params.Ip2[i]) * (vars.s0_y[i] - params.Vp2[i]) / (params.Vp3[i] - params.Vp2[i])) + params.VDL2_vdl2_segment3[i] * (params.Ip3[i] + (params.Ip4[i] - params.Ip3[i]) * (vars.s0_y[i] - params.Vp3[i]) / (params.Vp4[i] - params.Vp3[i]))
        Psel = vars.pfilt_y[i] * (1 - params.PFLAG[i]) + params.PFLAG[i] * vars.pfilt_y[i] * wg
        s5_y_rate = moving * (Psel - vars.s5_y[i]) / params.Tpord[i]
        # === Automaton Values ===
        values.PIV_unlimited_minus_Iqmax[i] = PIV_unlimited - vars.Iqmax[i]
        values.PIV_rate[i] = PIV_rate
        values.PIV_value[i] = vars.Iqmax[i] + PIV_unlimited
        values.PIQ_rate[i] = PIQ_rate
        values.PIQ_unlimited_minus_VMAX[i] = PIQ_unlimited - params.VMAX[i]
        values.PIV_side[i] = params.PIV_side[i]
        values.PIQ_unlimited_minus_VMIN[i] = PIQ_unlimited - params.VMIN[i]
        values.s5_y_minus_PMAX[i] = vars.s5_y[i] - params.PMAX[i]
        values.s5_y_rate[i] = s5_y_rate
        values.s5_y_minus_PMIN[i] = vars.s5_y[i] - params.PMIN[i]
        values.v_minus_Vup[i] = vars.v[i] - params.Vup[i]
        values.v_minus_Vdip[i] = vars.v[i] - params.Vdip[i]
        values.Thld[i] = params.Thld[i]
        values.t[i] = t
        values.injection_value[i] = t - params.Thld[i] - params.t_release[i]
        values.injection_value_2[i] = params.Thld[i] + t - params.t_release[i]
        values.Ipmax[i] = vars.Ipmax[i]
        values.ipmax_hold_value[i] = t - params.Thld2[i] - params.t_release2[i]
        values.vfloor_max_value[i] = vars.s0_y[i] - 0.01
        values.PIV_y_clamp_value[i] = vars.Iqmax[i] + vars.PIV_xi[i] + params.Kvp[i] * params.QFLAG[i] * (vars.Vsel[i] - vars.s0_y[i])
        values.PIV_y_clamp_value_2[i] = vars.PIV_xi[i] + params.Kvp[i] * params.QFLAG[i] * (vars.Vsel[i] - vars.s0_y[i]) - vars.Iqmax[i]
        values.PIQ_y_clamp_value[i] = vars.PIQ_xi[i] + params.Kqp[i] * vars.Qerr[i] * params.VFLAG[i] - params.VMAX[i]
        values.PIQ_y_clamp_value_2[i] = vars.PIQ_xi[i] + params.Kqp[i] * vars.Qerr[i] * params.VFLAG[i] - params.VMIN[i]
        values.Iqv_clamp_value[i] = Verr - params.dbd2[i]
        values.Iqv_clamp_value_2[i] = Verr - params.dbd1[i]
        values.VDL1_vdl1_value[i] = vars.s0_y[i] - params.Vq1[i]
        values.VDL1_vdl1_value_2[i] = vars.s0_y[i] - params.Vq2[i]
        values.VDL1_vdl1_value_3[i] = vars.s0_y[i] - params.Vq3[i]
        values.VDL1_vdl1_value_4[i] = vars.s0_y[i] - params.Vq4[i]
        values.VDL2_vdl2_value[i] = vars.s0_y[i] - params.Vp1[i]
        values.VDL2_vdl2_value_2[i] = vars.s0_y[i] - params.Vp2[i]
        values.VDL2_vdl2_value_3[i] = vars.s0_y[i] - params.Vp3[i]
        values.VDL2_vdl2_value_4[i] = vars.s0_y[i] - params.Vp4[i]
        values.Ipcmd_out_clamp_value[i] = vars.s5_y[i] / vfloor
        values.Ipcmd_out_clamp_value_2[i] = vars.s5_y[i] / vfloor - vars.Ipmax[i]
        values.Iqcmd_out_clamp_value[i] = vars.Iqinj[i] + vars.Iqmax[i] + PIV_y * params.QFLAG[i] + vars.s4_y[i] * (1 - params.QFLAG[i])
        values.Iqcmd_out_clamp_value_2[i] = vars.Iqinj[i] + PIV_y * params.QFLAG[i] + vars.s4_y[i] * (1 - params.QFLAG[i]) - vars.Iqmax[i]
        values.Qerr_clamp_value[i] = PFsel - params.QMax[i]
        values.Qerr_clamp_value_2[i] = PFsel - params.QMin[i]
        values.Vsel_clamp_value[i] = PIQ_y * params.VFLAG[i] + (1 - params.VFLAG[i]) * (PFsel + params.Vref1[i]) - params.VMAX[i]
        values.Vsel_clamp_value_2[i] = PIQ_y * params.VFLAG[i] + (1 - params.VFLAG[i]) * (PFsel + params.Vref1[i]) - params.VMIN[i]
        values.Iqinj_clamp_value[i] = Iqv - params.Iqh1[i]
        values.Iqinj_clamp_value_2[i] = Iqv - params.Iql1[i]
        values.Iqmax_max_value[i] = params.Imax[i] ^ 2 - vars.Ipcmd_out[i] ^ 2
        values.Iqmax_min_value[i] = sqrt(params.Iqmax_max_second[i]) * sqrt(max(0, params.Imax[i] ^ 2 - vars.Ipcmd_out[i] ^ 2)) - VDL1
        values.Iqmax_min_2_value[i] = VDL1 - params.Imax[i]
        values.Ipmax_min_value[i] = VDL2 - params.Imax[i]
        values.Ipmax_max_value[i] = params.Imax[i] ^ 2 - vars.Iqcmd_out[i] ^ 2
        values.Ipmax_min_2_value[i] = sqrt(params.Ipmax_max_second[i]) * sqrt(max(0, params.Imax[i] ^ 2 - vars.Iqcmd_out[i] ^ 2)) - VDL2
        values.pfilt_y_rate_clamp_value[i] = (vars.Pref[i] - vars.pfilt_y[i]) / tpfilt - params.dPmax[i]
        values.pfilt_y_rate_clamp_value_2[i] = (vars.Pref[i] - vars.pfilt_y[i]) / tpfilt - params.dPmin[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:REECA1},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    reeca1_automaton_values_kernel!(
        (PIV_unlimited_minus_Iqmax = ws.automata.values.PIV_unlimited_minus_Iqmax, PIV_rate = ws.automata.values.PIV_rate, PIV_value = ws.automata.values.PIV_value, PIQ_rate = ws.automata.values.PIQ_rate, PIQ_unlimited_minus_VMAX = ws.automata.values.PIQ_unlimited_minus_VMAX, PIV_side = ws.automata.values.PIV_side, PIQ_unlimited_minus_VMIN = ws.automata.values.PIQ_unlimited_minus_VMIN, s5_y_minus_PMAX = ws.automata.values.s5_y_minus_PMAX, s5_y_rate = ws.automata.values.s5_y_rate, s5_y_minus_PMIN = ws.automata.values.s5_y_minus_PMIN, v_minus_Vup = ws.automata.values.v_minus_Vup, v_minus_Vdip = ws.automata.values.v_minus_Vdip, Thld = ws.automata.values.Thld, t = ws.automata.values.t, injection_value = ws.automata.values.injection_value, injection_value_2 = ws.automata.values.injection_value_2, Ipmax = ws.automata.values.Ipmax, ipmax_hold_value = ws.automata.values.ipmax_hold_value, vfloor_max_value = ws.automata.values.vfloor_max_value, PIV_y_clamp_value = ws.automata.values.PIV_y_clamp_value, PIV_y_clamp_value_2 = ws.automata.values.PIV_y_clamp_value_2, PIQ_y_clamp_value = ws.automata.values.PIQ_y_clamp_value, PIQ_y_clamp_value_2 = ws.automata.values.PIQ_y_clamp_value_2, Iqv_clamp_value = ws.automata.values.Iqv_clamp_value, Iqv_clamp_value_2 = ws.automata.values.Iqv_clamp_value_2, VDL1_vdl1_value = ws.automata.values.VDL1_vdl1_value, VDL1_vdl1_value_2 = ws.automata.values.VDL1_vdl1_value_2, VDL1_vdl1_value_3 = ws.automata.values.VDL1_vdl1_value_3, VDL1_vdl1_value_4 = ws.automata.values.VDL1_vdl1_value_4, VDL2_vdl2_value = ws.automata.values.VDL2_vdl2_value, VDL2_vdl2_value_2 = ws.automata.values.VDL2_vdl2_value_2, VDL2_vdl2_value_3 = ws.automata.values.VDL2_vdl2_value_3, VDL2_vdl2_value_4 = ws.automata.values.VDL2_vdl2_value_4, Ipcmd_out_clamp_value = ws.automata.values.Ipcmd_out_clamp_value, Ipcmd_out_clamp_value_2 = ws.automata.values.Ipcmd_out_clamp_value_2, Iqcmd_out_clamp_value = ws.automata.values.Iqcmd_out_clamp_value, Iqcmd_out_clamp_value_2 = ws.automata.values.Iqcmd_out_clamp_value_2, Qerr_clamp_value = ws.automata.values.Qerr_clamp_value, Qerr_clamp_value_2 = ws.automata.values.Qerr_clamp_value_2, Vsel_clamp_value = ws.automata.values.Vsel_clamp_value, Vsel_clamp_value_2 = ws.automata.values.Vsel_clamp_value_2, Iqinj_clamp_value = ws.automata.values.Iqinj_clamp_value, Iqinj_clamp_value_2 = ws.automata.values.Iqinj_clamp_value_2, Iqmax_max_value = ws.automata.values.Iqmax_max_value, Iqmax_min_value = ws.automata.values.Iqmax_min_value, Iqmax_min_2_value = ws.automata.values.Iqmax_min_2_value, Ipmax_min_value = ws.automata.values.Ipmax_min_value, Ipmax_max_value = ws.automata.values.Ipmax_max_value, Ipmax_min_2_value = ws.automata.values.Ipmax_min_2_value, pfilt_y_rate_clamp_value = ws.automata.values.pfilt_y_rate_clamp_value, pfilt_y_rate_clamp_value_2 = ws.automata.values.pfilt_y_rate_clamp_value_2),
        (PFFLAG = model.PFFLAG, VFLAG = model.VFLAG, QFLAG = model.QFLAG, PFLAG = model.PFLAG, PQFLAG = model.PQFLAG, Vdip = model.Vdip, Vup = model.Vup, Trv = model.Trv, dbd1 = model.dbd1, dbd2 = model.dbd2, Kqv = model.Kqv, Iqh1 = model.Iqh1, Iql1 = model.Iql1, Vref0 = model.Vref0, Iqfrz = model.Iqfrz, Thld = model.Thld, Thld2 = model.Thld2, Tp = model.Tp, QMax = model.QMax, QMin = model.QMin, VMAX = model.VMAX, VMIN = model.VMIN, Kqp = model.Kqp, Kqi = model.Kqi, Kvp = model.Kvp, Kvi = model.Kvi, Vref1 = model.Vref1, Tiq = model.Tiq, dPmax = model.dPmax, dPmin = model.dPmin, PMAX = model.PMAX, PMIN = model.PMIN, Imax = model.Imax, Tpord = model.Tpord, Vq1 = model.Vq1, Iq1 = model.Iq1, Vq2 = model.Vq2, Iq2 = model.Iq2, Vq3 = model.Vq3, Iq3 = model.Iq3, Vq4 = model.Vq4, Iq4 = model.Iq4, Vp1 = model.Vp1, Ip1 = model.Ip1, Vp2 = model.Vp2, Ip2 = model.Ip2, Vp3 = model.Vp3, Ip3 = model.Ip3, Vp4 = model.Vp4, Ip4 = model.Ip4, Ipcmd0 = ws.setpoints.Ipcmd0, Iqcmd0 = ws.setpoints.Iqcmd0, v0 = ws.setpoints.v0, Pref0 = ws.setpoints.Pref0, pfaref0 = ws.setpoints.pfaref0, Qref0 = ws.setpoints.Qref0, Voltage_dip = ws.automata.discrete.Voltage_dip, in_state1 = ws.automata.discrete.in_state1, in_state2 = ws.automata.discrete.in_state2, Ipmax_hold = ws.automata.discrete.Ipmax_hold, PIV_gate = ws.automata.discrete.PIV_gate, PIV_side = ws.automata.discrete.PIV_side, PIQ_gate = ws.automata.discrete.PIQ_gate, PIQ_side = ws.automata.discrete.PIQ_side, s5_y_gate = ws.automata.discrete.s5_y_gate, vfloor_max_second = ws.automata.discrete.vfloor_max_second, vfloor_max_first = ws.automata.discrete.vfloor_max_first, PIV_y_clamp_high = ws.automata.discrete.PIV_y_clamp_high, PIV_y_clamp_low = ws.automata.discrete.PIV_y_clamp_low, PIV_y_clamp_mid = ws.automata.discrete.PIV_y_clamp_mid, PIQ_y_clamp_high = ws.automata.discrete.PIQ_y_clamp_high, PIQ_y_clamp_low = ws.automata.discrete.PIQ_y_clamp_low, PIQ_y_clamp_mid = ws.automata.discrete.PIQ_y_clamp_mid, Iqv_clamp_high = ws.automata.discrete.Iqv_clamp_high, Iqv_clamp_low = ws.automata.discrete.Iqv_clamp_low, Iqv_clamp_mid = ws.automata.discrete.Iqv_clamp_mid, VDL1_vdl1_below = ws.automata.discrete.VDL1_vdl1_below, VDL1_vdl1_segment1 = ws.automata.discrete.VDL1_vdl1_segment1, VDL1_vdl1_segment2 = ws.automata.discrete.VDL1_vdl1_segment2, VDL1_vdl1_segment3 = ws.automata.discrete.VDL1_vdl1_segment3, VDL1_vdl1_above = ws.automata.discrete.VDL1_vdl1_above, VDL2_vdl2_below = ws.automata.discrete.VDL2_vdl2_below, VDL2_vdl2_segment1 = ws.automata.discrete.VDL2_vdl2_segment1, VDL2_vdl2_segment2 = ws.automata.discrete.VDL2_vdl2_segment2, VDL2_vdl2_segment3 = ws.automata.discrete.VDL2_vdl2_segment3, VDL2_vdl2_above = ws.automata.discrete.VDL2_vdl2_above, Ipcmd_out_clamp_high = ws.automata.discrete.Ipcmd_out_clamp_high, Ipcmd_out_clamp_low = ws.automata.discrete.Ipcmd_out_clamp_low, Ipcmd_out_clamp_mid = ws.automata.discrete.Ipcmd_out_clamp_mid, Iqcmd_out_clamp_high = ws.automata.discrete.Iqcmd_out_clamp_high, Iqcmd_out_clamp_low = ws.automata.discrete.Iqcmd_out_clamp_low, Iqcmd_out_clamp_mid = ws.automata.discrete.Iqcmd_out_clamp_mid, Qerr_clamp_high = ws.automata.discrete.Qerr_clamp_high, Qerr_clamp_low = ws.automata.discrete.Qerr_clamp_low, Qerr_clamp_mid = ws.automata.discrete.Qerr_clamp_mid, Vsel_clamp_high = ws.automata.discrete.Vsel_clamp_high, Vsel_clamp_low = ws.automata.discrete.Vsel_clamp_low, Vsel_clamp_mid = ws.automata.discrete.Vsel_clamp_mid, Iqinj_clamp_high = ws.automata.discrete.Iqinj_clamp_high, Iqinj_clamp_low = ws.automata.discrete.Iqinj_clamp_low, Iqinj_clamp_mid = ws.automata.discrete.Iqinj_clamp_mid, Iqmax_max_second = ws.automata.discrete.Iqmax_max_second, Iqmax_max_first = ws.automata.discrete.Iqmax_max_first, Iqmax_min_second = ws.automata.discrete.Iqmax_min_second, Iqmax_min_first = ws.automata.discrete.Iqmax_min_first, Iqmax_min_2_second = ws.automata.discrete.Iqmax_min_2_second, Iqmax_min_2_first = ws.automata.discrete.Iqmax_min_2_first, Ipmax_min_second = ws.automata.discrete.Ipmax_min_second, Ipmax_min_first = ws.automata.discrete.Ipmax_min_first, Ipmax_max_second = ws.automata.discrete.Ipmax_max_second, Ipmax_max_first = ws.automata.discrete.Ipmax_max_first, Ipmax_min_2_second = ws.automata.discrete.Ipmax_min_2_second, Ipmax_min_2_first = ws.automata.discrete.Ipmax_min_2_first, pfilt_y_rate_clamp_high = ws.automata.discrete.pfilt_y_rate_clamp_high, pfilt_y_rate_clamp_low = ws.automata.discrete.pfilt_y_rate_clamp_low, pfilt_y_rate_clamp_mid = ws.automata.discrete.pfilt_y_rate_clamp_mid, t_release = ws.automata.memory.t_release, t_release2 = ws.automata.memory.t_release2, Ipmax_held = ws.automata.memory.Ipmax_held),
        (s4_y = ws.var_vals.s4_y, pfilt_y = ws.var_vals.pfilt_y, s5_y = ws.var_vals.s5_y, s0_y = ws.var_vals.s0_y, S1_y = ws.var_vals.S1_y, PIV_xi = ws.var_vals.PIV_xi, PIQ_xi = ws.var_vals.PIQ_xi, Pref = ws.var_vals.Pref, pfaref = ws.var_vals.pfaref, Qref = ws.var_vals.Qref, Ipcmd_out = ws.var_vals.Ipcmd_out, Iqcmd_out = ws.var_vals.Iqcmd_out, Qerr = ws.var_vals.Qerr, Vsel = ws.var_vals.Vsel, Iqinj = ws.var_vals.Iqinj, Iqmax = ws.var_vals.Iqmax, Ipmax = ws.var_vals.Ipmax, v = ws.var_vals.v, Pe = ws.var_vals.Pe, Qgen = ws.var_vals.Qgen),
        t,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{REECA1{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, REECA1)

    v = T.(data.v)
    Ipcmd_init = T.(data.Ipcmd_init)
    Iqcmd_init = T.(data.Iqcmd_init)
    Pe_init = T.(data.Pe_init)
    Qgen_init = T.(data.Qgen_init)

    # === Pre-allocate Result Arrays ===
    s4_y = Vector{T}(undef, n)
    pfilt_y = Vector{T}(undef, n)
    s5_y = Vector{T}(undef, n)
    s0_y = Vector{T}(undef, n)
    S1_y = Vector{T}(undef, n)
    PIV_xi = Vector{T}(undef, n)
    PIQ_xi = Vector{T}(undef, n)
    Pref = Vector{T}(undef, n)
    pfaref = Vector{T}(undef, n)
    Qref = Vector{T}(undef, n)
    Ipcmd_out = Vector{T}(undef, n)
    Iqcmd_out = Vector{T}(undef, n)
    Qerr = Vector{T}(undef, n)
    Vsel = Vector{T}(undef, n)
    Iqinj = Vector{T}(undef, n)
    Iqmax = Vector{T}(undef, n)
    Ipmax = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # injection: the mode of its first initial clause that holds
        Voltage_dip = ifelse(ifelse(v[i] > model.Vup[i], true, v[i] < model.Vdip[i]), 1.0, 0.0)
        in_state1 = ifelse(ifelse(v[i] > model.Vup[i], true, v[i] < model.Vdip[i]), 1.0, 0.0)
        in_state2 = 0.0

        # State initial values
        s0_y[i] = v[i]
        S1_y[i] = Pe_init[i]

        # Algebraic initial values
        Iqinj[i] = model.Iqfrz[i] * in_state2 + in_state1 * clamp(model.Kqv[i] * (ifelse(model.Vref0[i] == 0, v[i], model.Vref0[i]) - v[i] - clamp(ifelse(model.Vref0[i] == 0, v[i], model.Vref0[i]) - v[i], model.dbd1[i], model.dbd2[i])), model.Iql1[i], model.Iqh1[i])
        Qref[i] = model.QFLAG[i] * (Qgen_init[i] * model.VFLAG[i] + (1 - model.VFLAG[i]) * (v[i] - model.Vref1[i])) + (1 - model.QFLAG[i]) * (Qgen_init[i] - Iqinj[i] * max(0.01, v[i]))
        pfaref[i] = atan(Qref[i], S1_y[i])
        Pref[i] = Ipcmd_init[i] * max(0.01, v[i])
        Ipcmd_out[i] = Ipcmd_init[i]
        Iqmax[i] = model.PQFLAG[i] * min(sqrt(max(0, model.Imax[i] ^ 2 - Ipcmd_out[i] ^ 2)), model.Iq1[i] + (model.Iq2[i] - model.Iq1[i]) * (clamp(s0_y[i], model.Vq1[i], model.Vq2[i]) - model.Vq1[i]) / (model.Vq2[i] - model.Vq1[i]) + (model.Iq3[i] - model.Iq2[i]) * (clamp(s0_y[i], model.Vq2[i], model.Vq3[i]) - model.Vq2[i]) / (model.Vq3[i] - model.Vq2[i]) + (model.Iq4[i] - model.Iq3[i]) * (clamp(s0_y[i], model.Vq3[i], model.Vq4[i]) - model.Vq3[i]) / (model.Vq4[i] - model.Vq3[i])) + (1 - model.PQFLAG[i]) * min(model.Iq1[i] + (model.Iq2[i] - model.Iq1[i]) * (clamp(s0_y[i], model.Vq1[i], model.Vq2[i]) - model.Vq1[i]) / (model.Vq2[i] - model.Vq1[i]) + (model.Iq3[i] - model.Iq2[i]) * (clamp(s0_y[i], model.Vq2[i], model.Vq3[i]) - model.Vq2[i]) / (model.Vq3[i] - model.Vq2[i]) + (model.Iq4[i] - model.Iq3[i]) * (clamp(s0_y[i], model.Vq3[i], model.Vq4[i]) - model.Vq3[i]) / (model.Vq4[i] - model.Vq3[i]), model.Imax[i])
        Iqcmd_out[i] = -Iqcmd_init[i]
        Qerr[i] = clamp(Qref[i] * (1 - model.PFFLAG[i]) + model.PFFLAG[i] * S1_y[i] * tan(pfaref[i]), model.QMin[i], model.QMax[i]) - Qgen_init[i]

        # State initial values
        PIQ_xi[i] = ifelse(ifelse(model.QFLAG[i] == 1, ifelse(v[i] >= model.Vdip[i], v[i] <= model.Vup[i], false), false), v[i] - model.Kqp[i] * Qerr[i] * model.VFLAG[i], clamp(v[i] - model.Kqp[i] * Qerr[i] * model.VFLAG[i], model.VMIN[i], model.VMAX[i]))

        # Algebraic initial values
        Vsel[i] = clamp(model.VFLAG[i] * clamp(PIQ_xi[i] + model.Kqp[i] * Qerr[i] * model.VFLAG[i], model.VMIN[i], model.VMAX[i]) + (1 - model.VFLAG[i]) * (model.Vref1[i] + Qref[i] * (1 - model.PFFLAG[i]) + model.PFFLAG[i] * S1_y[i] * tan(pfaref[i])), model.VMIN[i], model.VMAX[i])

        # State initial values
        PIV_xi[i] = ifelse(ifelse(model.QFLAG[i] == 1, ifelse(v[i] >= model.Vdip[i], v[i] <= model.Vup[i], false), false), -Iqcmd_init[i] - Iqinj[i] - model.Kvp[i] * model.QFLAG[i] * (Vsel[i] - v[i]), clamp(model.QFLAG[i] * (-Iqcmd_init[i] - Iqinj[i] - model.Kvp[i] * model.QFLAG[i] * (Vsel[i] - v[i])), -Iqmax[i], Iqmax[i]))
        s4_y[i] = (Qref[i] * (1 - model.PFFLAG[i]) + model.PFFLAG[i] * S1_y[i] * tan(pfaref[i])) / max(0.01, v[i])
        pfilt_y[i] = Pref[i]
        s5_y[i] = Pref[i]

        # Algebraic initial values
        Ipmax[i] = model.PQFLAG[i] * min(model.Ip1[i] + (model.Ip2[i] - model.Ip1[i]) * (clamp(s0_y[i], model.Vp1[i], model.Vp2[i]) - model.Vp1[i]) / (model.Vp2[i] - model.Vp1[i]) + (model.Ip3[i] - model.Ip2[i]) * (clamp(s0_y[i], model.Vp2[i], model.Vp3[i]) - model.Vp2[i]) / (model.Vp3[i] - model.Vp2[i]) + (model.Ip4[i] - model.Ip3[i]) * (clamp(s0_y[i], model.Vp3[i], model.Vp4[i]) - model.Vp3[i]) / (model.Vp4[i] - model.Vp3[i]), model.Imax[i]) + (1 - model.PQFLAG[i]) * min(sqrt(max(0, model.Imax[i] ^ 2 - Iqcmd_out[i] ^ 2)), model.Ip1[i] + (model.Ip2[i] - model.Ip1[i]) * (clamp(s0_y[i], model.Vp1[i], model.Vp2[i]) - model.Vp1[i]) / (model.Vp2[i] - model.Vp1[i]) + (model.Ip3[i] - model.Ip2[i]) * (clamp(s0_y[i], model.Vp2[i], model.Vp3[i]) - model.Vp2[i]) / (model.Vp3[i] - model.Vp2[i]) + (model.Ip4[i] - model.Ip3[i]) * (clamp(s0_y[i], model.Vp3[i], model.Vp4[i]) - model.Vp3[i]) / (model.Vp4[i] - model.Vp3[i]))
    end

    # === Return Initial Values ===
    return (s4_y=s4_y, pfilt_y=pfilt_y, s5_y=s5_y, s0_y=s0_y, S1_y=S1_y, PIV_xi=PIV_xi, PIQ_xi=PIQ_xi, Pref=Pref, pfaref=pfaref, Qref=Qref, Ipcmd_out=Ipcmd_out, Iqcmd_out=Iqcmd_out, Qerr=Qerr, Vsel=Vsel, Iqinj=Iqinj, Iqmax=Iqmax, Ipmax=Ipmax)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{REECA1}, dyr_data, sys_pf, links) -> StructArray{REECA1}

Parse REECA1 model from DYR data.
"""
function parse_dynamic_model(
    ::Type{REECA1},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    ::Nothing=nothing
)
    n = length(dyr_data)
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{REECA1{T}}(undef, n)

    for idx in 1:n
        dyr_row = dyr_data[idx]
        instances[idx] = check_parameter_domains(REECA1(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            VREG = Int32(dyr_row.VREG),
            PFFLAG = T(dyr_row.PFFLAG),
            VFLAG = T(dyr_row.VFLAG),
            QFLAG = T(dyr_row.QFLAG),
            PFLAG = T(dyr_row.PFLAG),
            PQFLAG = T(dyr_row.PQFLAG),
            Vdip = T(dyr_row.Vdip),
            Vup = T(dyr_row.Vup),
            Trv = T(dyr_row.Trv),
            dbd1 = T(dyr_row.dbd1),
            dbd2 = T(dyr_row.dbd2),
            Kqv = T(dyr_row.Kqv),
            Iqh1 = T(dyr_row.Iqh1),
            Iql1 = T(dyr_row.Iql1),
            Vref0 = T(dyr_row.Vref0),
            Iqfrz = T(dyr_row.Iqfrz),
            Thld = T(dyr_row.Thld),
            Thld2 = T(dyr_row.Thld2),
            Tp = T(dyr_row.Tp),
            QMax = T(dyr_row.QMax),
            QMin = T(dyr_row.QMin),
            VMAX = T(dyr_row.VMAX),
            VMIN = T(dyr_row.VMIN),
            Kqp = T(dyr_row.Kqp),
            Kqi = T(dyr_row.Kqi),
            Kvp = T(dyr_row.Kvp),
            Kvi = T(dyr_row.Kvi),
            Vref1 = T(dyr_row.Vref1),
            Tiq = T(dyr_row.Tiq),
            dPmax = T(dyr_row.dPmax),
            dPmin = T(dyr_row.dPmin),
            PMAX = T(dyr_row.PMAX),
            PMIN = T(dyr_row.PMIN),
            Imax = T(dyr_row.Imax),
            Tpord = T(dyr_row.Tpord),
            Vq1 = T(dyr_row.Vq1),
            Iq1 = T(dyr_row.Iq1),
            Vq2 = T(dyr_row.Vq2),
            Iq2 = T(dyr_row.Iq2),
            Vq3 = T(dyr_row.Vq3),
            Iq3 = T(dyr_row.Iq3),
            Vq4 = T(dyr_row.Vq4),
            Iq4 = T(dyr_row.Iq4),
            Vp1 = T(dyr_row.Vp1),
            Ip1 = T(dyr_row.Ip1),
            Vp2 = T(dyr_row.Vp2),
            Ip2 = T(dyr_row.Ip2),
            Vp3 = T(dyr_row.Vp3),
            Ip3 = T(dyr_row.Ip3),
            Vp4 = T(dyr_row.Vp4),
            Ip4 = T(dyr_row.Ip4)
        ), "dyr record $(idx)")
    end

    return StructArray(instances)
end
