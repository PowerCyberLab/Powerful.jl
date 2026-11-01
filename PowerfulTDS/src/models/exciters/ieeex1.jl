"""IEEEX1: IEEE Type 1 Excitation System with a transient gain reduction

DO NOT EDIT - changes will be overwritten.
"""


using InlineStrings
using PowerfulCore: Automaton, Mode, Transition, InitialMode, Condition, Atom, AllOf

export IEEEX1

struct IEEEX1{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number"
    BUS::Int32
    "Generator ID"
    ID::String3

    # === Parameters ===
    TR::T  # Terminal voltage transducer time constant
    KA::T  # Regulator gain
    TA::T  # Regulator time constant
    TB::T  # Lag time constant of the transient gain reduction, zero only with TC zero
    TC::T  # Lead time constant of the transient gain reduction
    VRMAX::T  # Maximum regulator output
    VRMIN::T  # Minimum regulator output
    KE::T  # Exciter self-excitation constant
    TE::T  # Exciter time constant
    KF::T  # Rate feedback gain: the feedback is s KF / (1 + s TF), so KF carries the unit of time
    TF::T  # Rate feedback time constant
    Switch::T  # PSS/E record field the diagram does not draw; only 0 is accepted
    E1::T  # First exciter output voltage at which saturation is given
    SE1::T  # Saturation factor at E1
    E2::T  # Second exciter output voltage at which saturation is given
    SE2::T  # Saturation factor at E2; zero with SE1 zero is no saturation (the points must still be positive and different)

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function IEEEX1(; BUS::Int32, ID, TR, KA, TA, TB, TC, VRMAX, VRMIN, KE, TE, KF, TF, Switch, E1, SE1, E2, SE2, status::Bool=true)
    T = typeof(TR)
    check_parameter_domains(IEEEX1{T}(BUS, ID, T(TR), T(KA), T(TA), T(TB), T(TC), T(VRMAX), T(VRMIN), T(KE), T(TE), T(KF), T(TF), T(Switch), T(E1), T(SE1), T(E2), T(SE2), status),
                            "keyword construction")
end

# Parameter domains, checked where data enters
function check_parameter_domains(x::IEEEX1, from::AbstractString)
    (x.TR >= 0.0) || throw(ArgumentError("IEEEX1 ($from): TR = $(x.TR) is not >= 0"))
    (x.KA >= 0.0) || throw(ArgumentError("IEEEX1 ($from): KA = $(x.KA) is not >= 0"))
    (x.TA >= 0.0) || throw(ArgumentError("IEEEX1 ($from): TA = $(x.TA) is not >= 0"))
    (x.TB >= 0.0) || throw(ArgumentError("IEEEX1 ($from): TB = $(x.TB) is not >= 0"))
    (x.TC >= 0.0) || throw(ArgumentError("IEEEX1 ($from): TC = $(x.TC) is not >= 0"))
    (x.TE >= 0.0) || throw(ArgumentError("IEEEX1 ($from): TE = $(x.TE) is not >= 0"))
    (x.TF >= 0.0) || throw(ArgumentError("IEEEX1 ($from): TF = $(x.TF) is not >= 0"))
    (x.Switch >= 0.0 && x.Switch <= 0.0) || throw(ArgumentError("IEEEX1 ($from): Switch = $(x.Switch) is not >= 0 and <= 0"))
    (x.SE1 >= 0.0) || throw(ArgumentError("IEEEX1 ($from): SE1 = $(x.SE1) is not >= 0"))
    (x.SE2 >= 0.0) || throw(ArgumentError("IEEEX1 ($from): SE2 = $(x.SE2) is not >= 0"))
    x.SE1 >= 0.0 || throw(ArgumentError("IEEEX1 ($from) does not satisfy its constraint sat0_S1_min: QuadraticSaturation admits S1 >= 0"))
    x.SE2 >= 0.0 || throw(ArgumentError("IEEEX1 ($from) does not satisfy its constraint sat0_S2_min: QuadraticSaturation admits S2 >= 0"))
    ifelse(0 < x.E1, ifelse(0 < x.E2, x.E1 != x.E2, false), false) || throw(ArgumentError("IEEEX1 ($from) does not satisfy its constraint sat0_points_distinct: The points are at positive, different magnitudes"))
    ifelse(ifelse(x.E1 < x.E2, x.E2 * x.SE1 <= x.E1 * x.SE2, false), true, ifelse(x.E2 < x.E1, x.E1 * x.SE2 <= x.E2 * x.SE1, false)) || throw(ArgumentError("IEEEX1 ($from) does not satisfy its constraint sat0_knee_nonnegative: The saturation per unit of magnitude does not fall from the lower point to the upper, so the fitted knee A is at or above zero"))
    x.VRMIN <= x.VRMAX || throw(ArgumentError("IEEEX1 ($from) does not satisfy its constraint VR_limits_ordered: VR: the lower limit of VR.x is at most its upper limit"))
    (x.TF != 0) || throw(ArgumentError("IEEEX1 ($from): TF is zero, and the model divides by it (services.kfw)"))
    (ifelse(x.E1 < x.E2, x.SE2, x.SE1) <= 0 || max(x.E1, x.E2) != 0) || throw(ArgumentError("IEEEX1 ($from): max(E1, E2) is zero, and the model divides by it (services.sat0.services.r)"))
    (ifelse(x.E1 < x.E2, x.SE2, x.SE1) <= 0 || ifelse(ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0, ifelse(x.E1 < x.E2, x.SE2, x.SE1), 1) != 0) || throw(ArgumentError("IEEEX1 ($from): ifelse(ifelse(E1 < E2, SE2, SE1) > 0, ifelse(E1 < E2, SE2, SE1), 1) is zero, and the model divides by it (services.sat0.services.r)"))
    (1 - (ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0 ? sqrt(min(x.E1, x.E2) * ifelse(x.E1 < x.E2, x.SE1, x.SE2) / (max(x.E1, x.E2) * ifelse(ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0, ifelse(x.E1 < x.E2, x.SE2, x.SE1), 1))) : 0) != 0) || throw(ArgumentError("IEEEX1 ($from): 1 - (ifelse(E1 < E2, SE2, SE1) > 0 ? sqrt(min(E1, E2) * ifelse(E1 < E2, SE1, SE2) / (max(E1, E2) * ifelse(ifelse(E1 < E2, SE2, SE1) > 0, ifelse(E1 < E2, SE2, SE1), 1))) : 0) is zero, and the model divides by it (services.sat0.services.A)"))
    (max(x.E1, x.E2) - (min(x.E1, x.E2) - max(x.E1, x.E2) * (ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0 ? sqrt(min(x.E1, x.E2) * ifelse(x.E1 < x.E2, x.SE1, x.SE2) / (max(x.E1, x.E2) * ifelse(ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0, ifelse(x.E1 < x.E2, x.SE2, x.SE1), 1))) : 0)) / (1 - (ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0 ? sqrt(min(x.E1, x.E2) * ifelse(x.E1 < x.E2, x.SE1, x.SE2) / (max(x.E1, x.E2) * ifelse(ifelse(x.E1 < x.E2, x.SE2, x.SE1) > 0, ifelse(x.E1 < x.E2, x.SE2, x.SE1), 1))) : 0)) != 0) || throw(ArgumentError("IEEEX1 ($from): max(E1, E2) - (min(E1, E2) - max(E1, E2) * (ifelse(E1 < E2, SE2, SE1) > 0 ? sqrt(min(E1, E2) * ifelse(E1 < E2, SE1, SE2) / (max(E1, E2) * ifelse(ifelse(E1 < E2, SE2, SE1) > 0, ifelse(E1 < E2, SE2, SE1), 1))) : 0)) / (1 - (ifelse(E1 < E2, SE2, SE1) > 0 ? sqrt(min(E1, E2) * ifelse(E1 < E2, SE1, SE2) / (max(E1, E2) * ifelse(ifelse(E1 < E2, SE2, SE1) > 0, ifelse(E1 < E2, SE2, SE1), 1))) : 0)) is zero, and the model divides by it (services.sat0.services.B)"))
    (x.TE != 0) || throw(ArgumentError("IEEEX1 ($from): TE is zero, and the model divides by it (graph[6].efd_rate)"))
    (x.KA != 0) || throw(ArgumentError("IEEEX1 ($from): KA is zero, and the model divides by it (algebraics.vref.init)"))
    (ifelse(x.TB != 0, true, x.TC == 0) || x.TB != 0) || throw(ArgumentError("IEEEX1 ($from): TB is zero, and the model divides by it (graph[2].LL)"))
    return x
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:IEEEX1}) = true
get_dyr_model_name(::Type{<:IEEEX1}) = "IEEEX1"
get_controlled_group(::Type{<:IEEEX1}) = :generators


# ============================================================================
# Group Interface Implementation - :exciters
# ============================================================================

# Declare that IEEEX1 implements the :exciters group interface
implements_group(::Type{<:IEEEX1}) = :exciters

# Map group canonical variable names to IEEEX1's actual variable names
function group_var_mapping(::Type{<:IEEEX1})
    Dict(
        :E_fd => :E_fd, # Same name
        :vref => :vref, # Same name
    )
end

# Map group canonical residual names to IEEEX1's actual residual names
function group_res_mapping(::Type{<:IEEEX1})
    Dict(
        :g_vref => :g_vref, # Same name
    )
end

# Controllers don't replace anything - they control generators
get_replacement_map(::Type{<:IEEEX1}) = OrderedDict{Symbol, NamedTuple}()

# No data inheritance for controllers (all parameters from DYR)
function get_data_inheritance(::Type{<:IEEEX1})
    OrderedDict{Symbol, OrderedDict{Symbol, Symbol}}()
end

# ============================================================================
# Model Registration
# ============================================================================

supports_format(::Type{<:IEEEX1}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model IEEEX1

# Terminal declarations - where IEEEX1 connects
const IEEEX1_TERMINALS = (
    Terminal(:bus_group, :BUS),
    Terminal(:generators, [:BUS, :ID]),
)

get_terminals(::Type{<:IEEEX1}) = IEEEX1_TERMINALS

get_key_name(::Type{<:IEEEX1}) = [:BUS, :ID]

function model_metadata(::Type{IEEEX1})
    ModelMetadata(
        name = :IEEEX1,
        vars = [
            ModelVar(:LG_x, StateVar();
                description="Terminal voltage transducer"),
            ModelVar(:LL_x, StateVar();
                description="Transient gain reduction"),
            ModelVar(:VR_x, StateVar();
                description="Regulator, its state held within [VRMIN, VRMAX]"),
            ModelVar(:EFD_y, StateVar();
                description="Exciter output, with no limit"),
            ModelVar(:WF_x, StateVar();
                description="Rate feedback of the exciter output"),
            ModelVar(:vref, AlgebVar();
                description="Reference voltage setpoint (control input; the stabilizer adds to it)"),
            ModelVar(:yfb, AlgebVar();
                description="Rate feedback signal"),
            ModelVar(:E_fd, AlgebVar();
                description="Field voltage output, read before the integrator that holds it"),
            ModelVar(:v, AlgebVar(), :bus_group, :v, :BUS;
                description="Terminal voltage magnitude (not negative) (external)", units="pu")
        ],
        residuals = [
        ModelResidual(:dLG_x, StateRes(), PrivateRes();
            description="Terminal voltage transducer"),
        ModelResidual(:dLL_x, StateRes(), PrivateRes();
            description="Transient gain reduction"),
        ModelResidual(:dVR_x, StateRes(), PrivateRes();
            description="Regulator, its state held within [VRMIN, VRMAX]"),
        ModelResidual(:dEFD_y, StateRes(), PrivateRes();
            description="Exciter output, with no limit"),
        ModelResidual(:dWF_x, StateRes(), PrivateRes();
            description="Rate feedback of the exciter output"),
        ModelResidual(:g_vref, AlgebRes(), SharedRes(:add);
            description="Reference voltage setpoint (control input; the stabilizer adds to it)"),
        ModelResidual(:g_yfb, AlgebRes(), PrivateRes();
            description="Rate feedback signal"),
        ModelResidual(:g_E_fd, AlgebRes(), PrivateRes();
            description="Field voltage output, read before the integrator that holds it"),
        ModelResidual(:E_fd_out, AlgebRes(), SharedRes(:add), :generators, :g_E_fd, [:BUS, :ID];
            description="Field voltage contribution (delta from initial)")
    ],
        ext_params = ExtParam[],
        init_inputs = [
        InitInput(:E_fd0, :from_group_init, :generators, :E_fd;
            description="Field voltage from linked generator initialization"),
        InitInput(:v, :from_solution, :Bus, :v;
            description="Bus voltage from power flow solution")
    ],
        observable_vars = Symbol[],
        initialized_params = [
        InitParam(:E_fd0, :E_fd; description="Initial field voltage from generator"),
        InitParam(:vref0, :vref; description="Initial reference voltage setpoint")
    ],
        variable_limits = [VariableLimit(:VR_x, :VRMIN, :VRMAX)],
        automata = (
            Automaton(:VR;
                modes = [Mode(:free; holds = [:VR_gate => 1.0]),
                    Mode(:upper; holds = [:VR_gate => 0.0]),
                    Mode(:lower; holds = [:VR_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:VR_x_minus_VRMAX, :>=), Atom(:VR_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:VR_x_minus_VRMIN, :<=), Atom(:VR_rate, :<)])),
                    Transition(:upper, :free; when = Atom(:VR_rate, :<)),
                    Transition(:lower, :free; when = Atom(:VR_rate, :>))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:VR_x_minus_VRMAX, :>=), Atom(:VR_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:VR_x_minus_VRMIN, :<=), Atom(:VR_rate, :<)])),
                    InitialMode(:free)]),
        ),
        status_field = :status,
        evaluation_pass = 2,
    )
end

@inline function ieeex1_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t, ::Type{AT}
) where AT
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        kfw = params.KF[i] / params.TF[i]
        sat0_r = (ifelse(params.E1[i] < params.E2[i], params.SE2[i], params.SE1[i]) > 0 ? sqrt(min(params.E1[i], params.E2[i]) * ifelse(params.E1[i] < params.E2[i], params.SE1[i], params.SE2[i]) / (max(params.E1[i], params.E2[i]) * ifelse(ifelse(params.E1[i] < params.E2[i], params.SE2[i], params.SE1[i]) > 0, ifelse(params.E1[i] < params.E2[i], params.SE2[i], params.SE1[i]), 1))) : 0)
        sat0_A = (min(params.E1[i], params.E2[i]) - max(params.E1[i], params.E2[i]) * sat0_r) / (1 - sat0_r)
        sat0_B = max(params.E1[i], params.E2[i]) * ifelse(params.E1[i] < params.E2[i], params.SE2[i], params.SE1[i]) / (max(params.E1[i], params.E2[i]) - sat0_A) ^ 2

        # --- LG: LagBypass ---
        LG_y = ifelse((params.TR[i]) == 0, (1.0) * (vars.v[i]), vars.LG_x[i])
        LG_dx = ifelse((params.TR[i]) == 0, 0, ((1.0) * (vars.v[i]) - vars.LG_x[i]) / (params.TR[i]))

        err = vars.vref[i] - LG_y - vars.yfb[i]
        # --- LL: LeadLagBypass ---
        LL_y = ifelse((params.TB[i]) == 0, ifelse((params.TC[i]) == 0, (err), vars.LL_x[i] + ((params.TC[i]) / (params.TB[i])) * ((err) - vars.LL_x[i])), vars.LL_x[i] + ((params.TC[i]) / (params.TB[i])) * ((err) - vars.LL_x[i]))
        LL_dx = ifelse((params.TB[i]) == 0, 0, ((err) - vars.LL_x[i]) / (params.TB[i]))

        VR_rate = ifelse(params.TA[i] == 0, 0, (params.KA[i] * LL_y - vars.VR_x[i]) / params.TA[i])
        VR_y = ifelse(params.TA[i] == 0, clamp(params.KA[i] * LL_y, params.VRMIN[i], params.VRMAX[i]), vars.VR_x[i])
        SAT_excess = vars.E_fd[i] - sat0_A
        SAT_over = max(0, SAT_excess)
        SAT_us = sat0_B * SAT_over ^ 2
        vfe = params.KE[i] * vars.E_fd[i] + SAT_us
        efd_rate = (VR_y - vfe) / params.TE[i]
        # --- EFD: Integrator ---
        EFD_dy = (1.0) * (efd_rate) / (1.0)

        # --- WF: Washout ---
        WF_y = (kfw) * ((vars.E_fd[i]) - vars.WF_x[i])
        WF_dx = ((vars.E_fd[i]) - vars.WF_x[i]) / (params.TF[i])

        # === Residual Assignments ===
        outputs.dLG_x[i] = mask * (LG_dx)
        outputs.dLL_x[i] = mask * (LL_dx)
        outputs.dVR_x[i] = mask * (params.VR_gate[i] * VR_rate)
        outputs.dEFD_y[i] = mask * (EFD_dy)
        outputs.dWF_x[i] = mask * (WF_dx)
        outputs.g_vref[i] = mask * (params.vref0[i] - vars.vref[i])
        outputs.g_yfb[i] = mask * (WF_y - vars.yfb[i])
        outputs.g_E_fd[i] = mask * (vars.EFD_y[i] - vars.E_fd[i])
        outputs.E_fd_out[i] = mask * (vars.E_fd[i] - params.E_fd0[i])
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:IEEEX1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    ieeex1_residual_kernel!(
        (dLG_x = ws.res_vals.dLG_x, dLL_x = ws.res_vals.dLL_x, dVR_x = ws.res_vals.dVR_x, dEFD_y = ws.res_vals.dEFD_y, dWF_x = ws.res_vals.dWF_x, g_vref = ws.res_vals.g_vref, g_yfb = ws.res_vals.g_yfb, g_E_fd = ws.res_vals.g_E_fd, E_fd_out = ws.res_vals.E_fd_out),
        (TR = model.TR, KA = model.KA, TA = model.TA, TB = model.TB, TC = model.TC, VRMAX = model.VRMAX, VRMIN = model.VRMIN, KE = model.KE, TE = model.TE, KF = model.KF, TF = model.TF, Switch = model.Switch, E1 = model.E1, SE1 = model.SE1, E2 = model.E2, SE2 = model.SE2, E_fd0 = ws.setpoints.E_fd0, vref0 = ws.setpoints.vref0, VR_gate = ws.automata.discrete.VR_gate),
        (LG_x = ws.var_vals.LG_x, LL_x = ws.var_vals.LL_x, VR_x = ws.var_vals.VR_x, EFD_y = ws.var_vals.EFD_y, WF_x = ws.var_vals.WF_x, vref = ws.var_vals.vref, yfb = ws.var_vals.yfb, E_fd = ws.var_vals.E_fd, v = ws.var_vals.v),
        active_mask,
        t,
        AT
    )
end

@inline function ieeex1_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.VR_x_minus_VRMAX)
        # --- LG: LagBypass ---
        LG_y = ifelse((params.TR[i]) == 0, (1.0) * (vars.v[i]), vars.LG_x[i])
        LG_dx = ifelse((params.TR[i]) == 0, 0, ((1.0) * (vars.v[i]) - vars.LG_x[i]) / (params.TR[i]))

        err = vars.vref[i] - LG_y - vars.yfb[i]
        # --- LL: LeadLagBypass ---
        LL_y = ifelse((params.TB[i]) == 0, ifelse((params.TC[i]) == 0, (err), vars.LL_x[i] + ((params.TC[i]) / (params.TB[i])) * ((err) - vars.LL_x[i])), vars.LL_x[i] + ((params.TC[i]) / (params.TB[i])) * ((err) - vars.LL_x[i]))
        LL_dx = ifelse((params.TB[i]) == 0, 0, ((err) - vars.LL_x[i]) / (params.TB[i]))

        VR_rate = ifelse(params.TA[i] == 0, 0, (params.KA[i] * LL_y - vars.VR_x[i]) / params.TA[i])
        # === Automaton Values ===
        values.VR_x_minus_VRMAX[i] = vars.VR_x[i] - params.VRMAX[i]
        values.VR_rate[i] = VR_rate
        values.VR_x_minus_VRMIN[i] = vars.VR_x[i] - params.VRMIN[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:IEEEX1},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    ieeex1_automaton_values_kernel!(
        (VR_x_minus_VRMAX = ws.automata.values.VR_x_minus_VRMAX, VR_rate = ws.automata.values.VR_rate, VR_x_minus_VRMIN = ws.automata.values.VR_x_minus_VRMIN),
        (TR = model.TR, KA = model.KA, TA = model.TA, TB = model.TB, TC = model.TC, VRMAX = model.VRMAX, VRMIN = model.VRMIN, KE = model.KE, TE = model.TE, KF = model.KF, TF = model.TF, Switch = model.Switch, E1 = model.E1, SE1 = model.SE1, E2 = model.E2, SE2 = model.SE2, E_fd0 = ws.setpoints.E_fd0, vref0 = ws.setpoints.vref0, VR_gate = ws.automata.discrete.VR_gate),
        (LG_x = ws.var_vals.LG_x, LL_x = ws.var_vals.LL_x, VR_x = ws.var_vals.VR_x, EFD_y = ws.var_vals.EFD_y, WF_x = ws.var_vals.WF_x, vref = ws.var_vals.vref, yfb = ws.var_vals.yfb, E_fd = ws.var_vals.E_fd, v = ws.var_vals.v),
        t,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{IEEEX1{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, IEEEX1)

    E_fd0 = T.(data.E_fd0)
    v = T.(data.v)

    # === Pre-allocate Result Arrays ===
    LG_x = Vector{T}(undef, n)
    LL_x = Vector{T}(undef, n)
    VR_x = Vector{T}(undef, n)
    EFD_y = Vector{T}(undef, n)
    WF_x = Vector{T}(undef, n)
    vref = Vector{T}(undef, n)
    yfb = fill(zero(T), n)
    E_fd = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        sat0_r = (ifelse(model.E1[i] < model.E2[i], model.SE2[i], model.SE1[i]) > 0 ? sqrt(min(model.E1[i], model.E2[i]) * ifelse(model.E1[i] < model.E2[i], model.SE1[i], model.SE2[i]) / (max(model.E1[i], model.E2[i]) * ifelse(ifelse(model.E1[i] < model.E2[i], model.SE2[i], model.SE1[i]) > 0, ifelse(model.E1[i] < model.E2[i], model.SE2[i], model.SE1[i]), 1))) : 0)
        sat0_A = (min(model.E1[i], model.E2[i]) - max(model.E1[i], model.E2[i]) * sat0_r) / (1 - sat0_r)
        sat0_B = max(model.E1[i], model.E2[i]) * ifelse(model.E1[i] < model.E2[i], model.SE2[i], model.SE1[i]) / (max(model.E1[i], model.E2[i]) - sat0_A) ^ 2
        sat0_excess = E_fd0[i] - sat0_A
        sat0_over = max(sat0_excess, 0)
        sat0 = sat0_B * sat0_over ^ 2
        vr0 = model.KE[i] * E_fd0[i] + sat0

        # State initial values
        LG_x[i] = v[i]
        LL_x[i] = vr0 / model.KA[i]
        VR_x[i] = vr0
        EFD_y[i] = E_fd0[i]
        WF_x[i] = E_fd0[i]

        # Algebraic initial values
        E_fd[i] = E_fd0[i]
        vref[i] = v[i] + vr0 / model.KA[i]
    end

    # === Return Initial Values ===
    return (LG_x=LG_x, LL_x=LL_x, VR_x=VR_x, EFD_y=EFD_y, WF_x=WF_x, vref=vref, yfb=yfb, E_fd=E_fd)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{IEEEX1}, dyr_data, sys_pf, links) -> StructArray{IEEEX1}

Parse IEEEX1 model from DYR data.
"""
function parse_dynamic_model(
    ::Type{IEEEX1},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    ::Nothing=nothing
)
    n = length(dyr_data)
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{IEEEX1{T}}(undef, n)

    for idx in 1:n
        dyr_row = dyr_data[idx]
        instances[idx] = check_parameter_domains(IEEEX1(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            TR = T(dyr_row.TR),
            KA = T(dyr_row.KA),
            TA = T(dyr_row.TA),
            TB = T(dyr_row.TB),
            TC = T(dyr_row.TC),
            VRMAX = T(dyr_row.VRMAX),
            VRMIN = T(dyr_row.VRMIN),
            KE = T(dyr_row.KE),
            TE = T(dyr_row.TE),
            KF = T(dyr_row.KF),
            TF = T(dyr_row.TF),
            Switch = T(dyr_row.Switch),
            E1 = T(dyr_row.E1),
            SE1 = T(dyr_row.SE1),
            E2 = T(dyr_row.E2),
            SE2 = T(dyr_row.SE2)
        ), "dyr record $(idx)")
    end

    return StructArray(instances)
end
