"""TGOV1: Steam Turbine Governor with Lead-Lag Compensation

DO NOT EDIT - changes will be overwritten.
"""


using InlineStrings
using PowerfulCore: Automaton, Mode, Transition, InitialMode, Condition, Atom, AllOf

export TGOV1

struct TGOV1{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number"
    BUS::Int32
    "Generator ID"
    ID::String3

    # === Parameters ===
    R::T  # Speed regulation gain (droop)
    T1::T  # Governor lag time constant
    T2::T  # Governor lead time constant
    T3::T  # Turbine lag time constant
    Dt::T  # Turbine damping coefficient
    Vmax::T  # Maximum valve position
    Vmin::T  # Minimum valve position

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function TGOV1(; BUS::Int32, ID, R, T1, T2, T3, Dt, Vmax, Vmin, status::Bool=true)
    T = typeof(R)
    check_parameter_domains(TGOV1{T}(BUS, ID, T(R), T(T1), T(T2), T(T3), T(Dt), T(Vmax), T(Vmin), status),
                            "keyword construction")
end

# Parameter domains, checked where data enters
function check_parameter_domains(x::TGOV1, from::AbstractString)
    (x.R >= 0.0) || throw(ArgumentError("TGOV1 ($from): R = $(x.R) is not >= 0"))
    (x.T1 >= 0.0) || throw(ArgumentError("TGOV1 ($from): T1 = $(x.T1) is not >= 0"))
    (x.T2 >= 0.0) || throw(ArgumentError("TGOV1 ($from): T2 = $(x.T2) is not >= 0"))
    (x.T3 >= 0.0) || throw(ArgumentError("TGOV1 ($from): T3 = $(x.T3) is not >= 0"))
    x.Vmin <= x.Vmax || throw(ArgumentError("TGOV1 ($from) does not satisfy its constraint LAG_limits_ordered: LAG: the lower limit of LAG.y is at most its upper limit"))
    (x.R != 0) || throw(ArgumentError("TGOV1 ($from): R is zero, and the model divides by it (services.gain)"))
    (x.T1 != 0) || throw(ArgumentError("TGOV1 ($from): T1 is zero, and the model divides by it (graph[2].states.y.derivative)"))
    (x.T3 != 0) || throw(ArgumentError("TGOV1 ($from): T3 is zero, and the model divides by it (graph[3].LL)"))
    return x
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:TGOV1}) = true
get_dyr_model_name(::Type{<:TGOV1}) = "TGOV1"
get_controlled_group(::Type{<:TGOV1}) = :generators


# ============================================================================
# Group Interface Implementation - :governors
# ============================================================================

# Declare that TGOV1 implements the :governors group interface
implements_group(::Type{<:TGOV1}) = :governors

# Map group canonical variable names to TGOV1's actual variable names
function group_var_mapping(::Type{<:TGOV1})
    Dict(
        :pout => :pout, # Same name
        :pref => :pref, # Same name
    )
end

# Map group canonical residual names to TGOV1's actual residual names
function group_res_mapping(::Type{<:TGOV1})
    Dict(
        :g_pref => :g_pref, # Same name
    )
end

# Controllers don't replace anything - they control generators
get_replacement_map(::Type{<:TGOV1}) = OrderedDict{Symbol, NamedTuple}()

# No data inheritance for controllers (all parameters from DYR)
function get_data_inheritance(::Type{<:TGOV1})
    OrderedDict{Symbol, OrderedDict{Symbol, Symbol}}()
end

# ============================================================================
# Model Registration
# ============================================================================

supports_format(::Type{<:TGOV1}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model TGOV1

# Terminal declarations - where TGOV1 connects
const TGOV1_TERMINALS = (
    Terminal(:generators, [:BUS, :ID]),
)

get_terminals(::Type{<:TGOV1}) = TGOV1_TERMINALS

get_key_name(::Type{<:TGOV1}) = [:BUS, :ID]

function model_metadata(::Type{TGOV1})
    ModelMetadata(
        name = :TGOV1,
        vars = [
            ModelVar(:LAG_y, StateVar();
                description="Governor valve position lag with limits"),
            ModelVar(:LL_x, StateVar();
                description="Lead-lag compensator"),
            ModelVar(:pref, AlgebVar();
                description="Power reference setpoint (control input - receives AGC contribution)"),
            ModelVar(:pout, AlgebVar();
                description="Turbine output power"),
            ModelVar(:omega, StateVar(), :generators, :omega, [:BUS, :ID];
                description="Rotor speed (external)", units="pu")
        ],
        residuals = [
        ModelResidual(:dLAG_y, StateRes(), PrivateRes();
            description="Governor valve position lag with limits"),
        ModelResidual(:dLL_x, StateRes(), PrivateRes();
            description="Lead-lag compensator"),
        ModelResidual(:g_pref, AlgebRes(), SharedRes(:add);
            description="Power reference setpoint (control input - receives AGC contribution)"),
        ModelResidual(:g_pout, AlgebRes(), PrivateRes();
            description="Turbine output power"),
        ModelResidual(:tm, AlgebRes(), SharedRes(:add), :generators, :g_P_m, [:BUS, :ID];
            description="Mechanical power contribution (delta from initial)")
    ],
        ext_params = ExtParam[],
        init_inputs = [
        InitInput(:tm0, :from_group_init, :generators, :P_m;
            description="Mechanical power from linked generator initialization"),
        InitInput(:omega, :from_group_init, :generators, :omega;
            description="Rotor speed of the linked generator at the start")
    ],
        observable_vars = Symbol[],
        initialized_params = [
        InitParam(:tm0, :pout; description="Initial mechanical power"),
        InitParam(:pref0, :pref; description="Initial power reference setpoint")
    ],
        variable_limits = [VariableLimit(:LAG_y, :Vmin, :Vmax)],
        automata = (
            Automaton(:LAG;
                modes = [Mode(:free; holds = [:LAG_gate => 1.0]),
                    Mode(:upper; holds = [:LAG_gate => 0.0]),
                    Mode(:lower; holds = [:LAG_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:LAG_y_minus_Vmax, :>=), Atom(:LAG_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:LAG_y_minus_Vmin, :<=), Atom(:LAG_rate, :<)])),
                    Transition(:upper, :free; when = Atom(:LAG_rate, :<)),
                    Transition(:lower, :free; when = Atom(:LAG_rate, :>))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:LAG_y_minus_Vmax, :>=), Atom(:LAG_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:LAG_y_minus_Vmin, :<=), Atom(:LAG_rate, :<)])),
                    InitialMode(:free)],
                at_rest = [:LAG_rate]),
        ),
        status_field = :status,
        evaluation_pass = 2,
    )
end

@inline function tgov1_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t, ::Type{AT}
) where AT
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        gain = 1.0 / params.R[i]

        wd = vars.omega[i] - 1.0
        pref_in = gain * (vars.pref[i] - wd)
        LAG_rate = (1.0 * pref_in - vars.LAG_y[i]) / params.T1[i]
        # --- LL: LeadLag ---
        LL_y = vars.LL_x[i] + ((params.T2[i]) / (params.T3[i])) * ((vars.LAG_y[i]) - vars.LL_x[i])
        LL_dx = ((vars.LAG_y[i]) - vars.LL_x[i]) / (params.T3[i])

        pout_calc = LL_y - params.Dt[i] * wd
        # === Residual Assignments ===
        outputs.dLAG_y[i] = mask * (params.LAG_gate[i] * LAG_rate)
        outputs.dLL_x[i] = mask * (LL_dx)
        outputs.g_pref[i] = mask * (params.pref0[i] - vars.pref[i])
        outputs.g_pout[i] = mask * (pout_calc - vars.pout[i])
        outputs.tm[i] = mask * (vars.pout[i] - params.tm0[i])
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:TGOV1},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    tgov1_residual_kernel!(
        (dLAG_y = ws.res_vals.dLAG_y, dLL_x = ws.res_vals.dLL_x, g_pref = ws.res_vals.g_pref, g_pout = ws.res_vals.g_pout, tm = ws.res_vals.tm),
        (R = model.R, T1 = model.T1, T2 = model.T2, T3 = model.T3, Dt = model.Dt, Vmax = model.Vmax, Vmin = model.Vmin, tm0 = ws.setpoints.tm0, pref0 = ws.setpoints.pref0, LAG_gate = ws.automata.discrete.LAG_gate),
        (LAG_y = ws.var_vals.LAG_y, LL_x = ws.var_vals.LL_x, pref = ws.var_vals.pref, pout = ws.var_vals.pout, omega = ws.var_vals.omega),
        active_mask,
        t,
        AT
    )
end

@inline function tgov1_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.LAG_y_minus_Vmax)
        gain = 1.0 / params.R[i]
        wd = vars.omega[i] - 1.0
        pref_in = gain * (vars.pref[i] - wd)
        LAG_rate = (1.0 * pref_in - vars.LAG_y[i]) / params.T1[i]
        # === Automaton Values ===
        values.LAG_y_minus_Vmax[i] = vars.LAG_y[i] - params.Vmax[i]
        values.LAG_rate[i] = LAG_rate
        values.LAG_y_minus_Vmin[i] = vars.LAG_y[i] - params.Vmin[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:TGOV1},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    tgov1_automaton_values_kernel!(
        (LAG_y_minus_Vmax = ws.automata.values.LAG_y_minus_Vmax, LAG_rate = ws.automata.values.LAG_rate, LAG_y_minus_Vmin = ws.automata.values.LAG_y_minus_Vmin),
        (R = model.R, T1 = model.T1, T2 = model.T2, T3 = model.T3, Dt = model.Dt, Vmax = model.Vmax, Vmin = model.Vmin, tm0 = ws.setpoints.tm0, pref0 = ws.setpoints.pref0, LAG_gate = ws.automata.discrete.LAG_gate),
        (LAG_y = ws.var_vals.LAG_y, LL_x = ws.var_vals.LL_x, pref = ws.var_vals.pref, pout = ws.var_vals.pout, omega = ws.var_vals.omega),
        t,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{TGOV1{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, TGOV1)

    tm0 = T.(data.tm0)
    omega = T.(data.omega)

    # === Pre-allocate Result Arrays ===
    LAG_y = Vector{T}(undef, n)
    LL_x = Vector{T}(undef, n)
    pref = Vector{T}(undef, n)
    pout = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # State initial values
        LAG_y[i] = tm0[i] + model.Dt[i] * (omega[i] - 1.0)
        LL_x[i] = tm0[i] + model.Dt[i] * (omega[i] - 1.0)

        # Algebraic initial values
        pref[i] = omega[i] + model.R[i] * (tm0[i] + model.Dt[i] * (omega[i] - 1.0)) - 1.0
        pout[i] = tm0[i]
    end

    # === Return Initial Values ===
    return (LAG_y=LAG_y, LL_x=LL_x, pref=pref, pout=pout)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{TGOV1}, dyr_data, sys_pf, links) -> StructArray{TGOV1}

Parse TGOV1 model from DYR data.
"""
function parse_dynamic_model(
    ::Type{TGOV1},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    ::Nothing=nothing
)
    n = length(dyr_data)
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{TGOV1{T}}(undef, n)

    for idx in 1:n
        dyr_row = dyr_data[idx]
        instances[idx] = check_parameter_domains(TGOV1(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            R = T(dyr_row.R),
            T1 = T(dyr_row.T1),
            T2 = T(dyr_row.T2),
            T3 = T(dyr_row.T3),
            Dt = T(dyr_row.Dt),
            Vmax = T(dyr_row.Vmax),
            Vmin = T(dyr_row.Vmin)
        ), "dyr record $(idx)")
    end

    return StructArray(instances)
end
