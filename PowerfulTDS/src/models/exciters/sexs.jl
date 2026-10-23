"""SEXS: Simplified Excitation System

DO NOT EDIT - changes will be overwritten.
"""


using InlineStrings
using PowerfulCore: Automaton, Mode, Transition, InitialMode, Condition, Atom, AllOf

export SEXS

struct SEXS{T} <: AbstractModel
    # ===== Indexing Fields =====
    "Bus number"
    BUS::Int32
    "Generator ID"
    ID::String3

    # === Parameters ===
    TATB::T  # Time constant ratio TA/TB (lead/lag)
    TB::T  # Lag time constant of LeadLag block
    K::T  # Exciter gain
    TE::T  # Exciter field time constant
    EMIN::T  # Minimum excitation output
    EMAX::T  # Maximum excitation output

    # === Status ===
    "Model status: true=online, false=offline"
    status::Bool
end

# Keyword constructor with type promotion
function SEXS(; BUS::Int32, ID, TATB, TB, K, TE, EMIN, EMAX, status::Bool=true)
    T = typeof(TATB)
    check_parameter_domains(SEXS{T}(BUS, ID, T(TATB), T(TB), T(K), T(TE), T(EMIN), T(EMAX), status),
                            "keyword construction")
end

# Parameter domains, checked where data enters
function check_parameter_domains(x::SEXS, from::AbstractString)
    x.EMIN <= x.EMAX || throw(ArgumentError("SEXS ($from) does not satisfy its constraint LAW_limits_ordered: LAW: the lower limit of LAW.y is at most its upper limit"))
    (x.TE != 0) || throw(ArgumentError("SEXS ($from): TE is zero, and the model divides by it (graph[2].states.y.derivative)"))
    (x.K != 0) || throw(ArgumentError("SEXS ($from): K is zero, and the model divides by it (algebraics.vref.init)"))
    (x.TB != 0) || throw(ArgumentError("SEXS ($from): TB is zero, and the model divides by it (graph[1].LL)"))
    return x
end

# ============================================================================
# Trait Implementations
# ============================================================================

requires_initialization(::Type{<:SEXS}) = true
get_dyr_model_name(::Type{<:SEXS}) = "SEXS"
get_controlled_group(::Type{<:SEXS}) = :generators


# ============================================================================
# Group Interface Implementation - :exciters
# ============================================================================

# Declare that SEXS implements the :exciters group interface
implements_group(::Type{<:SEXS}) = :exciters

# Map group canonical variable names to SEXS's actual variable names
function group_var_mapping(::Type{<:SEXS})
    Dict(
        :E_fd => :E_fd, # Same name
        :vref => :vref, # Same name
    )
end

# Map group canonical residual names to SEXS's actual residual names
function group_res_mapping(::Type{<:SEXS})
    Dict(
        :g_vref => :g_vref, # Same name
    )
end

# Controllers don't replace anything - they control generators
get_replacement_map(::Type{<:SEXS}) = OrderedDict{Symbol, NamedTuple}()

# No data inheritance for controllers (all parameters from DYR)
function get_data_inheritance(::Type{<:SEXS})
    OrderedDict{Symbol, OrderedDict{Symbol, Symbol}}()
end

# ============================================================================
# Model Registration
# ============================================================================

supports_format(::Type{<:SEXS}, ::Type{PSSEDYR}) = FormatSupport{PSSEDYR}()

@register_model SEXS

# Terminal declarations - where SEXS connects
const SEXS_TERMINALS = (
    Terminal(:bus_group, :BUS),
    Terminal(:generators, [:BUS, :ID]),
)

get_terminals(::Type{<:SEXS}) = SEXS_TERMINALS

# External parameter declarations
const SEXS_EXT_PARAMS = [
    ExtParam(:bus_basekv, SEXS_TERMINALS[1], :basekv,
        BaseVoltage(), PrePerUnit(),
        required=false, description="Base voltage at exciter bus"),
]

get_key_name(::Type{<:SEXS}) = [:BUS, :ID]

function model_metadata(::Type{SEXS})
    ModelMetadata(
        name = :SEXS,
        vars = [
            ModelVar(:LL_x, StateVar();
                description="Lead-lag compensator"),
            ModelVar(:LAW_y, StateVar();
                description="Exciter field with gain, lag, and anti-windup limits"),
            ModelVar(:vref, AlgebVar();
                description="Reference voltage setpoint (control input)"),
            ModelVar(:E_fd, AlgebVar();
                description="Field voltage output"),
            ModelVar(:v, AlgebVar(), :bus_group, :v, :BUS;
                description="Bus voltage magnitude (external)", units="pu")
        ],
        residuals = [
        ModelResidual(:dLL_x, StateRes(), PrivateRes();
            description="Lead-lag compensator"),
        ModelResidual(:dLAW_y, StateRes(), PrivateRes();
            description="Exciter field with gain, lag, and anti-windup limits"),
        ModelResidual(:g_vref, AlgebRes(), SharedRes(:add);
            description="Reference voltage setpoint (control input)"),
        ModelResidual(:g_E_fd, AlgebRes(), PrivateRes();
            description="Field voltage output"),
        ModelResidual(:E_fd_out, AlgebRes(), SharedRes(:add), :generators, :g_E_fd, [:BUS, :ID];
            description="Field voltage contribution (delta from initial)")
    ],
        ext_params = SEXS_EXT_PARAMS,
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
        variable_limits = [VariableLimit(:LAW_y, :EMIN, :EMAX)],
        automata = (
            Automaton(:LAW;
                modes = [Mode(:free; holds = [:LAW_gate => 1.0]),
                    Mode(:upper; holds = [:LAW_gate => 0.0]),
                    Mode(:lower; holds = [:LAW_gate => 0.0])],
                transitions = [Transition(:free, :upper; when = AllOf(Condition[Atom(:LAW_y_minus_EMAX, :>=), Atom(:LAW_rate, :>)])),
                    Transition(:free, :lower; when = AllOf(Condition[Atom(:LAW_y_minus_EMIN, :<=), Atom(:LAW_rate, :<)])),
                    Transition(:upper, :free; when = Atom(:LAW_rate, :<)),
                    Transition(:lower, :free; when = Atom(:LAW_rate, :>))],
                initial = [InitialMode(:upper; when = AllOf(Condition[Atom(:LAW_y_minus_EMAX, :>=), Atom(:LAW_rate, :>)])),
                    InitialMode(:lower; when = AllOf(Condition[Atom(:LAW_y_minus_EMIN, :<=), Atom(:LAW_rate, :<)])),
                    InitialMode(:free)],
                at_rest = [:LAW_rate]),
        ),
        status_field = :status,
        evaluation_pass = 2,
    )
end

@inline function sexs_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t, ::Type{AT}
) where AT
    @inbounds @simd for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs[1])) : zero(eltype(outputs[1]))

        TA = params.TATB[i] * params.TB[i]

        vi = vars.vref[i] - vars.v[i]
        # --- LL: LeadLag ---
        LL_y = vars.LL_x[i] + ((TA) / (params.TB[i])) * ((vi) - vars.LL_x[i])
        LL_dx = ((vi) - vars.LL_x[i]) / (params.TB[i])

        LAW_rate = (params.K[i] * LL_y - vars.LAW_y[i]) / params.TE[i]
        E_fd_calc = vars.LAW_y[i]
        # === Residual Assignments ===
        outputs.dLL_x[i] = mask * (LL_dx)
        outputs.dLAW_y[i] = mask * (params.LAW_gate[i] * LAW_rate)
        outputs.g_vref[i] = mask * (params.vref0[i] - vars.vref[i])
        outputs.g_E_fd[i] = mask * (E_fd_calc - vars.E_fd[i])
        outputs.E_fd_out[i] = mask * (vars.E_fd[i] - params.E_fd0[i])
    end
    nothing
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:SEXS},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    sexs_residual_kernel!(
        (dLL_x = ws.res_vals.dLL_x, dLAW_y = ws.res_vals.dLAW_y, g_vref = ws.res_vals.g_vref, g_E_fd = ws.res_vals.g_E_fd, E_fd_out = ws.res_vals.E_fd_out),
        (TATB = model.TATB, TB = model.TB, K = model.K, TE = model.TE, EMIN = model.EMIN, EMAX = model.EMAX, E_fd0 = ws.setpoints.E_fd0, vref0 = ws.setpoints.vref0, LAW_gate = ws.automata.discrete.LAW_gate),
        (LL_x = ws.var_vals.LL_x, LAW_y = ws.var_vals.LAW_y, vref = ws.var_vals.vref, E_fd = ws.var_vals.E_fd, v = ws.var_vals.v),
        active_mask,
        t,
        AT
    )
end

@inline function sexs_automaton_values_kernel!(
    values::NamedTuple, params::NamedTuple, vars::NamedTuple,
    t, ::Type{AT}
) where AT
    @inbounds for i in eachindex(values.LAW_y_minus_EMAX)
        TA = params.TATB[i] * params.TB[i]
        vi = vars.vref[i] - vars.v[i]
        # --- LL: LeadLag ---
        LL_y = vars.LL_x[i] + ((TA) / (params.TB[i])) * ((vi) - vars.LL_x[i])
        LL_dx = ((vi) - vars.LL_x[i]) / (params.TB[i])

        LAW_rate = (params.K[i] * LL_y - vars.LAW_y[i]) / params.TE[i]
        # === Automaton Values ===
        values.LAW_y_minus_EMAX[i] = vars.LAW_y[i] - params.EMAX[i]
        values.LAW_rate[i] = LAW_rate
        values.LAW_y_minus_EMIN[i] = vars.LAW_y[i] - params.EMIN[i]
    end
    nothing
end

@inline function calc_automaton_values!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:SEXS},
    external::NamedTuple, ::PowerInjection, t
) where {T, AT}
    sexs_automaton_values_kernel!(
        (LAW_y_minus_EMAX = ws.automata.values.LAW_y_minus_EMAX, LAW_rate = ws.automata.values.LAW_rate, LAW_y_minus_EMIN = ws.automata.values.LAW_y_minus_EMIN),
        (TATB = model.TATB, TB = model.TB, K = model.K, TE = model.TE, EMIN = model.EMIN, EMAX = model.EMAX, E_fd0 = ws.setpoints.E_fd0, vref0 = ws.setpoints.vref0, LAW_gate = ws.automata.discrete.LAW_gate),
        (LL_x = ws.var_vals.LL_x, LAW_y = ws.var_vals.LAW_y, vref = ws.var_vals.vref, E_fd = ws.var_vals.E_fd, v = ws.var_vals.v),
        t,
        AT
    )
end

function initialize_dynamic_states(
    model::StructArray{SEXS{T}},
    ctx::InitContext
) where T
    n = length(model)

    # === Gather Init Inputs ===
    data = gather_init_data(ctx, SEXS)

    E_fd0 = T.(data.E_fd0)
    v = T.(data.v)

    # === Pre-allocate Result Arrays ===
    LL_x = Vector{T}(undef, n)
    LAW_y = Vector{T}(undef, n)
    vref = Vector{T}(undef, n)
    E_fd = Vector{T}(undef, n)

    # === Compute Initial Values ===
    @inbounds for i in 1:n
        # State initial values
        LL_x[i] = E_fd0[i] / model.K[i]
        LAW_y[i] = E_fd0[i]

        # Algebraic initial values
        E_fd[i] = E_fd0[i]
        vref[i] = v[i] + E_fd0[i] / model.K[i]
    end

    # === Return Initial Values ===
    return (LL_x=LL_x, LAW_y=LAW_y, vref=vref, E_fd=E_fd)
end

# ============================================================================
# DYR Parser - Trait Implementation
# ============================================================================

"""
    parse_dynamic_model(::Type{SEXS}, dyr_data, sys_pf, links) -> StructArray{SEXS}

Parse SEXS model from DYR data.
"""
function parse_dynamic_model(
    ::Type{SEXS},
    dyr_data::StructVector,
    sys_pf::SystemModel,
    ::Nothing=nothing
)
    n = length(dyr_data)
    T = eltype(sys_pf.models.Bus.vm)  # Infer float type from system

    instances = Vector{SEXS{T}}(undef, n)

    for idx in 1:n
        dyr_row = dyr_data[idx]
        instances[idx] = check_parameter_domains(SEXS(
            BUS = Int32(dyr_row.BUS),
            ID = InlineString3(dyr_row.ID),
            TATB = T(dyr_row.TATB),
            TB = T(dyr_row.TB),
            K = T(dyr_row.K),
            TE = T(dyr_row.TE),
            EMIN = T(dyr_row.EMIN),
            EMAX = T(dyr_row.EMAX)
        ), "dyr record $(idx)")
    end

    return StructArray(instances)
end
