"""
TDS Adapter - Evaluation machinery for time-domain simulation

This module provides the TDSAdapter type and evaluation methods that bridge
DiffEq.jl's interface with Powerful's internal residual evaluation.
"""

using SparseArrays: SparseMatrixCSC

# ============================================================================
# Unified TDS Adapter Architecture
# ============================================================================

# ----------------------------------------------------------------------------
# Interface Types (Strategy Pattern)
# ----------------------------------------------------------------------------

"""
    TDSInterface

Abstract type for TDS solver interface.
"""
abstract type TDSInterface end

"""
    ODEInterface <: TDSInterface

Evaluate the TDS system as a mass-matrix ODE `M du/dt = F(u, t)`, the form passed to
OrdinaryDiffEq. Selects the [`evaluate!`](@ref) method used by the ODE right-hand side.
"""
struct ODEInterface <: TDSInterface end

# ----------------------------------------------------------------------------
# Reinit Cache (Lazy-Initialized Jacobian Infrastructure)
# ----------------------------------------------------------------------------

"""
    ReinitCache

Cached infrastructure for algebraic reinitialization Newton solver.

Created lazily on first `reinitialize_algebraic!` call and reused for all
subsequent events. This saves ~50-200ms per event by avoiding recreation
of the closure, Jacobian prep, and sparsity pattern analysis.

# Fields
- `closure`: `g_only!(g_out, y_in)` closure using adapter.x_buf
- `jac_prep`: DifferentiationInterface prep for sparse Jacobian
- `jac_matrix`: Sparse Jacobian matrix (reused for pattern)
- `jac_backend`: AutoSparse(AutoForwardDiff()) or similar backend
- `linear_solver`: KLU factorization cache (numeric refactorization after the first solve,
  a full factorization again after `algebraic_equations_changed!`)
- `newton_bufs`: Preallocated Newton trial/step buffers
- `t`: The time the closure evaluates the residuals at, set by each reinitialization to the
  time of its state (`_reinitialize_algebraic`), as `JacobianClosureFunctor` holds `t_ref`
- `rows`: The bus power balance rows the solve measures as current mismatch (`BusPowerRows`)
- `check_buf`, `check_rows`: The DAE residual at the reinitialized state, and its rows without
  mass as the solve measures them, for the check of its result (`_reinit_mismatch`)

All fields use `Any` type to avoid closure type explosion that would cause
50+ second compilation times (see solve.jl comments for details).
"""
mutable struct ReinitCache
    closure::Any       # residual closure on the unknowns, using adapter.x_buf
    jac_prep::Any      # DifferentiationInterface prep for closure
    jac_matrix::Union{Nothing, SparseMatrixCSC{Float64, Int}}
    jac_backend::Any   # AutoSparse(AutoForwardDiff()) backend
    linear_solver::Any # KLULinearSolver, factorization reused across events whose equations are unchanged
    newton_bufs::Any   # NewtonBuffers (trial residual/point, step, Jacobian)
    massless::Vector{Int}       # States of x without mass, solved with y
    n_unknowns::Int             # length(massless) + n_alg
    unknowns::Vector{Float64}   # The unknowns [x without mass; y] (unused without such states)
    voltages::Vector{Int}       # The bus voltage magnitudes among the unknowns (kept positive)
    t::Base.RefValue{Float64}   # The time of the state being reinitialized, read by `closure`
    rows::Any                   # BusPowerRows or nothing (PowerInjection bus rows as current mismatch)
    check_buf::Vector{Float64}  # The DAE residual [f; g] at the result, for its check
    check_rows::Vector{Float64} # Its rows without mass, as the solve measures them
end

# ----------------------------------------------------------------------------
# Index Reduction of State-Only Constraints
# ----------------------------------------------------------------------------

"""
    ConstraintGate(modes, instance, mask, label)

When a replaced row of a `StateConstraintReduction` holds: the conjunction of terms, each the mode of an
automaton of one instance (`modes[t][instance[t]]`, the automaton's mode array) being one of those admitted
(`mask[t]` has bit `m` set for mode `m`). A gate without terms is always open. The rows of a mode-dependent
constraint are replaced while their modes last and are the original equations otherwise.
"""
struct ConstraintGate
    modes::Vector{Vector{Int}}
    instance::Vector{Int}
    mask::Vector{UInt64}
    label::String
end

"The gate of a constraint that holds in every mode"
const OPEN_GATE = ConstraintGate(Vector{Int}[], Int[], UInt64[], "")

@inline function gate_open(g::ConstraintGate)
    @inbounds for t in eachindex(g.instance)
        (g.mask[t] >>> g.modes[t][g.instance[t]]) & 0x1 == 0x1 || return false
    end
    return true
end

"""
    UnsupportedConstraint(gate, reason)

A state-only constraint found in some mode (`gate`) that the reduction does not make index 1 (`reason`
says why); the reduction does not replace its rows, and an integration that enters the mode fails with the
singular algebraic block (`SingularRatesError`), whose message names these.
"""
struct UnsupportedConstraint
    gate::ConstraintGate
    reason::String
end

"""
    StateConstraintReduction

Rows without mass of `M du/dt = F(u)` whose combination constrains differential states only,
`c(t)·x + b = 0`, each replaced by its time derivative (see [`install_state_reduction!`](@ref) and
[`reduction_algorithm`](@ref)).

A row without mass that depends on no variable without mass, e.g. the current balance of a
bus whose every connection is an inductor current, makes the algebraic Jacobian
singular: the DAE has index 2. So does a combination of rows: the balances of the buses of an inductor-only
cut-set joined to each other by resistors, whose sum is free of their voltages, or a clamp `Volim - v = 0`
beside a measurement `v - vm = 0`, whose sum `Volim - vm` is. The derivative `Σ_k (c_k(t)·f_k/M_k + ċ_k(t)·x_k)`
depends on the variables without mass through the state equations `f`, which restores index 1. The
coefficients are constant, or periodic at the system frequency, `c(t) = c₀ + c_c cos ωt + c_s sin ωt` (a
current injected in the synchronous frame at ωt).

# Fields
- `rows`: The replaced rows (indices into `u`, in the state space or among the algebraic rows): the pivot
  row of each combination
- `rowptr`, `colval`: compressed rows, one per replaced row, over the states its coefficients have
- `coef`, `coef_c`, `coef_s`: `c₀/M`, `c_c/M`, `c_s/M` at those entries (the replaced row evaluates to
  `Σ (coef + coef_c cos ωt + coef_s sin ωt)·F[colval] + ω Σ (-a_c sin ωt + a_s cos ωt)·x[colval]`,
  and for a periodic row also `+ ω Σ (a_0 + a_c cos ωt + a_s sin ωt)·x[colval]`, the stabilization)
- `a_0`, `a_c`, `a_s`: `c₀`, `c_c`, `c_s` at those entries (for `ċ(t)·x`, and the constraint itself,
  `constraint_violation`)
- `omega`: ω, the system frequency
- `periodic`: whether any coefficient varies with time
- `A0`, `Ac`, `As`: `c₀`, `c_c`, `c_s` as matrices (`length(rows) × n_diff`), the coefficients with
  every instance in service; `constraint_matrix(red, t)` is `c(t)`, used to restore the constraints
  after events
- `pattern`: Sparsity of the reduced DAE Jacobian (each replaced row takes the patterns of the
  state equations it sums, and those states, and of the rows of its combination)
- `owner_model`, `owner_index`: Model and instance owning each differential state
  (`:none`, 0 for states of models without a status), to hold the states of
  out-of-service instances when restoring the constraints
- `support_ptr`, `support_row`, `support_weight`: the combination of each replaced row, its rows of the
  unreduced residual and their weights (the row alone with weight one for a row that no variable enters);
  their weighted sum is the value of the constraint `c(t)·x + b`, read from the residual as evaluated
- `live`, `beta`: a row of a combination (`live`) adds the Baumgarte term `β · (c·x + b)` with that value
  read from the residual, so that it follows the parameters an event changes (`Volim`); a row alone
  (not `live`) adds the periodic stabilization of its coefficients above, and nothing if constant
- `gates`: when each row holds (`ConstraintGate`): always, or in the modes its combination exists in
- `reads`, `descriptions`: for every row, the model parameters its constraint reads (`(model, field, instance)`, found by
  tracing the kernels of the equations it replaces or combines: the shunt conductance `g` of a bus whose balance is
  reduced, REGCA1's `Volim` in the clamp), which an event cannot change (`check_constraint_modifiable`), and its
  description for messages
- `unsupported`: the constraints found in some mode that the reduction does not make index 1
"""
struct StateConstraintReduction
    rows::Vector{Int}
    rowptr::Vector{Int}
    colval::Vector{Int}
    coef::Vector{Float64}
    coef_c::Vector{Float64}
    coef_s::Vector{Float64}
    a_0::Vector{Float64}
    a_c::Vector{Float64}
    a_s::Vector{Float64}
    omega::Float64
    periodic::Bool
    A0::SparseMatrixCSC{Float64, Int}
    Ac::SparseMatrixCSC{Float64, Int}
    As::SparseMatrixCSC{Float64, Int}
    pattern::SparseMatrixCSC{Bool, Int}
    owner_model::Vector{Symbol}
    owner_index::Vector{Int}
    support_ptr::Vector{Int}
    support_row::Vector{Int}
    support_weight::Vector{Float64}
    live::Vector{Bool}
    beta::Vector{Float64}
    gates::Vector{ConstraintGate}
    reads::Vector{Vector{Tuple{Symbol, Symbol, Int}}}
    descriptions::Vector{String}
    unsupported::Vector{UnsupportedConstraint}
end

"""
    constraint_matrix(red::StateConstraintReduction, t) -> SparseMatrixCSC

The coefficients `c(t)` of the reduced constraints at time `t`.
"""
constraint_matrix(red::StateConstraintReduction, t::Real) =
    red.periodic ? red.A0 + red.Ac * cos(red.omega * t) + red.As * sin(red.omega * t) : red.A0

"""
    constraint_violation(red::StateConstraintReduction, u, t) -> Float64

The largest `|c(t)·x|` over the reduced constraints that are a row alone (homogeneous: `b = 0`) at the state `u` (its
states `x`) and time `t`. The value of a combination (`c·x + b`) is read from the residual, `constraint_values`.
"""
function constraint_violation(red::StateConstraintReduction, u::AbstractVector, t::Real)
    s, c = sincos(red.omega * t)
    worst = 0.0
    @inbounds for j in eachindex(red.rows)
        red.live[j] && continue
        acc = 0.0
        for p in red.rowptr[j]:(red.rowptr[j + 1] - 1)
            acc += (red.a_0[p] + red.a_c[p] * c + red.a_s[p] * s) * u[red.colval[p]]
        end
        worst = max(worst, abs(acc))
    end
    return worst
end

"""
    reduction_algorithm(solver_config) -> Union{Nothing, Function}

The state-reduction algorithm `solver_config` registers for detecting rows that constrain
differential states only (`StateConstraintReduction`), as `(adapter, prob) -> Union{Nothing,
StateConstraintReduction}`, or `nothing` if it registers none. Plain TDS registers none for
any of its own solver configs: a model with a zero-mass differential state and no registered
algorithm is left to the ordinary singular-Jacobian error. A package whose models can have
zero-mass states (e.g. an ideal inductor/capacitor, or an algebraic regulation of unbounded
bandwidth) adds a method of this function for its own solver config type.
"""
reduction_algorithm(::Any) = nothing

# The value of the constraint of the replaced row `j` from the unreduced residual `du`: the weighted sum of the rows of its
# combination
@inline function _constraint_value(red::StateConstraintReduction, du::AbstractVector{T}, j::Int) where T
    value = zero(T)
    @inbounds for q in red.support_ptr[j]:(red.support_ptr[j + 1] - 1)
        value += red.support_weight[q] * du[red.support_row[q]]
    end
    return value
end

# Replace each reduced row of `du` (evaluated [f; g] at the state `u`, time `t`) by the time derivative
# of its constraint, where its gate is open (the modes its combination exists in; the original equation
# otherwise). The summed rows are differential states (nonzero mass), never replaced rows, and the rows of
# a combination are not the replaced rows of another. A periodic row is replaced by ċ·x + c·ẋ + ω c·x
# (Baumgarte stabilization with rate ω, the system frequency): the same equation on the constraint c·x = 0,
# where the integration leaves it; but with coefficients that rotate, the Newton iterations on a Jacobian kept
# over many steps leave a residual in the row, which integrates into a drift of c·x, and the last term makes
# that drift decay at ω instead (round B6). A constant row needs none: its Jacobian does not age. A row of a
# combination (`live`) is replaced by ċ·x + c·ẋ + β (c·x + b), the constraint's value read from the residual
# it replaces (round 14a), whose constant b the parameters of the models set.
@inline function _apply_reduction!(du::AbstractVector{T}, red::StateConstraintReduction, u::AbstractVector, t) where T
    if red.periodic
        s, c = sincos(red.omega * t)
        w = red.omega
        @inbounds for j in eachindex(red.rows)
            gate_open(red.gates[j]) || continue
            acc = zero(T)
            if red.live[j]
                for p in red.rowptr[j]:(red.rowptr[j + 1] - 1)
                    k = red.colval[p]
                    acc += (red.coef[p] + red.coef_c[p] * c + red.coef_s[p] * s) * du[k] +
                           (w * (red.a_s[p] * c - red.a_c[p] * s)) * u[k]
                end
                acc += red.beta[j] * _constraint_value(red, du, j)
            else
                for p in red.rowptr[j]:(red.rowptr[j + 1] - 1)
                    k = red.colval[p]
                    acc += (red.coef[p] + red.coef_c[p] * c + red.coef_s[p] * s) * du[k] +
                           (w * (red.a_s[p] * c - red.a_c[p] * s) + w * (red.a_0[p] + red.a_c[p] * c + red.a_s[p] * s)) * u[k]
                end
            end
            du[red.rows[j]] = acc
        end
    else
        @inbounds for j in eachindex(red.rows)
            gate_open(red.gates[j]) || continue
            acc = zero(T)
            for p in red.rowptr[j]:(red.rowptr[j + 1] - 1)
                acc += red.coef[p] * du[red.colval[p]]
            end
            red.live[j] && (acc += red.beta[j] * _constraint_value(red, du, j))
            du[red.rows[j]] = acc
        end
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Unified TDS Adapter
# ----------------------------------------------------------------------------

"""
    TDSAdapter

Adapter bridging DiffEq.jl's combined state vector with internal buffers.

# Mass Matrix DAE Formulation
The TDS system is a DAE in mass matrix form:

    M * du/dt = [f(x,y); g(x,y)]

where:
- `u = [x; y]` - combined state vector (differential `x`, algebraic `y`)
- `M = diag(M_x, 0)` - singular mass matrix (zeros for algebraic equations)
- `f(x,y)` - differential residuals (dx/dt = f)
- `g(x,y) = 0` - algebraic constraints (network equations)

# Example
```julia
adapter = TDSAdapter(adcache, prob_tds)
ode_fn! = (du, u, p, t) -> evaluate!(du, adapter, u, t, ODEInterface())
ode_func = ODEFunction(ode_fn!; mass_matrix=adapter.mass_matrix)
sol = solve(ODEProblem(ode_func, u0, tspan), QNDF(linsolve=KLUFactorization()))
```
"""
mutable struct TDSAdapter{CACHE, WS, MT<:AbstractMatrix, REC}
    # Universal evaluator (shared with power flow!)
    adcache::CACHE

    # Direct reference to Float64 workspace for zero-allocation hot path
    # Note: This avoids accessing adcache.float_ws which allocates 16 bytes
    # due to ::Any fields in ADCache causing type instability
    sys_ws::WS

    # Problem dimensions
    n_diff::Int
    n_alg::Int

    # Mass matrix (problem metadata, not computation state)
    mass_matrix::MT

    # Unpacking/packing buffers (Float64 only; Dual buffers cached in ADCache)
    x_buf::Vector{Float64}
    y_buf::Vector{Float64}
    f_buf::Vector{Float64}
    g_buf::Vector{Float64}

    # Cached reinit infrastructure (lazy initialized on first reinitialize_algebraic! call)
    # See ReinitCache docstring for details on caching strategy
    reinit_cache::Union{Nothing, ReinitCache}

    # Optional injection recorder for tracking model contributions (typed for zero-alloc)
    recorder::REC  # Nothing or InjectionRecorder{T, K}

    # Status dict for each model (passed to evaluate_system_residuals!)
    # Built from system model, enables parallel ACCC with thread-local status copies
    status::Dict{Symbol, Vector{Bool}}

    # Algebraic addresses of the bus voltage magnitudes, for reinitialization guesses
    # (see `bus_voltage_addresses`)
    bus_v_addrs::Vector{UInt}

    # Zero-mass rows replaced by their time derivative (nothing when the DAE has none)
    reduction::Union{Nothing, StateConstraintReduction}

    # Rates of change of the differential states at the initial time (nonzero for
    # oscillating states of a sinusoidal steady state, zero otherwise)
    rate_buf::Vector{Float64}

    # Automata of the models (nothing when no model declares any), started at the
    # initial state; the simulation fires their transitions (see automata.jl)
    automata::Union{Nothing, AutomatonRuntime}
end

# Detect the state-only constraints via whichever algorithm, if any, a solver config of `prob.config`
# registers (`reduction_algorithm`); a DAE without zero-mass differential states costs nothing
function install_state_reduction!(adapter::TDSAdapter, prob::TDSProblem)
    mass = diag(adapter.mass_matrix)
    any(r -> iszero(mass[r]), 1:adapter.n_diff) || return nothing
    for solver_config in prob.config.solver
        algorithm = reduction_algorithm(solver_config)
        algorithm === nothing || return algorithm(adapter, prob)
    end
    return nothing
end

"""
    TDSAdapter(adcache::AbstractADCache, prob::TDSProblem; recorder=nothing)

Construct a TDS adapter from ADCache and problem specification.

Initializes buffers, populates setpoints, and builds mass matrix.
Jacobian infrastructure is lazily initialized in TDSAdapter when needed.
Services are computed inline in calc_residuals! to enable AD through parameters.

# Arguments
- `recorder`: Optional InjectionRecorder for tracking model contributions (default: nothing)
"""
function TDSAdapter(adcache::AbstractADCache, prob::TDSProblem; recorder::REC=nothing) where REC
    # Populate setpoints in ADCache workspaces (must be done before evaluation)
    populate_setpoints!(adcache, prob.init.setpoints)

    # Extract Float64 workspace reference for zero-allocation hot path
    # This avoids accessing adcache.float_ws in evaluate! which would allocate
    # due to ::Any fields in ADCache causing type instability
    sys_ws = get_float_ws(adcache)

    # Compute dimensions from AddressManager (convert to Int for struct)
    n_diff = Int(max_address_in_space(prob.address_manager, StateVar()))
    n_alg = Int(max_address_in_space(prob.address_manager, AlgebVar()))

    # Allocate Float64 buffers
    x_buf = zeros(Float64, n_diff)
    y_buf = zeros(Float64, n_alg)
    f_buf = zeros(Float64, n_diff)
    g_buf = zeros(Float64, n_alg)

    # Initialize buffers from initial conditions; out-of-service instances hold their
    # own algebraic variables and zero-mass states at the initial values
    init_tds_buffers!(x_buf, y_buf, prob)
    set_hold!(adcache, x_buf, y_buf)
    rate_buf = init_tds_rates(prob, n_diff)

    # Build mass matrix
    mass_matrix = build_tds_mass_matrix(prob, x_buf, y_buf)

    # Build status dict from system model for evaluate_system_residuals!
    status = build_status_dict(prob.sys)

    # Automata started at the initial state, so that every evaluation (the validation of the initialization included)
    # reads the discrete quantities they start with. The detection of the state-only constraints below reads them, and
    # the algebraic point is solved between the rounds of the start only once the constraints are known (the system
    # with them removed is the one the integration solves), so the start settles in rounds after the detection
    sync_memory! = let adcache = adcache
        () -> sync_automaton_memory!(adcache)
    end
    automata = automaton_runtime(sys_ws, prob.sys, status, n_diff, n_alg, sync_memory!, Vector(diag(mass_matrix));
                                 traces = prob.traces)

    # Note: Jacobian infrastructure is lazily initialized in TDSAdapter
    # when reinitialize_algebraic! is first called

    fired = automata === nothing ? false : start!(automata, x_buf, y_buf, prob.tspan[1])
    adapter = TDSAdapter(
        adcache,
        sys_ws,
        n_diff, n_alg,
        mass_matrix,
        x_buf, y_buf, f_buf, g_buf,
        nothing,  # reinit_cache (lazy initialized on first reinitialize_algebraic! call)
        recorder,  # typed via REC parameter for zero-allocation hot path
        status,
        bus_voltage_addresses(prob.sys, prob.address_manager),
        nothing,  # reduction, detected below from the unreduced residuals
        rate_buf,
        automata
    )
    adapter.reduction = install_state_reduction!(adapter, prob)
    check_constraint_modifiable(adapter, prob)
    automata === nothing || _start_automata!(adapter, automata, prob.tspan[1],
                                             get_solver_config(prob.config, :TDS, TDSConfig).abstol, fired)
    return adapter
end

# ----------------------------------------------------------------------------
# ODE Evaluation (Mass Matrix Formulation)
# ----------------------------------------------------------------------------

"""
    evaluate!(du_out, adapter, u, t, ::ODEInterface)

Evaluate TDS ODE right-hand side (explicit form for Rosenbrock/Rodas).

Implements: `du/dt = [f(x,y); g(x,y)]` with mass matrix handled by solver.

The Float64 method is the zero-allocation hot path and feeds the adapter's recorder.
Other element types (ForwardDiff Duals, e.g. when DiffEq's DAE initializer
differentiates the right-hand side) evaluate through the ADCache without recording.
"""
function evaluate!(
    du_out::AbstractVector{Float64},
    adapter::TDSAdapter,
    u::AbstractVector{Float64},
    t::Real,
    ::ODEInterface
)
    _evaluate_packed!(du_out, adapter, u, t, adapter.recorder)
end

function evaluate!(du_out::AbstractVector, adapter::TDSAdapter, u::AbstractVector, t::Real, ::ODEInterface)
    _evaluate_packed!(du_out, adapter, u, t, nothing)
end

"""
    _evaluate_packed!(du, adapter, u, t, recorder)

Unpack `u = [x; y]` into the adapter's buffers for element type `eltype(du)`,
evaluate `[f; g]` at time `t`, and pack the result into `du`. Float64 uses the
adapter's `sys_ws` (zero-allocation hot path); other types (Duals) use the ADCache
with cached typed buffers. `u` may have a different element type than `du` (e.g.
a Float64 state with a Dual time for ∂F/∂t); values are converted on unpacking.

This is the one evaluation routine behind `evaluate!` and the Jacobian/tgrad
functors. Only the ODE right-hand side passes a `recorder`. Rows of the adapter's
`reduction` are replaced by the time derivative of their constraint.
"""
@inline function _evaluate_packed!(du::AbstractVector, adapter::TDSAdapter,
                                   u::AbstractVector, t, recorder)
    _evaluate_unreduced!(du, adapter, u, t, recorder, adapter.status)
    red = adapter.reduction
    red === nothing || _apply_reduction!(du, red, u, t)
    nothing
end

"""
    _evaluate_unreduced!(du, adapter, u, t, recorder, status)

`[f; g]` at `u` and time `t` with the instance statuses `status`, without the
reduction of state-only constraints: the DAE as the models define it.
"""
@inline function _evaluate_unreduced!(du::AbstractVector{T}, adapter::TDSAdapter,
                                      u::AbstractVector, t, recorder, status) where T
    x_buf, y_buf, f_buf, g_buf = _get_typed_buffers(adapter, T)
    n_diff = adapter.n_diff
    n_alg = adapter.n_alg

    # Unpack: u = [x; y]
    @inbounds for i in 1:n_diff
        x_buf[i] = u[i]
    end
    @inbounds for i in 1:n_alg
        y_buf[i] = u[n_diff + i]
    end

    # Using adapter.sys_ws (not adcache.float_ws) avoids an allocation from ::Any
    # fields in ADCache on the Float64 hot path
    if T === Float64
        evaluate_system_residuals!(adapter.sys_ws, f_buf, g_buf, x_buf, y_buf, status, t; recorder)
    else
        adapter.adcache(f_buf, g_buf, x_buf, y_buf, status, t)
    end

    # Pack: du = [f; g]
    @inbounds for i in 1:n_diff
        du[i] = f_buf[i]
    end
    @inbounds for i in 1:n_alg
        du[n_diff + i] = g_buf[i]
    end
    nothing
end

# Type-stable buffer dispatch: the adapter's own buffers for Float64, cached Dual buffers otherwise
@inline function _get_typed_buffers(adapter::TDSAdapter, ::Type{Float64})
    (adapter.x_buf, adapter.y_buf, adapter.f_buf, adapter.g_buf)
end

@inline function _get_typed_buffers(adapter::TDSAdapter, ::Type{T}) where T
    buffers = get_tds_dual_buffers(adapter.adcache, T, adapter.n_diff, adapter.n_alg)
    # Type assertions help compiler (buffers.x_buf is already Vector{T})
    (buffers.x_buf::Vector{T}, buffers.y_buf::Vector{T},
     buffers.f_buf::Vector{T}, buffers.g_buf::Vector{T})
end

"""
    pack_state(adapter::TDSAdapter) -> Vector{Float64}

The adapter's current state as one vector `u = [x; y]` (differential then algebraic).
"""
pack_state(adapter::TDSAdapter) = vcat(adapter.x_buf, adapter.y_buf)

export TDSAdapter, ReinitCache, pack_state
export TDSInterface, ODEInterface
export evaluate!

# ============================================================================
# Service and Setpoint Population
# ============================================================================

"""
    populate_setpoints!(adcache::AbstractADCache, setpoints::NamedTuple)

Populate setpoints in ADCache workspaces from pre-computed initialization values.

This must be called after creating ADCache for TDS problems.

# Arguments
- `adcache`: ADCache with float_ws containing workspaces_by_name
- `setpoints`: From `prob.init.setpoints` (GENCLS=(P_m0=[], E_q0=[]), ...)
"""
function populate_setpoints!(adcache::AbstractADCache, setpoints::NamedTuple)
    workspaces_by_name = get_float_ws(adcache).workspaces_by_name

    for (model_name, ws) in pairs(workspaces_by_name)
        # Skip models without setpoints (network models don't need initialization)
        haskey(setpoints, model_name) || continue

        model_setpoints = setpoints[model_name]

        # Copy each setpoint to workspace setpoints buffer
        for (setpoint_name, setpoint_vals) in pairs(model_setpoints)
            dest_array = getfield(ws.setpoints, setpoint_name)
            copyto!(dest_array, setpoint_vals)
        end
    end

    nothing
end

# ============================================================================
# Buffer Initialization and Mass Matrix
# ============================================================================

"""
    init_tds_buffers!(x_buf, y_buf, prob)
    init_tds_buffers!(x_buf, y_buf, am, sys, init::TDSInitData)

Initialize TDS buffers from `prob.init` (network model values, dyn states, pf snapshot),
or from `init` with the addresses `am` of the system `sys`.

With `prob`, the state at the start of the problem, `prob.tspan[1]`: a sinusoidal steady state
(`SinusoidalSteadyState`, whose phasors are at `t = 0`) is taken there (`_steady_state_at!`); with
`init`, the state as the initialization builds it, at `t = 0`.

This is called once during adapter construction to populate initial conditions.
"""
function init_tds_buffers!(x_buf::AbstractVector, y_buf::AbstractVector, prob::TDSProblem)
    init_tds_buffers!(x_buf, y_buf, prob.address_manager, prob.sys, prob.init)
    _steady_state_at!(x_buf, nothing, prob)
    return nothing
end

function init_tds_buffers!(
    x_buf::AbstractVector,
    y_buf::AbstractVector,
    am::AddressManager,
    sys::SystemModel,
    init::TDSInitData
)
    dyn = init.dyn
    pf = init.pf

    # Populate the network models' own variables from init.network: states into x,
    # algebraic variables into y
    for (model_name, values_by_var) in pairs(init.network)
        haskey(sys.models, model_name) || continue
        metadata = model_metadata(get_base_model_type(sys.models[model_name]))

        for var in metadata.vars
            is_external(var) && continue
            haskey(values_by_var, var.name) || continue

            buf = is_state_var(var) ? x_buf : y_buf
            addrs = get_addresses(am, model_name, var.name)
            values = getfield(values_by_var, var.name)

            @inbounds for i in eachindex(addrs)
                buf[addrs[i]] = values[i]
            end
        end
    end

    # Populate differential states from init.dyn (dynamic models)
    for (model_name, states) in pairs(dyn)
        models = sys.models[model_name]
        model_type = get_base_model_type(models)
        metadata = model_metadata(model_type)

        for var in metadata.vars
            is_state_var(var) || continue

            var_name = var.name
            haskey(states, var_name) || continue

            addrs = get_addresses(am, model_name, var_name)
            values = getfield(states, var_name)

            @inbounds for i in eachindex(addrs)
                x_buf[addrs[i]] = values[i]
            end
        end
    end

    # Populate the phasor bus voltages (v, theta) from the PF snapshot, by bus number
    if haskey(sys.models, :Bus)
        pf_index = pf_snapshot_positions(pf)
        v_addrs = get_addresses(am, :Bus, :v)
        theta_addrs = get_addresses(am, :Bus, :theta)
        buses = sys.models.Bus.i
        @inbounds for i in eachindex(v_addrs)
            y_buf[v_addrs[i]] = pf.v[pf_index[buses[i]]]
            y_buf[theta_addrs[i]] = pf.theta[pf_index[buses[i]]]
        end
    end

    # Overwrite with dynamic model algebraic variables (like I_d, I_q from GENCLS, or the
    # power of a voltage playback, which has no states)
    for (model_name, states) in pairs(dyn)
        models = sys.models[model_name]
        model_type = get_base_model_type(models)
        metadata = model_metadata(model_type)

        for var in metadata.vars
            is_state_var(var) && continue
            is_external(var) && continue

            var_name = var.name
            haskey(states, var_name) || continue

            addrs = get_addresses(am, model_name, var_name)
            values = getfield(states, var_name)

            @inbounds for i in eachindex(addrs)
                y_buf[addrs[i]] = values[i]
            end
        end
    end

    nothing
end

"""
    init_tds_rates(prob, n_diff) -> Vector{Float64}

Rates of change of the differential states at the start of the problem, `prob.tspan[1]`, from
the steady state the problem starts in (`initial_rates(prob.init.steady_state)`, zero for states
not listed; a sinusoidal one taken at the start, `_steady_state_at!`): it satisfies `M·rate = f`,
which `validate_tds_initialization` checks.
"""
function init_tds_rates(prob::TDSProblem, n_diff::Int)
    rate_buf = _rates_at_zero(prob, n_diff)
    if prob.init.steady_state isa SinusoidalSteadyState && !iszero(prob.tspan[1])
        # (the rates at the start need the states at t = 0 too)
        x0 = zeros(Float64, n_diff)
        y0 = zeros(Float64, Int(max_address_in_space(prob.address_manager, AlgebVar())))
        init_tds_buffers!(x0, y0, prob.address_manager, prob.sys, prob.init)
        _steady_state_at!(x0, rate_buf, prob)
    end
    return rate_buf
end

# The rates of the states at t = 0, the time of the steady state's phasors (zero for states not listed)
function _rates_at_zero(prob::TDSProblem, n_diff::Int)
    rate_buf = zeros(Float64, n_diff)
    am = prob.address_manager
    for (model_name, rates) in pairs(initial_rates(prob.init.steady_state))
        haskey(prob.sys.models, model_name) || continue
        for (var_name, values) in pairs(rates)
            addrs = get_addresses(am, model_name, var_name)
            @inbounds for i in eachindex(addrs)
                rate_buf[addrs[i]] = values[i]
            end
        end
    end
    return rate_buf
end

"""
    _steady_state_at!(x, rate, prob)

Take the states `x` (and their rates `rate`, or `nothing`) of the steady state the problem starts
in from `t = 0`, where the initialization builds them, to the start of the problem `t₁ =
prob.tspan[1]`: for a `SinusoidalSteadyState`, each state it lists (the states that oscillate; those
at rest are not listed, `initial_rates`) oscillates at `ω`, the
system's `omega_base`, so `x(t₁) = x₀ cos ωt₁ + (ẋ₀/ω) sin ωt₁` and `ẋ(t₁) = ẋ₀ cos ωt₁ −
ω x₀ sin ωt₁`; the states at rest, and every state of an `Equilibrium`, are as at `t = 0`.
Without it, a problem starting at `t₁ ≠ 0` started from the steady state at `t = 0`: its
interfaces' angles off by `ω t₁` and DiffEq's initialization check failing.
"""
function _steady_state_at!(x::AbstractVector, rate, prob::TDSProblem)
    prob.init.steady_state isa SinusoidalSteadyState || return nothing
    t1 = prob.tspan[1]
    iszero(t1) && return nothing
    ω = prob.sys.settings.omega_base
    c, s = cos(ω * t1), sin(ω * t1)
    r0 = _rates_at_zero(prob, length(x))
    am = prob.address_manager
    for (model_name, rates) in pairs(initial_rates(prob.init.steady_state))
        haskey(prob.sys.models, model_name) || continue
        for var_name in keys(rates), a in get_addresses(am, model_name, var_name)
            x0, ẋ0 = x[a], r0[a]
            x[a] = x0 * c + ẋ0 / ω * s
            rate === nothing || (rate[a] = ẋ0 * c - ω * x0 * s)
        end
    end
    return nothing
end

"""
    build_tds_mass_matrix(prob::TDSProblem, x_buf::AbstractVector, y_buf::AbstractVector)
    build_tds_mass_matrix(sys, am, x_buf, y_buf)

Build sparse diagonal mass matrix `M` for the DAE: `M * du/dt = F(u)`, from the models of
`prob.sys` (or `sys`, with the addresses `am`).

Returns sparse matrix of size (n_diff + n_alg) x (n_diff + n_alg):

    M = diag(M_x, 0)

where `M_x` contains mass coefficients for differential states (e.g., inertia `2H/ω_s`)
and zeros for algebraic equations (making `M` singular).
"""
build_tds_mass_matrix(prob::TDSProblem, x_buf::AbstractVector, y_buf::AbstractVector) =
    build_tds_mass_matrix(prob.sys, prob.address_manager, x_buf, y_buf)

function build_tds_mass_matrix(sys::SystemModel, am::AddressManager, x_buf::AbstractVector,
                               y_buf::AbstractVector)
    n_diff = length(x_buf)
    n_alg = length(y_buf)
    n_total = n_diff + n_alg
    T = eltype(x_buf)

    # Initialize mass diagonal (1.0 for differential, 0.0 for algebraic)
    M_diag = Vector{T}(undef, n_total)
    M_diag[1:n_diff] .= one(T)
    M_diag[n_diff+1:end] .= zero(T)

    # Override with model-specific mass coefficients

    for (model_name, models) in pairs(sys.models)
        isempty(models) && continue

        model_type = get_base_model_type(models)
        metadata = model_metadata(model_type)

        for var_meta in metadata.vars
            is_state_var(var_meta) || continue

            mass_coeff = get_mass_coeff(var_meta)
            isnothing(mass_coeff) && continue

            # Get addresses and mass values
            addrs = get_addresses(am, model_name, var_meta.name)
            # Type-stable dispatch: Symbol (field lookup) or Function (computed)
            mass_values = eval_mass_coeff(mass_coeff, models)

            for (i, addr) in enumerate(addrs)
                M_diag[addr] = mass_values[i]
            end
        end
    end

    return spdiagm(M_diag)
end
