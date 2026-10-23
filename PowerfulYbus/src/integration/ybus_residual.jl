# ============================================================================
# Y-Bus DAE Residual
# ============================================================================
#
# The Y-bus TDS model as a mass-matrix DAE residual for OrdinaryDiffEq.jl, with
# two treatments of the network equations:
#
# - EliminatedNetwork: the network is solved inside every evaluation; the bus
#   voltage variables are pinned to the solved voltages. The Jacobian is dense in
#   the machine states (the solved voltages depend on every internal voltage).
# - SimultaneousNetwork: the bus voltages are DAE unknowns and the network
#   equations are DAE residuals. The Jacobian is sparse.
#
# The residual is a function of (u, t) only: every evaluation gathers its inputs
# from u, and the workspaces of each number type (Float64, ForwardDiff Duals) are
# built from the same initialization.
#
# ============================================================================

using PowerfulCore: StateVar, AlgebVar, max_address_in_space, get_addresses, PowerInjection

export NetworkTreatment, EliminatedNetwork, SimultaneousNetwork
export YBusResidual, consistent_algebraics!

"""
    NetworkTreatment

How the Y-bus DAE residual treats the network equations: [`EliminatedNetwork`](@ref)
or [`SimultaneousNetwork`](@ref).
"""
abstract type NetworkTreatment end

"""
    EliminatedNetwork()

The network equations are solved inside every residual evaluation
([`solve_network!`](@ref)) from the generator internal voltages; the bus voltage
variables (Bus `v`, `theta`) are pinned to the solved voltages V by the residuals
`v - |V|` and `theta - ∠V` at the Bus `q` and `p` rows. ForwardDiff derivatives of
the solved voltages come from the implicit function theorem. The Jacobian is dense
in the machine states: every solved voltage depends on every internal voltage.

With nonlinear network equations (constant-power or constant-current loads, static PV
generators) the network solve iterates from the previous solution until its residual
is below `abstol`, so the residual is a function of `(u, t)` up to that tolerance.

Current-source models (e.g. REGCA1) are not supported: their output current is
computed by the model from the bus voltage, outside the network solve.
"""
struct EliminatedNetwork <: NetworkTreatment end

"""
    SimultaneousNetwork()

The bus voltages V = v∠theta are DAE unknowns (the Bus `v`, `theta` variables) and
the network equations are DAE residuals at the Bus `p` and `q` rows: the real and
imaginary parts of the KCL mismatch `Y V - I_E - I_inj(V)`, or at a static generator
bus its constraints (see [`solve_network!`](@ref)). The Jacobian is sparse: machine
blocks and the pattern of `Y`.
"""
struct SimultaneousNetwork <: NetworkTreatment end

"""
    YBusWorkspace{T}

Model and network buffers of one number type `T`.
"""
struct YBusWorkspace{T, WS<:YBusSystemBuffers{T}, NC<:NetworkCache{T}}
    sys_ws::WS
    network::NC
end

"""
    YBusResidual(prob::YBusTDSProblem, treatment::NetworkTreatment; maxiters=50, abstol=1e-8)

The Y-bus TDS DAE residual `du = F(u, t)` (mass-matrix form `M du/dt = F(u, t)`)
with the network treatment `treatment`; callable as `r(du, u, p, t)` (the
OrdinaryDiffEq signature). `u` and `du` hold the differential states then the
algebraic variables, in the address layout of `prob.base_prob`.

The number type of the evaluation is `eltype(du)`. The Float64 buffers are built
at construction; the buffers of a ForwardDiff Dual type are built on its first
evaluation from the same initialization (setpoints included). The automata of the
models are those of `prob.sys_ws` (the buffers of the partitioned algebraic solve,
`consistent_algebraics!`): every buffer reads the discrete state the automaton runtime
sets there. The network solve
limits `maxiters`, `abstol` apply to `EliminatedNetwork` and to
[`consistent_algebraics!`](@ref).

Where the network does not converge in an `EliminatedNetwork` evaluation (a numerical
failure) the evaluation returns NaN residuals, which the integrator rejects as it does a step
that blows up, and records the failure (`failure`: the message naming the time, the iterations,
the residual and the worst rows). A right-hand side has no channel to terminate an integration;
`solve` reads the record when the integration ends without success (`network_failure_retcode`),
warns and returns `ConvergenceFailure`. A later evaluation that converges clears the record: a
failure at a rejected trial point is not one of the integration.
"""
struct YBusResidual{N<:NetworkTreatment, P<:YBusTDSProblem, W<:YBusWorkspace{Float64}}
    treatment::N
    prob::P
    float::W
    dual::Dict{DataType, YBusWorkspace}
    n_diff::Int
    n_alg::Int
    v_addrs::Vector{Int}      # Bus v (algebraic variable addresses)
    theta_addrs::Vector{Int}  # Bus theta
    p_rows::Vector{Int}       # Bus p (algebraic residual addresses)
    q_rows::Vector{Int}       # Bus q
    vmin::Float64             # constant-power load voltage threshold (SystemConfig)
    maxiters::Int
    abstol::Float64
    failure::Base.RefValue{Union{Nothing, String}}  # the last `EliminatedNetwork` evaluation's failure
end

function YBusResidual(prob::YBusTDSProblem, treatment::NetworkTreatment;
                      maxiters::Int = 50, abstol::Real = 1e-8)
    bp = prob.base_prob
    am, sys = bp.address_manager, bp.sys
    sys_ws = new_buffers(prob, Float64)
    float = YBusWorkspace(sys_ws, similar_network_cache(prob.network_cache, Float64))
    _check_supported(treatment, float)
    return YBusResidual(treatment, prob, float, Dict{DataType, YBusWorkspace}(),
                        Int(max_address_in_space(am, StateVar())),
                        Int(max_address_in_space(am, AlgebVar())),
                        Int.(get_addresses(am, :Bus, :v)), Int.(get_addresses(am, :Bus, :theta)),
                        Int.(get_addresses(am, :Bus, :p)), Int.(get_addresses(am, :Bus, :q)),
                        Float64(sys.settings.constant_power_vmin), maxiters, Float64(abstol),
                        Ref{Union{Nothing, String}}(nothing))
end

_check_supported(::SimultaneousNetwork, ::YBusWorkspace) = nothing
function _check_supported(::EliminatedNetwork, ws::YBusWorkspace)
    sources = first.(ws.network.cs_buf.src_model_map)
    isempty(sources) ||
        error("EliminatedNetwork does not support current-source models ($(join(sources, ", "))): " *
              "their output current is computed by the model from the bus voltage, outside " *
              "the network solve. Use SimultaneousNetwork().")
    nothing
end

function (r::YBusResidual)(du::AbstractVector{T}, u::AbstractVector, p, t) where T
    n = r.n_diff + r.n_alg
    length(u) == length(du) == n ||
        error("Y-bus residual: u and du must have $n entries, got $(length(u)) and $(length(du))")
    dae = (f = view(du, 1:r.n_diff), g = view(du, (r.n_diff + 1):n),
           x = view(u, 1:r.n_diff), y = view(u, (r.n_diff + 1):n))
    _evaluate!(r.treatment, r, _workspace(r, T), dae, t)
    return nothing
end

_workspace(r::YBusResidual, ::Type{Float64}) = r.float
function _workspace(r::YBusResidual, ::Type{T}) where T
    ws = get(r.dual, T, nothing)
    if ws === nothing
        sys_ws = new_buffers(r.prob, T)
        sync_automaton_memory!(sys_ws, r.float.sys_ws)
        ws = r.dual[T] = YBusWorkspace(sys_ws, similar_network_cache(r.float.network, T))
    end
    return ws::YBusWorkspace{T}
end

# The Y-bus DAE's residuals, for the automaton runtime's check of the algebraic equations
dae_residuals!(du, r::YBusResidual, u, t) = r(du, u, nothing, t)

# The Y-bus DAE is a phasor DAE that starts from the equilibrium of its power flow: every state
# at rest at the start (`PowerfulTDS._start_state_rates`, for the automata's start)
_start_state_rates(::YBusResidual, n_diff::Int) = zeros(n_diff)

# The rate of a driven start of the Y-bus DAE, for a value of an automaton on a surface there
# (`PowerfulTDS._start_algebraic_rate`): the states at `rates`, the algebraic variables at the rate
# the algebraic rows give (`algebraic_rate`, the simultaneous network's, which the eliminated one's
# solution satisfies too). Its algebraic solve has `y` alone for unknowns: a state without mass
# would need its row solved with them, and is refused by name
function _start_algebraic_rate(r::YBusResidual, x::AbstractVector, y::AbstractVector, t::Real,
                               rates::AbstractVector, h::Real)
    M = r.prob.M
    massless = [k for k in 1:r.n_diff if iszero(M[k, k])]
    isempty(massless) || return "the Y-bus DAE solves its algebraic variables alone, and states " *
        "$(join(massless[1:min(5, end)], ", ")) have no mass: their rate at a driven start is not available there"
    ydot = algebraic_rate(r.prob.algebraic, r.prob, x, y, Float64(t), rates, Float64(h))
    ydot isa AbstractString && return ydot
    return vcat(rates, ydot)
end

# Where the automata change the algebraic equations of the Y-bus DAE, nothing is recorded
# here: its reinitialization (`consistent_algebraics!`) solves with `AlgebraicNewton`, which
# watches the modes of the automata itself (`_stale`) and takes a full factorization with
# pivoting at the first solve after one changes, never refactorizing with the old pivots
algebraic_equations_changed!(::YBusResidual) = nothing

# The buffers of the Y-bus DAE of element type `T`, for the rates of the values of the
# automata (the automaton runtime's sliding check of located transitions)
dual_workspace(r::YBusResidual, ::Type{T}) where T = _workspace(r, T).sys_ws

# The memory of the automata, which the runtime sets in the Float64 buffers, in the
# buffers of every Dual type
function sync_automaton_memory!(r::YBusResidual)
    for ws in values(r.dual)
        sync_automaton_memory!(ws.sys_ws, r.float.sys_ws)
    end
    return nothing
end

# Network solve of the eliminated treatment: Float64, or Dual with implicit-function
# derivatives from the Float64 cache
_solve_network!(net::NetworkCache{Float64}, r::YBusResidual) =
    solve_network!(net, r.prob.base_prob.sys; maxiters=r.maxiters, abstol=r.abstol)
_solve_network!(net::NetworkCache, r::YBusResidual) =
    solve_network!(net, r.float.network, r.prob.base_prob.sys; maxiters=r.maxiters, abstol=r.abstol)

function _evaluate!(::EliminatedNetwork, r::YBusResidual, ws::YBusWorkspace, dae, t)
    net = ws.network
    gather_internal_voltages!(net.E_buf, ws.sys_ws, dae)
    sol = _solve_network!(net, r)
    if !sol.converged
        r.failure[] = network_failure_message(sol, t)
        fill!(dae.f, NaN)
        fill!(dae.g, NaN)
        return nothing
    end
    r.failure[] = nothing
    evaluate_ybus_residuals!(ws.sys_ws, dae, sol, net.V_bus, t)

    # Bus voltage variables pinned to the solved voltages. An angle is defined modulo 2π (the
    # equations are periodic in it, and a variable that follows a trajectory is not confined to
    # (-π, π]): the angle row is the difference of the angles, taken in (-π, π]
    g, y, V = dae.g, dae.y, net.V_bus
    @inbounds for k in eachindex(V)
        g[r.q_rows[k]] = y[r.v_addrs[k]] - abs(V[k])
        g[r.p_rows[k]] = -angle(V[k] * cis(-y[r.theta_addrs[k]]))
    end
    nothing
end

_evaluate!(::SimultaneousNetwork, r::YBusResidual, ws::YBusWorkspace, dae, t) =
    _simultaneous_residuals!(ws, r.v_addrs, r.theta_addrs, r.p_rows, r.q_rows, r.vmin, dae, t)

# The residuals of the Y-bus DAE with the network simultaneous at `dae` (`f`, `g`, `x`, `y`),
# into `dae.f` and `dae.g`: the bus voltages from the Bus `v` and `theta` variables of `y` (at
# `v_addrs`, `theta_addrs`), the models' residuals at them, and the network equations at the
# Bus `p` and `q` rows (`p_rows`, `q_rows`). The residual of the algebraic solve
# (`AlgebraicNewton`) and of the DAE integration, in the buffers `ws` of any number type.
function _simultaneous_residuals!(ws::YBusWorkspace{T}, v_addrs, theta_addrs, p_rows, q_rows, vmin, dae, t) where T
    net = ws.network
    g, y, V = dae.g, dae.y, net.V_bus
    @inbounds for k in eachindex(V)
        V[k] = y[v_addrs[k]] * cis(y[theta_addrs[k]])
    end
    gather_internal_voltages!(net.E_buf, ws.sys_ws, dae)
    sol = _network_solution!(net, true, 0, zero(T))  # generator outputs at V
    evaluate_ybus_residuals!(ws.sys_ws, dae, sol, V, t)

    # Network equations; current-source output currents from the evaluation above
    gather_source_currents!(net.cs_buf, ws.sys_ws, dae)
    _network_mismatch!(net, network_jacobian!(net), net.E_buf.E_internal, vmin)
    res = _equation_residual!(net)
    @inbounds for k in eachindex(res)
        g[p_rows[k]] = real(res[k])
        g[q_rows[k]] = imag(res[k])
    end
    return sol
end

"""
    network_failure_retcode(r::YBusResidual, retcode) -> retcode

The return code of an integration of `r` that ended with `retcode`: where it did not succeed
and the last `EliminatedNetwork` evaluation failed (`YBusResidual`), `ConvergenceFailure`, after
a warning with the failure's diagnostics; `retcode` otherwise (an integration already terminated
with `ConvergenceFailure` has reported its own, `consistent_algebraics!`).
"""
function network_failure_retcode(r::YBusResidual, retcode)
    failure = r.failure[]
    (failure === nothing || retcode in (ReturnCode.Success, ReturnCode.Terminated, ReturnCode.ConvergenceFailure)) &&
        return retcode
    @warn failure * "; stopping the simulation"
    return ReturnCode.ConvergenceFailure
end

"""
    consistent_algebraics!(u::Vector{Float64}, r::YBusResidual, t) -> Bool

Set the algebraic variables of `u` consistent with its differential states and the
current network topology, for the initial state and after a topology change.

Repeats the algebraic solution of the trapezoidal solver (`solve_algebraic!`: the Newton
iteration on the bus voltages and the models' algebraic variables) until the largest algebraic
residual of `r` is below `r.abstol`. Returns `true` when it is; where a solve does not converge
or the residual is not below `r.abstol` after `r.maxiters` passes (a numerical failure) it warns
(the time, the iterations, the residual and its worst rows) and returns `false`, the callers
stopping with a failing retcode. Errors in the models or the code throw.
"""
function consistent_algebraics!(u::Vector{Float64}, r::YBusResidual, t)
    du, work = similar(u), similar(u)
    g = view(du, (r.n_diff + 1):length(du))
    residual = Inf
    for pass in 0:r.maxiters
        r(du, u, nothing, t)
        residual = maximum(abs, g; init=0.0)
        residual < r.abstol && return true
        pass == r.maxiters && break
        sol = solve_algebraic!(u, work, r.prob, Float64(t), r.maxiters, r.abstol)
        if !sol.converged
            @warn network_failure_message(sol, t) * "; stopping the simulation"
            return false
        end
    end
    am = r.prob.base_prob.address_manager
    isnothing(am.reverse_map) && build_reverse_map!(am)
    worst = sortperm(abs.(g); rev = true)[1:min(5, length(g))]
    rows = join(("$(lookup_address(am, AlgebRes(), UInt(k))): $(abs(g[k]))" for k in worst), "; ")
    @warn "Y-bus algebraic variables not consistent at t = $t after $(r.maxiters) passes: " *
          "max algebraic residual $residual against $(r.abstol), worst rows (model, variable, instance) $rows; " *
          "stopping the simulation"
    return false
end
