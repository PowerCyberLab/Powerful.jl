# ============================================================================
# Coupled Newton Solve of the Algebraic Variables
# ============================================================================
#
# The algebraic variables `y` of the Y-bus DAE (the bus voltages, `v` and `theta` of every
# Bus, and the algebraic variables of the models) for the differential states `x` frozen
# are the solution of
#
#     g(x, y) = 0,
#
# the residuals of the DAE with the network simultaneous (`_simultaneous_residuals!`): at each
# Bus the real and imaginary parts of the KCL mismatch, at each model its algebraic rows. The
# network is in the equations, not solved apart from them, so that a block's variables and the
# bus voltages that its rows read move together in every Newton step, whatever the coupling
# between them (a converter current P/V, a clamp that holds a bus voltage): the Newton
# iteration
#
#     y <- y - J⁻¹ g(x, y),     J = ∂g/∂y,
#
# converges quadratically wherever it does at all, with no condition on the coupling that an
# alternation of the network solve and the models' steps needs.
#
# J is sparse, as the Jacobian of the DAE: the models' blocks, and the pattern of the
# admittance matrix at the Bus rows (a Bus row reads the voltages of the buses a branch
# connects it to, and the variables of the devices at it). It is taken by ForwardDiff with
# the coloring of the pattern from the model structure (`jacobian_sparsity_pattern`), and
# factored by KLU: no dense matrix is formed. The sparse system of a step is the one the DAE
# integration solves at every step, with `y` alone for unknowns.
#
# The iteration is the Newton engine of PowerfulCore (`newton_solve!`, shared with the
# reinitialization of the phasor simulation), with a policy of its own:
#  - J is kept as it is (a chord Newton) from one solve to the next, where the state that sets it
#    moves little (`reuse_jacobian`), and taken again
#     - where the equations differ from those it was taken at by something other than the state,
#       which is known structurally: a change of topology, of the data of a model (`invalidate!`)
#       or of a mode of an automaton (`_stale`), with a factorization afresh (pivoting); and
#     - where an iteration with it contracts the residual by less than `CHORD_CONTRACTION`: the
#       state has moved too far for it (the only criterion that is not structural: the
#       nonlinear dependence of the entries on the state is not known beforehand);
#  - the step is limited to 0.5 pu (`max_step`: at a bus whose voltage is near zero the angle
#    column of J is proportional to it, and a Newton step in the angle is a multiple of 2π long to
#    a root that is the same mod 2π), a trial point keeps every bus voltage magnitude above a
#    thousandth of its value (`positive`: a voltage that crosses zero reaches the other polar
#    solution), and the step is halved until it decreases the residual (down to 1/512 of it), the
#    Jacobian taken again where it finds no decrease (`strict_search`).
#
# ============================================================================

using ForwardDiff: ForwardDiff
using DifferentiationInterface: DifferentiationInterface
using LinearAlgebra: LinearAlgebra
using PowerfulCore: jacobian_sparsity_pattern, create_di_backend, ForwardDiffBackend, AlgebRes, AlgebVar
using PowerfulCore: build_reverse_map!, lookup_address
using PowerfulCore: NewtonPolicy, NewtonBuffers, BUS_VOLTAGE_STEP_LIMIT, KLULinearSolver, newton_solve!, reset_factorization!, factorize!,
    linear_solve!
using KLU: KLU

# An iteration that contracts the residual by less than this (new / old) with a kept
# Jacobian takes the Jacobian again at the next one
const CHORD_CONTRACTION = 0.25

# The iterations a solve may take by default, as many as the limit of the network solve
const MAX_ITERATIONS = 50

"""
    AlgebraicNewton

The Newton solve of the algebraic variables of a `YBusTDSProblem` (`solve_algebraic_system!`), with the
Jacobian it keeps from one solve to the next (`J`, factored in `solver`) and what it was taken
at (`topology`, `generation`, the modes of the automata `saved_modes`: `_stale`).

`history` is the largest residual at the start and after each step of the last solve, `refreshes` the
number of times a Jacobian was taken, `iterations` the number of steps taken, over all solves.
"""
mutable struct AlgebraicNewton
    n_diff::Int
    n_alg::Int
    v_addrs::Vector{Int}
    theta_addrs::Vector{Int}
    p_rows::Vector{Int}
    q_rows::Vector{Int}
    vmin::Float64
    pattern::SparseMatrixCSC{Bool, Int}
    J::SparseMatrixCSC{Float64, Int}
    solver::KLULinearSolver
    policy::NewtonPolicy
    bufs::NewtonBuffers{Vector{Float64}, Vector{Float64}, Vector{Float64}, SparseMatrixCSC{Float64, Int}}
    # The sparse ForwardDiff Jacobian (built when first taken)
    closure::Any
    prep::Any
    backend::Any
    workspaces::Dict{DataType, Any}
    float_ws::Any
    # What the Jacobian was taken at
    taken::Bool
    topology::Any
    generation::Int
    taken_generation::Int
    watched::Vector{Symbol}
    live_modes::Vector{Vector{Int}}
    saved_modes::Vector{Vector{Int}}
    # Whether the last step reused the Jacobian and contracted the residual too little
    refresh_due::Bool
    # The iterate, the differential rows of the residual (not solved for) and the output of the
    # Jacobian's evaluation
    y::Vector{Float64}
    f::Vector{Float64}
    scratch::Vector{Float64}
    history::Vector{Float64}
    refreshes::Int
    iterations::Int
end

# The live mode arrays of the automata of the models of the problem's buffers, with the
# model of each
function _mode_arrays(sys_ws)
    owners, arrays = Symbol[], Vector{Int}[]
    for name in sys_ws.automaton_models
        for modes in values(sys_ws.workspaces_by_name[name].automata.modes)
            push!(owners, name)
            push!(arrays, modes)
        end
    end
    return owners, arrays
end

function AlgebraicNewton(am, sys, sys_ws)
    n_diff = Int(max_address_in_space(am, StateVar()))
    n_alg = Int(max_address_in_space(am, AlgebVar()))
    pattern = jacobian_sparsity_pattern(am, sys, PowerInjection();
                                        models = filter(m -> haskey(am.assignments, m), keys(sys.models)))
    owners, modes = _mode_arrays(sys_ws)
    J = Float64.(pattern)
    v_addrs = Int.(get_addresses(am, :Bus, :v))
    policy = NewtonPolicy(; jacobian_every = typemax(Int), stall_ratio = CHORD_CONTRACTION, stall_reused_only = true, reuse_jacobian = true,
                          max_step = BUS_VOLTAGE_STEP_LIMIT, positive = v_addrs, max_tries = 10, merit_norm = 2.0, sufficient_decrease = 1e-4, strict_search = true)
    bufs = NewtonBuffers(zeros(n_alg), zeros(n_alg), zeros(n_alg), zeros(n_alg), J)
    return AlgebraicNewton(n_diff, n_alg, v_addrs, Int.(get_addresses(am, :Bus, :theta)),
                           Int.(get_addresses(am, :Bus, :p)), Int.(get_addresses(am, :Bus, :q)),
                           Float64(sys.settings.constant_power_vmin), pattern, J, KLULinearSolver(), policy, bufs,
                           nothing, nothing, nothing, Dict{DataType, Any}(), nothing,
                           false, nothing, 0, 0, owners, modes, map(copy, modes), false,
                           zeros(n_alg), zeros(n_diff), zeros(n_alg), Float64[], 0, 0)
end

"""
    invalidate!(an::AlgebraicNewton)

The data of a model changed (an event): the Jacobian is taken again at the next solve.
"""
invalidate!(an::AlgebraicNewton) = (an.generation += 1; nothing)

# Whether the equations differ from those the Jacobian was taken at, other than by the state: the
# topology or the data of a model changed (or there is no Jacobian yet), or an automaton changed a mode
function _stale(an::AlgebraicNewton, topology)
    (an.taken && an.topology === topology && an.generation == an.taken_generation) || return true
    @inbounds for k in eachindex(an.live_modes)
        an.live_modes[k] == an.saved_modes[k] || return true
    end
    return false
end

# The buffers of the models and the network of the number type `T` of the Jacobian's Dual
# evaluations (built on first use from the problem's initialization, sharing the automata
# of its buffers, with the setpoints it has now), with the full state they are evaluated at
function _dual_buffers(an::AlgebraicNewton, prob, ::Type{T}) where T
    return get!(an.workspaces, T) do
        sys_ws = new_buffers(prob, T)
        n = an.n_diff + an.n_alg
        (ws = YBusWorkspace(sys_ws, similar_network_cache(prob.network_cache, T)),
         u = zeros(T, n), du = zeros(T, n))
    end
end

# `y ↦ g(x, y)` at the differential states `x` and the time `t`, for DifferentiationInterface
struct JointResidual{A, P}
    an::A
    prob::P
    x::Vector{Float64}
    t::Base.RefValue{Float64}
end

function (r::JointResidual)(g::AbstractVector{T}, y::AbstractVector{T}) where T
    an = r.an
    w = _dual_buffers(an, r.prob, T)
    sync_automaton_memory!(w.ws.sys_ws, r.prob.sys_ws)
    n_diff = an.n_diff
    @inbounds for k in 1:n_diff
        w.u[k] = r.x[k]
    end
    dae = (f = view(w.du, 1:n_diff), g = g, x = view(w.u, 1:n_diff), y = y)
    _simultaneous_residuals!(w.ws, an.v_addrs, an.theta_addrs, an.p_rows, an.q_rows, an.vmin, dae, r.t[])
    return nothing
end

# The Jacobian ∂g/∂y at `y` (the engine's `jacobian!`), with what it is taken at recorded
struct AlgebraicJacobian{A, P}
    an::A
    prob::P
end

function (j::AlgebraicJacobian)(J::SparseMatrixCSC{Float64, Int}, y::Vector{Float64})
    an = j.an
    DifferentiationInterface.jacobian!(an.closure, an.scratch, J, an.prep, an.backend, y)
    an.refreshes += 1
    an.topology = j.prob.network_cache.topology_cache.current[]
    an.taken_generation = an.generation
    for k in eachindex(an.live_modes)
        copyto!(an.saved_modes[k], an.live_modes[k])
    end
    an.taken = true
    return nothing
end

# The algebraic rows `g` at `y` (the engine's `residual!`), with the differential states `x`
# and the time `t`: the network solution of the buffers is that of `y`
struct AlgebraicResidual{A, X}
    an::A
    x::X
    t::Float64
end

function (r::AlgebraicResidual)(g::Vector{Float64}, y::Vector{Float64})
    an = r.an
    dae = (f = an.f, g = g, x = r.x, y = y)
    return _simultaneous_residuals!(an.float_ws, an.v_addrs, an.theta_addrs, an.p_rows, an.q_rows, an.vmin, dae, r.t)
end

# The model, the instance and the variable at algebraic variable address `a`
function _variable_at(prob, a::Integer)
    am = prob.base_prob.address_manager
    isnothing(am.reverse_map) && build_reverse_map!(am)
    return lookup_address(am, AlgebVar(), UInt(a))
end

# The instance key of instance `i` of `model`
_instance_key(prob, model::Symbol, i::Int) = get_key(prob.sys_ws.models_by_name[model])[i]

# The automaton modes of an instance, as text (empty without automata)
function _modes_text(prob, model::Symbol, instance::Int)
    sys_ws = prob.sys_ws
    model in sys_ws.automaton_models || return ""
    meta = sys_ws.metadata_by_name[model]
    state = sys_ws.workspaces_by_name[model].automata.modes
    modes = ["$(a.name) $(a.modes[getfield(state, a.name)[instance]].name)" for a in meta.automata]
    return isempty(modes) ? "" : ", in $(join(modes, ", "))"
end

# A singular Jacobian is a model whose rows do not determine its variables (or a network that
# does not determine the voltages): the variable of the first column with no entry (or no
# row with entries), else the column KLU found singular
function _singular_message(prob, an::AlgebraicNewton, err)
    J = an.J
    empty_column = findfirst(c -> all(iszero, view(nonzeros(J), nzrange(J, c))), 1:size(J, 2))
    column = empty_column === nothing ? min(max(Int(err.info), 1), size(J, 2)) : empty_column
    model, variable, instance = _variable_at(prob, column)
    key = join(string.(_instance_key(prob, model, instance)), ", ")
    pairs_ = _algebraic_pairs(model_metadata(get_base_model_type(prob.base_prob.sys.models[model])))
    rows = [r for (r, v) in pairs_]
    vars = [v for (r, v) in pairs_]
    what = isempty(pairs_) ? "the Jacobian with respect to $variable" :
           "the Jacobian of $(join(rows, ", ")) with respect to $(join(vars, ", "))"
    return "The algebraic rows of $model ($key) do not determine their variables: $what is singular " *
           "in the system of the models' rows and the network equations" * _modes_text(prob, model, instance)
end

# The algebraic rows with the `n` largest residuals of `g` as text naming each (model, instance
# key, variable the row is paired with and its value)
function _worst_rows(prob, g::AbstractVector, n::Int = 5)
    am = prob.base_prob.address_manager
    isnothing(am.reverse_map) && build_reverse_map!(am)
    worst = sortperm(abs.(g); rev = true)[1:min(n, length(g))]
    lines = String[]
    for k in worst
        model, row, instance = lookup_address(am, AlgebRes(), UInt(k))
        push!(lines, "$model ($(join(string.(_instance_key(prob, model, instance)), ", "))) $row: $(abs(g[k]))")
    end
    return join(lines, "; ")
end

# The Float64 buffers and the Jacobian's closure of `an` (built at the first use), the closure at the differential
# states `x` and the time `t`
function _prepared!(an::AlgebraicNewton, prob, x::AbstractVector, t::Float64)
    an.float_ws === nothing && (an.float_ws = YBusWorkspace(prob.sys_ws, prob.network_cache))
    if an.closure === nothing
        an.closure = JointResidual(an, prob, Vector{Float64}(x), Ref(t))
        an.backend = create_di_backend(ForwardDiffBackend(), an.pattern)
        an.prep = DifferentiationInterface.prepare_jacobian(an.closure, similar(an.y), an.backend, an.y)
    end
    copyto!(an.closure.x, x)
    an.closure.t[] = t
    return nothing
end

"""
    algebraic_rate(an::AlgebraicNewton, prob, x, y, t, rates, h) -> Union{Vector{Float64}, String}

The rate of the algebraic variables `y` at `(x, y, t)` along the rate `rates` of the differential states `x` (all
with mass), `g_y ẏ = -(g_x ẋ + g_t)`: `g_x ẋ + g_t` a central difference of the algebraic rows over `±h` along
`(ẋ, 1)` at `y` (the residual `y` leaves there cancelling in it), `g_y` the Jacobian of the algebraic solve at the
start, taken and factorized with pivoting by its solver (and kept for its next solve, as a Jacobian it takes is).
Where `g_y` is singular (KLU finds it so, or estimates its condition number at or above `1 / eps()`), what failed, as
text. For the rate of a driven start on a surface of an automaton (`PowerfulTDS._start_algebraic_rate`).
"""
function algebraic_rate(an::AlgebraicNewton, prob, x::AbstractVector, y::AbstractVector, t::Float64,
                        rates::AbstractVector, h::Float64)
    copyto!(an.y, y)
    _prepared!(an, prob, x, t)
    at(s) = (g = zeros(an.n_alg); AlgebraicResidual(an, x .+ s .* rates, t + s)(g, copy(an.y)); g)
    driven = (at(h) .- at(-h)) ./ (2h)
    reset_factorization!(an.solver)
    try
        AlgebraicJacobian(an, prob)(an.J, an.y)
        factorize!(an.solver, an.J)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        return _singular_message(prob, an, err)
    end
    κ = KLU.condest(an.solver.fact)
    κ * eps() >= 1 && return "the Jacobian of the algebraic rows with respect to their variables is singular at " *
        "t = $t: KLU estimates its condition number at $κ"
    ydot = similar(driven)
    linear_solve!(an.solver, ydot, driven)
    return -ydot
end

"""
    algebraic_solve!(an::AlgebraicNewton, prob, u, du, t, abstol, limit) -> NetworkSolution

Solve `g(x, y) = 0` for `y` (the algebraic part of `u`) at the differential states `x` (its
first part, frozen), by the Newton engine of PowerfulCore with the policy of `an`, until the
largest residual of the algebraic rows is below `abstol` in at most `limit` steps, and leave the
residuals there in `du`. A Jacobian that is to be taken again (`_stale`) is taken at the first
solve after the change, converged or not, so that rows that determine nothing are reported where
the equations change, not where they first move; it throws (`ArgumentError`) where it is singular.

Returns the network solution at the result. Where the residual is not below `abstol` after
`limit` steps, is not finite, or finds no decrease along the step of a Jacobian taken at the
iterate, the solution is not `converged` and its `detail` names the time, the iterations, the
residual and its worst rows.
"""
function algebraic_solve!(an::AlgebraicNewton, prob, u::Vector{Float64}, du::Vector{Float64}, t::Float64,
                          abstol::Float64, limit::Int)
    n_diff, n_alg = an.n_diff, an.n_alg
    x = view(u, 1:n_diff)
    y = an.y
    copyto!(y, 1, u, n_diff + 1, n_alg)
    _prepared!(an, prob, x, t)
    jacobian! = AlgebraicJacobian(an, prob)
    result = try
        stale = _stale(an, prob.network_cache.topology_cache.current[])
        if stale
            reset_factorization!(an.solver)
            jacobian!(an.J, y)
            factorize!(an.solver, an.J)
        end
        empty!(an.history)
        newton_solve!(AlgebraicResidual(an, x, t), jacobian!, y, an.bufs, an.solver, an.policy;
                      abstol, maxiters = limit, history = an.history, fresh_jacobian = stale, refresh_due = an.refresh_due)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        throw(ArgumentError(_singular_message(prob, an, err)))
    end
    an.iterations += result.steps
    an.refresh_due = result.refresh_due
    copyto!(u, n_diff + 1, y, 1, n_alg)
    copyto!(du, 1, an.f, 1, n_diff)
    copyto!(du, n_diff + 1, an.bufs.g, 1, n_alg)
    detail = result.converged ? "" : _failure_detail(prob, an, result, t, abstol)
    return _network_solution!(prob.network_cache, result.converged, result.steps, result.residual_norm, detail)
end

# What a solve that did not converge reports: the time, the iterations, the residual and the
# worst rows
function _failure_detail(prob, an::AlgebraicNewton, result, t, abstol)
    why = result.status === :stalled ?
        "no decrease of the residual along the Newton step with a Jacobian taken at the iterate" :
        result.status === :diverged ? "the residual is not finite" :
        "the Newton iteration on the network and the models' rows"
    return "The algebraic rows of the models did not converge at t = $t in $(result.steps) iterations " *
           "($why): residual $(result.residual_norm) against $abstol, at " *
           _worst_rows(prob, an.bufs.g)
end
