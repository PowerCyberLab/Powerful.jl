# ============================================================================
# Y-Bus TDS Problem
# ============================================================================
#
# Wrapper for TDS problem that uses the Y-bus network formulation.
# Provides a fast network solution; every model takes its bus voltages from
# the network's bus voltages.
#
# ============================================================================

using PowerfulCore: Disturbance
using PowerfulCore: StateVar, AlgebVar, AlgebRes, max_address_in_space
using PowerfulCore: model_metadata, get_base_model_type, get_addresses
using PowerfulCore: is_state_var, is_algeb_var, is_algeb_res, get_mass_coeff, eval_mass_coeff
using PowerfulCore: is_external
using SparseArrays: spdiagm

# Note: YBusEventSchedule is from ybus_event_schedule.jl (available in module scope)

"""
    YBusTDSProblem{T, PROB, WS, CACHE, MT, S}

A time-domain simulation problem that represents the network by its sparse admittance
(Y-bus) matrix instead of solving the full power-flow-style network equations at every
step. It wraps a standard [`TDSProblem`](@ref) (same models, variables and events) and
replaces only how the network is solved: every bus's voltage still comes from the same
network, but a topology change (a fault, a line trip) switches to a different
pre-computed Y matrix instead of re-deriving the network equations on the fly.

Build one from an ordinary `TDSProblem`, optionally with its own events, and solve it the
same way:

```julia
tds_prob = TDSProblem(init; tspan=(0.0, 5.0))
ybus_prob = YBusTDSProblem(tds_prob)

# Or with events (a fault at bus 5 cleared by tripping line 5-1)
events = [BusFault(1.0, 1.1, 5; xf = 1e-4), BranchTrip(1.1, 1, 5, "1")]
ybus_prob = YBusTDSProblem(tds_prob; events=events)

sol = solve(ybus_prob)   # the time span comes from tds_prob
```

A fault or switching event gives the same simulation whether it is attached to the
original `TDSProblem` or to the `YBusTDSProblem` built from it.

# Extended help
Every topology the simulation can reach is factored once, when the problem is built, so
that a topology change during the solve is a cheap switch to an already-factored matrix
rather than a new factorization.

- `base_prob::PROB`: the original `TDSProblem` this was built from.
- `sys_ws::WS`: the buffers the Y-bus residual evaluation reads and writes.
- `network_cache::CACHE`: the pre-computed Y matrices and their factorizations, one per
  topology the events can reach.
- `event_schedule`: the pre-computed timeline of topology changes.
- `M`, `M_inv_diag`: the DAE mass matrix and its inverse diagonal.
- `algebraic`: the Newton solve of the models' own algebraic rows.
- `automata`, `memory_sync`: the running state of each model's limiter/controller mode
  logic, and how it is kept in sync across the differentiation types used internally.
"""
struct YBusTDSProblem{T, PROB, WS<:YBusSystemBuffers{T}, CACHE<:NetworkCache, MT<:AbstractMatrix{T}}
    base_prob::PROB
    sys_ws::WS
    network_cache::CACHE
    event_schedule::Union{Nothing, YBusEventSchedule{T}}
    config::SimConfig

    # Mass matrix for DAE formulation (pre-computed, immutable)
    M::MT                    # Sparse diagonal for DiffEq: M * dx/dt = f
    M_inv_diag::Vector{T}    # Inverse diagonal for trapezoidal: dx/dt = f * M_inv

    # The Newton solve of the models' algebraic rows (`solve_algebraic!`)
    algebraic::AlgebraicNewton
    automata::Union{Nothing, AutomatonRuntime}
    # What copies the memory of the automata into the workspaces of the Dual types, which the
    # integration that runs them sets (`sync_automaton_memory!` of its DAE residual)
    memory_sync::Base.RefValue{Function}
end

"""
    new_buffers(prob::YBusTDSProblem, ::Type{T}) -> YBusSystemBuffers{T}

Buffers of number type `T` for an evaluation of the equations of `prob`: built from its
initialization, sharing the automata of its buffers (`prob.sys_ws`), with the setpoints
`prob.sys_ws` has now (an event may have changed them since the initialization).
"""
function new_buffers(prob::YBusTDSProblem, ::Type{T}) where T
    bp = prob.base_prob
    sys_ws = YBusSystemBuffers(T, bp.address_manager, bp.sys, PowerInjection(), bp.init;
                               automata = automaton_states(prob.sys_ws))
    copy_setpoints!(sys_ws, prob.sys_ws)
    return sys_ws
end

"""
    _algebraic_pairs(metadata) -> Vector{Tuple{Symbol, Symbol}}

The algebraic rows of a model that `solve_algebraic!` updates, each with its variable:
`(residual, variable)` for each internal algebraic residual `g_<name>` whose model has an
internal algebraic variable `<name>`.
"""
function _algebraic_pairs(metadata)
    variables = Set(var.name for var in metadata.vars if is_algeb_var(var) && !is_external(var))
    found = Tuple{Symbol, Symbol}[]
    for res in metadata.residuals
        (is_algeb_res(res) && !is_external(res)) || continue
        name = String(res.name)
        startswith(name, "g_") || continue
        var = Symbol(name[3:end])
        var in variables && push!(found, (res.name, var))
    end
    return found
end

"""
    YBusTDSProblem(tds_prob; events=[])

Create a Y-bus TDS problem from a standard TDS problem.

# Arguments
- `tds_prob`: a `TDSProblem` with its dynamic models already initialized; its solver
  settings (from `setup_tds`) carry over
- `events`: disturbance events (default: empty), applied together with any already on
  `tds_prob`

Events are processed at construction time, and every Y-matrix topology the events can
reach is pre-computed, so a topology change during the simulation is a cheap switch to
an already-factored matrix rather than a new factorization.

# Supported Event Types
- `BusFault(time_on, time_off, bus, id)`: a bus fault
- `BranchTrip`, `BranchClose`: line switching
- `FieldModification`: a change of a field of a model (a status field, for a Branch or a
  Fault only)
- `SetpointChange`: a change of a controller's reference value; any limiter or
  operating-mode logic that depends on it is re-evaluated against the new value

Any other event type is refused by name before any topology is built. A fault attached
to `tds_prob` or to the Y-bus problem built from it gives the same simulation.

# Returns
YBusTDSProblem ready for simulation with the Y-bus formulation

# Example
```julia
# Simple simulation without events
ybus_prob = YBusTDSProblem(tds_prob)


# With a fault event (pre-fault, fault-on and post-fault topologies, all pre-computed)
events = [BusFault(1.0, 1.1, 5; xf = 1e-4), BranchTrip(1.1, 1, 5, "1")]  # Fault at bus 5, clear with line 1-5 trip
ybus_prob = YBusTDSProblem(tds_prob; events=events)
```
"""
function YBusTDSProblem(tds_prob; events::Vector{<:Disturbance}=Disturbance[])
    T = Float64
    # The problem owns its system: the fault models of the events go into a copy of
    # tds_prob's system, so they reach neither tds_prob nor other problems built from it
    tds_prob = working_copy(tds_prob)
    sys = tds_prob.sys
    am = tds_prob.address_manager
    # The events of the problem: those of tds_prob (expanded there: a BusFault is a fault model of
    # its system and the changes of its status) and `events`
    events = Disturbance[tds_prob.events; events]

    config = tds_prob.config

    # Create Y-bus system buffers (metadata-driven)
    # (setpoints E_q0 = E_q(t=0), P_m0 = P_m(t=0), etc. from the TDS initialization)
    sys_ws = YBusSystemBuffers(T, am, sys, PowerInjection(), tds_prob.init)

    # Bus voltages from the power flow (initial network iterate), by bus number
    pf = tds_prob.init.pf
    pf_index = pf_snapshot_positions(pf)
    buses = sys.models.Bus.i
    V_pf = Vector{Complex{T}}(undef, length(buses))
    for i in eachindex(buses)
        k = pf_index[buses[i]]
        V_pf[i] = pf.v[k] * cis(pf.theta[k])
    end

    # Create network cache (Y matrices, solver buffers); the power flow voltages
    # are the initial network iterate
    network_cache = NetworkCache(sys; T=T, V_pf=V_pf)

    # The internal voltage gatherers of sys_ws write the network's generator order
    for (k, (model_name, model_idx)) in enumerate(network_cache.E_buf.gen_model_map)
        sys_ws.model_to_gen_idx[model_name][model_idx] == k ||
            error("Generator $model_name[$model_idx] is network generator $k but has index " *
                  "$(sys_ws.model_to_gen_idx[model_name][model_idx]) in the system buffers")
    end
    n_machines = sum(v -> count(>(0), v), values(sys_ws.model_to_gen_idx); init=0)
    n_machines == length(network_cache.E_buf.gen_model_map) ||
        error("$n_machines in-service machines in the system buffers, " *
              "$(length(network_cache.E_buf.gen_model_map)) network generators")

    # Process events and pre-compute all required topologies
    event_schedule = nothing
    if !isempty(events)
        PowerfulTDS.check_events(sys, events)  # (an event that cannot take effect is refused before any work)
        # Process events → pre-compute all topologies from the cache's base Y
        topology_cache = network_cache.topology_cache
        event_schedule = process_ybus_events(events, sys, topology_cache.base, topology_cache, config)

        if !isnothing(event_schedule)
            @info "Pre-computed $(length(event_schedule.times)) event topologies"
            # A setpoint change names a model of the problem, an instance and a setpoint of it
            for change in event_schedule.changes
                change isa SetpointChange || continue
                get_index_by_key(sys, change.model, getfield(change, :key))
                setpoint_array(sys_ws, change.model, change.setpoint)
            end
        end
    end

    # Pre-compute mass matrix and its inverse (used by both DiffEq and trapezoidal solvers)
    M, M_inv_diag = create_mass_matrices(am, sys, T)

    # Pre-build algebraic var→res mapping for solve_algebraic!
    algebraic = AlgebraicNewton(am, sys, sys_ws)

    n_diff = Int(max_address_in_space(am, StateVar()))
    n_alg = Int(max_address_in_space(am, AlgebVar()))
    memory_sync = Ref{Function}(() -> nothing)
    automata = automaton_runtime(sys_ws, sys, sys_ws.status, n_diff, n_alg, () -> memory_sync[](), Vector(diag(M)))
    prob = YBusTDSProblem{T, typeof(tds_prob), typeof(sys_ws), typeof(network_cache), typeof(M)}(
        tds_prob,
        sys_ws,
        network_cache,
        event_schedule,
        config,
        M,
        M_inv_diag,
        algebraic,
        automata,
        memory_sync
    )
    if automata !== nothing
        u0 = zeros(T, n_diff + n_alg)
        init_tds_buffers!(view(u0, 1:n_diff), view(u0, (n_diff + 1):(n_diff + n_alg)), tds_prob)
        start!(automata, view(u0, 1:n_diff), view(u0, (n_diff + 1):(n_diff + n_alg)), T(tds_prob.tspan[1]))
    end
    return prob
end

"""
    create_mass_matrices(am, sys, ::Type{T}) -> (M, M_inv_diag)

Create mass matrix and its inverse diagonal for DAE formulation.

For DAE: M * dx/dt = f
- DiffEq.jl uses sparse M directly
- Trapezoidal uses M_inv_diag: dx/dt = f * M_inv_diag

# Returns
- `M`: Sparse diagonal matrix (n_total × n_total) for DiffEq
- `M_inv_diag`: Vector{T} of length n_diff for trapezoidal solver

# Safe Inverse
If M[i] = 0 for a differential state (algebraic constraint in differential slot),
sets M_inv[i] = 0 with warning (treats as frozen variable).
"""
function create_mass_matrices(am, sys, ::Type{T}) where T
    n_diff = Int(max_address_in_space(am, StateVar()))
    n_alg = Int(max_address_in_space(am, AlgebVar()))
    n_total = n_diff + n_alg

    # Initialize mass diagonal: 1.0 for differential, 0.0 for algebraic
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

            addrs = get_addresses(am, model_name, var_meta.name)
            mass_values = eval_mass_coeff(mass_coeff, models)
            for (i, addr) in enumerate(addrs)
                M_diag[addr] = mass_values[i]
            end
        end
    end

    # Create sparse diagonal matrix for DiffEq
    M = spdiagm(M_diag)

    # Create inverse diagonal for trapezoidal (differential states only)
    M_inv_diag = Vector{T}(undef, n_diff)
    for i in 1:n_diff
        if M_diag[i] == zero(T)
            @warn "Mass coefficient M[$i] = 0 for differential state - treating as frozen"
            M_inv_diag[i] = zero(T)
        else
            M_inv_diag[i] = one(T) / M_diag[i]
        end
    end

    return M, M_inv_diag
end

"""
    get_network_cache(prob::YBusTDSProblem)

Access the network cache for topology switching and inspection.
"""
get_network_cache(prob::YBusTDSProblem) = prob.network_cache

"""
    get_topology_cache(prob::YBusTDSProblem)

Access the topology cache for adding/switching topologies.
"""
get_topology_cache(prob::YBusTDSProblem) = prob.network_cache.topology_cache

export YBusTDSProblem
export get_network_cache, get_topology_cache
