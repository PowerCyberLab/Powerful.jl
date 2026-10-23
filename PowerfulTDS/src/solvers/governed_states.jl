# ============================================================================
# Governed States: A Mode That Holds a Differential State on a Moving Relation
# ============================================================================
#
# `Mode(:upper; states = [:y => :bound])` governs the differential state `y` while the model is in the
# mode: `y` is held on the right side `bound`, a number or a value the model computes
# (`calc_automaton_values!`), which may move with the states it reads. The generator lowers a bounded state to
# such modes (`bounded`), the bound being the right side. The model's rate for
# the state in the mode is the rate of the right side, `Σ_j (∂r/∂x_j) f_j`, gated by what the mode holds
# (the generator emits it; nothing here changes the rate), so the integration keeps the state on the
# relation: the relation is affine in differential states with constant coefficients, a linear invariant
# the multistep methods preserve, as the state-only constraints of index_reduction.jl.
#
# What is Powerful's:
#
# - Admissibility (`check_governed`), checked here and not left to the generator, from what the values
#   read (`ValueReads`, traced) and their degree (`PolyDegree`): a right side must be affine in the
#   variables (degree at most 1, so that its coefficients read no variable and no time), read only
#   differential states with mass of its own model, none of them governed by any mode of the model, nor
#   the governed state itself, and no discrete quantity, memory or time. Parameters are allowed, finite.
# - The entry (`_enter_governed!`): where a transition enters a mode that governs a state, the state is
#   written the right side at that instant, along that state alone (`project_states!`, `AlongStates`; the
#   states the right side reads do not move). A crossing located on `x - r` lands on the relation only to
#   the located root, and the write is made all the same.
# - What follows from the write (`_state_written!`): the values that read the state lose their sides
#   (they moved by a jump, as by a discrete quantity that a transition changed, `_quantity_changed!`); an
#   algebraic residual that reads it (`AlgebraicCheck` `switches`: traced with the governed states as
#   inputs, and the variables of other models that are the state) makes the settle of the instant
#   reinitialize the network before its next round, as a change of the algebraic equations by a mode does.
# - The start (`_check_start_relations!`): a mode that governs a state finds it on its right side or the
#   start is refused. A computed right side cannot move the start, as a computed bound of a
#   `VariableLimit` cannot.
# - The repair after an event (`_repair_writes`, round 19a): `check_governed`'s admissibility (above) proves a
#   model field (a "parameter", read from `components(model)`) is the only *continuous-state-free* channel that
#   can move a governed right side discontinuously -- a continuous input (another state, read through its own
#   stated rate) can only move it continuously, which the integration already tracks to its own tolerance
#   (`moving`, above). A mode a transition *enters* this instant is already on its relation, from
#   `_enter_governed!`'s own write; a mode that *stays active* through the instant is not touched by that write
#   at all, so an event that moves its right side leaves the state exactly where it was until this repair puts
#   it back. The order at an event (the code generator's own contract, not Powerful's alone): the event's data is
#   replaced first (`apply!`), then every guard is read on the retained state in the settle's first round
#   (ordinary transitions fire, including a release that already reads the new data, e.g. REGCA1's own position
#   atoms); in that *same* round, beside the entries transitions just made, every state whose mode stayed active
#   and whose right side the event touched is put back on its relation, through `_enter_governed!`'s own write
#   (`project_states!`, `_state_written!`, `_states_written!`); only then does the settle continue in ordinary
#   rounds. Which states need repair is found exactly (`g.param_reads[i]`, traced once when `Governed` is built,
#   `_value_param_reads`, tested against `touched`, the fields an event actually wrote): never by closeness to
#   the relation (`_off_relation` is the start's own check, and would also flag ordinary drift between events,
#   which is not what an event touched).
#
#   What writes `touched`, and so can be repaired: a `FieldModification` and a `SetpointChange` (`events.jl`'s
#   `to_callback` methods for each push to it; a setpoint is a column of the model like any other, which
#   `_value_param_reads` traces the same way). Simultaneous `FieldModification`s at one time each apply and
#   settle in sequence, and only the last of them may repair, reading every one of their fields together, not
#   the first alone (`AutomatonSettle.repair`, `PowerfulTDS/src/solvers/automata.jl`); a `SetpointChange`
#   simultaneous with another event at the same time is not grouped this way, so a relation that reads fields
#   from both would be repaired from the first alone -- narrower than the `FieldModification` guarantee, since
#   this is not a case any model in this codebase exercises today. A `Disturbance` subtype that resets a
#   differential state directly (bypassing `apply!`, e.g. a fault that reinitializes states through
#   `address_manager`) never writes `touched` at all: `check_governed` admits a right side that reads *another*
#   (non-governed) state through its own stated rate, and such a reset would move that state discontinuously
#   without the repair ever seeing it. No model in this codebase reads a governed right side off a state a
#   `Disturbance` resets, so this is a known, undemonstrated gap, not a guarantee this code provides.
#   `_value_param_reads` only ever traces a value's dependence on its *own* instance's fields (not another
#   instance's, as `ValueReads.cross`/`trace_value_reads` tracks for variable reads): it throws, rather than
#   silently under-reporting `param_reads`, if it ever finds one that does (no value in this codebase reads a
#   field this way today). It also only traces `AbstractFloat` fields, not a `Bool`/`Int` field such as a status
#   flag read as a branch condition whose own branches are each affine: no right side in this codebase reads a
#   non-`Float` field this way today, so this case is left untraced (not loudly checked, unlike the one above)
#   rather than paying to trace a kind of dependence that does not arise; revisit if a model's right side ever
#   does.
# ============================================================================

using SparseArrays: sparse
using PowerfulCore: is_algeb_var, is_external, governed_states, refresh_affine_governed!

"""
    GOVERNED_START_RTOL

How far a state governed in an initial mode may start from its right side, relative to the larger of 1 and
the right side: the tolerance of the initialization that computed it.
"""
const GOVERNED_START_RTOL = 1e-8

# The variables (indices among the keys of `ws.var_vals`) of the states the modes of `metadata` govern
function _governed_variables(ws::ResidualBuffers, metadata)
    metadata.governs_any_state || return Int[]
    names = collect(keys(ws.var_vals))
    return Int[let j = findfirst(==(s), names)
                   j === nothing && throw(ArgumentError(
                       "$(metadata.name) governs the state $s, which is not a variable of its workspace " *
                       "(variables: $(join(names, ", ")))"))
                   j
               end for s in governed_states(metadata.automata)]
end

# The variables `governed` of a traced workspace from `v`, after its first `offset` entries: a block of one
# entry per instance for each, in order
function _set_governed_inputs!(traced::ResidualBuffers, governed::Vector{Int}, v, offset::Int)
    n = traced.n_instances
    columns = values(traced.var_vals)
    for (q, j) in enumerate(governed)
        copyto!(columns[j], 1, v, offset + (q - 1) * n + 1, n)
    end
    return nothing
end

"""
    _state_readers(sys_ws) -> Union{Nothing, Dict{Int, Vector{Tuple{Symbol, Int, Int}}}}

For every address of a differential state, the variables of the workspaces of `sys_ws` gathered from it:
`(model name, variable index, instance)`. `nothing` where no model of the runtime governs a state, which
needs none of it.
"""
function _state_readers(sys_ws)
    any(name -> sys_ws.metadata_by_name[name].governs_any_state, sys_ws.automaton_models) || return nothing
    readers = Dict{Int, Vector{Tuple{Symbol, Int, Int}}}()
    for (name, ws) in pairs(sys_ws.workspaces_by_name)
        for (j, (addrs, _, is_state)) in enumerate(ws.var_gather_tuple)
            is_state || continue
            for (i, a) in enumerate(addrs)
                a == UNROUTED && continue
                push!(get!(() -> Tuple{Symbol, Int, Int}[], readers, Int(a)), (name, j, i))
            end
        end
    end
    return readers
end

# Per instance, whether a model other than its own gathers a state its modes govern (`nothing` where it governs none)
function _read_by_others(governed::Governed, n::Int)
    isempty(governed.states) && return nothing
    read = fill(false, n)
    for per_instance in governed.foreign, i in 1:n
        isempty(per_instance[i]) || (read[i] = true)
    end
    return read
end

"""
    Governed(name, ws, metadata, model, external, formulation, reads, mass, n_diff, readers) -> Governed

The states the modes of model `name` govern (`GovernedState`), checked for admissibility (`check_governed`):
`ws` its Float64 workspace, `metadata` its metadata, `reads` what its values read (`ValueReads`, with its
massless variables marked), `mass` the diagonal of the mass matrix and `readers` the variables of every
workspace by the state address they gather (`_state_readers`). Empty where no mode governs a state.
"""
function Governed(name::Symbol, ws::ResidualBuffers, metadata, model, external, formulation, reads::ValueReads,
                  mass::AbstractVector, n_diff::Int, readers)
    metadata.governs_any_state || return Governed()
    for a in metadata.automata
        isempty(a.held) && continue
        for h in a.held, mode in a.modes
            mode.name === h.mode && !isempty(mode.states) && throw(ArgumentError(
                "$name: automaton $(a.name) has a held clause for mode $(mode.name), which also governs " *
                "$(join(string.(first.(mode.states)), ", ")): a held mode (no transition of the automaton " *
                "fires there, round 19a's `_round_steps!` gives it no `Transition`, hence no `tr.from` to " *
                "read a farewell write from) combined with a governed state is not implemented -- no " *
                "installed model does this today (checked by exhaustive grep: `Held`/`held =` never appears " *
                "beside `states =` in any model file). Refused loudly rather than silently skipping the " *
                "farewell write or the entry write for this mode."))
        end
    end
    var_names = collect(keys(ws.var_vals))
    value_names = collect(keys(ws.automata.values))
    states = GovernedState[]
    by_mode = [[Int[] for _ in a.modes] for a in metadata.automata]
    foreign = Vector{Vector{Tuple{Symbol, Int, Int}}}[]
    refs = Union{Symbol, Nothing}[]   # the right side's value name (nothing for a number), for param_reads below
    for (ka, a) in enumerate(metadata.automata), (km, mode) in enumerate(a.modes), (s, ref) in mode.states
        j = findfirst(==(s), var_names)
        k = findfirst(v -> v.name === s, metadata.vars)
        (j === nothing || k === nothing) && throw(ArgumentError(
            "$name: mode $(mode.name) of automaton $(a.name) governs $s, which is not a variable of its workspace"))
        addrs, _, is_state = ws.var_gather_tuple[k]
        is_state || throw(ArgumentError(
            "$name: mode $(mode.name) of automaton $(a.name) governs $s, which is not a differential state"))
        addresses = Int[a == UNROUTED ? 0 : Int(a) for a in addrs]
        value, column = ref isa Symbol ? (findfirst(==(ref), value_names), getfield(ws.automata.values, ref)) :
                                         (ref, nothing)
        push!(by_mode[ka][km], length(states) + 1)
        push!(states, GovernedState(ka, km, j, s, addresses, value, column, Set{Symbol}[]))
        push!(refs, ref isa Symbol ? ref : nothing)
        push!(foreign, [Tuple{Symbol, Int, Int}[e for e in get(readers, addr, Tuple{Symbol, Int, Int}[])
                        if e[1] !== name] for addr in addresses])
    end
    recheck! = let name = name, metadata = metadata, ws = ws, reads = reads, states = states, mass = mass,
                   model = model, external = external, formulation = formulation, value_names = value_names
        function ()
            check_governed(name, metadata, ws, reads, states, mass,
                           trace_value_time_reads(ws, model, external, formulation))
            # Round 19a, Phase D: the residual path's affine cache is a function of the model's
            # parameters only (`check_governed`'s own admissibility proof); refresh it here, the same
            # retrace trigger that already re-checks admissibility (setup, and every event that changes
            # a parameter a relation reads), so the hot path never recomputes it itself. The read set
            # (`read_vars`) is recomputed from `reads.pattern` on every retrace too, not only the first
            # time: a retrace can change what `reads` itself holds (`check_governed`'s own docstring), and
            # this reuses that trace rather than keeping a possibly-stale one from the first call.
            read_vars = Dict{Symbol, Vector{Int}}(
                ref => _value_read_variables(ws, reads, findfirst(==(ref), value_names))
                for ref in keys(ws.automata.affine))
            refresh_affine_governed!(ws, model, external, formulation, read_vars)
            return nothing
        end
    end
    recheck!()
    # The right side is admissible (recheck! above did not throw): trace which of the model's own parameters each
    # value reads (`_value_param_reads`), once per distinct value (states that share a right side share its
    # trace), and rebuild each state with the result (`GovernedState` is immutable; `states` itself is not, and
    # `recheck!` above closed over this same vector, so a later recheck! still sees these same, now-completed,
    # states). Traced once here, not refreshed on every later `recheck!`: a right side's kernel branching on a
    # parameter (admissible, as long as each branch is itself affine in the variables) could in principle change
    # which fields it reads after an event retraces it, but no model in this codebase does this today, and
    # retracing it on every event a touched model receives would cost a repeated Jacobian-sparsity trace for a
    # case that does not arise; revisit this if a model's right side ever does branch on a parameter.
    cache = Dict{Symbol, Vector{Set{Symbol}}}()
    for (idx, g) in enumerate(states)
        ref = refs[idx]
        param_reads = ref === nothing ? [Set{Symbol}() for _ in 1:ws.n_instances] :
            get!(() -> _value_param_reads(ws, model, external, formulation, ref), cache, ref)
        states[idx] = GovernedState(g.automaton, g.mode, g.var, g.state, g.addresses, g.value, g.column, param_reads)
    end
    return Governed(states, by_mode, foreign, recheck!)
end

"""
    _governed_read_variables(ws::ResidualBuffers, reads::ValueReads, states::Vector{GovernedState}) -> Vector{Int}

The variables (indices among the keys of `ws.var_vals`) that the relation of some admissible governed state of
`states` reads (`reads.pattern`, restricted to the values those states' relations are, the same scan
`check_governed` makes to find what each relation reads): the union across every mode of the model that
governs a state, not only whichever mode `reads` happened to be traced in. A single trace of the values
(`trace_value_reads`) runs with the workspace in one mode -- typically the model's initial one, which governs
nothing -- so it alone would miss a variable a *different* mode's relation reads (REGCA1's `Vflt`, read only
by `Ip_upper_1`, the relation of modes `upper_0`/`upper_1`, never read while `:free`). Rather than retracing
the model once per mode of every governing automaton (which would need writing `km` into the shared, live
`ws.automata.modes` for the trace and restoring it after -- unsafe if the trace throws mid-way, corrupting the
live mode of a workspace other code already reads), this reuses exactly what `check_governed`'s own
admissibility proof already traced for every mode's relation, read once: nothing is retraced here.
"""
function _governed_read_variables(ws::ResidualBuffers, reads::ValueReads, states::Vector{GovernedState})
    wanted = Set(g.value for g in states if g.value isa Int)
    isempty(wanted) && return Int[]
    n = ws.n_instances
    n_values = length(ws.automata.values)
    cols = Set{Int}()
    rows = rowvals(reads.pattern)
    for c in 1:size(reads.pattern, 2), p in nzrange(reads.pattern, c)
        _, v = fldmod1(rows[p], n_values)
        v in wanted || continue
        j, _ = fldmod1(c, n)
        push!(cols, j)
    end
    return sort!(collect(cols))
end

"""
    _value_read_variables(ws::ResidualBuffers, reads::ValueReads, value_index::Int) -> Vector{Int}

The variables (indices among the keys of `ws.var_vals`) that a single value of `ws.automata.values`
(`value_index`, its position there) reads: the same trace `_governed_read_variables` reads (`check_governed`'s
own sparsity pattern, `reads.pattern`), restricted to one value instead of the union over every governed
state's value. What `refresh_affine_governed!` (`PowerfulCore`) needs per relation to store only its nonzero
coefficients (round 19a, Phase D's coordinator follow-up): reusing this trace, not a second one.
"""
function _value_read_variables(ws::ResidualBuffers, reads::ValueReads, value_index::Int)
    n = ws.n_instances
    n_values = length(ws.automata.values)
    cols = Set{Int}()
    rows = rowvals(reads.pattern)
    for c in 1:size(reads.pattern, 2), p in nzrange(reads.pattern, c)
        _, v = fldmod1(rows[p], n_values)
        v == value_index || continue
        j, _ = fldmod1(c, n)
        push!(cols, j)
    end
    return sort!(collect(cols))
end

"""
    trace_value_time_reads(ws::ResidualBuffers, model, external, formulation) -> Vector{Bool}

For each value of the automata of `model`, whether it reads the time (of any instance): traced from its
`calc_automaton_values!` with the time as the one input, as `trace_value_reads` traces the variables.
"""
function trace_value_time_reads(ws::ResidualBuffers, model, external, formulation)
    n, n_values = ws.n_instances, length(ws.automata.values)
    trace = function (at_values::Bool)
        values! = function (out, v)
            T = eltype(v)
            traced = ResidualBuffers(ws, T; setpoints_at_values = at_values)
            foreach(column -> fill!(column, zero(T)), values(traced.var_vals))
            calc_automaton_values!(traced, at_values ? model : with_float_columns(model, T),
                                   at_values ? external : _float_external(external, T), formulation, v[1])
            for (k, column) in enumerate(values(traced.automata.values)), i in 1:n
                out[(i - 1) * n_values + k] = column[i]
            end
            return nothing
        end
        return jacobian_sparsity(values!, zeros(n_values * n), [0.0], TracerSparsityDetector())
    end
    pattern = try
        trace(false)
    catch e
        _branches_on_tracer(e) || rethrow()
        trace(true)
    end
    return [any(pattern[(i - 1) * n_values + k, 1] for i in 1:n) for k in 1:n_values]
end

"""
    check_governed(name, metadata, ws, reads::ValueReads, states, mass, time_reads)

Throw an `ArgumentError` naming the mode, the state and the reason unless every right side of the governed
`states` of model `name` is admissible (see the header):

- the state is a differential state with mass at every instance (`mass` the diagonal of the mass matrix);
- a right side that is a value is a polynomial of degree at most 1 in the variables (`reads.degree`, traced
  with `PolyDegree`), so that it is affine with coefficients that read neither a variable nor the time (only
  the parameters, which are finite where it starts, `_check_start_relations!`);
- it reads no discrete quantity or memory (`reads.state`: a coefficient that was one would evolve), nor the
  time (`time_reads`, `trace_value_time_reads`), nor an algebraic variable, a variable without mass or a
  quantity of another model (`reads.massless`, the variables' kinds);
- it reads neither the governed state itself nor any state governed by a mode of the model: a state
  governed by another mode or automaton moves by the relation that governs it, not by its equation.

`reads` is shared with the model, so that a check after an event that changed the parameters on which the
kernel of the values branches (`retrace_automata!`) sees what they read then.
"""
function check_governed(name::Symbol, metadata, ws::ResidualBuffers, reads::ValueReads, states::Vector{GovernedState},
                        mass::AbstractVector, time_reads::AbstractVector{Bool})
    n = ws.n_instances
    var_names = collect(keys(ws.var_vals))
    value_names = collect(keys(ws.automata.values))
    n_values = length(value_names)
    quantities = (keys(ws.automata.discrete)..., keys(ws.automata.memory)...)
    governed_vars = Set(g.var for g in states)
    wanted = Set(g.value for g in states if g.value isa Int)
    read_variables = Dict(v => Set{Tuple{Int, Int}}() for v in wanted)
    read_quantities = Dict(v => Set{Int}() for v in wanted)
    for (matrix, into, per) in ((reads.pattern, read_variables, true), (reads.state, read_quantities, false))
        rows = rowvals(matrix)
        for c in 1:size(matrix, 2), p in nzrange(matrix, c)
            _, v = fldmod1(rows[p], n_values)
            v in wanted || continue
            if per
                push!(into[v], Tuple(fldmod1(c, n)))
            else
                push!(into[v], fldmod1(c, n)[1])
            end
        end
    end
    metadata_var(j) = metadata.vars[findfirst(v -> v.name === var_names[j], metadata.vars)]
    automata = metadata.automata
    for g in states
        a = automata[g.automaton]
        what = "$name: mode $(a.modes[g.mode].name) of automaton $(a.name) governs $(g.state)"
        for (i, address) in enumerate(g.addresses)
            (address != 0 && !iszero(mass[address])) || throw(ArgumentError(
                "$what, which has no mass at instance $i: a governed state is a differential state with mass " *
                "(a state without mass is a row of the DAE without mass, an algebraic equation)"))
        end
        g.value isa Int || continue
        by = "$what by $(value_names[g.value])"
        reasons = String[]
        degree = reads.degree[g.value]
        degree >= NOT_POLYNOMIAL ? push!(reasons, "is not a polynomial in the variables, so not affine") :
        degree > 1 ? push!(reasons, "is not affine in the variables (a polynomial of degree $degree: its " *
                                    "coefficients depend on variables or the time)") : nothing
        for (j, i) in sort!(collect(read_variables[g.value]))
            variable = metadata_var(j)
            reason = j == g.var ? "reads the governed state $(g.state) itself" :
                     j in governed_vars ? "reads $(variable.name), a state governed by a mode of the model " *
                                          "(governed by another mode or automaton)" :
                     is_external(variable) ? "reads $(variable.name), a quantity of another model" :
                     is_algeb_var(variable) ? "reads the algebraic variable $(variable.name)" :
                     reads.massless[(j - 1) * n + i] ? "reads $(variable.name), a variable without mass at instance " *
                                                      "$i (an algebraic variable)" : nothing
            reason === nothing || reason in reasons || push!(reasons, reason)
        end
        for q in sort!(collect(read_quantities[g.value]))
            push!(reasons, "reads the $(q <= length(ws.automata.discrete) ? "discrete quantity" : "memory") " *
                           "$(quantities[q]), which a transition changes (a coefficient that evolves)")
        end
        time_reads[g.value] && push!(reasons, "reads the time")
        isempty(reasons) || throw(ArgumentError(
            "$by, which $(join(reasons, "; ")). A governed state's right side must be affine in the differential " *
            "states with constant coefficients (parameters), and read no algebraic variable, quantity of another " *
            "model, discrete quantity, memory, the time, or state a mode of the model governs"))
    end
    return nothing
end

# ----------------------------------------------------------------------------
# The start
# ----------------------------------------------------------------------------

# The right side of the governed state `g` at instance `i`: its number, or its value there
_right_side(g::GovernedState, i::Int) = g.column === nothing ? g.value : g.column[i]

# Whether the state `x` is off its right side `r` by more than the start allows
_off_relation(x, r) = abs(x - r) > GOVERNED_START_RTOL * max(1.0, abs(r))

function _off_the_relation(rt::AutomatonRuntime, m::AutomatonModel, i::Int, g::GovernedState, x, r, t)
    a = m.automata[g.automaton].automaton
    throw(ErrorException(
        "Automaton start: $(m.name) ($(join(string.(m.keys[i]), ", "))), automaton $(a.name) in mode " *
        "$(a.modes[g.mode].name) governs $(g.state) by $(g.value isa Int ? keys(m.state.values)[g.value] : "the number $(g.value)"), " *
        "but $(g.state) = $x is off its relation $(g.state) = $r (residual $(x - r)) at t = $t. A start off its " *
        "relation is refused: the right side is computed from the model and cannot be moved to the start " *
        "(a computed bound of a `VariableLimit` is refused likewise)."))
end

"""
    _check_start_relations!(rt::AutomatonRuntime, x, t)

At the start of the automata, with the rounds of the start over (the settled point) and the values computed at
`(x, y, t)`: throw where an in-service instance starts in a mode that governs a state that is off its right side, naming the
state, the relation and the residual (`_off_the_relation`); and where a coefficient of a right side is not
finite (`x` shifted by one: every right side is affine, so a coefficient that is not finite makes it not so
there).
"""
function _check_start_relations!(rt::AutomatonRuntime, x, y, t)
    for (km, m) in enumerate(rt.models)
        isempty(m.governed.states) && continue
        for (ka, a) in enumerate(m.automata), i in 1:n_instances(m)
            m.active[i] || continue
            for k in m.governed.by_mode[ka][a.modes[i]]
                g = m.governed.states[k]
                r = _right_side(g, i)
                isfinite(r) || throw(ErrorException(
                    "Automaton start: the right side of $(g.state) governed by $(m.name) " *
                    "($(join(string.(m.keys[i]), ", "))) is $r at t = $t"))
                _off_relation(x[g.addresses[i]], r) && _off_the_relation(rt, m, i, g, x[g.addresses[i]], r, t)
            end
        end
        _check_coefficients!(rt, km, x, y, t)
    end
    return nothing
end

# The right sides of the governed states of model `km` are finite with every differential state shifted by
# one (the coefficients of an affine function that are not finite make it not finite there), and the values
# of the model are computed again at the state itself
function _check_coefficients!(rt::AutomatonRuntime, km::Int, x, y, t)
    m = rt.models[km]
    columns = [g.column for g in m.governed.states if g.column !== nothing]
    isempty(columns) && return nothing
    shifted = x .+ 1.0
    rt.evaluators[km](shifted, y, t)
    for g in m.governed.states, i in 1:n_instances(m)
        (g.column === nothing || !m.active[i]) && continue
        isfinite(g.column[i]) || throw(ErrorException(
            "Automaton start: a coefficient of the right side $(keys(m.state.values)[g.value]) of $(g.state), governed " *
            "by $(m.name) ($(join(string.(m.keys[i]), ", "))), is not finite (the right side is $(g.column[i]) with " *
            "every differential state shifted by one): a governed state's parameter-dependent coefficients must " *
            "be finite"))
    end
    rt.evaluators[km](x, y, t)
    return nothing
end

# ----------------------------------------------------------------------------
# The entry
# ----------------------------------------------------------------------------

"""
    _repair_writes(rt::AutomatonRuntime) -> Vector{Tuple{Int, Int, Int}}

The `(model, instance, governed state)` triples of the states whose mode stayed active through this round (every
in-service instance, in the mode its automaton is currently in, `by_mode[ka][a.modes[i]]`) and whose right side a
fixed-time event's parameters just changed: `g.param_reads[i]`, the fields the right side reads, intersects
`rt.touched`, the exact `(model, field, instance)` triples this instant's event(s) wrote (no tolerance: an event
that writes none of a state's `param_reads[i]` leaves it exactly where it was). `check_governed`'s admissibility
already proves a parameter is the only channel an event can move a governed right side through (no discrete
quantity, memory, algebraic variable or time; a state cannot itself be a `FieldModification` target), so this
is the relation's entire dependence on anything an event could touch, found exactly: `g.param_reads[i]` is traced
once, when the model's `Governed` is built (`_value_param_reads`), and this is a plain set-membership test, not a
numerical closeness check (`_off_relation`, the start's own check, would also catch ordinary drift between
events, which is not what an event touched). Empty where `rt.touched` is empty (most events touch no governed
relation at all). The caller gates this on `event` (a located instant never jumps a right side) and the first
round of the instant only (a later round's repair, if any were still needed, would already show up as an
ordinary transition: see the header of this file).
"""
function _repair_writes(rt::AutomatonRuntime)
    isempty(rt.touched) && return Tuple{Int, Int, Int}[]
    writes = Tuple{Int, Int, Int}[]
    for (km, m) in enumerate(rt.models)
        isempty(m.governed.states) && continue
        for (ka, a) in enumerate(m.automata), i in 1:n_instances(m)
            m.active[i] || continue
            for k in m.governed.by_mode[ka][a.modes[i]]
                g = m.governed.states[k]
                isempty(g.param_reads[i]) && continue
                any(f -> (m.name, f, i) in rt.touched, g.param_reads[i]) && push!(writes, (km, i, k))
            end
        end
    end
    return writes
end

"""
    _enter_governed!(rt, steps, stepped, x, y, t, starting, jump, integrator; event=false)
        -> Union{Nothing, NamedTuple{(:readers, :partial)}}

After a round of transitions (`steps` of the instances `stepped`, their modes entered, the discrete quantities
following and the values computed): write each state a mode entered governs its right side, at once and
along that state alone (`project_states!`, `AlongStates`), the others not moving. With `event` (the first round
of a fixed-time event's instant, from the caller, `settle!`), the states whose mode stayed active through the
round and whose right side the event just moved (`_repair_writes`) are added to the same write, through the
same path, beside the states a transition of this round entered: a transition that both entered and was
touched is written once (the repair's own write would recompute the identical value `_enter_governed!`'s own
entry already wrote, an idempotent no-op in effect, so there is nothing to gain, and a risk of reading
`before`/`at` twice at the same address, from writing it twice; see the header of this file). Returns `nothing`
where no state needed writing, else `readers`, the instances (`(model, instance)` pairs) whose values read
what was written (they take part in the next round), and `partial`, whether an algebraic equation reads it,
so that the settle of the instant must reinitialize the network before the next round
(`AutomatonRuntime` `settle_partial`).

What the write changes, in this order: the values that read a state lose their sides
(`_state_written!`: the state is exactly where its relation puts it, whether or not the write moved it, and the
values are decided by what they are, not by the side of a crossing that led there) and are computed again; an augmentation that carries derivatives is told, to set the
sensitivities of the written states by the reset map (`_states_written!`), whether or not the state moved: a
crossing located on the relation lands on it to the root's tolerance only, and a repair's write (no located
crossing, no switching time: see a trajectory-sensitivity augmentation's own handling of this, noted in the
design) the same way.
With `starting`, nothing is written (and no repair is attempted: `event` is never true at a start): the start
finds the state on its right side or refuses it, in the modes it settles in (`_check_start_relations!`, after
the rounds).
"""
function _enter_governed!(rt::AutomatonRuntime, steps, stepped, x, y, t, starting::Bool, jump, integrator;
                          event::Bool = false)
    writes = Tuple{Int, Int, Int}[]   # (model, instance, governed state)
    # Which (model, instance, address) a write already targets, so a later write (farewell, repair) for the
    # same address is dropped rather than added again: two `GovernedState` entries (e.g. a governing mode's
    # own and the mode it departed to/from) can name the *same* state, and `project_states!` cannot take two
    # different values for one address in one projection.
    written_addr = Set{Tuple{Int, Int, Int}}()
    _write_address(km, i, k) = (m = rt.models[km]; m.governed.states[k].addresses[i])
    for ((km, i), (m, _, ka, _, target)) in zip(stepped, steps)
        isempty(m.governed.states) && continue
        for k in m.governed.by_mode[ka][target]
            push!(writes, (km, i, k))
            push!(written_addr, (km, i, _write_address(km, i, k)))
        end
    end
    # The farewell write (round 19a): a mode a transition *leaves* (`tr.from`) may have governed a state
    # the entry write above does not touch (the entered mode governs a different state, or none, as every
    # transition out of REGCA1's `upper_0`/`upper_1` does, entering `:free`, which governs nothing). Land
    # it on the departed mode's relation, read now -- before this round's `evaluate_values!` below
    # recomputes anything, i.e. at what the relation actually was at the moment of release -- one instant
    # before the entered mode's own, unrelated rate takes over (file header, design `3). A held-clause
    # step has no transition (`tr === nothing`, `_round_steps!`): `Governed`'s own construction already
    # refuses a held mode that governs a state, so this is simply skipped, not silently wrong, for any
    # model this runtime could hold. The entry write above wins a conflict on the same *address* (added
    # first; a farewell write for an address a write already targets is dropped), the same tolerance
    # `_repair_writes`'s own addition below relies on.
    #
    # Never where `event` is this round's qualifying fixed-time event (the first round only,
    # `settle!`'s own `do_repair`, passed through as this parameter): the event's data is already
    # applied by the time this round's values were computed (the header of this file: "the event's
    # data is replaced first, then every guard is read"), so `_right_side` here would read the
    # relation *after* the jump, not "what the departed mode's relation actually was at the moment of
    # release" (the farewell write's whole point, `3 of the design) -- for PIReset, landing `xi` on
    # `upper - Kp * u_new` instead of leaving it at its own retained value is exactly the counterexample
    # the release is supposed to produce (`y` continuing from the release point, not re-pinned to the
    # new bound). An ordinary release, with no event this round, has no such "after the jump" reading to
    # avoid: the relation `_right_side` reads here is the same one this round's own values just used to
    # decide the release, which is what "at the moment of release" means there.
    for ((km, i), (m, _, ka, tr, _)) in zip(stepped, steps)
        (isempty(m.governed.states) || tr === nothing || event) && continue
        for k in m.governed.by_mode[ka][tr.from]
            addr = (km, i, _write_address(km, i, k))
            addr in written_addr && continue
            push!(writes, (km, i, k))
            push!(written_addr, addr)
        end
    end
    if event
        for triple in _repair_writes(rt)
            (km, i, k) = triple
            addr = (km, i, _write_address(km, i, k))
            addr in written_addr && continue
            push!(writes, triple)
            push!(written_addr, addr)
        end
    end
    isempty(writes) && return nothing
    n = length(writes)
    right, at = Vector{Float64}(undef, n), Vector{Int}(undef, n)
    for (q, (km, i, k)) in enumerate(writes)
        m = rt.models[km]
        g = m.governed.states[k]
        r = _right_side(g, i)
        isfinite(r) || throw(ErrorException(
            "Automaton settle at t = $t: the right side of $(g.state) governed by $(m.name) " *
            "($(join(string.(m.keys[i]), ", "))) is $r where its mode is entered or repaired after an event"))
        right[q], at[q] = r, g.addresses[i]
    end
    # (a start writes nothing, and checks the relations of the modes it settles in once the rounds are over:
    # `_check_start_relations!`)
    starting && return nothing
    before = x[at]
    project_states!(x, sparse(1:n, 1:n, ones(n), n, n), right, AlongStates(at))
    rt.projections[] += n
    readers = Tuple{Int, Int}[]
    partial = false
    for (q, (km, i, k)) in enumerate(writes)
        m = rt.models[km]
        # (the values that read the state are decided by their values from here on, which the relation puts
        # exactly where it says, moved or not: a side a crossing gave them, as it left them, is no longer theirs)
        _state_written!(rt, readers, km, i, m.governed.states[k], m.governed.foreign[k][i])
        partial |= x[at[q]] != before[q] && m.algebraic.switches[i]
    end
    jump === nothing || _states_written!(jump, integrator, rt, writes, t)
    # The values that read what moved, and the models of the instances that stepped, at the written state
    unique!(readers)
    evaluate_values!(rt, vcat([(km, i) for (km, i, _) in writes], readers), x, y, t)
    return (readers = readers, partial = partial)
end

"""
    _state_written!(rt, readers, km, i, g::GovernedState, foreign)

A transition wrote the governed state `g` of instance `i` of model `km`: the values that read it,
of the model (`ValueReads` `pattern`) and of the models that gather it (`foreign`: `(model name, variable,
instance)`), moved by that, a jump of the value and not its motion along the trajectory. Each that has a side
(`AutomatonModel` `side`: located at this instant, or left with one by an earlier) loses it, so that it is
decided by its value from here on (`atom_holds`), as where a transition changes a discrete quantity or memory
that a value reads (`_quantity_changed!`). Adds the instances of those values to `readers`.
"""
function _state_written!(rt::AutomatonRuntime, readers, km::Int, i::Int, g::GovernedState, foreign)
    lose!(k, m, instance, variable) =
        _lose_readers!(m, m.reads.pattern, (variable - 1) * n_instances(m) + instance, reader -> push!(readers, (k, reader)))
    lose!(km, rt.models[km], i, g.var)
    for (name, j, i′) in foreign
        k′ = findfirst(m -> m.name === name, rt.models)
        k′ === nothing && continue    # (a model without automata has no values to move)
        lose!(k′, rt.models[k′], i′, j)
    end
    return nothing
end

"""
    _states_written!(jump, integrator, rt::AutomatonRuntime, writes, t)

Extension point: what an augmentation carried through the integration beside the DAE
state (`jump`) does when a mode governing a state writes it directly (`writes`, the
states and the values written). The default, for `jump === nothing` (the plain
simulation), does nothing; an augmentation that must track what a governed state is set
to (e.g. its sensitivities to the parameters) adds a method for its own `jump` type.
"""
_states_written!(::Nothing, integrator, rt::AutomatonRuntime, writes, t) = nothing

"""
    _enforce_governed_sensitivity!(jump, ws, metadata, t)

Extension point (round 19a, design `4): what an augmentation carried through the integration beside the
DAE state (`jump`, the same object `_states_written!` dispatches on) does when [`enforce_governed!`](@ref)
substitutes a governed state's column on a workspace of `jump`'s own element type (a Dual/sensitivity
workspace). The default, for `jump === nothing` (the plain simulation) and for any `jump` type without its
own method (nothing needs doing beside what `enforce_governed!` itself already did to `ws.var_vals`), does
nothing: `enforce_governed!`, run on a Dual-typed workspace, already substitutes the *whole* Dual value
(primal and partials together, since `calc_automaton_values!`/`calc_residuals!` are the same
AD-compatible kernel regardless of element type), which already carries the relation's own correct
partials w.r.t. the parameters -- the natural, continuous generalization of `_states_written!`'s one-time
`G`-formula to "no switch, held the whole time". What is *not* covered by `enforce_governed!` alone, and
needs a `jump`-specific method, is a trajectory-sensitivity augmentation's own augmented state vector `u`
(`aj.S[a, :]`/`u[n*k+a]`, read directly from `u` by `_prepare_instant!`, not derived from a workspace
evaluation): a method for that `jump` type must substitute those rows the same way, called from the same
points `enforce_governed!` itself is (a trajectory-sensitivity augmentation provides this; not implemented
in this round -- see the round's own report).
"""
_enforce_governed_sensitivity!(jump, ws, metadata, t) = nothing
