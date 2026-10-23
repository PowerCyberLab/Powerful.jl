# ============================================================================
# Automata in the Simulation: Located Transitions and Settling After Events
# ============================================================================
#
# The automata of the models (automaton_runtime.jl) in an OrdinaryDiffEq integration:
#
# - One located `VectorContinuousCallback` with a component per instance for every
#   condition that can fire something (each condition transition and each held clause)
#   and every crossing guard. A component is zero (never an event) where it cannot fire
#   anything: for an out-of-service instance, for a transition from another mode, and
#   for every transition while a held clause holds. Which components are zero changes
#   at events only (modes, held clauses and the status change there), so within a step
#   DiffEq, which compares each point with the step's start, never sees a component
#   turn to or from zero. The components the modes and held clauses can fire are listed
#   by every settle (`AutomatonComponents`), so that an evaluation computes those alone.
# - A condition is reported by `signed_condition`: its atoms' `signed_truth` (sign =
#   truth of the comparison, never zero), the least over `AllOf`, the greatest over
#   `AnyOf`, negated through `Negation`. Its sign is the condition's truth, so the
#   callback locates exactly the instants the condition changes its truth, from an
#   exact zero of an atom too, and nothing else: an atom's zero is located only where
#   the atom decides the condition (in `AllOf` while the other terms hold, in `AnyOf`
#   while they fail), and there the condition is that atom. It is continuous, so two
#   atoms changing within one step (a rate turning outward, then a state reaching its
#   limit) change its sign too. A crossing guard is reported as it is, so that, as in
#   the core, a guard starting at its level and leaving it is not a crossing.
# - At a located event the values located (`_located_atoms!`) are those of the atoms of
#   the condition that changed their truth at the instant, between the two ends of the
#   root's bracket (in a condition of atoms combined once, those that decide it at one
#   end; in a nested one, also those that decide it at neither), and the guards of the
#   crossings, which fire in their direction. Every atom over a located value, of any
#   condition, is decided by the direction of the value's motion until the value moves
#   away from where the instant left it (`atom_holds`): through the instant (its settles and
#   reinitializations, which solve the algebraic variables to a tolerance) as the crossing
#   leaves it (`_locate_value!`; not by the value there, which is zero to rounding), and
#   after it by the value's displacement from where the instant left it (its anchor,
#   `_instant_left!`), as the crossing left it where that is zero. So a release that reads
#   the complement of a hold's atom over the same value reads it as the crossing left it,
#   and its next change of sign is a crossing where the value moves: where the discrete
#   state holds the value still, the integration's nonlinear solves still move it (by a
#   few ulps, or by their tolerance where the integrator keeps the Jacobian from before the
#   hold), and so does the integration's error over the short steps after an event: a
#   change of sign against the value's rate is no crossing (`_moving_crossings!`), and the
#   value is decided by its side until the next instant, or until a step starts with it
#   back on its side. The anchor is zero where the instant left the value on its side of
#   it (`_anchor`). A guard located at an instant is
#   reported as its displacement too, zero where the instant left it, so that, as at the
#   start, a guard starting at its level and leaving it is not a crossing. The sides of the
#   values an automaton reads are cleared where it changes mode (it reads them from another
#   mode then), unless they are located there; a value whose discrete quantities or memory
#   a transition changes loses its side at once, located or not (the change moves it, a
#   jump, not its motion: `_quantity_changed!`, from the values trace).
# - After every event of `events.jl` (and any other callback whose reinitialization is
#   bound by `_simulation_callbacks`), the DAE is reinitialized as before and the
#   automata then settle (`AutomatonSettle`).
# - Every settle (`settle_instant!`) keeps the algebraic variables consistent: where a
#   transition changes the algebraic equations (their residuals at the same state
#   differ before and after it), the network is reinitialized and the automata settle
#   again. Which instances' transitions can change them is traced from the residual
#   kernels (`AlgebraicCheck`); a settle evaluates nothing where none of those switches
#   (an anti-windup automaton), and otherwise the switched models' contributions where
#   they decide alone.
# - An integration that carries more than the DAE state passes a `jump`, which each
#   instant tells where it begins (with the values located
#   at their zero there) and ends (`_instant_begins!`, `_instant_ends!`), and the
#   runtime where memory starts and transitions fire (automaton_runtime.jl).
# ============================================================================

using DiffEqBase: isinplace

"""
    LocatedValue

A value of an automaton located at its zero by an event of the automaton callback: the
model (index into the runtime), the instance, the column of the values of the model (an
atom's or a crossing guard's), and the direction it crossed zero in (`+1` upward, `-1`
downward).
"""
const LocatedValue = Tuple{Int, Int, Vector{Float64}, Int8}

"""
    LocatedRates

The rates of the values of the automata at a located event, at the integrator's state and
time, in the discrete state of one instant's stage (before its transitions, or after a
round of them: `_rates_after`), shared by its checks: the motion of the values that crossed
where an earlier instant left them (`_moving_crossings!`) and the candidates for a sliding
mode (`_sliding_candidates`). `values[km]` are those of the values of model `km` (value,
instance), computed for the values asked for (`_located_value_rates!`): from the rate of the
rows with mass alone (`mass_rate`) where they read no variable without mass, the entries of
the others `NaN`; and from the rate of the whole state (`rate`, the rows without mass
completed, `complete` the models whose `values` rest on it) otherwise. `factor` keeps the
factorization of `F_ZZ`, and `reuse` whether the equations of the rows without mass are
those it was built for (it is then taken from it at the same state and time).
"""
mutable struct LocatedRates
    rate::Union{Nothing, Vector{Float64}}
    mass_rate::Union{Nothing, Vector{Float64}}
    values::Dict{Int, Matrix{Float64}}
    complete::Set{Int}
    factor::RateFactor
    reuse::Bool
end

LocatedRates(factor::RateFactor = RateFactor(), reuse::Bool = false) =
    LocatedRates(nothing, nothing, Dict{Int, Matrix{Float64}}(), Set{Int}(), factor, reuse)

_forget!(r::LocatedRates) = (r.rate = nothing; r.mass_rate = nothing; empty!(r.values); empty!(r.complete);
                             _forget!(r.factor); r.reuse = false; r)

"""
    automaton_conditions!(out, rt, evaluators, u, t)

The components of the automaton callback at `(u, t)` (see the header of automata.jl),
from the values of the automata, which it computes there with `evaluators` (the
runtime's, concretely typed where the callback holds them), and their atoms' signed
truths (`AutomatonModel` `signed`, measured from where located events left the values):
each block that can fire something (`AutomatonComponents` `conditions` and `crossings`,
the guards as `_guard_component`) for all its instances at once, the components that can
(`watch`) where their instance is in service, zero elsewhere.
"""
function automaton_conditions!(out, rt::AutomatonRuntime, evaluators::Tuple, u, t)
    nd = rt.n_diff
    _evaluate_values!(evaluators, view(u, 1:nd), view(u, (nd + 1):rt.n), t)
    comps = rt.components
    @inbounds for km in eachindex(rt.models)
        _signed_truths!(rt.models[km], comps.live[km])
    end
    fill!(out, 0.0)
    @inbounds for b in comps.conditions
        m = rt.models[comps.model[b]]
        n, offset, watch, active = n_instances(m), comps.offset[b], comps.watch[b], m.active
        signed = signed_block!(comps.scratch, comps.condition[b]::CompiledCondition, m.signed, n)
        for i in 1:n
            out[offset + i] = ifelse(watch[i] & active[i], signed[i], 0.0)
        end
    end
    @inbounds for b in comps.crossings
        m = rt.models[comps.model[b]]
        offset, watch, active = comps.offset[b], comps.watch[b], m.active
        tr = comps.transition[b]::CompiledTransition
        guard, v = tr.guard::Vector{Float64}, tr.guard_value
        for i in 1:n_instances(m)
            out[offset + i] = ifelse(watch[i] & active[i], _guard_component(m, guard, v, i), 0.0)
        end
    end
    return nothing
end

"""
    AtomTruths

The truths of the atoms of every instance, by their values, at each evaluation of the
automaton callback's components while DiffEq looks for events within one step (between
`tprev` and `t`; the record starts again with the next step): `times` the times of the
evaluations, `truths[km]` those of model `km`, atom `k` of instance `i` at evaluation `e`
at `((e - 1) * n_atoms + k - 1) * n_instances + i` (the signs of its `signed`, as they
lie). A located root's bracket ends at two adjacent times that DiffEq evaluated
(`_located_atoms!`); the record keeps the truths it saw there (those of the atoms the watched conditions read, `AutomatonComponents` `live`, are current; the rest are not read). Its vectors keep their
capacity, so a step allocates nothing once they have grown to the evaluations of the
busiest one.
"""
mutable struct AtomTruths
    tprev::Float64
    t::Float64
    times::Vector{Float64}
    truths::Vector{Vector{Bool}}
end

AtomTruths(rt::AutomatonRuntime) = AtomTruths(NaN, NaN, Float64[], [Bool[] for _ in rt.models])

# The components at `(u, t)` for a step that starts at `tprev` (`automaton_conditions!`); at the
# start of a step, the values pinned as they crossed against their rate that are back on their
# side are released first (`_returned!`), and the components are computed again
function _components!(out, rt::AutomatonRuntime, evaluators::Tuple, u, t, tprev)
    automaton_conditions!(out, rt, evaluators, u, t)
    (rt.returning[] > 0 && t == tprev && _returned!(rt)) && automaton_conditions!(out, rt, evaluators, u, t)
    return nothing
end

# The step an `AtomTruths` records: its two ends
mutable struct StepWindow
    tprev::Float64
    t::Float64
end

# Record the truths of the atoms at time `t` (their signed truths computed there) in the
# step of `integrator` (anything with `tprev` and `t`)
function _record_truths!(h::AtomTruths, rt::AutomatonRuntime, integrator, t)
    if integrator.tprev != h.tprev || integrator.t != h.t
        h.tprev, h.t = integrator.tprev, integrator.t
        empty!(h.times)
        foreach(empty!, h.truths)
    end
    push!(h.times, t)
    @inbounds for (km, m) in enumerate(rt.models)
        truths, signed = h.truths[km], m.signed
        base = length(truths)
        resize!(truths, base + length(signed))
        for j in eachindex(signed)
            truths[base + j] = signed[j] > 0
        end
    end
    return nothing
end

# The evaluation recorded at time `t` (the last, where there are several)
function _evaluation_at(h::AtomTruths, t::Float64)
    e = findlast(==(t), h.times)
    e === nothing && throw(ErrorException(
        "Automaton callback: no evaluation of its components was recorded at t = $t, an end of the bracket " *
        "of the root located there (recorded: $(length(h.times)) evaluations in the step from $(h.tprev) " *
        "to $(h.t)). DiffEq's root finding evaluates both ends of the bracket it returns."))
    return e
end

# The truth of atom `k` of instance `i` of model `km` at evaluation `e`
@inline function _truth(h::AtomTruths, rt::AutomatonRuntime, km::Int, e::Int, k::Int, i::Int)
    m = rt.models[km]
    return @inbounds h.truths[km][((e - 1) * length(m.atoms) + k - 1) * n_instances(m) + i]
end

"""
    _located_atoms!(located, prior, rt, h, km, i, c, t)

At an event of the automaton callback at time `t`, where condition `c` of instance `i` of
model `km` changed its truth: the values located at their zero there, between the ends of
the root's bracket (`prevfloat(t)` and `t`, both evaluated by DiffEq's bracketing root
finder, whose tolerance is zero; `t_before`, where a fixed-step integration brackets the
event, the near end it evaluated): the value of every atom of `c` that changed its truth
between the two, crossing zero in the direction that change gives (upward where a `>` or
`>=` turned true or a `<` or `<=` turned false; of its displacement from where an earlier
instant left it, for a value with a side, `atom_holds`). Each value is located by
`_locate_value!`: every atom over it holds as the crossing leaves it.

In a condition that combines atoms and their negations in one `AllOf` or `AnyOf`, every
atom that changed decides the condition at one end: at the end where an `AllOf` holds
(an `AnyOf` fails) every term does (fails), and flipping any one flips it. In a nested
condition an atom can change at the instant without deciding it at either end; its truth
is known from the bracket all the same.
"""
function _located_atoms!(located::Vector{LocatedValue}, prior::Vector{Int8}, rt::AutomatonRuntime, h::AtomTruths,
                         km::Int, i::Int, c::CompiledCondition, t::Float64, t_before::Float64 = prevfloat(t))
    m = rt.models[km]
    before, after = _evaluation_at(h, t_before), _evaluation_at(h, t)
    truth_before = k -> _truth(h, rt, km, before, k, i)
    truth_after = k -> _truth(h, rt, km, after, k, i)
    condition_holds(c, truth_before) != condition_holds(c, truth_after) || throw(ErrorException(
        "Automaton callback: an event of $(m.name) ($(join(string.(m.keys[i]), ", "))) was located at t = $t, " *
        "where the truth of its condition is the same at the two ends of the root's bracket."))
    for k in c.atoms
        truth_before(k) != truth_after(k) || continue
        upward = truth_after(k) == _is_upper(m.atom_relation[k])
        _locate_value!(located, prior, m, km, i, m.atom_column[k], upward ? Int8(1) : Int8(-1), t)
    end
    return nothing
end

"""
    _locate_value!(located, prior, m, km, i, column, direction, t)

At a located event at time `t`, value `column` of instance `i` of model `m` (index `km`)
crossed zero in `direction` (`+1` upward, `-1` downward; for a value with a side, its
displacement from where an earlier instant left it did): its side is that direction, and
it is at the instant (`AutomatonModel` `side`, `at_instant`), so that every atom over it
holds as the crossing leaves it (`>` and `>=` true and `<` and `<=` false after an upward
crossing, the reverse after a downward one) through the settles and reinitializations of
the instant, not by the value there, which is zero to rounding (`atom_holds`). The value is
pushed to `located` once, with its direction, and its side before the event to `prior`
(`0` for none). Throws where the value was already located at the instant crossing the
other way.
"""
function _locate_value!(located::Vector{LocatedValue}, prior::Vector{Int8}, m::AutomatonModel, km::Int, i::Int,
                        column::Vector{Float64}, direction::Int8, t::Float64)
    p = findfirst(l -> l[1] == km && l[2] == i && l[3] === column, located)
    (p === nothing || located[p][4] == direction) || throw(ErrorException(
        "Automaton callback: the value $(_value_name(m, column)) of $(m.name) " *
        "($(join(string.(m.keys[i]), ", "))) was located at t = $t crossing zero both upward and downward."))
    p === nothing || return nothing
    v = _value_index(m, column)
    push!(located, (km, i, column, direction))
    push!(prior, m.side[i, v])
    _set_side!(m, i, v, direction)
    _measure!(m, i, v, true)
    m.at_instant[i, v] = true
    return nothing
end

"""
    _instant_begins!(jump, integrator, located)
    _instant_ends!(jump, integrator, changed::Bool)

Extension points: an instant of transitions, bracketed for an augmentation carried
through the integration beside the DAE state (`jump`). `_instant_begins!` runs before the
automata settle at the integrator's state, with the values `located` at their zero there
(none at a start or after a fixed-time event); `_instant_ends!` runs after, with whether
any discrete state changed. The defaults, for `jump === nothing` (the plain simulation),
do nothing; an augmentation that must act once per instant (e.g. grazing corrections)
adds a method for its own `jump` type.
"""
_instant_begins!(::Nothing, integrator, located) = nothing
_instant_ends!(::Nothing, integrator, changed::Bool) = nothing

"""
    _memory_seeds(jump)
    _memory_seeds_changed(jump, seeds)

Extension points: the memory an augmentation carried through the integration beside the
DAE state (`jump`) is seeded with, and whether a start changed it (`seeds` from before the
start) — the start re-solves what `jump` carries only then. The defaults, for the plain
simulation, report no seeds (`_memory_seeds(jump) = nothing`) and treat that as changed
(`_memory_seeds_changed(jump, ::Nothing) = true`, so the start always re-solves); an
augmentation that tracks its own memory-like state adds methods for its own `jump` type.
"""
_memory_seeds(jump) = nothing
_memory_seeds_changed(jump, ::Nothing) = true

const NO_LOCATED_VALUES = LocatedValue[]

"""
    dae_residuals!(du, host, u, t)

The residuals `[f; g]` at `(u, t)` of the DAE `M u̇ = F(u, t)` the automata of a simulation
belong to, as its integration evaluates them: `host` is the simulation's `TDSAdapter`
(with the reduction of state-only constraints), or the DAE of a package that runs the
automata of its models with this runtime (the Y-bus DAE, PowerfulYbus), which adds its
method.
"""
dae_residuals!(du, adapter::TDSAdapter, u, t) = _evaluate_packed!(du, adapter, u, t, nothing)

"""
    HostReinit(reinit!) <: DAEInitializationAlgorithm

The reinitialization of a DAE that is not a `TDSAdapter`'s (see `dae_residuals!`), for the
settles of its automata: `reinit!(u, t)` makes the variables without mass of `u = [x; y]`
consistent at `t` with the states with mass held; the previous states of a multistep
method restart from it. A `reinit!` that returns `false` failed (a numerical failure it has
reported): the integration terminates with `ReturnCode.ConvergenceFailure`.
"""
struct HostReinit{F} <: DAEInitializationAlgorithm
    reinit!::F
end

function _initialize_dae!(integrator, prob, alg::HostReinit, ::Val{iip}) where {iip}
    _terminate_unless!(integrator, alg.reinit!(integrator.u, integrator.t) !== false) &&
        _restart_history!(integrator)
    return nothing
end

"""
    _terminate_unless!(integrator, converged::Bool) -> Bool

The one way every reinitialization of an integration ends on a numerical failure (the phasor
TDS's `PowerFlowInit`, the Y-bus DAE's `HostReinit`): where `converged` is false (the caller
has warned what failed) the integration terminates with `ReturnCode.ConvergenceFailure`, its
solution holding the trajectory up to the failed instant. Returns `converged`.
"""
function _terminate_unless!(integrator, converged::Bool)
    converged || DiffEq.terminate!(integrator, DiffEq.ReturnCode.ConvergenceFailure)
    return converged
end

# Whether the integration goes on: it has not been terminated (`terminate!`, with a failing or
# terminating return code); a running integration holds `Success` (or `Default`)
_running(integrator) = integrator.sol.retcode in (DiffEq.ReturnCode.Success, DiffEq.ReturnCode.Default)

# Reinitialize the DAE of `adapter` at the instant of `integrator` after the transitions of its
# automata changed the algebraic equations (their factorization taken afresh,
# `algebraic_equations_changed!`), with the simulation's `reinit_alg` started from the
# integrator's state (`at_current_state`; a network event's own reinitialization keeps the
# configured guess); whether the settle goes on (not where `reinit_alg` reinitializes nothing,
# or the reinitialization failed and terminated the integration, `_reinit_failed`). `jump`: see
# `_initialize_dae!`.
function _reinitialize!(integrator, adapter, reinit_alg, jump)
    algebraic_equations_changed!(adapter)
    prob = integrator.sol.prob
    alg = at_current_state(reinit_alg, jump)
    _initialize_dae!(integrator, prob, alg, jump, Val(isinplace(prob)))
    return !(alg isa DiffEq.NoInit) && _running(integrator)
end

# Whether the last reinitialization of an instant failed: it terminated the integration
# (`_terminate_unless!`, `ReturnCode.ConvergenceFailure`); for a fixed-step integration, its
# instant records it (`FixedInstant`)
_reinit_failed(integrator) = !_running(integrator)

"""
    settle_instant!(rt, host, integrator, reinit_alg; located = NO_LOCATED_VALUES,
                    jump = nothing, event = false, repair = event) -> Bool

Settle the automata at the integrator's state and time in rounds (`_settle_rounds!`), keeping the algebraic
variables consistent: every transition that holds, in every automaton of every instance, fires in a round, on the
point before it; where the round changed the residuals without mass at the same
state (the algebraic residuals and the state residuals whose mass is zero; their values
before and after differ, `_algebraic_before!` and `_algebraic_changed!`, on the DAE of
`host`, see `dae_residuals!`), `reinit_alg` reinitializes the DAE (the simulation's
`PowerFlowInit`) from the modes entered, and the next round reads that point, until no transition holds at the
point: the settled point, the only one the instant publishes and checks. `located`: at an event of the automaton callback, the values
located at their zero (`LocatedValue`s); its crossings fire in the first round, in which
the instances of those values alone can take a step (no other's condition changed its
truth there), and every instance in a round after a reinitialization; a transition
located there that enters a sliding mode throws (`_sliding_candidates`, `_check_sliding`).
A reinitialization that fails ends the instant at once, in any round: it terminated the integration with
`ReturnCode.ConvergenceFailure` and its warning, and the failed state decides nothing more
(no sliding check). The transitions' reinitializations start from the integrator's state
(`at_current_state`).
`jump`: see `_instant_begins!`; `rates`: those of the values at the event before its
transitions where the event computed them (`LocatedRates`). Returns whether any discrete
state changed. The located values are left where the instant leaves them, and the sides a
change of mode made irrelevant are cleared (`_instant_left!`).

`repair` (defaults to `event`, never true where `event` is not): whether the first round of this instant also
repairs a governed state whose right side a `FieldModification` just moved, in a mode that stays active
(`_repair_writes`, `governed_states.jl`); `false` only where `event` is `true` but this is not the last of several
`FieldModification`s applied at the same time (`AutomatonSettle` `repair`, `_simulation_callbacks`), so that a
relation that reads fields from two of them is read once, from all of them, not from the first alone.

The values the transitions hold at zero by continuity are decided by their rate in the
discrete state entered (`_held_at_zero!`), evaluated as the sliding check's are
(`AfterRates`, along `_state_rate` of the DAE of the integration, or the rate of the rows with
mass alone for a value that reads no variable without mass). The values
decided as zero since a start at rest are released where the instant changes the
equations (`_release_zeros!`): at an `event` (a fixed-time event of the simulation, before
the automata settle), and where a transition fired, after which every instance considers its
transitions again with them decided by their values.

Through the instant the atoms over the located values hold as the crossing leaves them,
whatever its rounds and reinitializations make of the values (`atom_holds`): the
reinitializations solve the algebraic variables to the tolerance of their Newton
iterations, which moves a value located at its zero by that much, either way.

The check evaluates what decides it alone: nothing where no instance that takes a
step has its discrete state read by a residual without mass (`AlgebraicCheck` `switches`;
an anti-windup automaton, whose gate holds a rate with mass), and of the others the
contributions of their models to the rows without mass, or the system's residuals where those do not decide alone
(`AlgebraicCheck` `alone`), before and after each round.
"""
function settle_instant!(rt::AutomatonRuntime, adapter, integrator, reinit_alg;
                         located::Vector{LocatedValue} = NO_LOCATED_VALUES, jump = nothing,
                         rates::Union{Nothing, LocatedRates} = nothing, event::Bool = false, repair::Bool = event)
    instances = isempty(located) ? nothing : unique!([(km, i) for (km, i, _) in located])
    _instant_begins!(jump, integrator, located)
    sliding = _sliding_candidates(rt, adapter, integrator, located, jump, rates)
    # (the rates of the values in the discrete state a round of transitions entered)
    held_rates = AfterRates(rt, adapter, integrator, jump, rates)
    event && _release_zeros!(rt)
    changed, failed = _settle_rounds!(rt, adapter, integrator, reinit_alg; located, instances, jump, rates = held_rates,
                                      event = repair)
    failed && return changed
    isempty(sliding) || _check_sliding(rt, adapter, integrator, sliding, jump, rates)
    _instant_ends!(jump, integrator, changed)
    u = integrator.u
    _instant_left!(rt, located, view(u, 1:rt.n_diff), view(u, (rt.n_diff + 1):rt.n), integrator.t)
    return changed
end

"""
    _settle_rounds!(rt, host, integrator, reinit_alg; located, instances, jump, rates, starting = false)
        -> (changed, failed)

The rounds of an instant at the state and time of `integrator` (a DiffEq integrator, or the `FixedInstant` of a
fixed-step integration). Every transition that holds, in every automaton of every instance, fires in a round
(`settle!`, which ends its round early where the round changed the equations without mass, `AlgebraicRounds`); the
algebraic point is then solved again from the modes entered (`_reinitialize!`) and the next round, every instance, reads
it. The point between rounds is provisional: nothing is published or checked there. Ends when no
transition holds at the point, or the equations are unchanged and nothing was released (`_release_zeros!`).
A reinitialization that fails fails the instant at once (`failed`). `starting`: see `settle!`; the equations
without mass in the initial modes are then solved before the first round where a residual of them at the
point is above `consistent_to` (`_equations_hold`). `event`: forwarded to `settle!` as its own `event` (so the
first round of the first call here also repairs a governed state the event moved, `_repair_writes`); `false`
unless `settle_instant!` was itself called with a `repair` that resolved true. A round that writes a governed
state this way without any instance switching sets `rt.settle_partial` with no `rt.settle_equations` (nothing
switched, so the ordinary algebraic-equations-changed check, keyed by what switched, cannot see it); this loop
reinitializes on either flag, not `settle_equations` alone, for exactly that reason (see the comment at the
call to `_reinitialize!`, below).
"""
function _settle_rounds!(rt::AutomatonRuntime, adapter, integrator, reinit_alg; located, instances, jump, rates,
                         starting::Bool = false, consistent_to::Float64 = Inf, event::Bool = false)
    u = integrator.u
    t = integrator.t
    nd, n = rt.n_diff, rt.n
    x = view(u, 1:nd)
    y = view(u, (nd + 1):n)
    changed = false
    # Round 19a: suspend `enforce_governed!` before `_evaluate_checked!` below is ever reached (the first
    # one, for a fixed-time event, reads the guards against the *retained* state and the event's new data,
    # the same comparison `settle!`'s own suspend protects -- but that call happens here, before `settle!`
    # itself is ever entered, so the suspend must start this early). `settle!` resumes it once its own
    # first round has decided its transitions.
    event && rt.governed && _set_suspended!(rt, true)
    if starting && !_equations_hold(rt, adapter, u, t, consistent_to)
        # The initial modes' own equations do not hold at the point the start was handed: solve it from them, as after
        # a round, before the first round reads it (the point is no longer the equilibrium the start was handed, so
        # the values decided as zero at rest are released)
        rt.reinitializations[] += 1
        if !_reinitialize!(integrator, adapter, reinit_alg, jump)
            _reinit_failed(integrator) && return (changed, true)
        end
        _release_zeros!(rt)
    end
    check = Ref(:none)
    rounds = AlgebraicRounds(taking -> (check[] = _algebraic_before!(rt, adapter, u, t, taking); nothing),
                             taking -> _algebraic_changed!(rt, adapter, u, t, check[], taking))
    while true
        candidates = changed ? nothing : instances
        # (at a located event, the event computed the values of the located models here)
        (changed || isempty(located)) && _evaluate_checked!(rt, candidates, x, y, t)
        starting || _hold_within_solve!(rt, adapter, integrator, reinit_alg, jump)
        settle!(rt, x, y, t; located = !isempty(located) && !changed, instances = candidates, jump,
                integrator, rates, starting, rounds, event) || break
        changed = true
        # The transitions changed the equations: what was at rest is released, and every
        # instance considers its transitions again with it decided by its value. A settle that ended its round
        # where a state a mode governs was written settles again whatever the check finds (it has
        # transitions left to consider)
        # (a start that fired transitions without changing the equations stays at rest: the values decided as zero are
        # released by an instant that changes the equations, and the start of a simulation is not one)
        released = (!starting || rt.settle_equations[]) && _release_zeros!(rt)
        if rt.settle_equations[] || rt.settle_partial[]
            # The round changed algebraic equations (`settle_equations`, from a switched instance's own rows, `rounds`),
            # or wrote a governed state an algebraic equation reads without any instance switching (`settle_partial`,
            # from `entered.partial`: a repair, `_repair_writes`, writes no discrete state, so `rounds.changed(taking)`
            # -- keyed by the instances that switched -- cannot see it; `entered.partial` is the only signal then):
            # reinitialize the network with the new point (factorized afresh where the equations changed,
            # `algebraic_equations_changed!`), and settle again at the consistent state
            rt.reinitializations[] += 1
            if !_reinitialize!(integrator, adapter, reinit_alg, jump)
                # A failed reinitialization terminated the integration with
                # `ReturnCode.ConvergenceFailure` and its warning (the worst rows): the state
                # does not satisfy the DAE, and nothing more is decided on it (no check for a
                # sliding mode, whose rates there would be those of no trajectory)
                _reinit_failed(integrator) && return (changed, true)
                break
            end
        elseif !released
            break
        end
    end
    return (changed, false)
end

"""
    StartResolve(host, instant, reinit_alg, abstol)

What solves the algebraic point at a start (`start!`): the DAE `host` (see `dae_residuals!`), the
`instant` whose state `u` the start's `x` and `y` are views of (a DiffEq integrator, or a `FixedInstant`), and the
reinitialization `reinit_alg` of the simulation (`nothing` for a `FixedInstant`, which solves by its own). The
initial modes are entered before the first round: where the equations without mass in them do not hold at the
initial point (a residual above `abstol`: the point is that of the initialization, which does not know the modes
the start chose), the point is solved from them first, and then again between the rounds as at any instant.
"""
struct StartResolve{H, I, A}
    host::H
    instant::I
    reinit_alg::A
    abstol::Float64
end

function _settle_start!(rt::AutomatonRuntime, r::StartResolve, x, y, t, jump, integrator)
    # (a reinitialization that fails has reported it: it terminated an integration, or marked its `FixedInstant`
    # `failed`, which the caller of the start reads)
    changed, _ = _settle_rounds!(rt, r.host, r.instant, r.reinit_alg; located = NO_LOCATED_VALUES, instances = nothing,
                                 jump, rates = nothing, starting = true, consistent_to = r.abstol)
    return changed
end

"""
    _equations_hold(rt, host, u, t, tol) -> Bool

Whether the rows without mass of the DAE of `host` (the algebraic rows, and the state rows whose mass is zero) hold at
`(u, t)`, in the modes the automata are in, to `tol` (one evaluation of the residuals of the DAE, not counted among the
evaluations of the algebraic check of a settle, `AutomatonRuntime` `system_evaluations`). Where no model has an instance whose
discrete state a residual without mass reads (`AlgebraicCheck` `switches`) the modes cannot have changed those rows, and
the answer is yes without evaluating: the point is the initialization's, whatever it satisfies.
"""
function _equations_hold(rt::AutomatonRuntime, host, u, t, tol::Float64)
    any(m -> any(m.algebraic.switches), rt.models) || return true
    nd, n = rt.n_diff, rt.n
    dae_residuals!(rt.after, host, view(u, 1:n), t)
    for a in (nd + 1):n
        abs(rt.after[a]) > tol && return false
    end
    for a in rt.massless
        abs(rt.after[a]) > tol && return false
    end
    return true
end

"""
    _start_automata!(adapter::TDSAdapter, rt::AutomatonRuntime, t, abstol, settled_without_network)

The start of the automata of `rt` on the initial buffers of `adapter` (the state the initialization of the DAE and DiffEq's
check at `t0` read), once the adapter exists to solve the algebraic point (`reinitialize_algebraic!`, as a fixed-step
integration solves its own), and its state-only constraints are detected. The adapter started the automata before
that, without a network (`start!`: the detection reads the discrete state), and `settled_without_network` is whether its
settle fired any transition. Where it fired none and the equations without mass of the initial modes hold at the initial
point (`_equations_hold`), that start was the start of the rule and nothing more is done. Otherwise the start is made again,
settling in rounds with the algebraic point solved between them from the modes entered (`StartResolve`), so that the
buffers the integration starts from satisfy the DAE in the modes the start settles in. A solve that fails is an error:
there is no integration to terminate.
"""
function _start_automata!(adapter::TDSAdapter, rt::AutomatonRuntime, t, abstol::Float64, settled_without_network::Bool)
    nd, n = adapter.n_diff, adapter.n_diff + adapter.n_alg
    u = vcat(adapter.x_buf, adapter.y_buf)
    !settled_without_network && _equations_hold(rt, adapter, u, t, abstol) && return nothing
    reinit! = (u, t) -> reinitialize_algebraic!(u, adapter, t)
    instant = FixedInstant(u, Float64(t), reinit!, () -> error("No rates are evaluated at the start of the automata"),
                           false)
    start!(rt, view(u, 1:nd), view(u, (nd + 1):n), t; resolve = StartResolve(adapter, instant, nothing, abstol))
    instant.failed && error("The start of the automata at t = $t: the algebraic point could not be solved again from " *
                            "the modes the first rounds entered")
    copyto!(adapter.x_buf, view(u, 1:nd))
    copyto!(adapter.y_buf, view(u, (nd + 1):n))
    return nothing
end

# ----------------------------------------------------------------------------
# Sliding at located transitions
# ----------------------------------------------------------------------------

"""
    SlidingCandidate

A transition of a located event that may enter a sliding mode (`_sliding_candidates`):
the value located at its zero (`LocatedValue`, with its direction `σ`), the automaton
(index) whose transition `transition` reads it and holds there, and its rate `ġ` before
the transitions of the instant.
"""
struct SlidingCandidate
    located::LocatedValue
    automaton::Int
    transition::CompiledTransition
    rate::Float64
end

const NO_SLIDING_CANDIDATES = SlidingCandidate[]

"""
    _sliding_candidates(rt, host, integrator, located, jump, rates) -> Vector{SlidingCandidate}

The transitions of the located event (`located`, before the automata settle) that may
enter a sliding mode, as `_check_sliding` decides after the settle: for each value `g`
located at its zero, each automaton of its instance (in service, no held clause holding)
whose first transition from its mode that holds (the step it takes) reads `g`, into a
mode with a transition that reads `g` too. The structure and the truths alone decide
this; for the candidates only, the rate of `g` in the discrete state before the
transitions, at the integrator's state (`_value_rates` along `_state_rate` of the DAE of
the integration, `_rate_dae`; from `rates` where the event computed them,
`LocatedRates`). Nothing is computed at an event without such a transition
(an anti-windup automaton reaching its limit: its release reads the rate, not the limited
state); its release, on the rate, which the transitions of `free` read, is a candidate.
"""
function _sliding_candidates(rt::AutomatonRuntime, host, integrator, located::Vector{LocatedValue}, jump,
                             rates::Union{Nothing, LocatedRates})
    isempty(located) && return NO_SLIDING_CANDIDATES
    found = Tuple{LocatedValue, Int, CompiledTransition}[]
    for l in located
        km, i, column = l
        m = rt.models[km]
        m.active[i] || continue
        for (ka, a) in enumerate(m.automata)
            any(c -> holds(c, m, i), a.held) && continue
            mode = a.modes[i]
            k = findfirst(tr -> tr.from == mode &&
                                (tr.crossing == 0 ? holds(tr.compiled::CompiledCondition, m, i) :
                                 m.crossed[tr.crossing, i] == tr.direction), a.transitions)
            k === nothing && continue
            tr = a.transitions[k]
            (_reads(m, tr, column) && any(back -> back.from == tr.to && _reads(m, back, column), a.transitions)) ||
                continue
            push!(found, (l, ka, tr))
        end
    end
    isempty(found) && return NO_SLIDING_CANDIDATES
    values = _located_value_rates!(rates === nothing ? LocatedRates() : rates, rt, host, integrator, jump,
                                   unique((f[1][1], _value_index(rt.models[f[1][1]], f[1][3])) for f in found))
    return [SlidingCandidate(l, ka, tr, values[l[1]][_value_index(rt.models[l[1]], l[3]), l[2]]) for (l, ka, tr) in found]
end


# The rates of the values of the automata asked for in `requests` (pairs of a model, an index
# into the runtime, and the index of one of its values) at the integrator's state and time, in the
# current discrete state: a matrix (value, instance) for each model, those not in `cache` computed
# (along one rate of the state) and kept there. A value that reads no variable without mass has
# its rate from the rows with mass alone, with no solve of `F_ZZ`; the rows without mass are
# completed (`_state_rate`) for a model with a value that does.
function _located_value_rates!(cache::LocatedRates, rt::AutomatonRuntime, host, integrator, jump, requests)
    all(((km, v),) -> _rates_known(cache, rt.models[km], km, v), requests) && return cache.values
    r = _rate_dae(jump, rt, host, integrator)
    u, t = integrator.u, integrator.t
    for (km, v) in requests
        m = rt.models[km]
        _rates_known(cache, m, km, v) && continue
        if cache.rate !== nothing || _reads_massless(m, v)
            if cache.rate === nothing
                built = cache.factor.built
                cache.rate = _state_rate(r, u, t; cache.factor, cache.reuse)
                rt.rate_completions[] += 1
                rt.rate_factorizations[] += cache.factor.built - built
            end
            cache.values[km] = _value_rates(r, m, u, cache.rate, t)[2]
            push!(cache.complete, km)
        else
            cache.mass_rate === nothing && (cache.mass_rate = _mass_rate(r, u, t))
            rates = _value_rates(r, m, u, cache.mass_rate, t)[2]
            for k in axes(rates, 1)
                _reads_massless(m, k) && (rates[k, :] .= NaN)
            end
            cache.values[km] = rates
        end
    end
    return cache.values
end

# Whether value `v` of model `m` reads a variable without mass
_reads_massless(m::AutomatonModel, v::Int) = reads_massless(m.reads, v, length(m.state.values))

# Whether the rate of value `v` of model `km` (`m`) is in `cache`: computed from the rate of the
# whole state, or from the rows with mass alone where the value reads none without mass
_rates_known(cache::LocatedRates, m::AutomatonModel, km::Int, v::Int) =
    haskey(cache.values, km) && (km in cache.complete || !_reads_massless(m, v))

"""
    AfterRates(rt, host, integrator, jump, before)

The rates of the values in the discrete state the transitions of the instant entered, as
`_held_at_zero!` asks for them, `before` those at the event (`LocatedRates`; `nothing`
where the event computed none). `_rates_after` makes the cache of one round.
"""
struct AfterRates{H, I, J}
    rt::AutomatonRuntime
    host::H
    integrator::I
    jump::J
    before::Union{Nothing, LocatedRates}
end

# The rates of value `w` of model `km` in the cache `cache` of a round (a matrix (value, instance))
(a::AfterRates)(cache::LocatedRates, km::Int, w::Int) =
    _located_value_rates!(cache, a.rt, a.host, a.integrator, a.jump, ((km, w),))[km]

# The tolerance of entry `j` of the integrator's state: what the integration leaves it uncertain by
@inline function _state_tolerance(opts, u, j::Int)
    abstol = opts.abstol isa AbstractArray ? opts.abstol[j] : opts.abstol
    reltol = opts.reltol isa AbstractArray ? opts.reltol[j] : opts.reltol
    return abstol + reltol * abs(u[j])
end

"""
    _resolved_rate(rates::AfterRates, computed::LocatedRates, km, w, i, rate) -> Float64

`rate` of value `w` of instance `i` of model `km` where its sign is resolved by the integration, `0.0` where it is not.

The rate is computed at a state the integration knows to its tolerance, and a rate is a sum of terms of the size of the
fastest motion of that state: a state oscillating at 10² to 10³ cycles per second has a rate of that order, so a state error of
10⁻⁶ moves the rate of a slow envelope (a measured magnitude, the value of a held clamp) by about 10⁻⁴. A rate smaller than that is
the integration's error, and its sign is a coin: a value held at zero by continuity and decided by it flips between runs of
different tolerance. Where the rate is exactly zero in the model (a state continuous across the transition that leaves the value's
first derivative unchanged: the release of a clamp that holds a measured state) the sign must not decide; the value then stays at its
zero (`zero`), a transition that reads it strictly does not fire, and a crossing that follows is located as any other.

The sign is resolved where `|rate|` exceeds the first-order bound of its error `Σ_j |∂ġ/∂u_j| tol_j`, `tol_j = abstol_j + reltol_j |u_j|` the
tolerance of entry `j` of the state (`_rate_error_bound`: the gradient of the rate with respect to the state, by Duals and one adjoint
solve of `F_ZZ`; see the header of that section of automaton_rates.jl). Nothing is evaluated for an exact zero, or without an integrator that has
a tolerance: a fixed-step integration has none, and decides by the raw rate.
"""
function _resolved_rate(a::AfterRates, computed::LocatedRates, km::Int, w::Int, i::Int, rate::Float64)
    (iszero(rate) || a.integrator === nothing) && return rate
    integrator = a.integrator
    hasproperty(integrator, :opts) || return rate
    r = _rate_dae(a.jump, a.rt, a.host, integrator)
    m = a.rt.models[km]
    state_rate = computed.rate !== nothing ? computed.rate : computed.mass_rate
    bound = _rate_error_bound(r, m, w, i, integrator.u, integrator.t, state_rate, integrator.opts)
    return abs(rate) > bound ? rate : 0.0
end

# The factor of each row without mass of the DAE between its entry of `F` and the residual the reinitialization's convergence test
# reads: the test is on `r = g(z)`, the rows without mass with the power balance of each in-service bus divided by its voltage (a current
# mismatch, `to_current_mismatch!`), so `F_k = s_k r_k` with `s_k` the bus voltage on those rows and 1 on the others. `nothing` for a host
# whose reinitialization this does not know (it has no rows to scale: `_hold_within_solve!` does nothing there)
_solve_row_scales(host, u, n_z::Int) = nothing
function _solve_row_scales(adapter::TDSAdapter, u, n_z::Int)
    cache = adapter.reinit_cache
    cache === nothing && return nothing
    n_massless = length(cache.massless)
    n_massless + adapter.n_alg == n_z || return nothing
    scales = ones(n_z)
    rows = cache.rows
    if rows !== nothing
        bus_status = adapter.status[:Bus]
        for k in eachindex(rows.v)
            bus_status[k] || continue
            v = u[adapter.n_diff + rows.v[k]]
            scales[n_massless + rows.p[k]] = v
            scales[n_massless + rows.q[k]] = v
        end
    end
    return scales
end

"""
    _hold_within_solve!(rt, host, integrator, reinit_alg, jump) -> Bool

Where the reinitialization has just left a value of an automaton within its solve tolerance of the surface it was left on at the
last instant, decide the value as at that surface (`zero`) through this instant. Returns whether any value was.

A value an instant left with a side (`side`, `anchor`: a located value, or one held by continuity) is decided from there by its
displacement from the anchor, and the reinitialization that begins the next event, or follows a round of transitions, solves the
algebraic variables to a tolerance, which moves the value by its sensitivity to them. The solve stops (`newton_solve!`) when the
∞-norm of its residual `r` is below `τ` (`reinit_alg.abstol`); `r` is the rows without mass of the DAE, the power balance of each in-service
bus divided by its voltage (`to_current_mismatch!`), so what remains of `F_Z` is `|F_k| <= τ s_k` (`_solve_row_scales`). With
`p = ∇_u g` and `λ = F_ZZ⁻ᵀ p_Z`, the value is off by `λᵀ δF_Z` at most

    B = τ Σ_k |λ_k| s_k      (the dual norm of the ∞-norm of the scaled residual: the weighted 1-norm of λ).

A displacement of the size of `B` is the solve's, and its sign decides nothing: the PI of REECA1 held at its limit by a transition read
-5e-11 after the reinitialization of the next event, where the solve tolerance of 1e-8 bounds the value by 2.0e-8 (`τ = 1e-8`, `Σ|λ| s = 2.04`),
and released, to be held again 4e-11 s later. Only displacements within the integration's tolerance are examined (the reinitialization's is a
hundredth of it); a displacement larger than the bound is the value's own motion. Nothing is done without an integrator that has a
tolerance, a reinitialization that names its own, or a host whose reinitialization `_solve_row_scales` does not know (the Y-bus
simulation).
"""
function _hold_within_solve!(rt::AutomatonRuntime, host, integrator, reinit_alg, jump)
    hasproperty(integrator, :opts) || return false
    τ = hasproperty(reinit_alg, :abstol) ? reinit_alg.abstol : nothing
    τ isa Real || return false
    limit = minimum(integrator.opts.abstol)
    nd, n = rt.n_diff, rt.n
    u, t = integrator.u, integrator.t
    held = false
    r = nothing
    K = nothing
    scales = nothing
    for (km, m) in enumerate(rt.models)
        pending = [(i, v) for v in axes(m.side, 2), i in 1:n_instances(m)
                   if m.active[i] && m.side[i, v] != 0 && !m.zero[i, v] && !m.pinned[i, v] && !m.at_instant[i, v]]
        isempty(pending) && continue
        evaluate_values!(rt, [(km, i) for i in unique(first.(pending))], view(u, 1:nd), view(u, (nd + 1):n), t)
        for (i, v) in pending
            d = m.value_columns[v][i] - m.anchor[i, v]
            (iszero(d) || abs(d) > limit) && continue
            if r === nothing
                r = _rate_dae(jump, rt, host, integrator)
                scales = _solve_row_scales(host, u, length(r.zero_mass))
                (scales === nothing || isempty(r.zero_mass)) && return false
                K = _zero_mass_factor(r, u, t)
            end
            p = _value_gradient(r, m, v, i, u, t)
            pZ = p[r.zero_mass]
            all(iszero, pZ) && continue
            bound = τ * sum(abs.(K' \ pZ) .* scales)
            @debug "value left within the reinitialization's tolerance of its anchor" model = m.name value = keys(m.state.values)[v] instance = i displacement = d solve_bound = bound tolerance = τ held = abs(d) <= bound
            abs(d) <= bound || continue
            _set_side!(m, i, v, Int8(0))
            _zero!(m, i, v, true)
            m.at_instant[i, v] = true
            push!(rt.held, (km, i, v))
            held = true
        end
    end
    return held
end

"""
    _rates_after(rt, before) -> LocatedRates

The cache for the rates at the state of the integrator after the transitions of the instant
so far (`rt.switched_automata`), `before` those at the event. It takes the factorization
of `F_ZZ` from `before` where no automaton that switched has its discrete state read by a
residual without mass (`AlgebraicCheck` `switches`): the equations of the rows without mass
are then those of the event, and the factorization is reused if the state and time are those
it was built at (a reinitialization moves the state only where the algebraic equations changed).
"""
function _rates_after(rt::AutomatonRuntime, before::Union{Nothing, LocatedRates})
    before === nothing && return LocatedRates()
    unchanged = !any(((km, i, _),) -> rt.models[km].algebraic.switches[i], rt.switched_automata)
    return LocatedRates(before.factor, unchanged)
end

"""
    _moving_crossings!(located, prior, rt, host, integrator, jump, rates) -> located

Of the values `located` at an event of the automaton callback, those that crossed from
where an earlier instant left them (`prior`, their sides before the event, nonzero; their
displacement from the anchor changed sign, `atom_holds`, `_guard_component`) cross where
they move across it: their rate `ġ` in the discrete state before the event, at the
integrator's state (`rates`, `_located_value_rates!`), has the direction `σ` of the
crossing, `σ ġ > 0`. The others are no crossing: the sign changed by what the integration
moves the value by besides its motion, which no step after an instant resolves (its
nonlinear solves, by their tolerance where the integrator keeps the Jacobian of the mode
before the instant, which it renews only where Newton struggles; its error, within the
tolerance, over the short steps after an event; the low order of the dense output after a
restart). Such a value keeps its side and its anchor, its atoms are decided by its side
and a guard over it is zero (`pinned`, `atom_holds`), it is no longer located, and a
crossing of a guard over it is forgotten; `rt.still_crossings` counts them. It stays
pinned:

- `ġ = 0` (the discrete state holds the value still, a held integrator, whose rate its
  gate zeroes): until the next instant measures it again (`_instant_left!`);
- `ġ` the other way (it moves back to its side): until a step starts with it back on its
  side of the anchor (`returning`, `_returned!`), or the next instant, whichever is first.

A value that moves crosses where its displacement from the anchor changes sign: where it
crosses its surface, or the anchor, which is off the surface by what the instant's
reinitialization (or the integration of a value held still) moved it across, `|g*|`; that
moves the crossing by `|g*| / |ġ|` in time, within the integration's tolerance where `g*`
is (the reinitialization's Newton tolerance, `|ġ|` that of a value crossing
transversally).
"""
function _moving_crossings!(located::Vector{LocatedValue}, prior::Vector{Int8}, rt::AutomatonRuntime, host,
                            integrator, jump, rates::LocatedRates)
    any(!iszero, prior) || return located
    values = _located_value_rates!(rates, rt, host, integrator, jump,
                                   unique((located[q][1], _value_index(rt.models[located[q][1]], located[q][3]))
                                          for q in eachindex(located) if prior[q] != 0))
    keep = fill(true, length(located))
    for (q, (km, i, column, σ)) in enumerate(located)
        prior[q] == 0 && continue
        m = rt.models[km]
        v = _value_index(m, column)
        ġ = values[km][v, i]
        isfinite(ġ) || throw(ErrorException(
            "Automaton callback: the rate of $(_value_name(m, column)) of $(m.name) " *
            "($(join(string.(m.keys[i]), ", "))) is $ġ where it crosses where an earlier instant left it, at " *
            "t = $(integrator.t)."))
        σ * ġ > 0 && continue
        keep[q] = false
        _set_side!(m, i, v, prior[q])
        m.at_instant[i, v] = false
        m.pinned[i, v] = true
        _measure!(m, i, v, true)
        m.returning[i, v] = !iszero(ġ)
        rt.returning[] += !iszero(ġ)
        for a in m.automata, tr in a.crossings
            tr.guard_value == v && (m.crossed[tr.crossing, i] = 0)
        end
        rt.still_crossings[] += 1
    end
    all(keep) && return located
    keepat!(prior, keep)
    return keepat!(located, keep)
end

"""
    _rate_dae(jump, rt, host, integrator) -> RateDAE

Extension point: the DAE of an integration whose rates the automata evaluate (a
`RateDAE`). The default, for the plain simulation (`jump === nothing`), or a `jump` that
only follows the instants and carries nothing more, is the rates of the DAE of `host`
(`StateRates`, built at the first use from the pattern of the integrator's Jacobian,
`jac_prototype`, or a full pattern where it has none, which DiffEq takes as a dense
Jacobian (the Y-bus DAE with the network eliminated), and kept in `rt.rates`); an
augmentation that integrates more than the DAE state (e.g. its sensitivities) adds a
method for its own `jump` type, returning its own `RateDAE`.
"""
function _rate_dae(jump, rt::AutomatonRuntime, host, integrator)
    rt.rates[] === nothing || return rt.rates[]
    prototype = integrator.f.jac_prototype
    sparsity = prototype === nothing ? sparse(trues(rt.n, rt.n)) : prototype
    sparsity isa AbstractSparseMatrix || throw(ErrorException(
        "The rates of the automata need the sparsity pattern of the Jacobian of the DAE, the integrator's " *
        "`jac_prototype`; this integration has $(typeof(sparsity))."))
    rt.rates[] = StateRates(host, rt.n_diff, rt.mass, sparsity)
    return rt.rates[]
end

"""
    _check_sliding(rt, host, integrator, candidates, jump, before)

Throw where a located transition entered a sliding mode (a Filippov sliding mode, which
the simulation does not support): the transition was located on a value `g` crossing
zero in direction `σ`, and after the instant (through which `g` is at its zero,
`atom_holds`) its automaton is in the mode `B` the transition entered, from which a transition would fire
at once on `g` crossing back (its condition holds with the atoms over `g` at their truths
before the crossing, the others as they hold; or a crossing back over the guard `g`),
and `g` moves back across zero there: `σ ġ_A > 0` before the transition (towards `B`)
and `σ ġ_B < 0` after it (towards `A`), `ġ_B` evaluated in the discrete state after the
instant at the integrator's state then. Names the instance, the automaton, `g`, both
modes and both rates, and every automaton that switched at the instant (`_switched_names`).

It does not throw where `g` stays on the surface in `B` (the anti-windup automaton,
whose gate makes the rate of the limited state exactly zero once held), where `B` moves
`g` further across (an ordinary reach or release), where nothing from `B` would fire
on `g`, or where the instant moved `g` off its surface by changing what it reads (it
lost its side there, `_quantity_changed!`).
"""
function _check_sliding(rt::AutomatonRuntime, host, integrator, candidates::Vector{SlidingCandidate}, jump,
                        before::Union{Nothing, LocatedRates})
    after = nothing
    after_cache = nothing
    for c in candidates
        km, i, column, σ = c.located
        m = rt.models[km]
        a = m.automata[c.automaton]
        B = c.transition.to
        (a.modes[i] == B && !a.held_now[i]) || continue
        # (a value the transitions moved by changing what it reads jumped off its surface)
        m.side[i, _value_index(m, column)] == 0 && continue
        flipped(k) = m.atom_column[k] === column ? (σ < 0) == _is_upper(m.atom_relation[k]) : atom_holds(m, k, i)
        k = findfirst(a.transitions) do back
            back.from == B || return false
            back.crossing == 0 || return back.guard === column && back.direction == -σ
            return _reads(m, back, column) && condition_holds(back.compiled::CompiledCondition, flipped)
        end
        k === nothing && continue
        σ * c.rate > 0 || continue
        if after === nothing
            after_cache = _rates_after(rt, before)
            after = _located_value_rates!(after_cache, rt, host, integrator, jump,
                                          unique((cc.located[1], _value_index(rt.models[cc.located[1]], cc.located[3]))
                                                 for cc in candidates))
        end
        w = _value_index(m, column)
        rate_B = after[km][w, i]
        # (a rate in B that the integration does not resolve, the rounding of a held integrator's rate, is no motion back across zero)
        σ * _resolved_rate(AfterRates(rt, host, integrator, jump, before), after_cache, km, w, i, rate_B) < 0 || continue
        modes = a.automaton.modes
        A, back = modes[c.transition.from].name, a.transitions[k]
        throw(ErrorException(
            "Sliding mode at t = $(integrator.t): $(m.name) ($(join(string.(m.keys[i]), ", "))), automaton " *
            "$(a.automaton.name): the transition $A -> $(modes[B].name) was located where its value " *
            "$(_value_name(m, column)) crosses zero $(σ > 0 ? "upward" : "downward") at the rate $(c.rate) in $A; " *
            "in $(modes[B].name) its rate is $rate_B, back across zero, where $(modes[B].name) -> " *
            "$(modes[back.to].name) fires at once. The trajectory slides on $(_value_name(m, column)) = 0 " *
            "(a Filippov sliding mode), which the simulation does not support. The automata that switched at this " *
            "instant: $(_switched_names(rt))."))
    end
    return nothing
end

# The automata whose mode the settles of the instant changed (`switched_automata`), each once, by model, instance and
# declaration order, as `MODEL (keys) automaton`: several automata may transition on one located event (a hold and a
# piece automaton on one surface), and which of them the event locates is decided by rounding
function _switched_names(rt::AutomatonRuntime)
    names = String[]
    for (km, i, ka) in sort!(unique(rt.switched_automata))
        m = rt.models[km]
        push!(names, "$(m.name) ($(join(string.(m.keys[i]), ", "))) $(m.automata[ka].automaton.name)")
    end
    return join(names, ", ")
end

"""
    _algebraic_before!(rt, host, u, t, candidates) -> Symbol

What the algebraic check of a settle evaluates at `(u, t)` before the transitions of
`candidates` (`(model, instance)` pairs; every instance where `nothing`), and how it
decides after them (`_algebraic_changed!`): `:none` where no candidate has its discrete
state read by a residual without mass (`AlgebraicCheck` `switches`), evaluating nothing;
`:models` where every model of those (`rt.checked`) decides alone, copying their
contributions to the rows without mass (`AlgebraicCheck` `before`); `:system` otherwise, evaluating the
system's residuals (`rt.before`).
"""
function _algebraic_before!(rt::AutomatonRuntime, adapter, u, t, candidates)
    checked = rt.checked
    empty!(checked)
    if candidates === nothing
        for (km, m) in enumerate(rt.models)
            any(m.algebraic.switches) && push!(checked, km)
        end
    else
        # (a model whose values read other instances' discrete state settles every instance)
        for (km, i) in _expand_cross(rt, candidates)
            (rt.models[km].algebraic.switches[i] && !(km in checked)) && push!(checked, km)
        end
    end
    isempty(checked) && return :none
    nd, n = rt.n_diff, rt.n
    if all(km -> rt.models[km].algebraic.alone, checked)
        x, y = view(u, 1:nd), view(u, (nd + 1):n)
        for km in checked
            a = rt.models[km].algebraic
            a.evaluate!(x, y, t)
            rt.model_evaluations[] += 1
            foreach(copyto!, a.before, a.columns)
        end
        return :models
    end
    dae_residuals!(rt.before, adapter, view(u, 1:n), t)
    rt.system_evaluations[] += 1
    return :system
end

"""
    _algebraic_changed!(rt, host, u, t, check::Symbol, switched) -> Bool

Whether the transitions of the last round of a settle (of the instances `switched`) changed the residuals without
mass at `(u, t)` (the algebraic rows, and the rows of `x` whose mass is zero), as
`_algebraic_before!` prepared it (`check`). Only an instance whose discrete state a
residual without mass reads can change them; where one switched, its model's contributions
to the rows without mass (`:models`) or the system's rows without mass (`:system`) are
evaluated again and compared, exactly, with those before.

The contributions decide alone (`AlgebraicCheck` `alone`): every other model's are
unchanged at the same state (they read their own discrete state, unchanged, and no
residual accumulated by another: none of them sets a shared residual where a switched
model adds), so the system's rows without mass, which are their sums, are unchanged
where the switched models' are. Where those changed, the system's differ unless the
changes cancel in a sum, or vanish in its rounding; the network is then reinitialized at
a state that was consistent.
"""
function _algebraic_changed!(rt::AutomatonRuntime, adapter, u, t, check::Symbol, switched)
    check === :none && return false
    nd, n = rt.n_diff, rt.n
    read = false
    for (km, i) in switched
        rt.models[km].algebraic.switches[i] || continue
        km in rt.checked || throw(ErrorException(
            "Automaton settle: $(rt.models[km].name) ($(join(string.(rt.models[km].keys[i]), ", "))) switched " *
            "without its algebraic contributions copied before; a settle switches only the instances it was given."))
        read = true
    end
    read || return false
    if check === :system
        dae_residuals!(rt.after, adapter, view(u, 1:n), t)
        rt.system_evaluations[] += 1
        view(rt.before, (nd + 1):n) == view(rt.after, (nd + 1):n) || return true
        return any(r -> rt.before[r] != rt.after[r], rt.massless)
    end
    x, y = view(u, 1:nd), view(u, (nd + 1):n)
    for km in rt.checked
        a = rt.models[km].algebraic
        any(((k, i),) -> k == km && a.switches[i], switched) || continue
        a.evaluate!(x, y, t)
        rt.model_evaluations[] += 1
        _contributions_changed(a) && return true
    end
    return false
end

"""
    AutomatonSettle(inner, runtime, adapter, jump) <: DAEInitializationAlgorithm

The reinitialization after an event of a simulation with automata: `inner` (the
simulation's reinitialization, e.g. a bound `PowerFlowInit`), then the automata settle at
the reinitialized state (`settle_instant!`, reinitializing with `inner` where a
transition changes the algebraic equations), telling `jump` (see `_instant_begins!`;
`nothing` for the plain simulation). Which instances' transitions can change the
algebraic equations is traced again first for the models of `retraced` (indices into the
runtime, `_event_models`), whose parameters the events may change, on which a residual
kernel may branch, and for those of `value_retraced` whose residual kernel branches on a
setpoint (`AlgebraicCheck` `retrace!`); and what the values of the automata read for the
models of `value_retraced`, whose parameters or setpoints the events may change, on which a
value kernel may branch (`retrace_value_reads!`). An augmentation that carries a `jump`
takes the dependencies of the values at its setup, and throws where a retrace finds that
they changed. `repair`: whether the settle considers the repair of a governed state after
this event (`event = repair`, `settle_instant!`); `true` unless the caller built this one specifically to skip it
(`_simulation_callbacks`, for a `FieldModification` that is not the last of several at the same time: see its own
header, and the header of `governed_states.jl`, for why only the last may).
"""
struct AutomatonSettle{A, H, J} <: DAEInitializationAlgorithm
    inner::A
    runtime::AutomatonRuntime
    adapter::H
    jump::J
    retraced::Vector{Int}
    value_retraced::Vector{Int}
    repair::Bool
end

AutomatonSettle(inner, runtime::AutomatonRuntime, adapter, jump, events::Vector, sys; repair::Bool = true) =
    AutomatonSettle(inner, runtime, adapter, jump, _event_models(runtime, events, sys),
                    _event_models(runtime, events, sys; setpoints = true), repair)

"""
    _with_repair(alg::AutomatonSettle, repair::Bool) -> AutomatonSettle

`alg` with its `repair` field set to `repair`, everything else shared (the same runtime, adapter and retraced
sets): one `AutomatonSettle` built per problem (`_simulation_callbacks`), a cheap copy of it per `FieldModification`
callback that must not repair (see `AutomatonSettle`'s header).
"""
_with_repair(alg::AutomatonSettle, repair::Bool) =
    AutomatonSettle(alg.inner, alg.runtime, alg.adapter, alg.jump, alg.retraced, alg.value_retraced, repair)

"""
    _event_models(rt::AutomatonRuntime, events, sys; setpoints = false) -> Vector{Int}

The models of `rt` (indices) whose residual kernels an event of `events` may make branch
otherwise: those whose fields a `FieldModification` changes (the fields enter the trace
at their values, `trace_algebraic_switches`), and those whose external parameters are copies
of the field it changes (`external_dependents` of the system `sys`). A `SetpointChange` changes none of the
fields; with `setpoints` (the traces that have the setpoints at their values where a
kernel branches on them, `trace_value_reads` and `trace_algebraic_switches`), the models
whose setpoints it changes as well.
Every model where an event of another kind is among them, whose changes are not declared.
"""
function _event_models(rt::AutomatonRuntime, events::Vector, sys; setpoints::Bool = false)
    names = Set{Symbol}()
    for event in events, e in (event isa Vector ? event : [event])
        if e isa SetpointChange
            setpoints && push!(names, e.model)
            continue
        end
        e isa FieldModification || return collect(eachindex(rt.models))
        push!(names, e.model)
        # (and the models whose external parameters are copies of the field)
        union!(names, external_dependents(sys, e.model, e.field))
    end
    return [km for (km, m) in enumerate(rt.models) if m.name in names]
end

"""
    retrace_automata!(rt, models_by_name, retraced, value_retraced; augmented=false)

Trace again, after events that changed the parameters or setpoints of models, what the
automata of the models of `rt` read: which transitions change the algebraic equations for
the models `retraced` (`AlgebraicCheck` `retrace!`), and for the models of `value_retraced`
(indices into `rt`, `_event_models`) whose residual kernel branches on a setpoint, and what
the values read for the models `value_retraced` (`retrace_value_reads!` of the model arrays
in `models_by_name`). The one place every integration with automata (the phasor
TDS's `AutomatonSettle`, and the Y-bus integrations) re-traces. With `augmented` set, an
augmentation that takes what the values read at its setup throws where the values read
otherwise.
"""
function retrace_automata!(rt::AutomatonRuntime, models_by_name, retraced, value_retraced; augmented::Bool = false)
    foreach(km -> rt.models[km].algebraic.retrace!(), retraced)
    for km in value_retraced
        m = rt.models[km]
        km in retraced || !m.algebraic.setpoints[] || m.algebraic.retrace!()
        changed = retrace_value_reads!(m.reads, models_by_name[m.name])
        m.governed.recheck!()   # (what a right side of a governed state reads may have changed)
        changed && augmented && throw(ErrorException(
            "The values of the automata of $(m.name) read other variables or discrete quantities after the event " *
            "that changed their parameters or setpoints, where such an augmentation takes what the values " *
            "read at their setup (`_value_dependencies`)"))
    end
    return nothing
end

function _initialize_dae!(integrator, prob, alg::AutomatonSettle, ::Val{iip}) where {iip}
    # DiffEq calls here with no parameter for `jump`: `alg` (the `AutomatonSettle` it is
    # holding) carries its own, read here and forwarded to the extension point
    _initialize_dae!(integrator, prob, alg.inner, alg.jump, Val(iip))
    _running(integrator) || return nothing
    # The events may have changed parameters on which residual kernels branch
    retrace_automata!(alg.runtime, alg.adapter.sys_ws.models_by_name, alg.retraced, alg.value_retraced;
                      augmented = alg.jump !== nothing)
    # Transitions here were located by no value at its zero; the event changed the equations. `repair`: see
    # `AutomatonSettle`'s header (every event releases the values held at zero, `event = true`; only the last of
    # several `FieldModification`s at the same time repairs a governed state, `repair = alg.repair`)
    settle_instant!(alg.runtime, alg.adapter, integrator, alg.inner; jump = alg.jump, event = true, repair = alg.repair)
    return nothing
end

"""
    _located_event!(rt, host, integrator, truths, fired, direction, t_before, jump, reinit_alg, located,
                    prior, rates) -> Bool

An event of the automaton callback, at the state and time of `integrator` (anything with `u`
and `t`: DiffEq's integrator, or the instant of a fixed-step integration, `FixedInstant`),
where the components `fired` (indices among `rt.components`, with `direction(j)` the
direction each crossed in) changed truth between the ends of the root's bracket, the near
one at `t_before` (as `truths` recorded them, `AtomTruths`). The values of the models of the
components at the event; the atoms located where a condition changed (`_located_atoms!`),
and the crossings that fired, with their direction; a value that crossed where an earlier
instant left it without moving is no crossing (`_moving_crossings!`); and the automata settle
(`settle_instant!`, reinitializing with `reinit_alg` where a transition changes the algebraic
equations). `located`, `prior` and `rates` are the work arrays of the caller. Returns whether
any discrete state changed (false where nothing was located). The callback of an adaptive
integration (`build_automaton_callback`) and the fixed-step one (`step_automata!`) both
handle their events here.
"""
function _located_event!(rt::AutomatonRuntime, host, integrator, truths::AtomTruths, fired::Vector{Int}, direction,
                         t_before::Float64, jump, reinit_alg, located::Vector{LocatedValue}, prior::Vector{Int8},
                         rates::LocatedRates)
    comps = rt.components
    u, t = integrator.u, integrator.t
    nd = rt.n_diff
    # The values of the models of the components at the event
    fill!(rt.marked, false)
    @inbounds for j in fired
        rt.marked[comps.model[searchsortedlast(comps.offset, j - 1)]] = true
    end
    _evaluate_marked!(rt, view(u, 1:nd), view(u, (nd + 1):rt.n), t)
    empty!(located)
    empty!(prior)
    rt.located_events[] += 1
    @inbounds for j in fired
        d = direction(j)
        b = searchsortedlast(comps.offset, j - 1)
        km = comps.model[b]
        i = j - comps.offset[b]
        c = comps.condition[b]
        if c === nothing
            tr = comps.transition[b]::CompiledTransition
            rt.models[km].crossed[tr.crossing, i] = d
            _locate_value!(located, prior, rt.models[km], km, i, tr.guard::Vector{Float64}, Int8(d), t)
        else
            _located_atoms!(located, prior, rt, truths, km, i, c, t, t_before)
        end
    end
    # Values that crossed where an earlier instant left them without moving: no event
    _moving_crossings!(located, prior, rt, host, integrator, jump, _forget!(rates))
    isempty(located) && return false
    return settle_instant!(rt, host, integrator, reinit_alg; located, jump, rates)
end

"""
    _start_state_rates(host, n_diff) -> Vector{Float64}

The rates of the states (with mass or not) at the start, in the steady state the integration
starts in: for a `TDSAdapter`, those of its problem's steady state (`rate_buf`,
`init_tds_rates`: zero at an equilibrium, those of the oscillation for a sinusoidal
steady state's quantities).
A host adds its method (the Y-bus DAE, PowerfulYbus).
"""
_start_state_rates(adapter::TDSAdapter, n_diff::Int) = (length(adapter.rate_buf) == n_diff ||
    throw(DimensionMismatch("start rates: $(length(adapter.rate_buf)) state rates for $n_diff states"));
    adapter.rate_buf)

"""
    build_automaton_callback(rt::AutomatonRuntime, host, reinit_alg; start_abstol, jump=nothing, locate_extrema=true)
        -> Union{VectorContinuousCallback, Nothing}

The located callback of the automata of `rt` (`nothing` when no condition or crossing
is watched): its components are `automaton_conditions!`, whose atom truths it records
(`AtomTruths`); at an event, the atoms located where a condition changed its truth hold
as they changed (`_located_atoms!`) and the crossings fire in their direction, and the
automata settle (`settle_instant!`, reinitializing with `reinit_alg` where a transition
changes the algebraic equations); a value that crossed where an earlier instant left it
without moving is no crossing (`_moving_crossings!`). An event that changes nothing is no
discontinuity.
The automata start at the integration's initial state, which settles in rounds like any instant (`start!`, `StartResolve`):
the algebraic point is solved from the initial modes where their equations are violated at it by more than
`start_abstol` (the tolerance the initialization is validated to, `TDSConfig` `abstol`). `jump` (`nothing` for the plain
simulation) is told of each instant (`_instant_begins!`), the start included. With `locate_extrema`
(the default), a condition that holds or fails only between two of the points DiffEq evaluates it at
(its atoms crossing their surfaces and coming back, or one after the other, inside one step) is located
too, at its first change (`_step_extrema!`, see the header of step_extremum.jl).
"""
function build_automaton_callback(rt::AutomatonRuntime, adapter, reinit_alg; start_abstol::Float64, jump = nothing,
                                  locate_extrema::Bool = true)
    comps = rt.components
    comps.n == 0 && return nothing

    truths = AtomTruths(rt)
    extrema = locate_extrema ? StepExtrema(rt) : nothing
    # (the evaluators captured with their concrete type)
    condition = let rt = rt, evaluators = rt.evaluators, truths = truths, extrema = extrema, host = adapter, jump = jump
        function (out, u, t, integrator)
            _components!(out, rt, evaluators, u, t, integrator.tprev)
            _record_truths!(truths, rt, integrator, t)
            _step_extrema!(extrema, out, rt, evaluators, host, jump, u, t, integrator)
            return nothing
        end
    end

    affect! = let rt = rt, comps = comps, truths = truths, adapter = adapter, reinit_alg = reinit_alg, jump = jump
        located = LocatedValue[]
        prior = Int8[]
        rates = LocatedRates()
        fired = Int[]
        function (integrator, events)
            empty!(fired)
            @inbounds for j in eachindex(events)
                iszero(events[j]) || push!(fired, j)
            end
            _located_event!(rt, adapter, integrator, truths, fired, j -> events[j], prevfloat(integrator.t), jump,
                            reinit_alg, located, prior, rates) ||
                DiffEq.derivative_discontinuity!(integrator, false)
            return nothing
        end
    end

    initialize = let rt = rt, jump = jump, host = adapter, reinit_alg = reinit_alg, start_abstol = start_abstol
        function (cb, u, t, integrator)
            nd, n = rt.n_diff, rt.n
            # The start sets the discrete state; it moves the integrated state only where
            # `jump` re-solves what it carries (the sensitivities after memory starts). The
            # initialization of the DAE made those consistent with the discrete state the
            # simulation started with at its setup, which the start sets again: they are
            # re-solved only where the start changed the memory they are seeded with. Only
            # a moved state replaces the saved start
            before = copy(integrator.u)
            _instant_begins!(jump, integrator, NO_LOCATED_VALUES)
            seeds = _memory_seeds(jump)
            start!(rt, view(u, 1:nd), view(u, (nd + 1):n), t; jump, integrator, host,
                   resolve = StartResolve(host, integrator, reinit_alg, start_abstol))
            _instant_ends!(jump, integrator, _memory_seeds_changed(jump, seeds))
            # DiffEq saved t0 before the callbacks' initialization. A start that moved the state replaces that saved point
            # (the state before the start) with the settled one, and the start is no discontinuity either way: saving t0
            # again would leave two points at t0 with an empty interval between them, which no dense output evaluates
            integrator.u == before || isempty(integrator.sol.u) || copyto!(integrator.sol.u[1], integrator.u)
            DiffEq.derivative_discontinuity!(integrator, false)
            return nothing
        end
    end

    return erase_callback(DiffEq.VectorContinuousCallback(
        condition, affect!, comps.n;
        initialize = initialize,
        rootfind = DiffEq.RightRootFind,
        save_positions = (true, true),
        initializealg = DiffEq.NoInit()
    ))
end

# ============================================================================
# Automata in a fixed-step integration
# ============================================================================

"""
    FixedInstant

The instant of a fixed-step integration at which its automata settle, in the place of DiffEq's
integrator for the code the two share (`settle_instant!`, `_located_event!`): the state `u`
and time `t`, `reinit!(u, t)` (makes the variables without mass consistent, returning whether
it converged; `failed` records that it did not) and `sparsity()`, the pattern of the Jacobian of
the DAE, from which the rates of the values are evaluated (`_rate_dae`; taken when first needed).
"""
mutable struct FixedInstant
    u::Vector{Float64}
    t::Float64
    reinit!::Function
    sparsity::Function
    failed::Bool
end

# The rates of the plain simulation's values at an instant of a fixed-step integration
function _rate_dae(::Nothing, rt::AutomatonRuntime, host, instant::FixedInstant)
    rt.rates[] === nothing || return rt.rates[]
    rt.rates[] = StateRates(host, rt.n_diff, rt.mass, instant.sparsity())
    return rt.rates[]
end

# (the fixed-step integration reinitializes by its own algebraic solve, whose Jacobians are
# kept by the host; it carries no augmentation, `jump` always `nothing`)
function _reinitialize!(instant::FixedInstant, host, ::Nothing, ::Nothing)
    algebraic_equations_changed!(host)
    converged = instant.reinit!(instant.u, instant.t)
    instant.failed = !converged
    return converged
end
_reinit_failed(instant::FixedInstant) = instant.failed

"""
    FixedStepAutomata(rt::AutomatonRuntime, host, sparsity, reinit!)

The automata of the runtime `rt` in a fixed-step integration of the DAE of `host` (see
`dae_residuals!`; `sparsity()` the pattern of its Jacobian), whose algebraic solve is
`reinit!(u, t)`: the work arrays of `step_automata!` (the state at the start of a step and
at the near and far ends of the bracket of its first change of truth, the components of the
automaton callback at them and the truths of the atoms, `AtomTruths`), and the instant at
which they settle (`FixedInstant`). The rates of the values that the transitions hold at zero
by continuity (`Transition.continuity`) are evaluated on the DAE of `host`, as in an adaptive
integration.
"""
struct FixedStepAutomata{H}
    runtime::AutomatonRuntime
    host::H
    u0::Vector{Float64}
    u_hi::Vector{Float64}
    u_lo::Vector{Float64}
    c0::Vector{Float64}
    c_lo::Vector{Float64}
    c_hi::Vector{Float64}
    c_mid::Vector{Float64}
    truths::AtomTruths
    window::StepWindow
    instant::FixedInstant
    fired::Vector{Int}
    located::Vector{LocatedValue}
    prior::Vector{Int8}
    rates::LocatedRates
end

function FixedStepAutomata(rt::AutomatonRuntime, host, sparsity::Function, reinit!)
    nc = rt.components.n
    return FixedStepAutomata(rt, host, zeros(rt.n), zeros(rt.n), zeros(rt.n), zeros(nc), zeros(nc), zeros(nc), zeros(nc),
                             AtomTruths(rt), StepWindow(NaN, NaN),
                             FixedInstant(zeros(rt.n), NaN, reinit!, sparsity, false), Int[], LocatedValue[], Int8[],
                             LocatedRates())
end

# The time a fraction `s` of the step from `t` to `t_end` (of length `H`) is at: the end of the
# step exactly where `s` is 1
_at(t, H, t_end, s) = s == 1.0 ? t_end : t + s * H

# Whether component j changed sign from `a` to `b` (the components are never zero)
_changed(a, b, j) = signbit(a[j]) != signbit(b[j])

"""
    step_automata!(work::FixedStepAutomata, step!, u, t, h; time_tol) -> Bool

Advance `u = [x; y]` from `t` to `t + h` by an integration that locates nothing itself,
firing the transitions of the automata (`work.runtime`) where their conditions change
truth within the step, as the automaton callback does for an adaptive integration:
`step!(u, t, h)` takes one step of the integration (returning whether it succeeded). The
components of the automaton callback (`automaton_conditions!`) are evaluated at the ends of
the step; where one changes sign, the first change is bracketed to within `time_tol` by
steps from `t` of other lengths (false position on the component that changes first,
safeguarded towards bisection), the state taken at the far end of the bracket (where the
change has happened), the atoms and crossings that changed in the bracket located there
(`_located_event!`, the code of the adaptive callback, with the rates of the values on the
DAE of the host), and the automata settle there (`settle_fixed!`); the rest of the step
follows, with the same search. At the end of the step the automata settle again. Returns
`false` where `step!` or the algebraic solve failed (the state is then that of the failed
call).
"""
function step_automata!(work::FixedStepAutomata, step!, u, t, h; time_tol::Real)
    rt, window, truths = work.runtime, work.window, work.truths
    t_end = t + h
    while true
        H = t_end - t
        copyto!(work.u0, u)
        _components!(work.c0, rt, rt.evaluators, u, t, t)
        step!(u, t, H) || return false
        _components!(work.c_hi, rt, rt.evaluators, u, t_end, t)
        any(j -> _changed(work.c0, work.c_hi, j), eachindex(work.c0)) || break
        # The bracket (lo, hi] of the first change, as fractions of H
        lo, hi = 0.0, 1.0
        copyto!(work.c_lo, work.c0)
        copyto!(work.u_lo, work.u0)
        copyto!(work.u_hi, u)
        while (hi - lo) * H > time_tol
            # False position on the component that changes first between lo and hi, kept
            # within the middle 7/8 of the bracket
            s = hi
            for j in eachindex(work.c0)
                _changed(work.c_lo, work.c_hi, j) || continue
                a, b = work.c_lo[j], work.c_hi[j]
                s = min(s, lo + (hi - lo) * a / (a - b))
            end
            margin = (hi - lo) / 16
            s = clamp(s, lo + margin, hi - margin)
            copyto!(u, work.u0)
            step!(u, t, s * H) || return false
            _components!(work.c_mid, rt, rt.evaluators, u, _at(t, H, t_end, s), t)
            if any(j -> _changed(work.c0, work.c_mid, j), eachindex(work.c0))
                hi = s
                copyto!(work.c_hi, work.c_mid)
                copyto!(work.u_hi, u)
            else
                lo = s
                copyto!(work.c_lo, work.c_mid)
                copyto!(work.u_lo, u)
            end
        end
        # The event: the state where the change has happened, at the time of the evaluation
        # that found it; the components that changed in the bracket, and the truths of the
        # atoms at the two ends of the bracket
        copyto!(u, work.u_hi)
        t_event = _at(t, H, t_end, hi)
        window.tprev, window.t = _at(t, H, t_end, lo), t_event
        _components!(work.c_mid, rt, rt.evaluators, work.u_lo, _at(t, H, t_end, lo), NaN)
        _record_truths!(truths, rt, window, _at(t, H, t_end, lo))
        _components!(work.c_mid, rt, rt.evaluators, u, t_event, NaN)
        _record_truths!(truths, rt, window, t_event)
        empty!(work.fired)
        for j in eachindex(work.c0)
            _changed(work.c_lo, work.c_hi, j) && push!(work.fired, j)
        end
        instant = work.instant
        instant.u, instant.t, instant.failed = u, t_event, false
        _located_event!(rt, work.host, instant, truths, work.fired, j -> signbit(work.c_hi[j]) ? -1 : 1, _at(t, H, t_end, lo), nothing,
                        nothing, work.located, work.prior, work.rates)
        instant.failed && return false
        t = t_event
        t_end - t > time_tol || return true
    end
    return settle_fixed!(work, u, t_end)
end

"""
    settle_fixed!(work::FixedStepAutomata, u, t; event = false) -> Bool

Settle the automata of `work` at `(u, t)` for an integration that locates nothing itself, as
`settle_instant!` does for an adaptive one (the code is the same: the algebraic equations
the transitions change are solved again, `work.instant.reinit!`, and the values they hold at
zero by continuity are decided by their rate). The values decided as zero since a start at
rest are released at an `event`. Returns `false` where the algebraic solve failed.
"""
function settle_fixed!(work::FixedStepAutomata, u, t; event::Bool = false)
    instant = work.instant
    instant.u, instant.t, instant.failed = u, t, false
    settle_instant!(work.runtime, work.host, instant, nothing; event)
    return !instant.failed
end
