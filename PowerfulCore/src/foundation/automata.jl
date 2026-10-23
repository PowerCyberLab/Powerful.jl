# ============================================================================
# Automata: Discrete Behaviour of a Model
# ============================================================================
#
# The Julia mirror of the semantic core's automaton (Mode, Transition, Held,
# InitialMode, Automaton, Memory).
# A model declares its automata and memory in `ModelMetadata.automata` and
# `ModelMetadata.memory`; the simulation keeps each instance's mode, discrete
# quantities and memory in the model's workspace (`AutomatonState`).
#
# Guards, conditions, holds, sets and memory starts are expressions of the model.
# The model computes each one it names, for every instance, in its
# `calc_automaton_values!` method (the kernel pattern of `calc_residuals!`); the
# metadata refers to them by name. A condition is a boolean combination of atoms,
# each a named value compared with zero, so that a simulation can locate the first
# instant a condition holds at the zero of the atom whose comparison changed.
# ============================================================================

"""
    ValueRef

What the metadata of an automaton reads from the model: the name of a value the model's
`calc_automaton_values!` computes for every instance, or a number.
"""
const ValueRef = Union{Symbol, Float64}

_value_ref(v::Symbol) = v
_value_ref(v::Real) = Float64(v)

"""
    Condition

A boolean combination of [`Atom`](@ref)s: `a & b`, `a | b` and `!a` combine them.
"""
abstract type Condition end

"""
    Atom(value::Symbol, relation::Symbol) <: Condition

The comparison of the value `value` of the model (computed by its
`calc_automaton_values!` as `a - b` for a comparison of `a` with `b`) with zero:
`relation` is `:>`, `:>=`, `:<` or `:<=`.

# Example
```julia
Atom(:y_minus_upper, :>=) & Atom(:rate, :>)     # y >= upper and rate > 0
```
"""
struct Atom <: Condition
    value::Symbol
    relation::Symbol
    function Atom(value::Symbol, relation::Symbol)
        relation in (:>, :>=, :<, :<=) ||
            throw(ArgumentError("Atom $value: relation must be :>, :>=, :< or :<=, got :$relation"))
        new(value, relation)
    end
end

"""
    AllOf(terms) <: Condition

Holds where every one of `terms` holds (`a & b`).
"""
struct AllOf <: Condition
    terms::Vector{Condition}
end

"""
    AnyOf(terms) <: Condition

Holds where any of `terms` holds (`a | b`).
"""
struct AnyOf <: Condition
    terms::Vector{Condition}
end

"""
    Negation(term) <: Condition

Holds where `term` does not (`!a`).
"""
struct Negation <: Condition
    term::Condition
end

Base.:&(a::Condition, b::Condition) = AllOf(Condition[a, b])
Base.:|(a::Condition, b::Condition) = AnyOf(Condition[a, b])
Base.:!(a::Condition) = Negation(a)

"""
    atoms(condition) -> Vector{Atom}

The atoms `condition` combines, in order of appearance (with repeats).
"""
atoms(a::Atom) = Atom[a]
atoms(c::Union{AllOf, AnyOf}) = reduce(vcat, (atoms(t) for t in c.terms); init = Atom[])
atoms(c::Negation) = atoms(c.term)
atoms(::Nothing) = Atom[]

"""
    Mode(name; holds = [], relations = [], states = [])

One mode of an automaton and what each discrete quantity holds while the model is in it:
`holds` maps each discrete quantity to a `ValueRef`. The simulation realizes numbers only
and refuses a computed value (`check_automata`). `relations` (branch current => the residual that replaces its relation in
this mode) mirrors the core's `Mode.relations`; Powerful has no multi-port branch
contract, and its simulation refuses a mode that names any.

`states` (differential state => `ValueRef`) governs states: while the model is in the mode, each
state it names is held on the right side, a number or a value the model computes in
`calc_automaton_values!` (the core's `Mode.relations` of a state, solved for the state: in
the semantic core's `bounded`, the bound). The transition that enters the mode writes the right side to
the state at that instant (a projection along the governed state alone, `project_states!`),
and the model's rate for the state, in the mode, is the rate of the right side
(the model's, gated by what the mode holds). The right side must be affine in differential states
with constant coefficients and read nothing else of the model's evolving quantities
(`check_governed` of PowerfulTDS refuses the others, naming the mode and the state).

# Example
```julia
Mode(:upper; holds = [:gate => 0], states = [:y => :upper_bound])   # y = upper_bound while here
```
"""
struct Mode
    name::Symbol
    holds::Vector{Pair{Symbol, ValueRef}}
    relations::Vector{Pair{Symbol, Symbol}}
    states::Vector{Pair{Symbol, ValueRef}}
end

# The states of a mode that governs none: one empty vector for all of them (a model's metadata is built at every query,
# and a mode that governs nothing costs nothing for it); no mode's `states` is ever changed
const NO_GOVERNED_STATES = Pair{Symbol, ValueRef}[]

Mode(name::Symbol; holds = Pair{Symbol, ValueRef}[], relations = Pair{Symbol, Symbol}[], states = NO_GOVERNED_STATES) =
    Mode(name, Pair{Symbol, ValueRef}[q => _value_ref(v) for (q, v) in holds],
         Pair{Symbol, Symbol}[r for r in relations],
         isempty(states) ? NO_GOVERNED_STATES : Pair{Symbol, ValueRef}[s => _value_ref(v) for (s, v) in states])

"""
    Transition(from, to; when = nothing, guard = nothing, direction = nothing,
               period = nothing, sets = [], continuity = [], crossing_continuity = [],
               continuity_writers = [])

`from` to `to` either where the condition `when` holds (a condition transition: it fires
at the first instant `when` holds while the model is in `from`, at the start, right
after another transition or where it becomes true), or where the value `guard` crosses
zero in `direction` (`:rising` from below, `:falling` from above; a crossing
transition). `period` (a periodic crossing) mirrors the core's `Transition.period`; the
simulation refuses it.

`sets` maps memory to the `ValueRef` it takes when the transition fires. Every right
side reads the quantities as they were before the transition (the assignments are
simultaneous) and the time of the transition. A transition may return to the mode it
leaves only where it sets memory.

`continuity` (`surface => values`) names, for a value `surface` the transition compares
(an atom of its condition, or its guard), the values continuity holds at zero where the
transition fires on the zero of `surface`: each vanishes wherever the model, in `from`, is
on that surface, and the transition leaves the algebraic state where it is, so each is
zero after it too (derived from the semantic core's continuity rules: REGCA1's `dq` on
entering its clamp at `v = Volim`). A simulation decides the atoms over them at that
instant by their rate in the discrete state the transition enters, not by their value,
which is zero there to the rounding of the reinitialization.

`crossing_continuity` has the form of `continuity` and the same meaning, for facts that hold
only where the transition fires at a crossing located while the simulation was integrating
its source mode `from` (the first round of an instant): the source mode was settled before
the instant, so none of its other exits held. REGCA1's high-voltage clamp, released where
`dq` reaches zero, holds `Vt - Volim` at zero only from a settled clamp: its exit to `floor`
not holding there is what makes the clamp's reactive power and the clear mode's agree. A
clamp entered earlier in the same instant has no such guarantee. The values are decided by
their rate at that crossing, as with `continuity`, and never when the transition fires in a
later round of the instant, at a start, at an event's reinitialization or at any other
instant. `continuity` holds however the source mode was entered (derived the same way).

`continuity_writers` (`value => automata`) names, for a value held by `continuity` or
`crossing_continuity` on this transition, the automata of the same instance whose firing
voids that value's hold: the value reads a quantity one of them writes (REECA1's PI clip
holds values that read the limit `Iqmax`, which two other automata can write), so that the
fact derived without their firing no longer holds where one fires. A simulation holds the
value only where none of the listed automata fires in the round of the transition (the
algebraic point is solved again between rounds, so an earlier round's writer is in the point the
round reads), and decides it by its own value otherwise, as it does a value no marker names.
A value not listed has no writers.

A marker of any kind for a value that reads a quantity the network solves (a variable without
mass) is not applied in a round in which another transition, of any instance, changes a quantity
delivered to the network: the derivation reads the inputs of the model as unchanged across the
transition, and what other devices do reaches them through the network. Each listed value is one the transition holds, each automaton is one of the
model, and the transition's own automaton is not listed (it fires at most once per round, so
it is never a value's writer).
"""
struct Transition
    from::Symbol
    to::Symbol
    condition::Union{Nothing, Condition}
    guard::Union{Nothing, Symbol}
    direction::Union{Nothing, Symbol}
    period::Union{Nothing, Symbol}
    sets::Vector{Pair{Symbol, ValueRef}}
    continuity::Vector{Pair{Symbol, Vector{Symbol}}}
    crossing_continuity::Vector{Pair{Symbol, Vector{Symbol}}}
    continuity_writers::Vector{Pair{Symbol, Vector{Symbol}}}
end

function Transition(from::Symbol, to::Symbol; when::Union{Nothing, Condition} = nothing,
                    guard::Union{Nothing, Symbol} = nothing, direction::Union{Nothing, Symbol} = nothing,
                    period::Union{Nothing, Symbol} = nothing, sets = Pair{Symbol, ValueRef}[],
                    continuity = Pair{Symbol, Vector{Symbol}}[],
                    crossing_continuity = Pair{Symbol, Vector{Symbol}}[],
                    continuity_writers = Pair{Symbol, Vector{Symbol}}[])
    Transition(from, to, when, guard, direction, period,
               Pair{Symbol, ValueRef}[m => _value_ref(v) for (m, v) in sets],
               Pair{Symbol, Vector{Symbol}}[s => collect(Symbol, vs) for (s, vs) in continuity],
               Pair{Symbol, Vector{Symbol}}[s => collect(Symbol, vs) for (s, vs) in crossing_continuity],
               Pair{Symbol, Vector{Symbol}}[v => collect(Symbol, as) for (v, as) in continuity_writers])
end

is_condition(tr::Transition) = tr.condition !== nothing

"""
    InitialMode(mode; when = nothing)

One clause of where an automaton starts: `mode` where the condition `when` holds at the
start, or unconditionally where it is `nothing` (the last clause).
"""
struct InitialMode
    condition::Union{Nothing, Condition}
    mode::Symbol
end

InitialMode(mode::Symbol; when::Union{Nothing, Condition} = nothing) = InitialMode(when, mode)

"""
    Held(mode; when)

Where the condition `when` holds, the automaton is in `mode` and no transition of it
fires.
"""
struct Held
    condition::Condition
    mode::Symbol
end

Held(mode::Symbol; when::Condition) = Held(when, mode)

"""
    Automaton(name; modes, transitions, initial, held = [], at_rest = [])

Modes, the transitions between them, where the model starts (`initial`: the mode of the
first `InitialMode` clause whose condition holds, the last clause having none) and the
`Held` clauses. The discrete quantities of an automaton are those its modes hold.

`at_rest` names the values its conditions compare that a start from equilibrium holds at
zero: the initialization makes them zero, whatever the data, flags and start modes (the
rates of states it puts at rest; derived from the semantic core's continuity rules). A
simulation decides the atoms over them at the start as the exact zero they are, not by
the rounding of what is computed there.

# Example: anti-windup lag
```julia
reach_upper = Atom(:y_minus_upper, :>=) & Atom(:rate, :>)
reach_lower = Atom(:y_minus_lower, :<=) & Atom(:rate, :<)
Automaton(:LAG;
    modes = [Mode(:free; holds = [:gate => 1]), Mode(:upper; holds = [:gate => 0]),
             Mode(:lower; holds = [:gate => 0])],
    transitions = [Transition(:free, :upper; when = reach_upper),
                   Transition(:free, :lower; when = reach_lower),
                   Transition(:upper, :free; when = Atom(:rate, :<)),
                   Transition(:lower, :free; when = Atom(:rate, :>))],
    initial = [InitialMode(:upper; when = reach_upper), InitialMode(:lower; when = reach_lower),
               InitialMode(:free)])
```
"""
struct Automaton
    name::Symbol
    modes::Vector{Mode}
    transitions::Vector{Transition}
    initial::Vector{InitialMode}
    held::Vector{Held}
    at_rest::Vector{Symbol}
end

Automaton(name::Symbol; modes, transitions, initial, held = Held[], at_rest = Symbol[]) =
    Automaton(name, collect(Mode, modes), collect(Transition, transitions),
              collect(InitialMode, initial), collect(Held, held), collect(Symbol, at_rest))

"""
    MemoryQuantity(name; initial)

A value the model keeps between transitions: set by the `sets` of a transition, held
until another sets it, and starting at `initial` (a `ValueRef`, which reads no memory).
"""
struct MemoryQuantity
    name::Symbol
    initial::ValueRef
end

MemoryQuantity(name::Symbol; initial) = MemoryQuantity(name, _value_ref(initial))

"""
    discrete_quantities(automaton::Automaton) -> Vector{Symbol}

The discrete quantities `automaton` speaks for: those its first mode holds (every mode
holds the same ones, see `check_automata`).
"""
discrete_quantities(a::Automaton) = isempty(a.modes) ? Symbol[] : Symbol[first(h) for h in first(a.modes).holds]

"""
    governed_states(automata) -> Vector{Symbol}

The differential states the modes of `automata` govern (`Mode.states`), each once, in order of
appearance.
"""
governed_states(automata) = unique!(Symbol[s for a in automata for m in a.modes for (s, _) in m.states])

"""
    governs_states(automata) -> Bool

Whether any mode of `automata` governs a state (`Mode.states`), without building the list (`governed_states`).
"""
governs_states(automata) = any(a -> any(m -> !isempty(m.states), a.modes), automata)

"""
    automaton_values(automata, memory) -> Vector{Symbol}

For the `Automaton`s `automata` and the `MemoryQuantity`s `memory` of a model, the names of every value the model computes for its automata (atoms, crossing guards,
holds, the right sides of governed states, sets and memory starts that are not numbers), each once, in order of appearance:
the values its `calc_automaton_values!` writes.
"""
function automaton_values(automata, memory)
    names = Symbol[]
    add!(v::Symbol) = (v in names || push!(names, v); nothing)
    add!(::Float64) = nothing
    for a in automata
        for m in a.modes, (_, v) in m.holds
            add!(v)
        end
        for m in a.modes, (_, v) in m.states
            add!(v)
        end
        for tr in a.transitions
            foreach(at -> add!(at.value), atoms(tr.condition))
            tr.guard === nothing || add!(tr.guard)
            for (_, v) in tr.sets
                add!(v)
            end
        end
        for c in a.initial
            foreach(at -> add!(at.value), atoms(c.condition))
        end
        for h in a.held
            foreach(at -> add!(at.value), atoms(h.condition))
        end
    end
    for m in memory
        add!(m.initial)
    end
    return names
end

"""
    check_automata(model::Symbol, automata, memory; states = Symbol[])

Throw an `ArgumentError` unless the automata and memory of `model` (collections of
`Automaton` and `MemoryQuantity`, as `ModelMetadata` holds them) are well formed and
realizable by Powerful's simulation, stating what is wrong and where:

- what the semantic core requires: distinct modes that each hold the same discrete
  quantities; transitions between
  declared modes, returning to the mode they leave only where they set memory, with a
  condition or a crossing in a known direction; initial clauses naming declared modes,
  the last and only the last without a condition; held clauses naming declared modes;
  memory that starts at a value, and sets that set memory (memory no transition sets warns,
  being a parameter that stays what it starts at);
- what Powerful needs: at least one initial clause (a start mode); names that are
  distinct across automata, discrete quantities and memory;
- that the automata a value's `continuity_writers` name are automata of the model other than the transition's
  own, and that the value is one the transition holds (`continuity` or `crossing_continuity`);
- that the values continuity holds at zero name what the automata compare: a
  transition's `continuity` and `crossing_continuity` surfaces are values it compares (an atom of its condition, or
  its guard), and the values held there and an automaton's `at_rest` are values the
  automata of the model compare;
- and refuses what it does not realize: holds that are computed values rather than
  numbers (a discrete quantity changes only when a transition fires, and the core does
  not yet define a hold that reads the model), periodic crossings (`Transition.period`)
  and mode-dependent branch relations (`Mode.relations`);
- that the states its modes govern (`Mode.states`) are differential states of the model (`states`: the names of
  its internal differential states), each governed once by a mode and by the modes of one automaton only, and
  that the right side of each is a finite number or a value (`ValueRef`). What the right side may read of the
  model is checked where the values are traced (PowerfulTDS `check_governed`).
"""
const NO_STATE_NAMES = Symbol[]

function check_automata(model::Symbol, automata, memory; states = NO_STATE_NAMES)
    where_(a) = "automaton $(a.name) of $model"
    fail(msg) = throw(ArgumentError(msg))
    governing = nothing   # (the automaton that governs each state, a Dict made by the first mode that governs one)

    names = [a.name for a in automata]
    allunique(names) || fail("$model declares automata with the same name: $(join(names, ", "))")
    memory_names = [m.name for m in memory]
    allunique(memory_names) || fail("$model declares memory with the same name: $(join(memory_names, ", "))")
    discrete = Symbol[]
    set_somewhere = Set{Symbol}()

    for a in automata
        mode_names = [m.name for m in a.modes]
        isempty(mode_names) && fail("$(where_(a)) has no modes")
        allunique(mode_names) || fail("$(where_(a)) declares a mode twice: $(join(mode_names, ", "))")
        spoken = discrete_quantities(a)
        for m in a.modes
            held_here = [first(h) for h in m.holds]
            allunique(held_here) || fail("mode $(m.name) of $(where_(a)) holds a discrete quantity twice")
            Set(held_here) == Set(spoken) ||
                fail("mode $(m.name) of $(where_(a)) holds $(join(held_here, ", ")), and mode " *
                     "$(first(a.modes).name) holds $(join(spoken, ", ")); every mode says what every " *
                     "discrete quantity is")
            computed = [q for (q, v) in m.holds if v isa Symbol]
            isempty(computed) ||
                fail("mode $(m.name) of $(where_(a)) holds a computed value for " *
                     "$(join(computed, ", ")): a discrete quantity changes only when a transition " *
                     "fires, and a hold that reads the model would change between transitions, so " *
                     "Powerful's simulation realizes holds that are numbers only")
            isempty(m.relations) ||
                fail("mode $(m.name) of $(where_(a)) replaces the relations of " *
                     "$(join(first.(m.relations), ", ")) (Mode.relations): Powerful has no multi-port " *
                     "branch contract, so a mode-dependent branch relation is not realized")
            isempty(m.states) && continue
            allunique(first.(m.states)) || fail("mode $(m.name) of $(where_(a)) governs a state twice")
            for (s, v) in m.states
                s in states ||
                    fail("mode $(m.name) of $(where_(a)) governs $s, which is not a differential state of $model " *
                         "(Mode.states governs the differential states of the model: " *
                         "$(isempty(states) ? "it has none" : join(states, ", ")))")
                v isa Float64 && !isfinite(v) &&
                    fail("mode $(m.name) of $(where_(a)) governs $s at $v, which is not a finite number")
                governing === nothing && (governing = Dict{Symbol, Symbol}())
                owner = get!(governing, s, a.name)
                owner == a.name ||
                    fail("the state $s of $model is governed by the modes of two automata, $owner and $(a.name): " *
                         "two modes active together would govern it")
            end
        end
        for q in spoken
            (q in discrete || q in memory_names) &&
                fail("discrete quantity $q of $(where_(a)) is also a discrete quantity of another " *
                     "automaton or memory of $model")
            push!(discrete, q)
        end
        declared(mode, what) = mode in mode_names ||
            fail("$what names mode $mode, which is not a mode of $(where_(a)) ($(join(mode_names, ", ")))")

        for tr in a.transitions
            what = "the transition $(tr.from) -> $(tr.to) of $(where_(a))"
            declared(tr.from, what)
            declared(tr.to, what)
            tr.from == tr.to && isempty(tr.sets) &&
                fail("$what changes nothing; a transition returns to the mode it leaves only where it sets memory")
            if is_condition(tr)
                (tr.guard === nothing && tr.direction === nothing && tr.period === nothing) ||
                    fail("$what has both a condition and a crossing guard; it is one or the other")
            else
                tr.guard === nothing && fail("$what has neither a condition nor a crossing guard")
                tr.direction in (:rising, :falling) ||
                    fail("$what crosses its guard $(tr.guard) in direction $(repr(tr.direction)); a " *
                         "guard is crossed :rising or :falling")
                tr.period === nothing ||
                    fail("$what is a periodic crossing (guard $(tr.guard), period $(tr.period)): " *
                         "periodic crossings are not realized by Powerful's simulation")
            end
            for (m, _) in tr.sets
                m in memory_names ||
                    fail("$what sets $m, which is not memory of $model; a transition sets memory")
                push!(set_somewhere, m)
            end
            allunique(first.(tr.sets)) || fail("$what sets a memory twice")
        end

        isempty(a.initial) &&
            fail("$(where_(a)) states no initial mode, and Powerful's simulation needs one to start")
        for (k, c) in enumerate(a.initial)
            declared(c.mode, "initial clause $k of $(where_(a))")
            last = k == length(a.initial)
            if last && c.condition !== nothing
                fail("the last initial clause of $(where_(a)) has a condition; a start where no " *
                     "condition holds needs a mode, so the list ends with a clause without one")
            elseif !last && c.condition === nothing
                fail("initial clause $k ($(c.mode)) of $(where_(a)) has no condition; only the last " *
                     "clause may, since the ones after it could never be reached")
            end
        end
        for h in a.held
            declared(h.mode, "a held clause of $(where_(a))")
        end
    end

    for m in memory
        m.name in set_somewhere ||
            @warn "Memory $(m.name) of $model is set by no transition: a value nothing changes is a parameter, " *
                  "and it stays what it starts at"
    end

    compared = Set{Symbol}()
    for a in automata
        for tr in a.transitions
            foreach(at -> push!(compared, at.value), atoms(tr.condition))
            tr.guard === nothing || push!(compared, tr.guard)
        end
        foreach(c -> foreach(at -> push!(compared, at.value), atoms(c.condition)), a.initial)
        foreach(h -> foreach(at -> push!(compared, at.value), atoms(h.condition)), a.held)
    end
    automaton_names = Set{Symbol}(a.name for a in automata)
    for a in automata
        for tr in a.transitions
            what = "the transition $(tr.from) -> $(tr.to) of $(where_(a))"
            held_here = Set{Symbol}(v for pairs in (tr.continuity, tr.crossing_continuity) for (_, vs) in pairs for v in vs)
            for (v, writers) in tr.continuity_writers
                v in held_here ||
                    fail("$what lists the writers of $v (continuity_writers), which it does not hold at zero " *
                         "(neither its continuity nor its crossing_continuity names it)")
                for w in writers
                    w in automaton_names ||
                        fail("$what lists $w as a writer of $v (continuity_writers), and $model has no automaton $w")
                    w == a.name &&
                        fail("$what lists its own automaton $w as a writer of $v (continuity_writers): an " *
                             "automaton fires at most once per round, so it is never a value's writer")
                end
            end
            own = Set{Symbol}(at.value for at in atoms(tr.condition))
            tr.guard === nothing || push!(own, tr.guard)
            for (field, pairs) in ((:continuity, tr.continuity), (:crossing_continuity, tr.crossing_continuity)),
                (surface, held_values) in pairs
                surface in own ||
                    fail("$what holds values at zero ($field) where it fires on the zero of $surface, which it does not compare")
                for v in held_values
                    v in compared ||
                        fail("$what holds $v at zero ($field) where it fires on the zero of $surface, and no automaton of $model compares $v")
                end
            end
        end
        for v in a.at_rest
            v in compared ||
                fail("$(where_(a)) names $v as held at zero at rest, and no automaton of $model compares $v")
        end
    end
    return nothing
end

"""
    has_automata(metadata) -> Bool

Whether the model of `metadata` declares automata.
"""
has_automata(metadata) = !isempty(metadata.automata)
