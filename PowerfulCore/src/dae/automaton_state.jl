# ============================================================================
# Automaton State - The Discrete State of a Model's Instances
# ============================================================================
#
# Where the modes, discrete quantities and memory of a model's automata live during a
# simulation (see Foundation/automata.jl for the metadata). One `AutomatonState` per
# model, held by its Float64 workspace. The workspaces of every other element type share
# its modes and discrete quantities: Jacobians see the
# discrete state the simulation is in. They hold their memory and values in their own
# element type, so that memory and the values of the automata can carry derivatives when
# the simulation is evaluated with dual numbers; the simulation copies the memory it
# sets into them (`sync_automaton_memory!`). The simulation (PowerfulTDS,
# solvers/automata.jl) starts the state and fires the transitions; the residual kernels
# read the discrete quantities and memory.
#
# The contract of the kernels (`calc_residuals!` and `calc_automaton_values!` of a model
# with automata): they compute what belongs to instance `i` from the variables, the
# parameters, the setpoints, the external parameters and instance `i`'s own discrete
# quantities and memory only, never from another instance's. A settle relies on it to
# consider only the instances that switched (PowerfulTDS, solvers/automaton_runtime.jl,
# `settle!`), and so does an augmentation that integrates derivatives through the
# switches of several instances at one instant (`_split_jump!`). The simulation
# warns of values found to break it when it is set up (`trace_value_reads`) and settles
# every instance of that model in every round, and such an augmentation refuses rates
# found to break it at simultaneous switches.
# ============================================================================

"""
    AffineRelation{T<:Number}

The affine form of a governed state's right side (round 19a, Phase D): `constant[i] +
sum_k coefficients[k, i] * var_vals[vars[k]][i]`, SPARSE over `vars`, the variables (indices among
the keys of `ws.var_vals`) the relation actually reads -- `check_governed`'s own admissibility
trace (`reads.pattern`, PowerfulTDS `governed_states.jl`) already knows this exactly, so the hot
path (`enforce_governed_affine!`) need not scan every variable of the model's workspace to find
the handful (often one or two) a relation's right side actually names. Computed once per
parameter set, not per residual evaluation (`refresh_affine_governed!`), in the workspace's own
element type `T`, so that a parameter carried as a Dual number leaves its partials in the
coefficients themselves, not only in the variables they multiply -- the same requirement
`AutomatonState`'s own `memory`/`values` fields hold for per-`T` recomputation (see its own
docstring). `vars` and `coefficients` are reassigned (not merely mutated elementwise) by
`refresh_affine_governed!`, since a retrace can in principle change which variables are read (an
admissible right side's kernel branching on a parameter, see its own docstring) -- hence a
`mutable struct`, so the same `AffineRelation` object stays referenced from `ws.automata.affine`
across a resize rather than needing the immutable `NamedTuple` that holds it rebuilt.
"""
mutable struct AffineRelation{T<:Number}
    constant::Vector{T}        # per instance
    vars::Vector{Int}          # the variables (indices among keys(ws.var_vals)) this relation reads, ascending
    coefficients::Matrix{T}    # length(vars) x n_instances (column i: the coefficients at instance i)
end

# `AffineRelation`s of `affine` promoted (elementwise `T.(...)`) to a workspace of element type `T`: the
# coefficients are a function of the model's PARAMETERS only (`check_governed`'s own admissibility
# proof), never of the variable `T` perturbs for an ordinary state/time Jacobian (ForwardDiff seeds the
# DAE's own x/y there, with every parameter still plain `Float64`) -- so for that workspace, which is
# every Jacobian `_gather_and_calc!` is actually called on today, the coefficients are the SAME numbers
# `affine` already holds, promoted with a zero partial, exactly as `AutomatonState(shared, T)` already
# promotes `memory` (`T.(shared.memory)`) rather than recomputing it. A workspace whose PARAMETERS are
# themselves Dual-seeded (a parameter-sensitivity Jacobian, not built anywhere in this codebase today)
# would need the coefficients genuinely recomputed on that workspace instead of promoted (promoting
# would silently zero a partial the coefficient should carry); no such workspace reads `ws.automata.affine`
# today, so this is a known, narrower gap left for whoever adds one, not silently miscomputed for the
# case that exists.
_fresh_affine(affine::NamedTuple, ::Type{T}) where T =
    NamedTuple{keys(affine)}(Tuple(AffineRelation(T.(r.constant), r.vars, T.(r.coefficients)) for r in values(affine)))

# The distinct value names (`ref`, a Symbol into `ws.automata.values`) that govern a state of some mode
# of some automaton of `automata`: the right sides that need an affine cache. A right side that is a
# plain number (not a Symbol) is excluded: `enforce_governed!`/`enforce_governed_affine!` write it
# directly, with no kernel or cache needed.
function _affine_relation_names(automata)
    names = Symbol[]
    for a in automata, mode in a.modes, (_, ref) in mode.states
        (ref isa Symbol && ref ∉ names) && push!(names, ref)
    end
    return names
end

# Fresh (NaN, uncomputed) affine relations for every value name `metadata.automata` governs a state by,
# of element type `T`, dense over all `n_vars` variables and `n_instances` instances until the first
# real `refresh_affine_governed!` narrows each to the variables it actually reads (empty for a model
# without a governed state, or one every one of whose governed states is a plain number). Never read
# in this dense, unrefreshed shape: `Governed`'s own constructor (PowerfulTDS `governed_states.jl`)
# calls `refresh_affine_governed!` before returning, so nothing reads `ws.automata.affine` before then.
function _affine_relations(::Type{T}, automata, n_instances::Int, n_vars::Int) where T
    names = _affine_relation_names(automata)
    isempty(names) && return NamedTuple()
    return NamedTuple{Tuple(names)}(ntuple(_ -> AffineRelation(fill(T(NaN), n_instances), collect(1:n_vars),
                                                               fill(T(NaN), n_vars, n_instances)),
                                           length(names)))
end

"""
    AutomatonState{Modes, Discrete, Memory, Values, Affine}

The discrete state of every instance of a model with automata:

- `modes`: automaton name => the mode of each instance, an index into the automaton's
  `modes`;
- `discrete`: discrete quantity => its value for each instance, what the mode holds;
- `memory`: memory name => its value for each instance;
- `values`: value name => its value for each instance, as the model's
  [`calc_automaton_values!`](@ref) last computed it (atoms, guards, holds, sets and
  memory starts, see `automaton_values`);
- `affine`: value name => the affine form of a governed-state relation by that value
  ([`AffineRelation`](@ref), round 19a Phase D), empty for a model that governs no state.

The simulation sets the modes, discrete quantities and memory at its start and where
transitions fire; before a start they are `0` (modes) and `NaN` (values), so that a
residual evaluated before one is not silently computed with a default. The residual
kernels read `discrete` and `memory`; nothing else writes them.

**The contract of the kernels.** `calc_residuals!` computes the rates and residuals of
instance `i`, and `calc_automaton_values!` its values, from the variables, parameters,
setpoints, external parameters and instance `i`'s own entries of `discrete` and `memory`
only: never from another instance's (a rate `r[1] * r[2]` of two instances' discrete
quantities breaks it). With the state fixed, a switch of one instance then changes that
instance's rates alone, and switches of several instances at one instant combine without
interacting, which an augmentation that integrates sensitivities through them relies on;
and at the same state a switch of one instance changes that instance's values alone, so
that a settle considers only the instances that switched. Every simulation checks the
contract for `calc_automaton_values!` when it is set up, from the trace of the values with
the discrete quantities and memory of each instance as inputs of their own (PowerfulTDS,
`trace_value_reads`), and, for a model whose values break it, warns and settles every
instance of the model in every round (correct, slower). Where switches of several
instances fire at one instant, a parameter-sensitivity augmentation checks the contract on
the rows with mass (each switch alone from the discrete state before the instant, and
each switch left out from the state after it, changes its own instance's rates only) and
refuses a kernel that breaks it there; an interaction that shows at neither end remains
the kernel author's responsibility.
"""
struct AutomatonState{Modes<:NamedTuple, Discrete<:NamedTuple, Memory<:NamedTuple, Values<:NamedTuple,
                      Affine<:NamedTuple}
    modes::Modes
    discrete::Discrete
    memory::Memory
    values::Values
    # The affine form of each governed-state relation (round 19a, Phase D), one per distinct value name a
    # mode governs a state by (`_affine_relation_names`): recomputed by `refresh_affine_governed!`,
    # written by PowerfulTDS whenever it retraces a model's governed states (setup, and after an event
    # that changes the parameters a relation reads -- the same trigger that already re-checks
    # admissibility, `Governed.recheck!`). Read by `enforce_governed_affine!` in the residual path in
    # place of a second kernel call. Of this workspace's own element type `T` (never shared/converted
    # across element types, see `AffineRelation`'s own docstring): empty for a model that governs no
    # state (the common case), costing nothing.
    affine::Affine
    # Round 19a: true while a fixed-time event's first settle round is deciding its transitions
    # (`settle!`, PowerfulTDS automaton_runtime.jl) -- `enforce_governed!` is a no-op while this holds,
    # so a guard reads the *retained* state against the event's new data, as the header of
    # `governed_states.jl` requires ("guards read first, then entered and surviving modes put their
    # states on their relations"): substituting a governed state with its relation *before* that
    # guard decides would make a guard whose own formula is the relation minus the state (REGCA1's
    # position atoms, and PIReset's release condition alike) always read exactly zero, so a release an
    # event should cause could never fire. Shared across every element-type workspace beside the
    # Float64 one (the same `Ref`, like `modes`/`discrete`): PowerfulTDS sets it once per settle, not
    # once per workspace.
    suspended::Base.RefValue{Bool}
end

"""
    AutomatonState()

The empty state of a model without automata.
"""
AutomatonState() = AutomatonState(NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), NamedTuple(), Ref(false))

_nan_columns(names, n::Int) = NamedTuple{Tuple(names)}(ntuple(_ -> fill(NaN, n), length(names)))

"""
    AutomatonState(metadata::ModelMetadata, n_instances, array_type) -> AutomatonState

New state (not started) for the automata of the model of `metadata`, after checking them
(`check_automata`). Automata are realized on CPU workspaces only: errors for another
`array_type`.
"""
function AutomatonState(metadata::ModelMetadata, n_instances::Int, array_type::Type{<:AbstractVector})
    if metadata.governs_any_state
        check_automata(metadata.name, metadata.automata, metadata.memory;
                       states = Symbol[v.name for v in metadata.vars if is_state_var(v) && is_internal(v)])
    else
        check_automata(metadata.name, metadata.automata, metadata.memory)
    end
    array_type === Vector || throw(ArgumentError(
        "$(metadata.name) declares automata, which are realized on CPU workspaces only; " *
        "its workspace has array type $array_type"))
    names = Tuple(a.name for a in metadata.automata)
    modes = NamedTuple{names}(ntuple(_ -> zeros(Int, n_instances), length(names)))
    discrete = _nan_columns(reduce(vcat, (discrete_quantities(a) for a in metadata.automata); init = Symbol[]),
                            n_instances)
    memory = _nan_columns([m.name for m in metadata.memory], n_instances)
    values = _nan_columns(automaton_values(metadata.automata, metadata.memory), n_instances)
    affine = _affine_relations(Float64, metadata.automata, n_instances, length(metadata.vars))
    return AutomatonState(modes, discrete, memory, values, affine, Ref(false))
end

"""
    AutomatonState(shared::AutomatonState, ::Type{T}) -> AutomatonState

The state of the automata of a workspace of element type `T` built beside the Float64
workspace that holds `shared`: the modes and discrete quantities of `shared` (the same
arrays), and memory, values and governed-state affine relations of element type `T` of its
own, the memory at the values of `shared`, the values and affine relations `NaN` until
computed (`refresh_affine_governed!` recomputes the latter in `T`, never a conversion of
`shared`'s own Float64 coefficients, which would discard a parameter Dual's partials).
`shared` itself for `T = Float64`.
"""
function AutomatonState(shared::AutomatonState, ::Type{T}) where T
    T === Float64 && return shared
    return AutomatonState(shared.modes, shared.discrete, map(m -> T.(m), shared.memory),
                          map(v -> fill(T(NaN), length(v)), shared.values),
                          _fresh_affine(shared.affine, T), shared.suspended)
end

"""
    calc_automaton_values!(ws::ResidualBuffers{T, AT}, model::StructArray, external::NamedTuple,
                           formulation, t)

Compute, for every instance of `model`, each value its automata name (see
`automaton_values`: the atoms `a - b` of its conditions, its crossing guards, and the
holds, sets and memory starts that are not numbers) into `ws.automata.values`, from the
variables gathered into `ws.var_vals` at time `t`, its parameters, setpoints, external
parameters and the discrete quantities and memory of `ws.automata`: the inputs its
`calc_residuals!` reads.

A model that declares automata implements it in the kernel pattern of
`calc_residuals!` (a data adapter calling a kernel that dispatches on `AT`). Out-of-service
instances are computed as the others; the simulation does not fire their transitions.
Which variables it reads must not depend on the variables, the time, the discrete
quantities and memory (no branch on them): the simulation traces the kernel for the
variables it reads, and gathers only those (`evaluate_automaton_values!`). It may branch
on a real parameter or a setpoint: it is traced at their values then, and again where an
event changes them, gathering every variable (PowerfulTDS `trace_value_reads`).

Both this kernel and the model's `calc_residuals!` follow the contract of
[`AutomatonState`](@ref): what belongs to instance `i` reads instance `i`'s own discrete
quantities and memory only.
"""
function calc_automaton_values!(ws, model::StructArray, external::NamedTuple, formulation, t)
    error("$(eltype(model)) declares automata but has no calc_automaton_values! method for " *
          "$(typeof(formulation)) to compute the values they name ($(join(keys(ws.automata.values), ", ")))")
end

"""
    enforce_governed!(ws, metadata::ModelMetadata, t)

A governed state is its relation's value at every point a reader decides on, not whatever the
integration carries between the discrete writes that land it there exactly (`_enter_governed!`,
PowerfulTDS `governed_states.jl`): for every state a mode of `metadata.automata` governs
(`Mode.states`), at every instance whose automaton `ws.automata.modes` currently holds in that mode,
overwrite the state's column of `ws.var_vals` with the mode's right side -- a number, or the matching,
already-computed column of `ws.automata.values`. Every reader of `ws.var_vals` downstream of this call
(a value, a guard, a residual) then sees the relation, not the drifted state.

Zero cost where no automaton of `metadata` governs a state (`governs_states`): nothing is read or
written. The branch this function makes, per instance, is on `ws.automata.modes[a.name][i]`, a
concrete `Int` -- automaton-state data, never a Dual or a tracer quantity -- so tracing through it
(PowerfulTDS `trace_value_reads`) never refuses it as a branch on a variable.

This is a *view*-level correction only: it never writes `x`/`y`, the DAE vectors the integrator
advances, so it never itself makes the automata settle again; what is actually *stored* there is
corrected only at the discrete instants `_enter_governed!` already writes (entry, event repair, and
the farewell write at a mode's exit), not here.
"""
function enforce_governed!(ws, metadata, t)::Nothing
    (metadata.governs_any_state && !ws.automata.suspended[]) || return nothing
    columns = ws.var_vals
    names = keys(columns)
    # `ws.automata.modes[a.name]` (a runtime Symbol index into a multi-field NamedTuple) allocates --
    # confirmed by isolating it (round 19a Phase D, coordinator's third follow-up): 80 bytes/call for
    # REGCA1's 8 automata, zero for the `findfirst`-into-the-key-tuple-then-index-the-value-tuple form
    # below, which gives the identical `Vector{Int}`. Computed once per call, not once per automaton.
    mode_names = keys(ws.automata.modes)
    mode_values = values(ws.automata.modes)
    value_names = keys(ws.automata.values)
    value_values = values(ws.automata.values)
    @inbounds for a in metadata.automata
        any(m -> !isempty(m.states), a.modes) || continue
        modes = mode_values[findfirst(==(a.name), mode_names)]
        for (km, mode) in enumerate(a.modes)
            isempty(mode.states) && continue
            for (state_name, ref) in mode.states
                j = findfirst(==(state_name), names)
                j === nothing && throw(ArgumentError(
                    "$(metadata.name): mode $(mode.name) of automaton $(a.name) governs $state_name, " *
                    "which is not a variable of its workspace"))
                col = columns[j]
                if ref isa Symbol
                    vcol = value_values[findfirst(==(ref), value_names)]
                    for i in eachindex(modes)
                        modes[i] == km && (col[i] = vcol[i])
                    end
                else
                    r = oftype(col[1], ref)
                    for i in eachindex(modes)
                        modes[i] == km && (col[i] = r)
                    end
                end
            end
        end
    end
    return nothing
end

"""
    refresh_affine_governed!(ws::ResidualBuffers{T}, model, external, formulation, read_vars) where T

Recompute the affine form (`AffineRelation`) of every governed-state relation of `metadata.automata`,
in `ws`'s own element type `T`: once at the current (zeroed) point for the constant term, then once
more per variable any relation reads (the union of `read_vars`'s values), shifted by `one(T)` with
every other variable held at zero, for that variable's coefficient in every relation it affects --
exact, not a finite-difference approximation, because `check_governed` has already proven every such
relation affine (degree at most 1) in the variables it reads, with parameter-dependent coefficients
only (`PowerfulTDS governed_states.jl`). `ws.var_vals` is saved before and restored after (this never
leaves the workspace's variables at the zeroed probe point); `t` is not needed, since the same
admissibility proof refuses a relation that reads the time.

`read_vars[ref]` (a `Dict{Symbol, Vector{Int}}`, one entry per key of `ws.automata.affine`): the
variables (indices among the keys of `ws.var_vals`) relation `ref` actually reads, from
`check_governed`'s own trace of the model's values (`reads.pattern`, `PowerfulTDS governed_states.jl`
`_value_read_variables`) -- reused, not retraced here. Each relation's `vars`/`coefficients` are
REASSIGNED (not merely overwritten elementwise) to this sparsity, since `AffineRelation` is a
`mutable struct` precisely so an earlier, differently-sized relation can be replaced without
rebuilding the `NamedTuple` (`ws.automata.affine`) that holds it: the hot path
(`enforce_governed_affine!`) then sums over `length(vars)` terms, typically one or two, never over
every variable of the model's workspace, regardless of how many other variables the model gathers.

Zero cost (and a no-op) where `ws.automata.affine` is empty (a model that governs no state, or one
whose every right side is a plain number). Called once when `Governed` is built, and again whenever a
fixed-time event retraces a model's governed states (`Governed.recheck!`, the same trigger that already
re-checks admissibility) -- never in the residual hot path itself (`enforce_governed_affine!` reads the
cache there instead of recomputing it).
"""
function refresh_affine_governed!(ws, model, external, formulation, read_vars)
    isempty(ws.automata.affine) && return nothing
    T = eltype(first(values(ws.var_vals)))
    cols = values(ws.var_vals)
    colvec = collect(cols)
    n = ws.n_instances
    saved = map(copy, cols)
    foreach(col -> fill!(col, zero(T)), cols)
    calc_automaton_values!(ws, model, external, formulation, zero(T))
    for (ref, relation) in pairs(ws.automata.affine)
        copyto!(relation.constant, getfield(ws.automata.values, ref))
    end
    # Resize every relation to its own sparsity before perturbing anything, so the loop below can fill
    # each relation's (now correctly shaped) coefficients column by column as it goes.
    for (ref, relation) in pairs(ws.automata.affine)
        vars = get(read_vars, ref, Int[])
        relation.vars = vars
        relation.coefficients = fill(zero(T), length(vars), n)
    end
    wanted = sort!(unique(j for vars in values(read_vars) for j in vars))
    for j in wanted
        col = colvec[j]
        fill!(col, one(T))
        calc_automaton_values!(ws, model, external, formulation, zero(T))
        for (ref, relation) in pairs(ws.automata.affine)
            k = findfirst(==(j), relation.vars)
            k === nothing && continue
            vcol = getfield(ws.automata.values, ref)
            @inbounds for i in 1:n
                relation.coefficients[k, i] = vcol[i] - relation.constant[i]
            end
        end
        fill!(col, zero(T))
    end
    foreach((col, s) -> copyto!(col, s), cols, saved)
    return nothing
end

"""
    enforce_governed_affine!(ws, metadata)

The residual path's substitution (round 19a, Phase D): for every state a mode of `metadata.automata`
governs, at every instance whose automaton `ws.automata.modes` currently holds in that mode, overwrite
the state's column of `ws.var_vals` with its relation's *affine form* (`AffineRelation`,
[`refresh_affine_governed!`](@ref)) -- a dot product over `relation.vars`, the SPARSE set of variables
the relation actually reads (often one or two, never every variable the model gathers, regardless of
how large the model's own variable set is) -- instead of calling `calc_automaton_values!` to get a
fresh value. `calc_residuals!` (what this feeds, in `_gather_and_calc!`) never reads
`ws.automata.values`, only `ws.var_vals`, so there is no self-reference to resolve the way
[`evaluate_automaton_values!`](@ref)'s two-pass shape must for the values/condition path: one pass over
the governed states suffices here, and it costs one multiply-add per variable the relation reads
instead of a whole kernel invocation.

Zero cost where no automaton of `metadata` governs a state, or while the automata are suspended
(`ws.automata.suspended[]`, the same guard `enforce_governed!` itself makes, for the same reason: a
guard deciding a release during a fixed-time event's first settle round must read the *retained* state,
not a value this substitution would already have put on the relation). Where used, this REPLACES the
`calc_automaton_values!` call the residual path used to pay (see the private round 19a design notes):
the values/condition path (`evaluate_automaton_values!`) is
untouched and keeps calling the kernel twice, where exactness (not speed) is required.
"""
function enforce_governed_affine!(ws, metadata)::Nothing
    (metadata.governs_any_state && !ws.automata.suspended[]) || return nothing
    columns = ws.var_vals
    cols = values(columns)
    names = keys(columns)
    n = ws.n_instances
    # `ws.automata.modes[a.name]` (a runtime Symbol index into a multi-field NamedTuple) allocates --
    # found and fixed in round 19a Phase D's coordinator follow-up (see the private design notes):
    # isolating every piece of this function in its own real
    # calling context (not a standalone microbenchmark, which gave misleading numbers for this very
    # line) found this one line responsible for the entire ~80 bytes/call this function cost REGCA1,
    # explaining a real share of its measured ~60 ns/call gap over main's un-enforced residual path.
    # `findfirst` into the key tuple, then indexing the value tuple, gives the identical `Vector{Int}`
    # with zero allocation; computed once per call, not once per automaton.
    mode_names = keys(ws.automata.modes)
    mode_values = values(ws.automata.modes)
    # Same fix, same reason, for the affine relations: `getfield(ws.automata.affine, ref)` (also a
    # runtime Symbol index into a NamedTuple) accounted for the rest of this function's allocation once
    # the modes lookup above was fixed -- `ws.automata.affine` has only one field for REGCA1 today, but
    # nothing guarantees that in general, and the fix is the same regardless of field count.
    affine_names = keys(ws.automata.affine)
    affine_values = values(ws.automata.affine)
    @inbounds for a in metadata.automata
        any(m -> !isempty(m.states), a.modes) || continue
        modes = mode_values[findfirst(==(a.name), mode_names)]
        for (km, mode) in enumerate(a.modes)
            isempty(mode.states) && continue
            for (state_name, ref) in mode.states
                j = findfirst(==(state_name), names)
                j === nothing && throw(ArgumentError(
                    "$(metadata.name): mode $(mode.name) of automaton $(a.name) governs $state_name, " *
                    "which is not a variable of its workspace"))
                col = cols[j]
                if ref isa Symbol
                    relation = affine_values[findfirst(==(ref), affine_names)]
                    c0 = relation.constant
                    vars = relation.vars
                    C = relation.coefficients
                    for i in 1:n
                        modes[i] == km || continue
                        v = c0[i]
                        for (k, jj) in enumerate(vars)
                            v += C[k, i] * cols[jj][i]
                        end
                        col[i] = v
                    end
                else
                    r = oftype(col[1], ref)
                    for i in 1:n
                        modes[i] == km && (col[i] = r)
                    end
                end
            end
        end
    end
    return nothing
end

"""
    evaluate_automaton_values!(ws, gather::Tuple, model, external, x, y, formulation, t, metadata)

Gather the variables of `gather` from `x` and `y` into `ws` and compute the values of
the automata of `model` at time `t` into `ws.automata.values` (`calc_automaton_values!`).
`gather` holds entries of `ws.var_gather_tuple` (`variable_gather`): every variable the
model's `calc_automaton_values!` reads, or more. The variables outside it keep whatever
`ws.var_vals` held, which the kernel then must not read; the simulation gathers exactly
those its trace of the kernel found it reads (PowerfulTDS, `trace_value_reads`).

Where `metadata.automata` governs a state (`governs_states`), the kernel is run **twice**: once to
get a fresh relation for every value that is not a mirror of the governed state itself (this first
pass never reads the governed state, `check_governed`'s own admissibility proof), then
[`enforce_governed!`](@ref) substitutes the governed state's column with that relation, then the
kernel is run again so that every value that *is* a function of the governed state (e.g. a position
atom comparing it with its own bound) reads the enforced value, not the first pass's stale one. This
is `_enter_governed!`'s own "write, then re-evaluate the values" shape, applied on every call instead
of only at a discrete write. For a model that governs no state this is the original, single call.
"""
@inline function evaluate_automaton_values!(ws, gather::Tuple, model, external, x, y, formulation, t, metadata)
    gather_variables!(gather, x, y)
    calc_automaton_values!(ws, model, external, formulation, t)
    if metadata.governs_any_state
        enforce_governed!(ws, metadata, t)
        calc_automaton_values!(ws, model, external, formulation, t)
    end
    return nothing
end
