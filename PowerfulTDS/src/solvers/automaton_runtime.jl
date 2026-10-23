# ============================================================================
# Automaton Runtime: Starting Automata and Firing Their Transitions
# ============================================================================
#
# The semantics of the automata of the models of a simulation (PowerfulCore
# Foundation/automata.jl, mirroring the semantic core), on the
# `AutomatonState` of each model's Float64 workspace (shared by its Dual workspaces):
#
# - the start: memory at its initial values, each automaton in the mode of its first
#   initial clause that holds, the discrete quantities at what the modes hold, then a
#   settle;
# - a settle, at the start and at every instant something may have changed (an event
#   located by the automaton callback, a fixed-time event, a transition), in rounds: in a
#   round every transition that holds, of every automaton of every in-service instance, fires
#   on the point before the round (at most one per automaton; two that hold at once from one
#   mode are refused), among
#     * a held clause that holds: the automaton is in its mode, and none of its
#       transitions fires;
#     * a condition transition whose condition holds;
#     * a crossing transition whose guard was located crossing zero in its direction
#       (first round of a located event only; a crossing is an instant);
#   the sets of a transition read the values from before it; after each round the
#   discrete quantities follow the modes, and where the round changed the equations
#   without mass the host solves the algebraic point again from the modes entered
#   (settle_instant!, automata.jl) before the next round reads it. Rounds repeat until no
#   transition holds at the point (the settled point). An instance that takes part in more
#   than AUTOMATON_MAX_ROUNDS rounds at one instant (rounds within ZENO_ULPS ulps of each other,
#   which located events that advance by rounding are) is Zeno and throws.
#
# An integration that carries more than the DAE state passes a `jump` and its
# integrator to `start!` and `settle!`, which tell it where memory
# starts (`_memory_started!`) and where a round of transitions fires
# (`_transitions_firing!`); `nothing` for the plain simulation.
#
# The values of the automata of a model are computed from the variables its
# `calc_automaton_values!` reads, traced once when the simulation is set up
# (`trace_value_reads`): the simulation gathers those alone.
#
# The DiffEq side (locating atoms and crossings, settling after events, algebraic
# consistency) is in automata.jl.
# ============================================================================

using PowerfulCore: Atom, AllOf, AnyOf, Negation, Condition, Transition, Automaton, ValueRef,
    MemoryQuantity, atoms, discrete_quantities, is_condition, has_automata, AutomatonState,
    evaluate_automaton_values!, sync_automaton_memory!, automaton_memory_synced, calc_automaton_values!, ResidualBuffers,
    variable_gather, with_float_columns, gather_variables!, is_algeb_res, is_state_res, is_shared_res, UNROUTED,
    enforce_governed!
using PowerfulCore.DAEInfra: value_buffers, _beside
using SparseConnectivityTracer: TracerSparsityDetector, TracerLocalSparsityDetector, MissingPrimalError,
    AbstractTracer
using ADTypes: jacobian_sparsity

"""
    AUTOMATON_MAX_ROUNDS

The default budget of rounds (see `settle!`) one instance takes part in at one instant, the
`max_rounds` of the runtime; one more needed is Zeno: an automaton that switches forever
without time advancing, which throws.
"""
const AUTOMATON_MAX_ROUNDS = 100

"""
    ZENO_ULPS

Switches of one instance at times within this many ulps of each other are at one
instant (see `_count_rounds!`): the rounds of a settle, which share their time, and the
switches of located events that advance by rounding only, as a relay does whose switch
drives its own input back across its threshold. Switches a step of the integration
apart are never within it.
"""
const ZENO_ULPS = 8

# Relations of atoms, as codes
const _GT, _GE, _LT, _LE = Int8(1), Int8(2), Int8(3), Int8(4)
_relation_code(r::Symbol) = r === :> ? _GT : r === :>= ? _GE : r === :< ? _LT : _LE

"""
    compare(v, relation::Int8) -> Bool

Whether the comparison `v relation 0` holds.
"""
@inline compare(v, rel::Int8) = rel === _GT ? v > 0 : rel === _GE ? v >= 0 : rel === _LT ? v < 0 : v <= 0

"""
    signed_truth(v, relation::Int8)

A number whose sign is the truth of the comparison `v relation 0`: positive exactly where
it holds, negative where it does not, never zero. It is `v` (for `>` and `>=`) or `-v`
(for `<` and `<=`), except at `v = 0`, where it is the smallest subnormal on the side the
comparison puts zero (true for `>=` and `<=`, false for `>` and `<`).

This is what the automaton callback locates: a comparison changes its truth exactly where
this changes sign, and DiffEq, which never reports a crossing from an exact zero, reports
it also where the value starts at zero and leaves it (a start at a limit).
"""
@inline function signed_truth(v, rel::Int8)
    s = ifelse(rel === _GT || rel === _GE, v, -v)
    tiny = nextfloat(zero(s))
    return ifelse(iszero(s), ifelse(rel === _GE || rel === _LE, tiny, -tiny), s)
end

"""
    CompiledCondition

A condition over the atom table of a model, compiled for evaluation without allocation
(`signed_condition`): `code`, a program in postfix order over a stack of numbers, each
instruction `op | (k << 3)` with `k` an atom (its index in the atom table) or, for
`_CONSTANT`, `1` for true and `0` for false:

- `_LEAF k`: push atom `k`;
- `_MIN_LEAF k`, `_MAX_LEAF k`: replace the top by its least (greatest) with atom `k`;
- `_MIN`, `_MAX`: replace the top two by their least (greatest);
- `_NEGATE`: negate the top;
- `_CONSTANT`: push `Inf` (true: an empty `AllOf`) or `-Inf` (false: an empty `AnyOf`).

An `AllOf` is its first term, then the least with each further term (`_MIN_LEAF` for an
atom), an `AnyOf` the greatest, a `Negation` its term negated. `atoms` holds the atoms it
reads, each once; `stack` is the scratch of the depth the program reaches.
"""
struct CompiledCondition
    code::Vector{Int}
    atoms::Vector{Int}
    stack::Vector{Float64}
end

const _LEAF, _MIN_LEAF, _MAX_LEAF, _MIN, _MAX, _NEGATE, _CONSTANT = 0, 1, 2, 3, 4, 5, 6

function CompiledCondition(c::Condition, atom_index::Dict{Atom, Int})
    code = Int[]
    depth = 0
    # Push the value of `c` onto a stack `level - 1` deep
    function emit!(c, level)
        depth = max(depth, level)
        if c isa Atom
            push!(code, _LEAF | (atom_index[c] << 3))
        elseif c isa Negation
            emit!(c.term, level)
            push!(code, _NEGATE)
        else
            all = c isa AllOf
            isempty(c.terms) && return push!(code, _CONSTANT | (Int(all) << 3))
            for (q, term) in enumerate(c.terms)
                if q == 1
                    emit!(term, level)
                elseif term isa Atom
                    push!(code, (all ? _MIN_LEAF : _MAX_LEAF) | (atom_index[term] << 3))
                else
                    emit!(term, level + 1)
                    push!(code, all ? _MIN : _MAX)
                end
            end
        end
    end
    emit!(c, 1)
    return CompiledCondition(code, unique(atom_index[at] for at in atoms(c)), zeros(depth))
end

"""
    signed_condition(c::CompiledCondition, leaf) -> Float64

A number whose sign is the truth of `c`, `leaf(k)` that of atom `k` (never zero, as
`signed_truth`): the least of its terms for `AllOf`, the greatest for `AnyOf`, negated
for `Negation` (its `code`). It is continuous where the atoms are, and changes sign
exactly where the condition changes its truth, at the zero of an atom whose truth
decides it there (see the header of automata.jl).
"""
@inline function signed_condition(c::CompiledCondition, leaf::F) where {F}
    stack = c.stack
    top = 0
    @inbounds for instruction in c.code
        op, k = instruction & 7, instruction >> 3
        if op == _LEAF
            top += 1
            stack[top] = leaf(k)::Float64
        elseif op == _MIN_LEAF
            stack[top] = min(stack[top], leaf(k)::Float64)
        elseif op == _MAX_LEAF
            stack[top] = max(stack[top], leaf(k)::Float64)
        elseif op == _MIN
            top -= 1
            stack[top] = min(stack[top], stack[top + 1])
        elseif op == _MAX
            top -= 1
            stack[top] = max(stack[top], stack[top + 1])
        elseif op == _NEGATE
            stack[top] = -stack[top]
        else
            top += 1
            stack[top] = k == 1 ? Inf : -Inf
        end
    end
    return @inbounds stack[1]
end

"""
    condition_holds(c::CompiledCondition, truth) -> Bool

Whether `c` holds with atom `k` true exactly where `truth(k)` is: the sign of
`signed_condition` with the atoms at `±1`.
"""
condition_holds(c::CompiledCondition, truth::F) where {F} =
    signed_condition(c, k -> truth(k)::Bool ? 1.0 : -1.0) > 0

"""
    signed_block!(scratch, c::CompiledCondition, signed::Matrix{Float64}, n) -> column

`signed_condition` of `c` for each of `n` instances at once, atom `k` of instance `i`
at `signed[i, k]`: its `code` run with a stack of columns (`scratch[:, level]`), each
instruction a loop over the instances, a leaf with the atom term after it in one. Returns the
result's column: a view of the first `n` rows of the first column of `scratch`, or, for a
condition that is one atom, of that atom's column of `signed` itself.
"""
function signed_block!(scratch::Matrix{Float64}, c::CompiledCondition, signed::Matrix{Float64}, n::Int)
    code = c.code
    # (a lone atom is its column of `signed`, which nothing writes while the callback reads it)
    length(code) == 1 && code[1] & 7 == _LEAF && return view(signed, 1:n, code[1] >> 3)
    top = 0
    fused = false   # the atom term after a leaf, taken with it
    @inbounds for (j, instruction) in enumerate(code)
        fused && (fused = false; continue)
        op, k = instruction & 7, instruction >> 3
        if op == _LEAF
            top += 1
            follower = j < length(code) ? code[j + 1] : _CONSTANT
            if follower & 7 == _MIN_LEAF
                k2 = follower >> 3
                @simd for i in 1:n
                    scratch[i, top] = min(signed[i, k], signed[i, k2])
                end
                fused = true
            elseif follower & 7 == _MAX_LEAF
                k2 = follower >> 3
                @simd for i in 1:n
                    scratch[i, top] = max(signed[i, k], signed[i, k2])
                end
                fused = true
            else
                @simd for i in 1:n
                    scratch[i, top] = signed[i, k]
                end
            end
        elseif op == _MIN_LEAF
            @simd for i in 1:n
                scratch[i, top] = min(scratch[i, top], signed[i, k])
            end
        elseif op == _MAX_LEAF
            @simd for i in 1:n
                scratch[i, top] = max(scratch[i, top], signed[i, k])
            end
        elseif op == _MIN
            top -= 1
            @simd for i in 1:n
                scratch[i, top] = min(scratch[i, top], scratch[i, top + 1])
            end
        elseif op == _MAX
            top -= 1
            @simd for i in 1:n
                scratch[i, top] = max(scratch[i, top], scratch[i, top + 1])
            end
        elseif op == _NEGATE
            @simd for i in 1:n
                scratch[i, top] = -scratch[i, top]
            end
        else
            top += 1
            v = k == 1 ? Inf : -Inf
            @simd for i in 1:n
                scratch[i, top] = v
            end
        end
    end
    return view(scratch, 1:n, 1)
end

# A value an automaton reads: a column of the model's values, or a number
@inline _ref_value(r::Float64, i::Int) = r
@inline _ref_value(r::Vector{Float64}, i::Int) = @inbounds r[i]
_resolve(ref::Float64, values::NamedTuple) = ref
_resolve(ref::Symbol, values::NamedTuple) = getfield(values, ref)

# Whether `e` is what SparseConnectivityTracer raises where a kernel branches on a tracer: a
# `TypeError` for a tracer used as a `Bool` (`if`, `?:`, `&&`, a comparison used as a
# condition) and its `MissingPrimalError` for a predicate that needs a value (`iszero`,
# `isnan`). Any other exception is the kernel's own and is not for the trace to take.
_branches_on_tracer(e) = e isa MissingPrimalError || (e isa TypeError && e.got isa AbstractTracer)

"""
    ValueReads

What the values of the automata of a model (its `calc_automaton_values!`) read of its
variables and of its discrete state, traced from the kernel (`trace_value_reads`):
`variables`, those the simulation gathers to compute the values (the ones any value of any
instance reads, in declaration order; every variable where `branches`); `pattern`, by
value and instance: entry `((i - 1) * n_values + k, (j - 1) * n_instances + i′)` is whether
value `k` (in the order of the `AutomatonState`'s `values`) of instance `i` reads variable
`j` (in declaration order) of instance `i′`; and `state`, the same of the discrete
quantities and then the memory (in the order of the `AutomatonState`): entry
`((i - 1) * n_values + k, (q - 1) * n_instances + i′)` is whether value `k` of instance `i`
reads quantity `q` of instance `i′`.

`massless`, once the simulation knows its mass matrix (`mark_massless!`): for each column
of `pattern`, whether the variable instance it stands for is a row without mass (an
algebraic variable, or a state whose mass is zero), so that `reads_massless` tells whether
a value reads one. `lines` are what the trace found of values reading another instance's
discrete state (the lines of the warning, `value_reads`).

`cross` is whether any value reads a discrete quantity or memory of another instance than
its own, outside the contract of the kernels (`AutomatonState`): a settle then considers
every instance of the model in every round (`_expand_cross`), not those that switched.
`branches` is whether the kernel branches on a real parameter or a setpoint, so that the
trace holds for their values at the moment of tracing: the pattern is traced again where
an event changes them (`retrace_value_reads!`) and every variable is gathered, whichever
the values read.

The simulation gathers `variables` alone to compute the values
(`evaluate_automaton_values!`); grazing reads from `pattern` the entries of the state each
value reads (`_value_dependencies`); a settle reads from `state` the values a transition
moves by changing what they read (`_quantity_changed!`).
"""
mutable struct ValueReads
    variables::Vector{Symbol}
    pattern::SparseMatrixCSC{Bool, Int}
    state::SparseMatrixCSC{Bool, Int}
    cross::Bool
    branches::Bool
    trace::Function
    massless::BitVector
    lines::Vector{String}
    degree::Vector{Int}
end

"""
    trace_value_degrees(ws::ResidualBuffers, model, external, formulation, at_values::Bool) -> Vector{Int}

For each value of the automata of `model`, the degree of the polynomial it is in the variables and the time, for
every instance (`PolyDegree`, value_degree.jl: `NOT_POLYNOMIAL` where it is none), traced from its
`calc_automaton_values!` with the discrete quantities and the memory (constant within a step) as constants.
`at_values` as in `trace_value_reads`: where a kernel compares parameters the trace is made with them at their
values. A kernel the tracer does not know an operation of is taken as not polynomial, which only costs its values
the exact path of `step_extremum.jl`.
"""
function trace_value_degrees(ws::ResidualBuffers, model, external, formulation, at_values::Bool)
    n, names = ws.n_instances, collect(keys(ws.var_vals))
    n_vars, n_values = length(names), length(ws.automata.values)
    n_state = length(ws.automata.discrete) + length(ws.automata.memory)
    T = PolyDegree
    v = [PolyDegree(j <= n_vars * n ? 1 : 0) for j in 1:(n_vars + n_state) * n]
    degrees = try
        traced = ResidualBuffers(ws, T; setpoints_at_values = at_values)
        for (j, column) in enumerate(values(traced.var_vals))
            copyto!(column, 1, v, (j - 1) * n + 1, n)
        end
        _set_discrete_inputs!(traced, v, n_vars * n)
        calc_automaton_values!(traced, at_values ? model : with_float_columns(model, T),
                               at_values ? external : _float_external(external, T), formulation, PolyDegree(1))
        [maximum(x -> x isa PolyDegree ? x.d : 0, column; init = 0) for column in values(traced.automata.values)]
    catch e
        e isa PolyBranch && !at_values && return trace_value_degrees(ws, model, external, formulation, true)
        e isa InterruptException && rethrow()
        fill(NOT_POLYNOMIAL, n_values)
    end
    return degrees
end

"""
    trace_value_reads(ws::ResidualBuffers, model, external, formulation) -> ValueReads

What the values of the automata of `model` read of its variables (`ValueReads`), traced
from its `calc_automaton_values!` in `formulation` on a copy of its Float64 workspace
`ws` of the tracer's element type: the sparsity pattern of the values with respect to the
gathered variables and to the discrete quantities and memory (transitions set them), each
an indexed input per instance, with SparseConnectivityTracer's global detector. The time
enters as a tracer too, not indexed.

The real parameters (events modify them), the external parameters (copies of fields of other
models, which events modify) and the setpoints first enter as tracers as well:
a global pattern holds at every value of them, so the variables it lists are those the
values read whatever the state, the time, the events and the transitions. A kernel that
branches on a parameter, an external parameter or a setpoint (`if`, `?:`, `&&`, a comparison;
an `ifelse` reads both of its branches) cannot be traced so. It is traced with them
at their values then, as the residuals are (`trace_algebraic_switches`), which holds until
an event changes them (`ValueReads` `branches`, `retrace_value_reads!`). A kernel that
branches on a variable, the time, a discrete quantity or memory cannot be traced either
way, and throws. Only SparseConnectivityTracer's refusal of a branch takes that fallback
(`_branches_on_tracer`); any other error of the kernel surfaces as it is.

Where a value of an instance reads a discrete quantity or memory of another, which breaks
the contract of the kernels (`AutomatonState`) on which a settle considers only the
instances that switched, the settle takes every instance of the model in every
round (`ValueReads` `cross`) and `value_reads` warns, from `lines`. The pattern with respect to the discrete quantities and memory
is kept (`ValueReads` `state`): where a transition changes what a value reads, the value
moves at the instant by that change, not by the trajectory (`_quantity_changed!`).
"""
function trace_value_reads(ws::ResidualBuffers, model, external, formulation)
    n, names = ws.n_instances, collect(keys(ws.var_vals))
    n_vars, n_values = length(names), length(ws.automata.values)
    n_state = length(ws.automata.discrete) + length(ws.automata.memory)
    trace = function (at_values::Bool)
        values! = function (out, v)
            T = eltype(v)
            traced = ResidualBuffers(ws, T; setpoints_at_values = at_values)
            for (j, column) in enumerate(values(traced.var_vals))
                copyto!(column, 1, v, (j - 1) * n + 1, n)
            end
            _set_discrete_inputs!(traced, v, n_vars * n)
            calc_automaton_values!(traced, at_values ? model : with_float_columns(model, T),
                                   at_values ? external : _float_external(external, T), formulation, zero(T))
            for (k, column) in enumerate(values(traced.automata.values)), i in 1:n
                out[(i - 1) * n_values + k] = column[i]
            end
            return nothing
        end
        return jacobian_sparsity(values!, zeros(n_values * n), zeros((n_vars + n_state) * n),
                                 TracerSparsityDetector())
    end
    branches = false
    pattern = try
        trace(false)
    catch e
        _branches_on_tracer(e) || rethrow()
        branches = true
        try
            trace(true)
        catch e
            _branches_on_tracer(e) || rethrow()
            throw(ErrorException(
                "The values of the automata of $(eltype(model)) could not be traced for the variables they read " *
                "(`calc_automaton_values!` with SparseConnectivityTracer's global detector, which refuses branches on " *
                "a variable, the time, a discrete quantity or memory): $(sprint(showerror, e))"))
        end
    end
    outputs = keys(ws.automata.values)
    split = function (pattern)
        variables = pattern[:, 1:(n_vars * n)]
        read = [any(view(variables, :, ((j - 1) * n + 1):(j * n))) for j in 1:n_vars]
        return (variables = branches ? names : [names[j] for j in 1:n_vars if read[j]],
                pattern = SparseMatrixCSC{Bool, Int}(variables),
                state = SparseMatrixCSC{Bool, Int}(pattern[:, (n_vars * n + 1):end]),
                lines = _instance_reads(pattern, n_vars * n, ws, model, "value", outputs))
    end
    read = split(pattern)
    degree = () -> trace_value_degrees(ws, model, external, formulation, branches)
    again = () -> merge(split(trace(true)), (; degree = trace_value_degrees(ws, model, external, formulation, true)))
    return ValueReads(read.variables, read.pattern, read.state, !isempty(read.lines), branches, again, BitVector(),
                      read.lines, degree())
end

"""
    _value_param_reads(ws, model, external, formulation, value_name::Symbol) -> Vector{Set{Symbol}}

The numeric fields of `model` (its own parameters) that the value `value_name` of its automata reads, for every
instance: traced from its `calc_automaton_values!`, perturbing the variables and the fields together, at two
points distinct from their real values and from each other (as a three-phase residual path's own analogous
`_parameter_reads` does for a constraint's residual, which this mirrors), so that a coefficient zero at one of them -- a variable the
right side reads at its rest value, as much as a parameter whose own effect the kernel happens to cancel there --
is not missed at the other. A parallel trace to, not a reuse of, `ValueReads` (`reads.pattern`): that one traces
variable reads for a different purpose (refusing a right side that is not affine, or reads a variable it must
not) and is not indexed by field, which is what the repair of a governed state after an event needs (`check_governed`'s
admissibility already proves a parameter is the only channel an event can move a governed right side through, so this
is its entire dependence on anything a `FieldModification` could touch). `TracerLocalSparsityDetector` resolves a
kernel that branches on a field or a variable at the concrete point of each trace (unlike the global detector
`trace_value_reads` tries first), so no fallback is needed here.
"""
function _value_param_reads(ws::ResidualBuffers, model, external, formulation, value_name::Symbol)
    n = ws.n_instances
    columns = StructArrays.components(model)
    fields = [f for (f, col) in pairs(columns) if eltype(col) <: AbstractFloat]
    read = [Set{Symbol}() for _ in 1:n]
    isempty(fields) && return read
    var_names = collect(keys(ws.var_vals))
    n_vars = length(var_names)
    n_in = (n_vars + length(fields)) * n
    for scale in (1.0, 7.0)
        values! = function (out, v)
            T = eltype(v)
            traced = ResidualBuffers(ws, T)
            for (j, column) in enumerate(values(traced.var_vals))
                copyto!(column, 1, v, (j - 1) * n + 1, n)
            end
            traced_columns = map(c -> eltype(c) <: AbstractFloat ? copyto!(similar(c, T), c) : c, columns)
            for (q, f) in enumerate(fields)
                copyto!(getproperty(traced_columns, f), 1, v, (n_vars + q - 1) * n + 1, n)
            end
            # (the model's own element type, not `get_base_model_type(model){T}`: a model with a type parameter of
            # its own besides the numeric one -- as `BadRight{K, T}` in the admissibility tests -- is not
            # `get_base_model_type(model){T}`, which only ever substitutes the last slot; `with_float_columns`
            # avoids the same trap the same way, by keeping the model's own element type)
            traced_model = StructArray{eltype(model)}(traced_columns)
            calc_automaton_values!(traced, traced_model, external, formulation, zero(T))
            copyto!(out, getproperty(traced.automata.values, value_name))
            return nothing
        end
        point = Vector{Float64}(undef, n_in)
        for (j, column) in enumerate(values(ws.var_vals)), i in 1:n
            point[(j - 1) * n + i] = column[i] + scale * (abs(column[i]) + 1.0)
        end
        for (q, f) in enumerate(fields)
            col = getproperty(columns, f)
            for i in 1:n
                point[(n_vars + q - 1) * n + i] = col[i] + scale * (abs(col[i]) + 1.0)
            end
        end
        pattern = jacobian_sparsity(values!, zeros(n), point, TracerLocalSparsityDetector())
        for (q, f) in enumerate(fields), i in 1:n
            pattern[i, (n_vars + q - 1) * n + i] && push!(read[i], f)
            # This trace only ever reads off the diagonal block (instance i's own fields): a value that read
            # another instance's field (as `ValueReads.cross`/`trace_value_reads` allows for a variable read)
            # would have that dependence missed here, under-reporting `param_reads` and letting the repair skip
            # a relation an event did move. No value in this codebase reads a field this way; throw loudly
            # rather than silently miss it if one ever does.
            for i2 in 1:n
                i2 == i && continue
                pattern[i, (n_vars + q - 1) * n + i2] && throw(ErrorException(
                    "$value_name of $(nameof(eltype(model))) reads $(f) of another instance ($i2, not its own " *
                    "$i): the repair of a governed state after an event (`_value_param_reads`) only traces a " *
                    "value's dependence on its own instance's fields, which this value does not have"))
            end
        end
    end
    return read
end

# What a trace that does not branch has for a retrace: nothing to trace again
_no_retrace() = throw(ErrorException(
    "The values of these automata do not branch on parameters or setpoints: there is nothing to trace again"))

# `reads` with the arrays it shares and its own `massless` (set by `mark_massless!`), and no
# closure over the workspace it was traced on: for a trace that does not branch
_shared_reads(reads::ValueReads) =
    ValueReads(reads.variables, reads.pattern, reads.state, reads.cross, reads.branches, _no_retrace, BitVector(),
               reads.lines, reads.degree)

"""
    value_reads(traces::StructuralTraces, name, ws, model, external, formulation) -> ValueReads

What the values of the automata of model `name` read (`trace_value_reads`), traced by the
first call for the problem of `traces` and taken from it by the next ones where the trace
does not branch on parameters or setpoints (`ValueReads` `branches`): the pattern then
holds whatever their values are, so it is the problem's, not the solve's. A trace that
branches depends on the values at the trace and is made in every solve. Warns, from the
lines of the trace, where a value reads another instance's discrete state, as every solve
of the problem would.
"""
function value_reads(traces::StructuralTraces, name::Symbol, ws::ResidualBuffers, model, external, formulation)
    stored = lock(() -> get(traces.reads, name, nothing), traces.lock)
    reads = if stored === nothing
        traced = trace_value_reads(ws, model, external, formulation)
        traced.branches || lock(() -> traces.reads[name] = _shared_reads(traced), traces.lock)
        traced
    else
        _shared_reads(stored)
    end
    isempty(reads.lines) || _warn_instance_reads(model, reads.lines)
    return reads
end

"""
    retrace_value_reads!(reads::ValueReads, model) -> Bool

Trace the values of the automata of `model` again for what they read, with the parameters
and setpoints at the values they have now, where `reads` is of a kernel that branches on
them (`ValueReads` `branches`; nothing otherwise): update `reads` and return whether what
the values read of the variables or of the discrete state changed. Warns where a value
newly reads another instance's discrete state (`ValueReads` `cross`).
"""
function retrace_value_reads!(reads::ValueReads, model)
    reads.branches || return false
    read = reads.trace()
    changed = read.pattern != reads.pattern || read.state != reads.state
    if !isempty(read.lines) && !reads.cross
        _warn_instance_reads(model, read.lines)
    end
    reads.pattern, reads.state, reads.cross = read.pattern, read.state, !isempty(read.lines)
    reads.degree = read.degree
    return changed
end

"""
    mark_massless!(reads::ValueReads, ws::ResidualBuffers, mass, n_diff) -> reads

Record in `reads` which variable instances of the model of the workspace `ws` are rows
without mass of the DAE whose mass diagonal is `mass` (`n_diff` the length of `x`): the
entries of `[x; y]` its variables are gathered from (`ws.var_gather_tuple`). Structural: it
holds whatever the pattern of `reads` becomes.
"""
function mark_massless!(reads::ValueReads, ws::ResidualBuffers, mass::AbstractVector, n_diff::Int)
    n = ws.n_instances
    reads.massless = BitVector(iszero(mass[(is_state ? 0 : n_diff) + Int(addrs[i])])
                               for (addrs, _, is_state) in ws.var_gather_tuple for i in 1:n)
    return reads
end

"""
    reads_massless(reads::ValueReads, k::Int, n_values::Int) -> Bool

Whether value `k` (of `n_values`) of any instance reads a variable without mass
(`mark_massless!`): its rate then needs the rows without mass completed, and is that of
the rows with mass alone (`F/M`) where it does not.
"""
function reads_massless(reads::ValueReads, k::Int, n_values::Int)
    rows = rowvals(reads.pattern)
    for c in findall(reads.massless), p in nzrange(reads.pattern, c)
        mod1(rows[p], n_values) == k && return true
    end
    return false
end

# The external parameters `external` of a model with the arrays of floats in element type `T`
# (the others as they are), for a trace that takes them as the parameters are: a kernel
# that branches on one is traced at their values then
_float_external(external::NamedTuple, ::Type{T}) where T =
    map(v -> eltype(v) <: AbstractFloat ? copyto!(similar(v, T), v) : v, external)

# The discrete quantities and then the memory of the automata of `traced` (a workspace of
# a tracer's element type) from `v`, after its first `offset` entries: a block of one
# entry per instance for each, in the order of the `AutomatonState`
function _set_discrete_inputs!(traced::ResidualBuffers, v, offset::Int)
    n = traced.n_instances
    columns = (values(traced.automata.discrete)..., values(traced.automata.memory)...)
    for (q, column) in enumerate(columns)
        copyto!(column, 1, v, offset + (q - 1) * n + 1, n)
    end
    return nothing
end

"""
    _instance_reads(pattern, offset, ws, model, what, outputs) -> Vector{String}

Where an output of an instance of `model` reads a discrete quantity or memory of another,
a line saying so: `pattern` is that of the outputs (`outputs[k]` of instance `i` at row
`(i - 1) * length(outputs) + k`) with respect to inputs whose columns after the first
`offset` are the discrete quantities and then the memory of the `AutomatonState` of its
workspace `ws` (quantity `q` of instance `i′` at column `offset + (q - 1) * n + i′`).
Each names the output (`what` and its name), the quantity and both instances.
"""
function _instance_reads(pattern::SparseMatrixCSC, offset::Int, ws::ResidualBuffers, model,
                         what::AbstractString, outputs)
    n, n_out = ws.n_instances, length(outputs)
    n_discrete = length(ws.automata.discrete)
    rows = rowvals(pattern)
    lines = String[]
    quantities = (keys(ws.automata.discrete)..., keys(ws.automata.memory)...)
    keys_ = get_key(model)
    for col in (offset + 1):size(pattern, 2), p in nzrange(pattern, col)
        row = rows[p]
        i, k = divrem(row - 1, n_out) .+ 1
        q, j = divrem(col - offset - 1, n) .+ 1
        i == j && continue
        kind = q <= n_discrete ? "discrete quantity" : "memory"
        push!(lines, "$what $(outputs[k]) of instance $i ($(join(string.(keys_[i]), ", "))) reads $kind " *
                     "$(quantities[q]) of instance $j ($(join(string.(keys_[j]), ", ")))")
    end
    return lines
end

# The warning for values that read another instance's discrete state: the first of `lines`
const _INSTANCE_READS_LISTED = 6
function _warn_instance_reads(model, lines::Vector{String})
    listed = join(("  " * l for l in first(lines, _INSTANCE_READS_LISTED)), "\n")
    more = length(lines) > _INSTANCE_READS_LISTED ? "\n  ... and $(length(lines) - _INSTANCE_READS_LISTED) more" : ""
    @warn "The automata of $(eltype(model)) break the contract of the kernels (`AutomatonState`): what " *
          "belongs to an instance reads its own discrete quantities and memory only, never another's. " *
          "Every instance of $(eltype(model)) is settled in every round instead of those that switched " *
          "(correct, slower). Traced from its kernel:\n$listed$more"
    return nothing
end

"""
    MasslessResiduals

The residuals of a model that are rows without mass of the DAE `M u̇ = F(u)`, which
constrain the state as algebraic equations: `names[k]` a residual, `instances[k]` the
instances at which its row has zero mass. Every algebraic residual at every instance; a
state residual at the instances whose row's mass is zero (a mass coefficient that is
zero for them, e.g. a zero inertia), where any instance's is.
"""
struct MasslessResiduals
    names::Vector{Symbol}
    instances::Vector{Vector{Int}}
end

"""
    MasslessResiduals(ws::ResidualBuffers, metadata, mass) -> MasslessResiduals

The residuals of the model of `metadata` (its Float64 workspace `ws`) without mass, from
the diagonal `mass` of the simulation's mass matrix (over `[x; y]`).
"""
function MasslessResiduals(ws::ResidualBuffers, metadata, mass::AbstractVector)
    n = ws.n_instances
    names, instances = Symbol[], Vector{Int}[]
    for r in metadata.residuals
        if is_algeb_res(r)
            push!(names, r.name)
            push!(instances, collect(1:n))
        elseif is_state_res(r)
            addrs = getfield(ws.res_addrs, r.name)
            zero_mass = [i for i in 1:n if iszero(mass[Int(addrs[i])])]
            isempty(zero_mass) && continue
            push!(names, r.name)
            push!(instances, zero_mass)
        end
    end
    return MasslessResiduals(names, instances)
end

# The rows (indices into `[x; y]`) of the residuals without mass of a model, its residual
# addresses `res_addrs`: of a state residual its address, of an algebraic one `n_diff`
# after it
function _massless_rows(massless::MasslessResiduals, metadata, res_addrs, n_diff::Int)
    rows = Int[]
    for (name, instances) in zip(massless.names, massless.instances)
        offset = is_algeb_res(_residual(metadata, name)) ? n_diff : 0
        addrs = getfield(res_addrs, name)
        # an instance not delivered through the residual has no row in it
        foreach(i -> addrs[i] == UNROUTED || push!(rows, offset + Int(addrs[i])), instances)
    end
    return rows
end

_residual(metadata, name::Symbol) = metadata.residuals[findfirst(r -> r.name === name, metadata.residuals)]

"""
    trace_algebraic_switches(ws::ResidualBuffers, model, external, metadata, formulation,
                             massless::MasslessResiduals) -> NamedTuple

For each instance of `model`, whether its discrete quantities or memory are read by a
residual without mass (`massless`: an algebraic residual, or a state residual at an
instance whose row has zero mass; of any instance), so that a transition of that instance
can change the algebraic equations: traced from the model's `calc_residuals!` in
`formulation`, with `metadata` its metadata there, on a copy of its Float64 workspace `ws`
of the tracer's element type, with SparseConnectivityTracer's global detector. The modes act on the
residuals through the discrete quantities they hold (`AutomatonState`: the kernels read
`discrete` and `memory`), and through the states they govern, which a transition into the mode
writes (`Mode.states`). Returns `(switches, inputs, setpoints)`: the vector of instances, `inputs[q, i]`
whether input `q` (a discrete quantity, a memory, or a governed state, in this order) of instance `i` is read by one, and
whether the trace holds for the setpoints at their values.

The discrete quantities and memory are the indexed inputs, one per quantity and instance, and
after them each governed state, one per instance.
The variables, the time and the setpoints enter as tracers, not indexed: what the trace
finds holds whatever the state, the time and the setpoints, and a kernel that branches on
them (`if`, `?:`, `&&`; an `ifelse` reads both of its branches) cannot be traced so. A
kernel that branches on a setpoint is traced with the setpoints at their values then, as
the values of the automata are (`trace_value_reads`), which holds until an event changes
them (`AlgebraicCheck` `retrace!`; `setpoints` is true). The parameters (the model's
columns) and the external parameters enter at their values, on which kernels branch (a
time constant that may be zero): the
pattern holds at the parameters of the trace, and the simulation traces it again after the
events that change the model's fields (`AutomatonSettle`). Every instance is traced in
service (an event may switch it on). A kernel that branches on a variable, the time, a
discrete quantity or memory cannot be traced and throws; so does any other error of the
kernel, as it is (`_branches_on_tracer`).
"""
function trace_algebraic_switches(ws::ResidualBuffers, model, external, metadata, formulation,
                                  massless::MasslessResiduals)
    n = ws.n_instances
    n_state = length(ws.automata.discrete) + length(ws.automata.memory)
    governed = _governed_variables(ws, metadata)
    n_inputs = n_state + length(governed)
    names = massless.names
    n_rows = length(names)
    (n_rows == 0 || n_inputs == 0) && return (switches = fill(false, n), inputs = falses(n_inputs, n), setpoints = false)
    # Where each residual is a row without mass (a state residual's other instances are
    # rows with mass, which the trace leaves out)
    without_mass = [(m = fill(false, n); m[instances] .= true; m) for instances in massless.instances]
    in_service = fill(true, n)
    trace = function (at_values::Bool)
        residuals! = function (out, v)
            T = eltype(v)
            traced = ResidualBuffers(ws, T; setpoints_at_values = at_values)
            _set_discrete_inputs!(traced, v, 0)
            _set_governed_inputs!(traced, governed, v, n_state * n)
            calc_residuals!(traced, model, external, formulation, zero(T); active_mask = in_service)
            for (k, name) in enumerate(names), i in 1:n
                out[(i - 1) * n_rows + k] = without_mass[k][i] ? getfield(traced.res_vals, name)[i] : zero(T)
            end
            return nothing
        end
        return jacobian_sparsity(residuals!, zeros(n_rows * n), zeros(n_inputs * n), TracerSparsityDetector())
    end
    setpoints = false
    pattern = try
        trace(false)
    catch e
        _branches_on_tracer(e) || rethrow()
        setpoints = true
        try
            trace(true)
        catch e
            _branches_on_tracer(e) || rethrow()
            throw(ErrorException(
                "The residuals of $(eltype(model)) could not be traced for the algebraic equations its automata " *
                "change (`calc_residuals!` with SparseConnectivityTracer's global detector, the parameters and " *
                "setpoints at their values, which refuses branches on a variable, the time, a discrete quantity or " *
                "memory): " * sprint(showerror, e)))
        end
    end
    # input q of instance i at (q - 1) n + i
    inputs = BitMatrix(any(view(pattern, :, (q - 1) * n + i)) for q in 1:n_inputs, i in 1:n)
    return (switches = [any(view(pattern, :, i:n:(n_inputs * n))) for i in 1:n], inputs, setpoints)
end

"""
    rows_without_their_variable(ws::ResidualBuffers, model, external, metadata, formulation,
                                pairs, x, y) -> Vector{NamedTuple}

The algebraic rows of a model with automata that do not read the variable they are paired
with, in some combination of the modes of its automata: `pairs` holds `(residual,
variable)` names (an algebraic residual and the algebraic variable an update that solves
each row through its own variable pairs it with). For each combination of one mode per
automaton, every instance's discrete quantities at what the modes hold, the residuals of
`model` (its `calc_residuals!` in `formulation`, on a copy of its Float64 workspace `ws` of
the tracer's element type, every instance in service) are traced at the variables gathered
from `x` and `y` with SparseConnectivityTracer's local detector, which drops a product
with an exact zero (a mode's discrete quantity that is zero) and the branch of a piecewise
function not taken there: the pattern each mode gives its rows. Returns, for each row and
instance whose pattern lacks its own variable, `(instance, residual, variable, modes)`, the
modes as `automaton => mode` names. The discrete quantities of `ws` are restored.
"""
function rows_without_their_variable(ws::ResidualBuffers, model, external, metadata, formulation, pairs, x, y)
    found = NamedTuple{(:instance, :residual, :variable, :modes), Tuple{Int, Symbol, Symbol, Vector{Pair{Symbol, Symbol}}}}[]
    isempty(pairs) && return found
    n = ws.n_instances
    names = collect(keys(ws.var_vals))
    gather_variables!(ws, x, y)
    point = reduce(vcat, (copy(c) for c in values(ws.var_vals)); init = Float64[])
    in_service = fill(true, n)
    rows! = function (out, v)
        traced = ResidualBuffers(ws, eltype(v))
        for (j, column) in enumerate(values(traced.var_vals))
            copyto!(column, 1, v, (j - 1) * n + 1, n)
        end
        calc_residuals!(traced, model, external, formulation, zero(eltype(v)); active_mask = in_service)
        for (p, (residual, _)) in enumerate(pairs), i in 1:n
            out[(p - 1) * n + i] = getfield(traced.res_vals, residual)[i]
        end
        return nothing
    end
    automata = metadata.automata
    discrete = ws.automata.discrete
    saved = map(copy, values(discrete))
    try
        for combo in Iterators.product((eachindex(a.modes) for a in automata)...)
            for (a, k) in zip(automata, combo), (q, v) in a.modes[k].holds
                fill!(getfield(discrete, q), v)
            end
            pattern = jacobian_sparsity(rows!, zeros(length(pairs) * n), point, TracerLocalSparsityDetector())
            modes = Pair{Symbol, Symbol}[a.name => a.modes[k].name for (a, k) in zip(automata, combo)]
            for (p, (residual, variable)) in enumerate(pairs), i in 1:n
                j = findfirst(==(variable), names)
                j === nothing && throw(ArgumentError(
                    "$(metadata.name): $variable, paired with the row $residual, is not a variable of its workspace"))
                pattern[(p - 1) * n + i, (j - 1) * n + i] || push!(found, (instance = i, residual, variable, modes))
            end
        end
    finally
        foreach(copyto!, values(discrete), saved)
    end
    return found
end

"""
    AlgebraicCheck

Whether the transitions of the instances of a model change the algebraic equations (the
rows without mass of the DAE: its algebraic residuals, and its state residuals where their
mass is zero, `MasslessResiduals`), as a settle checks it (`settle_instant!`):
`switches[i]`, whether instance `i`'s discrete quantities or memory, or a state its modes
govern (a transition into the mode writes it, `Mode.states`), are read by a residual
without mass (`trace_algebraic_switches`), or the state it governs is gathered by another
model (`Governed` `foreign`); `inputs[q, i]`, which of them: the discrete quantity or memory `q` (the order of
`AutomatonModel` `quantity_columns`) of instance `i`, or after them a governed state (`trace_algebraic_switches`). Then, for a model with any such instance, how a
settle decides whether its transitions changed the residuals without mass at the same
state:

- `alone`: by the model's own contributions to the rows without mass before and after
  (its kernel alone, `evaluate!` at `(x, y, t)` into `columns`, its residual buffers
  without mass, compared at `instances`, the instances where each is without mass;
  `before` a copy of them). They decide it where the model sets no shared residual
  without mass (`SharedRes(:set)`, which reads what the others accumulated there) and no
  model sets one in a row without mass where a residual of the model adds (which would
  replace what it adds; every row without mass of the model, so that this holds whichever
  rows a trace finds reading the discrete state): the system's rows without mass are
  then the sums of every model's contributions, and only the switched models' change;
- otherwise by the system's rows without mass before and after.

`setpoints` is whether the residual kernel branches on a setpoint, so that `switches`
holds for their values at the trace (`trace_algebraic_switches`); `retrace!()` traces
`switches` again (after an event that changes the model's fields, or its setpoints where
`setpoints`, on which its kernel may branch).
"""
struct AlgebraicCheck
    switches::Vector{Bool}
    inputs::BitMatrix
    alone::Bool
    evaluate!::Function
    columns::Vector{Vector{Float64}}
    instances::Vector{Vector{Int}}
    before::Vector{Vector{Float64}}
    setpoints::Base.RefValue{Bool}
    retrace!::Function
end

"""
    AlgebraicCheck(ws, model, external, metadata, formulation, active, massless, set_rows, n_diff)
        -> AlgebraicCheck

The `AlgebraicCheck` of `model`, its Float64 workspace `ws` and in-service mask `active`,
its residuals without mass `massless`; `set_rows` the rows without mass (indices into
`[x; y]`, `n_diff` the length of `x`) any model of the simulation sets (`SharedRes(:set)`).
"""
function AlgebraicCheck(ws::ResidualBuffers, model, external, metadata, formulation, active::Vector{Bool},
                        massless::MasslessResiduals, set_rows::BitSet, n_diff::Int;
                        foreign::Union{Nothing, Vector{Bool}} = nothing)
    traced = trace_algebraic_switches(ws, model, external, metadata, formulation, massless)
    switches, inputs, setpoints = traced.switches, copy(traced.inputs), Ref(traced.setpoints)
    foreign === nothing || (switches .|= foreign)
    columns = Vector{Float64}[getfield(ws.res_vals, name) for name in massless.names]
    # Whether the model's contributions decide alone: it sets no shared residual without
    # mass, and no model sets a row without mass its residuals add to or own, and no other
    # model reads a state its modes govern (a write of it changes that model's contributions)
    sets = any(name -> is_shared_res(:set)(_residual(metadata, name)), massless.names)
    alone = !sets && (foreign === nothing || !any(foreign)) && !any(in(set_rows), _massless_rows(massless, metadata, ws.res_addrs, n_diff))
    evaluate! = let ws = ws, model = model, external = external, metadata = metadata, formulation = formulation,
                    active = active
        function (x, y, t)
            gather_variables!(ws, x, y)
            calc_residuals!(ws, model, external, formulation, t; active_mask = active)
            return nothing
        end
    end
    retrace! = let ws = ws, model = model, external = external, metadata = metadata, formulation = formulation,
                   switches = switches, inputs = inputs, massless = massless, setpoints = setpoints, foreign = foreign
        function ()
            traced = trace_algebraic_switches(ws, model, external, metadata, formulation, massless)
            copyto!(switches, traced.switches)
            copyto!(inputs, traced.inputs)
            foreign === nothing || (switches .|= foreign)
            setpoints[] = traced.setpoints
            return nothing
        end
    end
    return AlgebraicCheck(switches, inputs, alone, evaluate!, columns, massless.instances, [similar(c) for c in columns],
                          setpoints, retrace!)
end

# Whether the model's contributions to its rows without mass differ from those copied
# before (exactly, at the instances where each residual is without mass)
function _contributions_changed(a::AlgebraicCheck)
    for (column, before, instances) in zip(a.columns, a.before, a.instances)
        for i in instances
            column[i] == before[i] || return true
        end
    end
    return false
end

# The rows without mass (indices into `[x; y]`, `n_diff` the length of `x`) set by a
# `SharedRes(:set)` residual of any model of `sys_ws`: its algebraic rows, and its state
# rows whose mass (the diagonal `mass` of the mass matrix) is zero
function _set_massless_rows(sys_ws, mass::AbstractVector, n_diff::Int)
    rows = BitSet()
    for (name, ws) in sys_ws.workspaces_by_name, r in sys_ws.metadata_by_name[name].residuals
        is_shared_res(:set)(r) || continue
        addrs = getfield(ws.res_addrs, r.name)
        if is_algeb_res(r)
            foreach(a -> push!(rows, n_diff + Int(a)), addrs)
        elseif is_state_res(r)
            foreach(a -> iszero(mass[Int(a)]) && push!(rows, Int(a)), addrs)
        end
    end
    return rows
end

"""
    CompiledTransition

A transition of an automaton with its modes as indices, its condition as declared and
compiled (`compiled`, for the automaton callback), its guard and sets as the columns
they read (`sets`, with the names of the memory and values in `declared`, the
transition's own; `guard_value` the guard's index among the values, 0 for a condition),
and `crossing` its index among the crossings of the automaton (0 for a condition
transition), and `continuity` its `Transition.continuity` as indices among the values:
`(surface, held)`, the values continuity holds at zero where it fires on the zero of
`surface` (`_held_at_zero!`); `crossing_continuity` its `Transition.crossing_continuity`
in the same form, applied only where it fires in the first round of a located instant, and `writers`
its `Transition.continuity_writers` as `(value, automata)`, the indices of the value among the values and of the
automata among those of the model: the value is not held at the instant where one of those fired (`_held_at_zero!`).
"""
struct CompiledTransition
    from::Int
    to::Int
    condition::Union{Nothing, Condition}
    compiled::Union{Nothing, CompiledCondition}
    guard::Union{Nothing, Vector{Float64}}
    guard_value::Int
    direction::Int8                     # +1 rising, -1 falling, 0 for a condition
    crossing::Int
    sets::Vector{Tuple{Vector{Float64}, Union{Float64, Vector{Float64}}}}
    continuity::Vector{Tuple{Int, Vector{Int}}}
    crossing_continuity::Vector{Tuple{Int, Vector{Int}}}
    writers::Vector{Tuple{Int, Vector{Int}}}
    declared::Transition
end

"""
    CompiledAutomaton

An automaton of a model at run time: its declaration, its mode array (one index per
instance, in the `AutomatonState`), what each mode holds, its transitions, its initial
clauses (`initial`: each condition compiled, `nothing` for an unconditional clause, with
its mode), its held clauses' conditions compiled (`held`, in declaration order) with their
modes (`held_mode`), its crossing transitions (`crossings`, in declaration order), and
`held_now`: per instance, whether a held clause holds (then no transition fires).
"""
struct CompiledAutomaton
    automaton::Automaton
    modes::Vector{Int}
    holds::Vector{Vector{Tuple{Vector{Float64}, Union{Float64, Vector{Float64}}}}}
    transitions::Vector{CompiledTransition}
    initial::Vector{Tuple{Union{Nothing, CompiledCondition}, Int}}
    held::Vector{CompiledCondition}
    held_mode::Vector{Int}
    crossings::Vector{CompiledTransition}
    held_now::Vector{Bool}
end

"""
    GovernedState

A state a mode of an automaton governs (`Mode.states`), in a model of the runtime: the indices of its
`automaton` and its `mode`, its variable `var` among those of the workspace and its name `state`,
`addresses` its position in `x` for each instance (`0` where the instance is not delivered a variable),
and `value`, its right side: a number, or the index among the model's values of the column that
computes it (`column`; `nothing` for a number). `param_reads[i]`, for each instance, the model's own
parameter fields the right side reads there (traced, `_value_param_reads`; empty for a number): an
event that writes none of them cannot move this instance's right side (`check_governed`'s admissibility
already proves a parameter is the only channel an event can use to move it), the exact test the repair
after an event uses to find which governed states it must put back on their relation.
"""
struct GovernedState
    automaton::Int
    mode::Int
    var::Int
    state::Symbol
    addresses::Vector{Int}
    value::Union{Float64, Int}
    column::Union{Nothing, Vector{Float64}}
    param_reads::Vector{Set{Symbol}}
end

"""
    Governed

The states the modes of a model govern: `states` (`GovernedState`), `by_mode[ka][km]` the indices
among them that mode `km` of automaton `ka` governs, and, for each, `foreign[k][i]` the variables of
other models that are the state of instance `i` (`(model name, variable index, instance)`; their
values and equations read it). `recheck!()` checks the right sides again for what they read
(`check_governed`), after events that change parameters on which the kernel of the values branches.
Empty (no state governed) for a model whose modes govern none.
"""
struct Governed
    states::Vector{GovernedState}
    by_mode::Vector{Vector{Vector{Int}}}
    foreign::Vector{Vector{Vector{Tuple{Symbol, Int, Int}}}}
    recheck!::Function
end

_no_recheck() = nothing
const NO_GOVERNED = Governed(GovernedState[], Vector{Vector{Int}}[], Vector{Vector{Tuple{Symbol, Int, Int}}}[], _no_recheck)
Governed() = NO_GOVERNED

"""
    AutomatonModel

The automata of one model of a simulation: its name and the keys of its instances (for
messages), its `AutomatonState`, its in-service mask, its external parameters (gathered
once at setup, as its kernels read them), what its values read of its variables (`reads`), whether its transitions change the algebraic equations
(`algebraic`, `AlgebraicCheck`), its atoms (`atom_column` the value each compares,
`atom_value` its index among the values, `atom_relation` how) and compiled automata, the
values each automaton reads (`automaton_values[ka]`: of its atoms and guards, by index),
its memory (`memory`, as declared, and `memory_initial`, the columns each starts from),
`quantity_columns`, the columns of its discrete quantities and then of its memory (the
order of `ValueReads` `state`), `signed[i, k]` the signed truth of atom `k` of instance `i` (`signed_truth`, measured as
below) at the last evaluation of the automaton callback's components, and the scratch of
a settle:

- `side[i, v]`, `anchor[i, v]`, `at_instant[i, v]`, `pinned[i, v]`, `returning[i, v]`: value `v` (index among
  the values, `value_columns[v]` its column) of instance `i` as the last located event
  that located it left it (`LocatedValue`, the header of automata.jl): `side` the
  direction it crossed zero in there (`+1` upward, `-1` downward, `0` where no event
  located it, or a change of mode of the automata that read it made that irrelevant,
  `_instant_left!`), `anchor` where it is measured from after the last instant (`_anchor`:
  where the instant left it if that is across zero from its side, else zero; `0.0`
  without a side), `at_instant` whether it was located at the instant being settled, and
  `pinned` whether a crossing of it was refused since as no motion (`_moving_crossings!`):
  the discrete state holds it still, or (`returning`) it crossed against its rate and is
  pinned until a step starts with it back on its side. The atoms over a value with a side
  are decided by the direction of its motion: at the instant and while it is pinned by
  its side, otherwise by its
  displacement `g - anchor` from where the instant left it, by its side where that is
  zero (`atom_holds`); a guard with a side is reported as its displacement, zero while it
  is pinned (`_guard_component`). A value whose discrete quantities or memory a
  transition changes moves by that change, a jump, and loses its side there
  (`_quantity_changed!`): it is decided by its value. `zero[i, v]`: the value is decided as
  the exact zero it is, whatever is computed for it (`atom_holds`: `>=` and `<=` hold, `>`
  and `<` do not; a guard over it is zero): from a start from equilibrium where the
  initialization holds it at zero (`at_rest`, the values of `Automaton.at_rest`), and where
  a transition holds it at zero by continuity and its rate there is zero too
  (`_held_at_zero!`), until an instant changes the equations (`_release_zeros!`);
  `n_zero[v]` counts the instances that are (`_zero!`).
  `measured[i, v]` is whether the atoms over it need
  more than its value (`_measure!`: it is pinned or zero, its anchor is off zero, or the instant
  left it exactly at zero, where its side decides; a value the instant left strictly on its
  side of zero holds there as its value does), and `n_measured[v]` how many instances do,
  so that the atoms over a value with none are evaluated by their value alone;
- `crossed[c, i]`: the direction in which crossing `c` (over all automata) of instance
  `i` was located, `0` if not;
- `last_switch[i]`, `rounds[i]`: the time instance `i` last switched, and the rounds in
  which it switched at that instant (see `_count_rounds!`).
- `governed`: the states its modes govern (`Governed`).
"""
struct AutomatonModel
    name::Symbol
    keys::Vector{Any}
    state::AutomatonState
    active::Vector{Bool}
    external::NamedTuple
    reads::ValueReads
    algebraic::AlgebraicCheck
    atoms::Vector{Atom}
    atom_index::Dict{Atom, Int}
    atom_column::Vector{Vector{Float64}}
    atom_value::Vector{Int}
    atom_relation::Vector{Int8}
    signed::Matrix{Float64}
    memory::Vector{MemoryQuantity}
    memory_initial::Vector{Tuple{Vector{Float64}, Union{Float64, Vector{Float64}}}}
    quantity_columns::Vector{Vector{Float64}}
    automata::Vector{CompiledAutomaton}
    automaton_values::Vector{Vector{Int}}
    value_columns::Vector{Vector{Float64}}
    side::Matrix{Int8}
    measured::Matrix{Bool}
    n_measured::Vector{Int}
    anchor::Matrix{Float64}
    at_instant::Matrix{Bool}
    pinned::Matrix{Bool}
    returning::Matrix{Bool}
    zero::Matrix{Bool}
    n_zero::Vector{Int}
    at_rest::Vector{Int}
    crossed::Matrix{Int8}
    last_switch::Vector{Float64}
    rounds::Vector{Int}
    governed::Governed
end

n_instances(m::AutomatonModel) = length(m.active)

function _compile_automaton(a::Automaton, state::AutomatonState, atom_index::Dict{Atom, Int},
                            crossing_offset::Int, n::Int, automaton_index::Dict{Symbol, Int})
    values = state.values
    mode_index = Dict(mode.name => k for (k, mode) in enumerate(a.modes))
    holds = [[(getfield(state.discrete, q), _resolve(v, values)) for (q, v) in mode.holds] for mode in a.modes]
    transitions = CompiledTransition[]
    n_crossings = 0
    value_index(name) = findfirst(==(name), keys(values))
    for tr in a.transitions
        crossing = is_condition(tr) ? 0 : (n_crossings += 1)
        push!(transitions, CompiledTransition(
            mode_index[tr.from], mode_index[tr.to], tr.condition,
            is_condition(tr) ? CompiledCondition(tr.condition, atom_index) : nothing,
            is_condition(tr) ? nothing : getfield(values, tr.guard),
            is_condition(tr) ? 0 : findfirst(==(tr.guard), keys(values)),
            is_condition(tr) ? Int8(0) : tr.direction === :rising ? Int8(1) : Int8(-1),
            crossing == 0 ? 0 : crossing_offset + crossing,
            [(getfield(state.memory, m), _resolve(v, values)) for (m, v) in tr.sets],
            [(value_index(surface), value_index.(held)) for (surface, held) in tr.continuity],
            [(value_index(surface), value_index.(held)) for (surface, held) in tr.crossing_continuity],
            [(value_index(v), [automaton_index[w] for w in writers]) for (v, writers) in tr.continuity_writers], tr))
    end
    initial = Tuple{Union{Nothing, CompiledCondition}, Int}[
        (c.condition === nothing ? nothing : CompiledCondition(c.condition, atom_index), mode_index[c.mode])
        for c in a.initial]
    held = [CompiledCondition(h.condition, atom_index) for h in a.held]
    crossings = [tr for tr in transitions if tr.crossing != 0]
    return CompiledAutomaton(a, getfield(state.modes, a.name), holds, transitions, initial, held,
                             [mode_index[h.mode] for h in a.held], crossings, fill(false, n))
end

# The conditions of an automaton: its transitions' (none for a crossing), its initial modes' and its
# held clauses' (nothing for the initial mode that holds otherwise), one iterator over the three
# (`vcat` of them dispatches at run time on their abstract element types)
_conditions(a::Automaton) = Iterators.flatten(((tr.condition for tr in a.transitions),
                                               (c.condition for c in a.initial), (h.condition for h in a.held)))

"""
    AutomatonModel(name, state, metadata, active, keys, external, reads::ValueReads, algebraic::AlgebraicCheck,
                   governed::Governed = Governed())

The run-time automata of model `name` of `metadata`, whose Float64 workspace holds its
`AutomatonState` `state`; `active` is its in-service mask (the simulation's status),
`external` its external parameters, `reads` what its values read (`trace_value_reads`), `algebraic` whether its
transitions change the algebraic equations and `governed` the states its modes govern.
"""
function AutomatonModel(name::Symbol, state::AutomatonState, metadata, active::Vector{Bool}, keys::AbstractVector,
                        external::NamedTuple, reads::ValueReads, algebraic::AlgebraicCheck,
                        governed::Governed = Governed())
    all_atoms = Atom[]
    for a in metadata.automata
        for c in _conditions(a), at in atoms(c)
            at in all_atoms || push!(all_atoms, at)
        end
    end
    atom_index = Dict(at => k for (k, at) in enumerate(all_atoms))
    n = length(active)
    automata = CompiledAutomaton[]
    offset = 0
    automaton_index = Dict(a.name => k for (k, a) in enumerate(metadata.automata))
    for a in metadata.automata
        compiled = _compile_automaton(a, state, atom_index, offset, n, automaton_index)
        offset += length(compiled.crossings)
        push!(automata, compiled)
    end
    memory_initial = [(getfield(state.memory, m.name), _resolve(m.initial, state.values)) for m in metadata.memory]
    value_names = Base.keys(state.values)
    atom_value = [findfirst(==(at.value), value_names) for at in all_atoms]
    # The values each automaton reads: those of the atoms of its conditions, and its guards
    automaton_values = Vector{Int}[]
    for a in metadata.automata
        read = Int[atom_value[atom_index[at]] for c in _conditions(a) for at in atoms(c)]
        append!(read, (findfirst(==(tr.guard), value_names) for tr in a.transitions if !is_condition(tr)))
        push!(automaton_values, sort!(unique!(read)))
    end
    n_values = length(value_names)
    at_rest = sort!(unique!(Int[findfirst(==(v), value_names) for a in metadata.automata for v in a.at_rest]))
    return AutomatonModel(name, collect(Any, keys), state, active, external, reads, algebraic, all_atoms, atom_index,
                          [getfield(state.values, at.value) for at in all_atoms], atom_value,
                          Int8[_relation_code(at.relation) for at in all_atoms], zeros(n, length(all_atoms)),
                          collect(MemoryQuantity, metadata.memory), memory_initial,
                          Vector{Float64}[values(state.discrete)..., values(state.memory)...], automata, automaton_values,
                          collect(Vector{Float64}, values(state.values)), zeros(Int8, n, n_values), fill(false, n, n_values), zeros(Int, n_values),
                          zeros(Float64, n, n_values), fill(false, n, n_values), fill(false, n, n_values),
                          fill(false, n, n_values), fill(false, n, n_values), zeros(Int, n_values), at_rest,
                          zeros(Int8, offset, n), fill(NaN, n), zeros(Int, n), governed)
end

# ----------------------------------------------------------------------------
# Conditions
# ----------------------------------------------------------------------------

"""
    atom_holds(m::AutomatonModel, k, i) -> Bool

Whether atom `k` of instance `i` holds. Over a value decided as zero (`AutomatonModel`
`zero`), as the exact zero it is: `>=` and `<=` hold, `>` and `<` do not. Over a value that
a located event left with a side (`AutomatonModel` `side`), by the direction of its motion: at the instant that
located it, as the crossing leaves it (`>` and `>=` true after an upward crossing, `<` and
`<=` false, the reverse after a downward one), whatever the settles and
reinitializations there make of the value, and so while it is `pinned` (a crossing of it
was refused as no motion since); otherwise by its displacement from where the last instant left it
(`anchor`), and as the crossing left it where that is zero. Over any other value, by the
value.
"""
@inline function atom_holds(m::AutomatonModel, k::Int, i::Int)
    v, rel = m.atom_value[k], m.atom_relation[k]
    m.zero[i, v] && return _holds_at_zero(rel)
    side = m.side[i, v]
    (m.at_instant[i, v] || m.pinned[i, v]) && return (side > 0) == _is_upper(rel)
    s = m.atom_column[k][i] - m.anchor[i, v]
    (iszero(s) && side != 0) && return (side > 0) == _is_upper(rel)
    return compare(s, rel)
end

# Whether an atom of relation `rel` holds above zero (`>`, `>=`) rather than below
@inline _is_upper(rel::Int8) = rel === _GT || rel === _GE

# Whether an atom of relation `rel` holds at an exact zero (`>=`, `<=`)
@inline _holds_at_zero(rel::Int8) = rel === _GE || rel === _LE

# Set the side of value `v` of instance `i` of `m` (`AutomatonModel` `side`)
@inline function _set_side!(m::AutomatonModel, i::Int, v::Int, side::Int8)
    m.side[i, v] = side
    return nothing
end

# Whether value `v` of instance `i` of `m` is decided as zero (`AutomatonModel` `zero`),
# counting the instances that are (`n_zero`)
@inline function _zero!(m::AutomatonModel, i::Int, v::Int, zero::Bool)
    m.n_zero[v] += zero - m.zero[i, v]
    m.zero[i, v] = zero
    return nothing
end

# Whether the atoms over value `v` of instance `i` of `m` need more than its value
# (`AutomatonModel` `measured`), counting the instances that do (`n_measured`)
@inline function _measure!(m::AutomatonModel, i::Int, v::Int, needed::Bool)
    m.n_measured[v] += needed - m.measured[i, v]
    m.measured[i, v] = needed
    return nothing
end

"""
    anchored_signed_truth(g, anchor, side, pinned, rel)

The signed truth (`signed_truth`) of an atom of relation `rel` over a value `g` that a
located event left at `anchor` with `side` (see `atom_holds`): that of its displacement
`g - anchor`, and where that is zero or the value is `pinned` the smallest subnormal on
the side the crossing left it. Without a side (`side = 0`, `anchor = 0`) it is
`signed_truth(g, rel)`.
"""
@inline function anchored_signed_truth(g::Float64, anchor::Float64, side::Int8, pinned::Bool, rel::Int8)
    s = pinned ? 0.0 : g - anchor
    (iszero(s) && side != 0) || return signed_truth(s, rel)
    tiny = nextfloat(0.0)
    return (side > 0) == _is_upper(rel) ? tiny : -tiny
end

"""
    _guard_component(m::AutomatonModel, guard, v, i)

The component of the automaton callback for a crossing over `guard` (value `v` of `m`) of
instance `i`: the guard, or, where a located event left it with a side, its displacement
from where the instant left it (`anchor`), zero there and while the discrete state holds
it pinned (`pinned`) or it is decided as zero (`zero`), so that, as at the start, a guard starting at its level and leaving it
is not a crossing.
"""
@inline function _guard_component(m::AutomatonModel, guard::Vector{Float64}, v::Int, i::Int)
    return @inbounds ifelse(m.pinned[i, v] | m.zero[i, v], 0.0, guard[i] - m.anchor[i, v])
end

# Whether condition `c` of instance `i` of `m` holds, its atoms by `atom_holds`
holds(c::CompiledCondition, m::AutomatonModel, i::Int) = condition_holds(c, k -> atom_holds(m, k, i))

# ----------------------------------------------------------------------------
# The runtime of a simulation
# ----------------------------------------------------------------------------

"""
    AutomatonComponents

The components of the automaton callback (automata.jl), in blocks of one component per
instance of a model: for each block, the model and automaton (indices into the runtime),
the transition (a condition or a crossing; `nothing` for a held clause), the condition it
reports (`nothing` for a crossing: its guard), and the offset of the block; `n` in all.

`watch[b][i]` is whether component `i` of block `b` can fire something in the current
modes and held clauses (a held clause always; a transition from the mode of an instance
whose automaton no held clause holds), `watched[b]` how many components of block `b` can,
and `conditions` and `crossings` the blocks with any such component; `blocks[km]` lists
the blocks of model `km`, and `live[km]` the atoms of model `km` that the conditions of those
blocks read: the only ones whose signed truths the callback computes and records. The modes and held clauses change in `settle!` only (and
`start!` through it), which sets them after every settle for the instances it considered
(`_watch!`); the callback reports those components where their instance is in service,
every other component zero. `scratch` holds the stack of `signed_block!` (one column per
level, one row per instance).
"""
struct AutomatonComponents
    model::Vector{Int}
    automaton::Vector{Int}
    transition::Vector{Union{Nothing, CompiledTransition}}
    condition::Vector{Union{Nothing, CompiledCondition}}
    offset::Vector{Int}
    n::Int
    watch::Vector{Vector{Bool}}
    watched::Vector{Int}
    blocks::Vector{Vector{Int}}
    conditions::Vector{Int}
    crossings::Vector{Int}
    live::Vector{Vector{Int}}
    scratch::Matrix{Float64}
end

function AutomatonComponents(models::Vector{AutomatonModel})
    model, automaton, offset = Int[], Int[], Int[]
    transition = Union{Nothing, CompiledTransition}[]
    condition = Union{Nothing, CompiledCondition}[]
    n = 0
    add!(km, ka, tr, c) = (push!(model, km); push!(automaton, ka); push!(transition, tr);
                           push!(condition, c); push!(offset, n); n += n_instances(models[km]))
    for (km, m) in enumerate(models), (ka, a) in enumerate(m.automata)
        foreach(tr -> add!(km, ka, tr, tr.compiled), a.transitions)
        foreach(c -> add!(km, ka, nothing, c), a.held)
    end
    depth = maximum((length(c.stack) for c in condition if c !== nothing); init = 0)
    return AutomatonComponents(model, automaton, transition, condition, offset, n,
                               [fill(false, n_instances(models[km])) for km in model], zeros(Int, length(model)),
                               [findall(==(km), model) for km in eachindex(models)], Int[], Int[],
                               [Int[] for _ in models], zeros(maximum(n_instances, models), depth))
end

"""
    ExtremumCounts

What the search inside each step of the automaton callback for conditions that change between the
points it is evaluated at did in a simulation (`StepExtrema`): `steps` the steps it scanned,
`candidate_steps` those with a component whose sampled values turned though its sign did not (the
prefilter), `candidates` those components (each turn examined: the rates of the values of its atoms and
the extremum), `located` those whose condition changes its truth inside the window (reported folded: the
callback locates its first change), and `nanoseconds` the time of the scans (the sampled values are kept at each evaluation of the
components, which costs a copy of them).
"""
mutable struct ExtremumCounts
    steps::Int
    candidate_steps::Int
    candidates::Int
    located::Int
    nanoseconds::Int
    dual_windows::Int
end

ExtremumCounts() = ExtremumCounts(0, 0, 0, 0, 0, 0)

"""
    AutomatonRuntime

The automata of the models of a simulation (`AutomatonModel`s), the lengths `n_diff` of
the differential state and `n` of the state `[x; y]` (an integration may carry more
after it), `evaluators`, one function `(x, y, t)` per model computing its values
(`evaluate_automaton_values!` on its Float64 workspace, gathering the variables the
values read alone, `ValueReads`), `bounds`, one function `(t)` per model computing the running error bounds
of its values at the state of its last evaluation (`RoundingBound`, read where a start decides which values are on a surface, `_on_surfaces`), the
components of the automaton callback (`AutomatonComponents`), and `sync_memory!`, which copies the memory
the runtime sets from the Float64 workspace into the workspaces of other element types
(`sync_automaton_memory!` on the simulation's AD cache, or the Y-bus DAE's).

The scratch of a settle: `switched`, the instances (`(model, instance)` pairs) whose
discrete state the last `settle!` changed, `held`, the values (`(model, instance, value)`)
the transitions of the instant held at zero by continuity and gave a side by their rate
(`_held_at_zero!`, left where the instant leaves them by `_instant_left!`),
`switched_automata`, the automata
(`(model, instance, automaton)`) whose mode the settles since the start of the instant
changed (read where it ends, `_instant_left!`), `marked`, the models whose values are to be
computed (`evaluate_values!` of some instances); and of its algebraic check
(`settle_instant!`): `checked`, the models whose contributions to the rows without mass it
copied before the transitions (`AlgebraicCheck`), and `before` and `after`, the system's
residuals where it evaluates them, compared on the algebraic rows and on `massless`, the
rows of `x` whose mass is zero. `system_evaluations` and `model_evaluations` count
the evaluations of the system's residuals and of a model's by the algebraic checks, and
`located_events` and `still_crossings` the events of the automaton callback and the
crossings among them that were no motion (`_moving_crossings!`; read by tests), and
`returning` the values pinned until a step starts with them back on their side
(`_returned!`). `mass` is the diagonal of the mass matrix, and `rates` the rates of the
plain simulation's DAE (`StateRates`, built at the first located transition that needs
them, `_rate_dae`; `nothing` before), and `rate_completions` and `rate_factorizations` count
the rates of the whole state they evaluated (the rows without mass completed) and the
factorizations of `F_ZZ` among them (read by tests).

`evaluators` is a `Tuple`, concretely typed where a function holds it: the automaton
callback calls it at every evaluation of its components, without dispatch or
allocation; elsewhere (`evaluate_values!`, at events) the field's type is a barrier.

The buffers of the variables the values do not read (`ws.var_vals`) keep what the last
evaluation of the residuals gathered there, at another state: only the model's kernels
read those buffers, and every evaluation of them gathers what they read first (the
residuals all variables, the values those of `ValueReads`).

`max_rounds` is the budget of transitions one instance takes at one instant (`AUTOMATON_MAX_ROUNDS` unless
set: a transition still needed after the budget is spent throws, `_count_rounds!`). `governed` is whether any
model governs a state (`Mode.states`): a runtime without governed states does no work for them. Where a
transition entered a mode that governs a state and wrote it, an algebraic equation that reads it changed:
`settle!` ends its round there (`settle_partial`, set while transitions remain to be considered), so that the
settle of the instant reinitializes the network before the next round; `projections` counts the writes.

`touched` is the exact `(model, field, instance)` triples a fixed-time event's `FieldModification`s have written
since it was last cleared (`events.jl`'s `affect!`, which pushes to it right before the DAE reinitializes and the
automata settle; `_repair_writes`/`governed_states.jl` reads it and clears it at the end of the first round of the
instant). Simultaneous events at one time each run their own apply-then-settle in sequence (DiffEq applies one
`FieldModification`'s callback, including its reinitialization, before the next), so `touched` at any one settle
holds only that callback's own write -- but the model's fields it reads are already cumulative by then (every
earlier callback at this same time already applied), so a relation that reads fields from two simultaneous
events is still repaired correctly, at the last of them to touch it.
"""
# Work counters of the rounds (read by tests and benchmarks): `rounds_taken`, the rounds in which any instance took a step;
# `multi_round`, `(t, model, instance)` for each instance that took part in a second round at one instant (the instants
# that do not settle in one round); `reinitializations`, the solves of the algebraic point between rounds and those
# of a start (`_settle_rounds!`).
struct AutomatonRuntime
    models::Vector{AutomatonModel}
    n_diff::Int
    n::Int
    evaluators::Tuple
    bounds::Tuple
    components::AutomatonComponents
    sync_memory!::Function
    switched::Vector{Tuple{Int, Int}}
    held::Vector{Tuple{Int, Int, Int}}
    switched_automata::Vector{Tuple{Int, Int, Int}}
    marked::Vector{Bool}
    checked::Vector{Int}
    before::Vector{Float64}
    after::Vector{Float64}
    massless::Vector{Int}
    system_evaluations::Base.RefValue{Int}
    model_evaluations::Base.RefValue{Int}
    located_events::Base.RefValue{Int}
    still_crossings::Base.RefValue{Int}
    returning::Base.RefValue{Int}
    mass::Vector{Float64}
    rates::Base.RefValue{Any}
    rate_completions::Base.RefValue{Int}
    rate_factorizations::Base.RefValue{Int}
    extrema::ExtremumCounts
    max_rounds::Base.RefValue{Int}
    governed::Bool
    settle_partial::Base.RefValue{Bool}
    settle_equations::Base.RefValue{Bool}
    projections::Base.RefValue{Int}
    rounds_taken::Base.RefValue{Int}
    multi_round::Vector{Tuple{Float64, Symbol, Int}}
    reinitializations::Base.RefValue{Int}
    touched::Set{Tuple{Symbol, Symbol, Int}}
end

"""
    automaton_runtime(sys_ws, sys, status, n_diff, n_alg, sync_memory!, mass; traces=StructuralTraces())
        -> Union{AutomatonRuntime, Nothing}

The automata of the models of `sys` evaluated with the Float64 workspace `sys_ws` (its
`automaton_models`, `workspaces_by_name`, `metadata_by_name`, `models_by_name` and
`formulation`: a `SystemBuffers`, or the Y-bus DAE's), `status` the simulation's
in-service masks, `n_diff` and `n_alg` the lengths of `x` and `y`, `sync_memory!()`
copying the memory of the automata from `sys_ws` into the workspaces of other element
types, `mass` the diagonal of the mass matrix; or `nothing` when no model declares any. What the values of each
model read is traced here (`value_reads`, kept in `traces`, the problem's, for its next solve), and its
evaluator gathers that alone;
and which instances' transitions can change the algebraic equations (`AlgebraicCheck`).
"""
function automaton_runtime(sys_ws, sys, status::AbstractDict, n_diff::Int, n_alg::Int, sync_memory!::Function,
                           mass::AbstractVector; traces::StructuralTraces = StructuralTraces())
    isempty(sys_ws.automaton_models) && return nothing
    models = AutomatonModel[]
    evaluators = Any[]
    bounds = Any[]
    set_rows = _set_massless_rows(sys_ws, mass, n_diff)
    readers = _state_readers(sys_ws)   # nothing where no model governs a state
    for name in sys_ws.automaton_models
        ws = sys_ws.workspaces_by_name[name]
        meta = sys_ws.metadata_by_name[name]
        model = sys_ws.models_by_name[name]
        external = getfield(sys.external, name)
        reads = mark_massless!(value_reads(traces, name, ws, model, external, sys_ws.formulation), ws, mass, n_diff)
        governed = Governed(name, ws, meta, model, external, sys_ws.formulation, reads, mass, n_diff, readers)
        if !isempty(governed.states)
            # Widen the gathered variables with whatever a governing mode's relation reads, which a single
            # trace of `reads` may have missed (`_governed_read_variables`, correction for the per-mode
            # retrace this round's design rejected as unsafe on the live, shared mode arrays)
            var_names = collect(keys(ws.var_vals))
            for j in _governed_read_variables(ws, reads, governed.states)
                nm = var_names[j]
                nm in reads.variables || push!(reads.variables, nm)
            end
        end
        algebraic = AlgebraicCheck(ws, model, external, meta, sys_ws.formulation, status[name],
                                   MasslessResiduals(ws, meta, mass), set_rows, n_diff;
                                   foreign = _read_by_others(governed, ws.n_instances))
        push!(models, AutomatonModel(name, ws.automata, meta, status[name], get_key(model), external, reads, algebraic,
                                     governed))
        gather = variable_gather(ws, meta, reads.variables)
        push!(evaluators, _values_evaluator(ws, gather, model, external, sys_ws.formulation, meta))
        push!(bounds, _bounds_evaluator(ws, gather, model, external, sys_ws.formulation, meta))
    end
    _wire_foreign_governed!(sys_ws, sys, models)
    n = n_diff + n_alg
    return AutomatonRuntime(models, n_diff, n, Tuple(evaluators), Tuple(bounds), AutomatonComponents(models), sync_memory!,
                            Tuple{Int, Int}[], Tuple{Int, Int, Int}[], Tuple{Int, Int, Int}[], fill(false, length(models)),
                            Int[], zeros(n),
                            zeros(n), [r for r in 1:n_diff if iszero(mass[r])], Ref(0), Ref(0), Ref(0), Ref(0), Ref(0),
                            collect(Float64, mass), Ref{Any}(nothing), Ref(0), Ref(0), ExtremumCounts(),
                            Ref(AUTOMATON_MAX_ROUNDS), any(m -> !isempty(m.governed.states), models), Ref(false), Ref(false), Ref(0),
                            Ref(0), Tuple{Float64, Symbol, Int}[], Ref(0), Set{Tuple{Symbol, Symbol, Int}}())
end

"""
    _wire_foreign_governed!(sys_ws, sys, models::Vector{AutomatonModel})

Where a model gathers *another* model's governed state as a foreign reference (`GovernedState.foreign`,
found from every workspace's `var_gather_tuple` by address, `_state_readers`): point the reading model's
box in `sys_ws.foreign_refresh` (read by `_gather_and_calc!`, PowerfulCore) at a closure that, before the
reading model's own residual kernel runs, refreshes the owning model's relation at the current `(x, y, t)`
(`gather_variables!` + `calc_automaton_values!` on the owner's own workspace, then its relation's value,
`_right_side`) and overwrites the reading model's own `var_vals` column with it -- the same correction
[`enforce_governed!`](@ref) makes for a model's own governed state, carried to a model that reads someone
else's. A reading model's own `gather_variables!` call (inside its own `_gather_and_calc!`) runs *before*
this hook and would otherwise leave the raw, possibly drifted address value in that column; this hook runs
right after it and overwrites it, so ordering between the owner's and the reader's own evaluator in
`sys_ws.evaluators` does not matter: the owner is always freshly recomputed here, not merely read.

No model in this codebase does this today (`Governed.foreign` is empty everywhere); where it is
empty, nothing is set and `sys_ws.foreign_refresh[name][]` stays `nothing`, its default, at zero cost.

Round 19a, Phase D: `sys_ws` without a `foreign_refresh` field (`YBusSystemBuffers`, which has no such
mechanism) is refused loudly here when the detection above finds a cross-model reader, rather than
silently returning as it did through Phase C (which left this as a known, documented gap rather than a
silent one): a Y-bus reader would otherwise see the raw, possibly stale address instead of the owner's
relation, exactly the bug [`enforce_governed!`](@ref) exists to fix for a model's own governed state.
Implementing the mechanism for Y-bus is out of scope here (no installed model exercises it); the refusal
names every owner/state/reader triple found, so a future model that does can be dealt with on its own
schedule instead of silently computing a wrong answer.
"""
function _wire_foreign_governed!(sys_ws, sys, models::Vector{AutomatonModel})
    by_reader = Dict{Symbol, Vector{NTuple{5, Int}}}()  # reader name => (owner km, state k, owner i, foreign j, foreign i')
    for (km, m) in enumerate(models), (k, _) in enumerate(m.governed.states), i in eachindex(m.governed.foreign[k])
        for (reader_name, j, i′) in m.governed.foreign[k][i]
            push!(get!(() -> NTuple{5, Int}[], by_reader, reader_name), (km, k, i, j, i′))
        end
    end
    isempty(by_reader) && return nothing
    if !hasproperty(sys_ws, :foreign_refresh)
        triples = sort!(unique(
            "$(models[km].name).$(models[km].governed.states[k].state) (read by $reader_name)"
            for (reader_name, entries) in by_reader for (km, k, i, j, i′) in entries))
        throw(ArgumentError(
            "$(nameof(typeof(sys_ws))) has no foreign_refresh mechanism for a governed state read across " *
            "models, found here: $(join(triples, "; ")). A Y-bus reader would see the raw, possibly stale " *
            "address instead of the owner's relation; refused rather than silently proceeding with a wrong " *
            "value (round 19a, Phase D). Cross-model governed-state reads are not implemented for Y-bus."))
    end
    for (reader_name, entries) in by_reader
        ref = get(sys_ws.foreign_refresh, reader_name, nothing)
        ref === nothing && continue   # the reader has no workspace of this kind (should not happen)
        ws_reader = sys_ws.workspaces_by_name[reader_name]
        reader_columns = values(ws_reader.var_vals)
        ref[] = let models = models, sys_ws = sys_ws, sys = sys, entries = entries, reader_columns = reader_columns
            function (x, y, t)
                refreshed = Int[]
                for (km, k, i, j, i′) in entries
                    m = models[km]
                    if km ∉ refreshed
                        ws_owner = sys_ws.workspaces_by_name[m.name]
                        model_owner = sys_ws.models_by_name[m.name]
                        external_owner = getfield(sys.external, m.name)
                        gather_variables!(ws_owner, x, y)
                        calc_automaton_values!(ws_owner, model_owner, external_owner, sys_ws.formulation, t)
                        push!(refreshed, km)
                    end
                    g = m.governed.states[k]
                    # Only where the owner's instance is *currently* in the mode this state's relation belongs to
                    # (`g.mode`): the owner may be in a different mode of the same automaton (e.g. REGCA1's `:free`,
                    # governing nothing) at this instant, in which case the reader must see the owner's own, freely
                    # evolving address value (whatever its own `gather_variables!` just read), not this relation.
                    # `CompiledAutomaton.modes[i]` is instance `i`'s current mode (the same array `ws.automata.modes`
                    # shares; `m.automata[g.automaton]` is this state's automaton).
                    m.automata[g.automaton].modes[i] == g.mode || continue
                    reader_columns[j][i′] = _right_side(g, i)
                end
                return nothing
            end
        end
    end
    return nothing
end

# The values of the automata of a model at (x, y, t), on its Float64 workspace `ws`,
# gathering the variables of `gather` (those the values read); `metadata` is used only for
# `enforce_governed!` where the model governs a state (`evaluate_automaton_values!`)
_values_evaluator(ws, gather::Tuple, model, external, formulation, metadata) =
    let ws = ws, gather = gather, model = model, external = external, formulation = formulation, metadata = metadata
        (x, y, t) -> evaluate_automaton_values!(ws, gather, model, external, x, y, formulation, t, metadata)
    end

# The running error bounds of the values of the automata of a model at the time `t` (`RoundingBound`), on buffers of
# that element type of its own (`value_buffers`), at the variables of `gather` as the last evaluation of its values
# gathered them into the Float64 workspace `ws` (`evaluate_values!`), copied there, each its own rounding; its kernel of
# the values evaluated on them, reading the parameters, setpoints, discrete quantities and memory of `ws`. Returns the
# values' bounds (their `values`). Where the model governs a state, the bound is also substituted with its
# relation's bound (`enforce_governed!`, run on `bounded`, the `RoundingBound` workspace: the governed state's
# own uncertainty is not meaningful once it is read as the relation, whose bound is already a function of what
# the relation reads), with the same two-pass shape `evaluate_automaton_values!` uses.
function _bounds_evaluator(ws, gather::Tuple, model, external, formulation, metadata)
    bounded = value_buffers(ws, RoundingBound)
    pairs = map(entry -> (entry[2], _beside(entry[2], ws.var_vals, bounded.var_vals)), gather)
    return let bounded = bounded, pairs = pairs, model = model, external = external, formulation = formulation,
               metadata = metadata
        function (t)
            for (from, to) in pairs
                to .= RoundingBound.(from)
            end
            calc_automaton_values!(bounded, model, external, formulation, RoundingBound(t))
            if metadata.governs_any_state
                enforce_governed!(bounded, metadata, RoundingBound(t))
                calc_automaton_values!(bounded, model, external, formulation, RoundingBound(t))
            end
            return bounded.automata.values
        end
    end
end

"""
    evaluate_values!(rt::AutomatonRuntime, x, y, t)
    evaluate_values!(rt::AutomatonRuntime, instances, x, y, t)

Compute the values of the automata at `(x, y, t)` of every model, or of the models of
`instances` (`(model, instance)` pairs; every model where `nothing`): a model's kernel
computes all its instances at once.
"""
evaluate_values!(rt::AutomatonRuntime, x, y, t) = _evaluate_values!(rt.evaluators, x, y, t)

function evaluate_values!(rt::AutomatonRuntime, instances, x, y, t)
    instances === nothing && return evaluate_values!(rt, x, y, t)
    fill!(rt.marked, false)
    for (km, _) in instances
        rt.marked[km] = true
    end
    return _evaluate_marked!(rt, x, y, t)
end

# The values of the models marked in `rt.marked` (the evaluator of each called through
# the tuple's type, at events only)
function _evaluate_marked!(rt::AutomatonRuntime, x, y, t)
    for km in eachindex(rt.marked)
        rt.marked[km] && rt.evaluators[km](x, y, t)
    end
    return nothing
end

# Each evaluator of a tuple of them, unrolled
@inline _evaluate_values!(::Tuple{}, x, y, t) = nothing
@inline function _evaluate_values!(evaluators::Tuple, x, y, t)
    first(evaluators)(x, y, t)
    return _evaluate_values!(Base.tail(evaluators), x, y, t)
end

# The signed truths of the atoms `live` of `m` (`m.signed`) from its values, measured from where
# located events left them (`anchored_signed_truth`), those over a value decided as zero
# as at the exact zero (`signed_truth(0.0, rel)`, which never changes sign). The others are
# not read: the conditions that read them cannot fire (`AutomatonComponents` `live`).
function _signed_truths!(m::AutomatonModel, live::Vector{Int})
    signed = m.signed
    @inbounds for k in live
        column, relation, v = m.atom_column[k], m.atom_relation[k], m.atom_value[k]
        if m.n_measured[v] == 0
            for i in axes(signed, 1)
                signed[i, k] = signed_truth(column[i], relation)
            end
        elseif m.n_zero[v] == size(signed, 1)
            # (every instance at the exact zero, as at rest: one signed truth for all)
            truth = signed_truth(0.0, relation)
            for i in axes(signed, 1)
                signed[i, k] = truth
            end
        else
            _measured_truths!(m, k, column, relation, v)
        end
    end
    return nothing
end

# The signed truths of atom `k` (over value `v`, `column`) of the instances of `m` where
# some instance's value needs more than its value (`measured`) and not every one is decided
# as zero: at the exact zero where it is, measured from its anchor where it has a side, by
# its value elsewhere
@noinline function _measured_truths!(m::AutomatonModel, k::Int, column::Vector{Float64}, relation::Int8, v::Int)
    signed, measured, zero, side, anchor, pinned = m.signed, m.measured, m.zero, m.side, m.anchor, m.pinned
    @inbounds for i in axes(signed, 1)
        signed[i, k] = !measured[i, v] ? signed_truth(column[i], relation) :
            zero[i, v] ? signed_truth(0.0, relation) :
            anchored_signed_truth(column[i], anchor[i, v], side[i, v], pinned[i, v], relation)
    end
    return nothing
end

"""
    _watch!(rt::AutomatonRuntime)
    _watch!(rt::AutomatonRuntime, instances)

Mark the components of the automaton callback that can fire something in the current
modes and held clauses, of every instance or of `instances` (`(model, instance)` pairs,
the others' unchanged), and list the blocks with any (`AutomatonComponents` `watch`,
`watched`, `conditions` and `crossings`).
"""
function _watch!(rt::AutomatonRuntime)
    comps = rt.components
    for b in eachindex(comps.offset)
        watch = comps.watch[b]
        for i in eachindex(watch)
            watch[i] = _can_fire(rt, b, i)
        end
        comps.watched[b] = count(watch)
    end
    return _list_watched!(comps)
end

function _watch!(rt::AutomatonRuntime, instances)
    comps = rt.components
    for (km, i) in instances, b in comps.blocks[km]
        watch, can = comps.watch[b], _can_fire(rt, b, i)
        watch[i] == can && continue
        watch[i] = can
        comps.watched[b] += can ? 1 : -1
    end
    return _list_watched!(comps)
end

# Whether component `i` of block `b` can fire something: a held clause always; a
# transition from the mode, unless a held clause holds
function _can_fire(rt::AutomatonRuntime, b::Int, i::Int)
    comps = rt.components
    tr = comps.transition[b]
    tr === nothing && return true
    a = rt.models[comps.model[b]].automata[comps.automaton[b]]
    return !a.held_now[i] && a.modes[i] == tr.from
end

# The blocks with any component that can fire, by kind, and the atoms their conditions read
function _list_watched!(comps::AutomatonComponents)
    empty!(comps.conditions)
    empty!(comps.crossings)
    foreach(empty!, comps.live)
    for b in eachindex(comps.offset)
        comps.watched[b] > 0 || continue
        c = comps.condition[b]
        push!(c === nothing ? comps.crossings : comps.conditions, b)
        c === nothing && continue
        live = comps.live[comps.model[b]]
        for k in c.atoms
            k in live || push!(live, k)
        end
    end
    return nothing
end

# Every instance of every model, as `(model, instance)` pairs
_all_instances(rt::AutomatonRuntime) = [(km, i) for (km, m) in enumerate(rt.models) for i in 1:n_instances(m)]

# The discrete quantities of instance `i` at what the modes hold; whether any changed (the
# values that read a changed one lose their sides, `_quantity_changed!`)
function _follow_holds!(m::AutomatonModel, i::Int)
    changed = false
    for a in m.automata
        for (column, ref) in a.holds[a.modes[i]]
            v = _ref_value(ref, i)
            column[i] === v && continue
            column[i] = v
            _quantity_changed!(m, i, column)
            changed = true
        end
    end
    return changed
end

# Whether the values of any model of `rt` read another instance's discrete state or memory
_has_cross(rt::AutomatonRuntime) = any(m -> m.reads.cross, rt.models)

"""
    _expand_cross(rt::AutomatonRuntime, instances) -> Vector{Tuple{Int, Int}}

The `instances` (`(model, instance)` pairs) with every instance of a model whose values read
another instance's discrete state or memory (`ValueReads` `cross`) present once an instance
of it is (a new vector; `instances` itself where no model has any, `_has_cross`): what such
an instance's values read changes with another's transitions, so a settle considers all of
them.
"""
function _expand_cross(rt::AutomatonRuntime, instances)
    _has_cross(rt) || return instances
    out = Tuple{Int, Int}[]
    expanded = Set{Int}()
    for (km, i) in instances
        if rt.models[km].reads.cross
            km in expanded && continue
            push!(expanded, km)
            append!(out, ((km, j) for j in 1:n_instances(rt.models[km])))
        else
            push!(out, (km, i))
        end
    end
    return out
end

"""
    _quantity_changed!(m::AutomatonModel, i, column)

A transition changed the discrete quantity or memory `column` of instance `i` of `m` (a
mode it entered holds another number, or a set wrote another value). The values of `i`
that read it (`ValueReads` `state`, traced from the kernel; those of other instances too
where a value reads another's, `ValueReads` `cross`) move at the instant by that
change: a jump of the value, not its motion along the trajectory. Each that has a side
(`AutomatonModel` `side`: located at this instant, or left with one by an earlier) loses
it, so that it is decided by its value from here on (`atom_holds`), in the rounds of this
settle and after the instant (no anchor: `_instant_left!`). A side it kept would decide
its atoms as the crossing left them, though the value is now elsewhere: a timer whose
transition resets its own start would fire again at once.
"""
function _quantity_changed!(m::AutomatonModel, i::Int, column::Vector{Float64})
    q = findfirst(c -> c === column, m.quantity_columns)
    q === nothing && error("A quantity set by the automata of $(m.name) is not one of its discrete quantities or memory")
    _lose_readers!(m, m.reads.state, (q - 1) * n_instances(m) + i)
    return nothing
end

# Each value that reads column `column` of `matrix` (`ValueReads` `state` or `pattern`: the quantity or variable
# instance whose entries it reads) loses its side (decided by its value from here on, `_lose_side!`) unless it has
# none; `visit(instance)` is told the instance of each reader
function _lose_readers!(m::AutomatonModel, matrix::SparseMatrixCSC, column::Int, visit::F = _ -> nothing) where {F}
    n_values = length(m.value_columns)
    rows = rowvals(matrix)
    for p in nzrange(matrix, column)
        reader, v = fldmod1(rows[p], n_values)
        (m.side[reader, v] == 0 && !m.zero[reader, v]) || _lose_side!(m, reader, v)
        visit(reader)
    end
    return nothing
end

# Value `v` of instance `i` of `m` without a side: decided by its value, measured from zero
function _lose_side!(m::AutomatonModel, i::Int, v::Int)
    _set_side!(m, i, v, Int8(0))
    _zero!(m, i, v, false)
    _measure!(m, i, v, false)
    m.anchor[i, v] = 0.0
    m.at_instant[i, v] = false
    m.pinned[i, v] = false
    m.returning[i, v] = false
    return nothing
end

# The same of every instance of `rt`, or of `instances` (`(model, instance)` pairs)
_follow_holds!(rt::AutomatonRuntime) = _follow_holds!(rt, _all_instances(rt))

function _follow_holds!(rt::AutomatonRuntime, instances)
    changed = false
    for (km, i) in instances
        changed |= _follow_holds!(rt.models[km], i)
    end
    return changed
end

# The rounds of the start: with a network to solve (`StartResolve`, automata.jl) the algebraic point is solved again
# between them; without, they read the point they start from
_settle_start!(rt::AutomatonRuntime, ::Nothing, x, y, t, jump, integrator) =
    settle!(rt, x, y, t; jump, integrator, starting = true)

"""
    _memory_started!(jump, integrator, rt::AutomatonRuntime, t)
    _transitions_firing!(jump, integrator, rt::AutomatonRuntime, steps, t)

Extension points: what an augmentation carried through the integration beside the DAE
state (`jump`, with its `integrator`) does at two points of the settle: `_memory_started!`
once every memory is at its initial value (before the modes start), and
`_transitions_firing!` with the steps of each round of a settle, before their sets and
mode changes apply (the values are those the sets read). The defaults, for
`jump === nothing` (the plain simulation), do nothing; an augmentation that must track a
start or a transition (e.g. its sensitivities' jump at the switching time) adds a method
for its own `jump` type.
"""
_memory_started!(::Nothing, integrator, rt::AutomatonRuntime, t) = nothing
_transitions_firing!(::Nothing, integrator, rt::AutomatonRuntime, steps, t) = nothing

"""
    start!(rt::AutomatonRuntime, x, y, t; jump=nothing, integrator=nothing, host=nothing, resolve=nothing)

Start the automata at `(x, y, t)`: every memory at its initial value, each automaton of
each instance in the mode of its first initial clause that holds, the discrete
quantities at what the modes hold, and then a settle in rounds as at any instant
(`settle!`), in which the conditions that hold at the start fire, and the algebraic point is solved again
between the rounds from the modes entered. `resolve` (a `StartResolve`: the DAE `host`, an integrator or `FixedInstant`
whose state `u` holds `x` and `y` as views, and the reinitialization) is what solves it; without it
(a runtime run alone, with no network) the rounds read the point they start from. Clears what located events recorded (the sides
and anchors of the values, the crossings) and the switches counted for Zeno. `jump` and `integrator`: see `_memory_started!`.

A value zero at the start to its rounding (exactly, or by a few ulps of the terms it subtracts,
`within_rounding`) is on a surface of its automata (an anti-windup state at its limit, the
argument of a piece at its breakpoint: a start on a kink, e.g. where the start moved a limit to
the state, `Qmax = Qinv`, or computed one onto it), and the start is no crossing of it: the
integration's start (`host`, the DAE it integrates) decides it before the initial clauses
(`_start_on_surfaces!`), as a transition decides a value it holds at zero by continuity
(`_held_at_zero!`). From an equilibrium every rate is zero and it is decided as the exact zero
it is (`zero`, as a value of `Automaton.at_rest`) until an instant changes the equations;
elsewhere a nonzero rate along the start gives it the side it leaves the zero to, so
the initial clauses and the settle take the piece the rate leaves the start into. A start
without a host (the adapter's, which only sets the discrete state its evaluations read, and
the fixed-step Y-bus integration's) decides such a value by its value, as any other.

The start is from equilibrium, the initialization's: the values it holds at zero
(`Automaton.at_rest`) are decided as the exact zero they are (`AutomatonModel` `zero`),
in the initial clauses, the settle and after, until an instant changes the equations
(`_release_zeros!`). Their rate is zero too, the whole state being at rest, so the
initial clauses decide by the zero itself: a comparison `>=` or `<=` holds and `>` or `<`
does not, as it does for the exact values the initialization computes, and the clauses
choose, in their order, the mode the model was initialized in. A mode whose transitions
compare such a value strictly (an anti-windup lag at its limit with its rate at rest:
`free` to `upper` needs `rate > 0`) is not left, and its dynamics leave the value at rest.

A mode that governs a state (`Mode.states`) must find it on its right side at the start, in an initial mode
or a mode the first settle enters: the start refuses a state off its relation (`_check_start_relations!`),
which it cannot move (the right side is computed from the model, as a computed bound of a `VariableLimit`).

Returns whether the settle of the start fired any transition.
"""
function start!(rt::AutomatonRuntime, x, y, t; jump = nothing, integrator = nothing, host = nothing,
                resolve = nothing)
    evaluate_values!(rt, x, y, t)
    empty!(rt.switched_automata)
    empty!(rt.held)
    for m in rt.models
        fill!(m.side, 0)
        fill!(m.measured, false)
        fill!(m.n_measured, 0)
        fill!(m.anchor, 0.0)
        fill!(m.at_instant, false)
        fill!(m.pinned, false)
        fill!(m.returning, false)
        fill!(m.zero, false)
        fill!(m.n_zero, 0)
        for v in m.at_rest, i in 1:n_instances(m)
            _zero!(m, i, v, true)
            _measure!(m, i, v, true)
        end
        fill!(m.crossed, 0)
        fill!(m.last_switch, NaN)
        fill!(m.rounds, 0)
        for (column, ref) in m.memory_initial, i in 1:n_instances(m)
            column[i] = _ref_value(ref, i)
        end
    end
    _memory_started!(jump, integrator, rt, t)
    rt.sync_memory!()
    evaluate_values!(rt, x, y, t)  # values may read memory
    host === nothing || _start_on_surfaces!(rt, host, x, y, t, minimum(integrator.opts.abstol))
    for m in rt.models, a in m.automata, i in 1:n_instances(m)
        clause = findfirst(((c, _),) -> c === nothing || holds(c, m, i), a.initial)
        a.modes[i] = a.initial[clause][2]
        a.held_now[i] = false
    end
    _follow_holds!(rt)
    evaluate_values!(rt, x, y, t)  # values may read discrete quantities
    changed = _settle_start!(rt, resolve, x, y, t, jump, integrator)
    # (the relation of a state a mode governs is checked at the settled point only: a mode of an earlier round
    # may be left in a later one, and the point between rounds is provisional)
    rt.governed && _check_start_relations!(rt, x, y, t)
    _instant_left!(rt, (), x, y, t)  # (no value is located at the start)
    return changed
end

"""
    _start_on_surfaces!(rt::AutomatonRuntime, host, x, y, t, tol)

At the integration's start `(x, y, t)` of the DAE of `host`: the values of the in-service
instances on a surface (zero to their rounding: a start on a kink, `_on_surfaces`), other than
those decided as zero already (`Automaton.at_rest`), decided by their rate at the start, as
`_held_at_zero!` decides a value a transition holds at zero. The start is in the steady state the
host declares (`_start_state_rates`) where its residuals agree with it, `|F - M·rate| <= tol` on
the rows with mass (`tol` the integration's absolute tolerance; `validate_tds_initialization`
checks the same):

- **From an equilibrium** (in a steady state with every state at rest, and the residuals at the
  start the same at a later time, `_start_at_equilibrium`), every rate is zero, of every order:
  such a value is decided as the zero it is (`zero`), until an instant changes the equations
  (`_release_zeros!`). Its computed rates, the rounding of the initialization, decide nothing.
- **Elsewhere**, its rate along the start (a central difference of the values over `±h`, at
  `(x ± h ẋ, t ± h)`) gives it the side it leaves the zero to where it is nonzero (`side`,
  through the start as at an instant, `at_instant`), and the start leaves it there (`rt.held`,
  `_instant_left!`). A zero rate there (a value at rest in a network that oscillates, or one a
  drive moves at a higher order) leaves it decided by its value. The states move at the steady
  state's rates; at a start in none (a model driven from its start) at their rates `F/M` on the
  rows with mass, and the variables without mass `z` (the algebraic variables, and the states
  whose mass is zero) at the rate their equations give, `F_z ż = -(F_x ẋ + F_t)`
  (`_start_algebraic_rate`): a value that reads one that moves (a bus voltage without
  capacitance, an interface angle) moves with it. In a steady state the algebraic variables are at
  rest (the phase quantities of a sinusoidal steady state are states, whose rates it declares) and
  are held there. Where `F_z` is singular the start stops, naming the values whose rate needed it:
  it never takes the rate of the variables held at rest instead.

Before this, such a value was decided by its value alone: the initial clauses took the piece of
a kink by their order (the order of the arguments of the piece), or by the sign of the rounding
of a value computed onto it, and the first step, moving the value off its zero by its rounding
alone, was located as a crossing of it at the start (a grazing, which an augmentation carrying
derivatives through it refuses).
"""
function _start_on_surfaces!(rt::AutomatonRuntime, host, x, y, t, tol::Real)
    on = _on_surfaces(rt, x, y, t)
    isempty(on) && return nothing
    nd = length(x)
    declared = _start_state_rates(host, nd)
    F = zeros(nd + length(y))
    dae_residuals!(F, host, vcat(x, y), t)
    steady = all(k -> iszero(rt.mass[k]) || abs(F[k] - rt.mass[k] * declared[k]) <= tol, 1:nd)
    if steady && _start_at_equilibrium(host, declared, x, y, t)
        for (km, i, v) in on
            m = rt.models[km]
            _zero!(m, i, v, true)
            _measure!(m, i, v, true)
        end
        return nothing
    end
    # The rate of the states at the start: the steady state's where the start is in it; else F/M
    # on the rows with mass and the declared rates on the states without (which, in a driven start,
    # `_start_algebraic_rate` replaces by the rate their equations give)
    rates = steady ? declared : [iszero(rt.mass[k]) ? declared[k] : F[k] / rt.mass[k] for k in 1:nd]
    # The rate of the whole state: in a steady state the variables without mass at rest; in a driven start at the
    # rate their equations give (`_start_algebraic_rate`)
    h = START_RATE_STEP
    udot = steady ? vcat(rates, zeros(length(y))) : _start_algebraic_rate(host, x, y, t, rates, h)
    udot isa AbstractString && throw(ErrorException(
        "Automaton start at t = $t: the rate of " * join((_value_name(rt.models[km], i, v) for (km, i, v) in on), "; ") *
        ", on its surface at the start, needs the rate of the variables without mass (the algebraic variables), " *
        "F_z ż = -(F_x ẋ + F_t), and $udot. The start does not hold them at rest instead, whose rate is not the start's."))
    # The rates of the values along it: a central difference over ±h
    plus, minus = _values_along(rt, x, y, t, udot, h), _values_along(rt, x, y, t, udot, -h)
    evaluate_values!(rt, x, y, t)
    offsets = cumsum([0; [length(m.value_columns) for m in rt.models]])
    for (km, i, v) in on
        m = rt.models[km]
        k = offsets[km] + v
        rate = (plus[k][i] - minus[k][i]) / (2h)
        isfinite(rate) || throw(ErrorException(
            "Automaton start at t = $t: the rate of $(_value_name(m, i, v)), on its surface at the start, is $rate."))
        iszero(rate) && continue
        _set_side!(m, i, v, rate > 0 ? Int8(1) : Int8(-1))
        m.at_instant[i, v] = true
        push!(rt.held, (km, i, v))
        _measure!(m, i, v, true)
    end
    return nothing
end

# Value `v` of instance `i` of `m`, for messages: its name, the model and the instance's key
_value_name(m::AutomatonModel, i::Int, v::Int) = "$(keys(m.state.values)[v]) of $(m.name) ($(join(string.(m.keys[i]), ", ")))"

"""
    _on_surfaces(rt::AutomatonRuntime, x, y, t) -> Vector{Tuple{Int, Int, Int}}

The values (`(model, instance, value)`) of the in-service instances at `(x, y, t)` on a surface
of their automata: zero to their rounding (`within_rounding`, from the running error bound of
each, `RoundingBound`, which the model's kernel of the values computes from the terms the value
adds and subtracts, `bounds`), an exact zero among them; other than those decided as zero
already (`AutomatonModel` `zero`). The values of `rt` are those at `(x, y, t)`: its last
evaluation (`evaluate_values!`) was there, and the bounds read the variables it gathered.
"""
function _on_surfaces(rt::AutomatonRuntime, x, y, t)
    on = Tuple{Int, Int, Int}[]
    _on_surfaces!(on, rt.bounds, rt.models, 1, t)
    return on
end

# Each model's bounds (a tuple of evaluators, unrolled)
_on_surfaces!(on, ::Tuple{}, models, km, t) = nothing
function _on_surfaces!(on, bounds::Tuple, models, km, t)
    m = models[km]
    columns = values(first(bounds)(t))
    for v in eachindex(m.value_columns)
        value, bound = m.value_columns[v], columns[v]
        for i in 1:n_instances(m)
            (m.active[i] && !m.zero[i, v] && within_rounding(value[i], bound[i].e)) && push!(on, (km, i, v))
        end
    end
    return _on_surfaces!(on, Base.tail(bounds), models, km + 1, t)
end

# The values of every model (copies of their columns) at the start `(x, y, t)` moved by `s` along `udot`, the rate of
# the whole state `[x; y]`, and the time by `s`
function _values_along(rt::AutomatonRuntime, x, y, t, udot, s)
    nd = length(x)
    u = vcat(x, y) .+ s .* udot
    evaluate_values!(rt, view(u, 1:nd), view(u, (nd + 1):length(u)), t + s)
    return [copy(c) for m in rt.models for c in m.value_columns]
end

"""
    _start_algebraic_rate(host, x, y, t, rates, h) -> Union{Vector{Float64}, String}

The rate of the whole state `[x; y]` at a driven start `(x, y, t)` of the DAE of `host`, for a
value on a surface there (`_start_on_surfaces!`): the states with mass at `rates`, the variables
without mass `z` (the algebraic variables, and the states whose mass is zero) at the rate their
equations give, `F_z ż = -(F_x ẋ + F_t)` on the rows without mass. `F_x ẋ + F_t` is a central
difference of those rows over `±h` along `(ẋ, 1)` at the start's `z`, the residual the start
leaves there cancelling in it; `F_z` the Jacobian of the host's reinitialization, factorized
with pivoting by its solver. Where `F_z` is singular (KLU finds it so, or estimates its condition
number at or above `1 / eps()`), what failed, as text. A host adds its method: the
`TDSAdapter`'s (reinit.jl), the Y-bus DAE's (PowerfulYbus).
"""
function _start_algebraic_rate end

# The step of the central difference of the values' rates at a start that is no equilibrium
# (`_start_on_surfaces!`): small against the period of a sinusoidal steady state's oscillation
# (the error, (ω h)²/6 of the rate, is 2e-8 at 60 Hz) and large against the rounding of the values
const START_RATE_STEP = 1e-6

"""
    _start_at_equilibrium(host, rates, x, y, t) -> Bool

Whether the start `(x, y, t)` of the DAE of `host`, in the steady state it declares (see
`_start_on_surfaces!`), is an equilibrium: every state at rest in it (`rates`,
`_start_state_rates`, all zero) and the residuals at the start the same at a later time
(bitwise: no residual reads the time there, as a drive or a model whose residuals depend
on time directly would). Then every rate of the start is zero, of every order.
"""
function _start_at_equilibrium(host, rates::AbstractVector, x, y, t)
    all(iszero, rates) || return false
    u = vcat(x, y)
    now, later = similar(u), similar(u)
    dae_residuals!(now, host, u, t)
    dae_residuals!(later, host, u, t + START_TIME_PROBE)
    return now == later
end

# The later time `_start_at_equilibrium` evaluates the residuals at: a time no period of a
# drive in the tests divides (an irrational multiple of a millisecond)
const START_TIME_PROBE = sqrt(2.0) * 1e-3

# The steps instance `i` of `m` (index `km`) takes in this round (see `settle!`), one per automaton that has one: the
# automaton (index), the transition (or nothing for a held clause) and the mode it goes to. Every automaton reads the
# values of the point before the round, so the steps are independent of one another and of the order of the automata.
# Pushes each onto `steps` (`stepped` gets `(km, i)` for each) and updates `held_now`; returns whether it pushed any.
function _round_steps!(steps, stepped, m::AutomatonModel, km::Int, i::Int, first_round::Bool, t)
    any_step = false
    for (ka, a) in enumerate(m.automata)
        mode = a.modes[i]
        held = findfirst(c -> holds(c, m, i), a.held)
        if held !== nothing
            a.held_now[i] = true
            target = a.held_mode[held]
            target == mode && continue
            push!(steps, (m, i, ka, nothing, target))
            push!(stepped, (km, i))
            any_step = true
            continue
        end
        a.held_now[i] = false
        for tr in a.transitions
            tr.from == mode || continue
            fires = tr.crossing == 0 ? holds(tr.compiled::CompiledCondition, m, i) :
                    first_round && m.crossed[tr.crossing, i] == tr.direction
            fires || continue
            _refuse_second!(m, i, a, tr, first_round, t)
            push!(steps, (m, i, ka, tr, tr.to))
            push!(stepped, (km, i))
            any_step = true
            break
        end
    end
    return any_step
end

# Throw if a transition after `first`, the one a round of a settle takes for automaton `a` of instance `i` of `m`, from
# the same mode, holds too: two that hold at once are refused, never ordered by their declaration, as the semantic
# core's reference settle has it. Decided for the transition that fires only, so that a settle in which
# nothing fires does no more than before
@noinline function _refuse_second!(m::AutomatonModel, i::Int, a::CompiledAutomaton, first::CompiledTransition,
                                   first_round::Bool, t::Real)
    mode = a.modes[i]
    transitions = a.transitions
    for j in (findfirst(tr -> tr === first, transitions) + 1):length(transitions)
        tr = transitions[j]
        tr.from == mode || continue
        fires = tr.crossing == 0 ? holds(tr.compiled::CompiledCondition, m, i) :
                first_round && m.crossed[tr.crossing, i] == tr.direction
        fires && _simultaneous(m, i, a, first, tr, t)
    end
    return nothing
end

# Throw for the transitions `first` and `second` of automaton `a` of instance `i` of `m`, which hold at once
function _simultaneous(m::AutomatonModel, i::Int, a::CompiledAutomaton, first::CompiledTransition,
                       second::CompiledTransition, t)
    modes = a.automaton.modes
    from = modes[first.from].name
    throw(ErrorException(
        "Automaton settle at t = $t: $(m.name) ($(join(string.(m.keys[i]), ", "))), automaton $(a.automaton.name): " *
        "the transitions $from -> $(modes[first.to].name) and $from -> $(modes[second.to].name) hold at once. " *
        "Two transitions holding at once from one mode are refused, never ordered (the semantic core's reference " *
        "settle): the conditions of the transitions from a mode must exclude one another."))
end

"""
    AutomatonStep

One step of a round of a settle: instance `i` of `model`, whose automaton `automaton`
(index) goes to mode `target` by `transition` (`nothing` for a held clause).
"""
const AutomatonStep = Tuple{AutomatonModel, Int, Int, Union{Nothing, CompiledTransition}, Int}

# The values of the automata of `instances` (`nothing`: every instance) at `(x, y, t)`, which must be finite for every
# in-service instance: a nonfinite value (or a value outside its domain, which the kernel returns as NaN) fails at once, in
# any round of an instant, never decides a guard silently
function _evaluate_checked!(rt::AutomatonRuntime, instances, x, y, t)
    evaluate_values!(rt, instances, x, y, t)
    for (km, i) in (instances === nothing ? _all_instances(rt) : instances)
        m = rt.models[km]
        m.active[i] || continue
        for (v, column) in enumerate(m.value_columns)
            isfinite(column[i]) && continue
            throw(ErrorException(
                "Automaton settle at t = $t: the value $(keys(m.state.values)[v]) of $(m.name) " *
                "($(join(string.(m.keys[i]), ", "))) is $(column[i]), not finite: a value outside its domain fails " *
                "at once, in any round of an instant"))
        end
    end
    return nothing
end

"""
    AlgebraicRounds(before!, changed)

What a host gives a settle (`settle!`) to decide, round by round, whether the transitions of the round changed the
equations without mass, so that the algebraic point must be solved again before the next round (the rule: an instant
settles in rounds, and the algebraic point of a round is solved from the modes it entered). `before!(instances)` is
called with the instances (`(model, instance)` pairs) about to take steps, before the round applies them, and prepares
the comparison; `changed(instances)` is called after the round, with the same instances, and tells whether the
equations without mass differ (`_algebraic_before!`, `_algebraic_changed!`).
"""
struct AlgebraicRounds{B, C}
    before!::B
    changed::C
end

"""
    settle!(rt::AutomatonRuntime, x, y, t; located=false, instances=nothing, jump=nothing,
            integrator=nothing, rates=nothing, starting=false, rounds=nothing) -> Bool

Fire the transitions of the automata at `(x, y, t)` until nothing changes, and return
whether anything did (a mode, a discrete quantity or a memory). The discrete quantities
first follow what the modes hold (they already do after a start or a settle: holds are
numbers, `check_automata`); then, in rounds: in a round every transition that holds, of every
automaton of every in-service instance, fires (`_round_steps!`: at most one per automaton, two that
hold at once from one mode are refused), all on the values of the point before the round, the sets of each
transition reading the values from before it; the discrete quantities follow the new modes and the values are
recomputed. Crossings fire in the first round only, and only where a located event
(`located`) recorded them. Throws where an instance switches without end at one instant
(Zeno, `_count_rounds!`, which counts rounds). Lists the components of the automaton callback that the modes
and held clauses it leaves can fire (`_watch!`), and the instances whose discrete
state changed (a mode, a discrete quantity or a memory) in `rt.switched`. `jump` and
`integrator`: see `_transitions_firing!`.

`rounds` (`AlgebraicRounds`): the host's test of whether a round changed the equations without mass. A round that did
ends the settle (`rt.settle_equations`), without recomputing the values: the caller solves the algebraic point
again from the modes entered (`settle_instant!`) and settles again, every instance, on the new point. A round that did
not leaves the point as it is (it is still a solution), and the next round reads it. Without `rounds` (a runtime
run alone, with no network to solve) no round changes the point. A governed state written where an algebraic equation
reads it (`_enter_governed!`) ends the settle the same way (`rt.settle_partial`).

After each round, the values the transitions that fired hold at zero by continuity
(`Transition.continuity`, and `Transition.crossing_continuity` where the round is the first
of a located instant; except the values a writer listed in `Transition.continuity_writers` fired with in this
round voided, and those that read a quantity the network solves in a round in which another transition changes a quantity
delivered to it, `_held_at_zero!`)
are decided by their rate in the discrete state entered
(`_held_at_zero!`, from `rates`, an `AfterRates`: the rates of the values in the current
discrete state; `nothing` where no rate can be evaluated).

`starting` (the start of the automata, `start!`): a mode a transition enters that governs a state
finds it on its right side or refuses the start, writing nothing (`_enter_governed!`).

`instances` (`(model, instance)` pairs) are those that can take a step in the first
round, the others' conditions unchanged since the last settle (at a located event, the
instances of the values located); every instance where `nothing`. After a round that left the
algebraic point as it was, only the instances that took a step can take another: at the same
`(x, y, t)` the values of an instance read its own discrete quantities and memory only (the contract
of the kernels, `AutomatonState`), which its own steps alone change. The work of a settle is then that of
the instances that switch.
"""
# Set the `enforce_governed!` suspend flag (round 19a, `AutomatonState.suspended`, shared by every
# element-type workspace of a model) of every model of `rt` that governs a state
function _set_suspended!(rt::AutomatonRuntime, v::Bool)
    for m in rt.models
        isempty(m.governed.states) || (m.state.suspended[] = v)
    end
    return nothing
end

function settle!(rt::AutomatonRuntime, x, y, t; located::Bool = false,
                 instances::Union{Nothing, Vector{Tuple{Int, Int}}} = nothing, jump = nothing, integrator = nothing,
                 rates = nothing, starting::Bool = false, rounds::Union{Nothing, AlgebraicRounds} = nothing,
                 event::Bool = false)
    rt.settle_partial[] = false
    rt.settle_equations[] = false
    # A fixed-time event's first round decides its transitions by comparing the *retained* state
    # against the event's new data (the header of governed_states.jl): suspend `enforce_governed!`
    # (round 19a) for exactly that comparison, so a guard that is itself the relation minus the
    # governed state (REGCA1's position atoms, PIReset's release condition) reads the real,
    # un-enforced state, not an identical zero that could never release. Resumed once round 1's
    # steps are decided, below (`_set_suspended!(rt, false)`), before `_enter_governed!` -- which
    # itself writes the state onto its relation where a mode enters or the event's own repair fires
    # -- and before any later round, so every other evaluation of this settle sees the enforced view.
    suspend_for_event = event && rt.governed
    suspend_for_event && _set_suspended!(rt, true)
    candidates = instances === nothing ? _all_instances(rt) : _expand_cross(rt, instances)
    switched = rt.switched
    empty!(switched)
    for (km, i) in candidates
        _follow_holds!(rt.models[km], i) && push!(switched, (km, i))
    end
    changed = !isempty(switched)
    changed && evaluate_values!(rt, instances, x, y, t)
    first_round = true
    memory_set = false
    steps = AutomatonStep[]
    round = candidates
    stepped = Tuple{Int, Int}[]
    delivers = Bool[]
    while true
        empty!(steps)
        empty!(stepped)
        # (the first round of a settle prepares the comparison for the instances it may take steps of, as the settle
        # before it did; a later round, and any round of a start, for the instances that took steps, once it has any)
        first_round && !starting && rounds !== nothing && rounds.before!(instances)
        for (km, i) in round
            m = rt.models[km]
            m.active[i] || continue
            _round_steps!(steps, stepped, m, km, i, first_round && located, t)
        end
        # (what tells a crossing located while integrating the source mode from a later round)
        at_crossing = first_round && located
        later_round = !first_round
        first_round = false
        # At a fixed-time event, the first round also repairs a governed state whose right side the event moved in
        # a mode that stayed active (`_repair_writes`, through the same write as an entry, `_enter_governed!`): this
        # must be considered even where no transition fires this round at all (the ordinary case: the mode holding
        # the state was already active before the event), so the "nothing fired" exit below is deferred until it
        # has had its chance. Never at a later round (the header of this file) or at a start (`event` is never true
        # there).
        do_repair = event && !later_round && rt.governed
        # Round 1's comparison against the retained state is made (immediately above); resume
        # enforcement before `_enter_governed!` (below) and every later round of this settle
        suspend_for_event && !later_round && _set_suspended!(rt, false)
        (isempty(steps) && !do_repair) && break
        taking = Tuple{Int, Int}[]
        if !isempty(steps)
            changed = true
            taking = unique(stepped)
            append!(switched, taking)
            _count_rounds!(rt, steps, t, rt.max_rounds[])
            (later_round || starting) && rounds !== nothing && rounds.before!(taking)
            _transitions_firing!(jump, integrator, rt, steps, t)
            _deliveries!(delivers, steps)
            # Every set reads the values from before the transitions of this round (a memory
            # set to another value moves the values that read it, `_quantity_changed!`)
            for (m, i, ka, tr, target) in steps
                tr === nothing && continue
                for (column, ref) in tr.sets
                    v = _ref_value(ref, i)
                    column[i] === v || _quantity_changed!(m, i, column)
                    column[i] = v
                    memory_set = true
                end
            end
            for ((km, i), (m, _, ka, _, target)) in zip(stepped, steps)
                m.automata[ka].modes[i] = target
                push!(rt.switched_automata, (km, i, ka))
            end
            # The discrete quantities follow the new modes (holds are numbers, which read no
            # value); the values are computed again below, at the point the next round reads
            _follow_holds!(rt, taking)
        end
        # A mode entered puts the states it governs on their right sides (governed_states.jl), and (`do_repair`)
        # a state whose right side the event moved in a mode that stayed active is put back on its relation the
        # same way: the values that read them move by that, and an algebraic equation that reads one ends the
        # round here, to be reinitialized by the settle of the instant before the next round
        entered = (rt.governed && (!isempty(steps) || do_repair)) ?
            _enter_governed!(rt, steps, stepped, x, y, t, starting, jump, integrator; event = do_repair) : nothing
        entered !== nothing && (changed = true)
        # `_repair_writes` (above) has now read `rt.touched` for this instant's first round, whether or not it
        # found anything: clear it so a later round, or the next event, does not see it again
        do_repair && empty!(rt.touched)
        _held_at_zero!(rt, steps, stepped, delivers, rates, t, at_crossing)
        # The equations without mass the modes entered changed: the algebraic point is solved again by the caller,
        # from them, before the next round reads it
        if rounds !== nothing && rounds.changed(taking)
            rt.settle_equations[] = true
            entered !== nothing && entered.partial && (rt.settle_partial[] = true)
            break
        end
        round = _has_cross(rt) ? _expand_cross(rt, taking) : taking
        if entered !== nothing
            round = copy(round)
            append!(round, entered.readers)
            unique!(round)
            if entered.partial
                rt.settle_partial[] = true
                break
            end
        end
        _evaluate_checked!(rt, taking, x, y, t)
    end
    memory_set && rt.sync_memory!()
    unique!(switched)
    instances === nothing ? _watch!(rt) : _watch!(rt, candidates)
    return changed
end

# Whether each step of a round, before it applies, changes a quantity delivered to the network (`delivers[s]` for
# step `s`): a discrete quantity its mode entered holds, or a memory it sets, to another number, that a residual
# without mass reads (`AlgebraicCheck` `inputs`), or a state the mode entered governs where a residual without mass
# reads one (`AlgebraicCheck` `switches`)
function _deliveries!(delivers::Vector{Bool}, steps)
    resize!(delivers, length(steps))
    for (s, (m, i, ka, tr, target)) in enumerate(steps)
        delivers[s] = _delivers_change(m, i, ka, tr, target)
    end
    return delivers
end

function _delivers_change(m::AutomatonModel, i::Int, ka::Int, tr::Union{Nothing, CompiledTransition}, target::Int)
    inputs = m.algebraic.inputs
    reads(column) = (q = findfirst(c -> c === column, m.quantity_columns); q !== nothing && q <= size(inputs, 1) && inputs[q, i])
    for (column, ref) in m.automata[ka].holds[target]
        column[i] === _ref_value(ref, i) || reads(column) && return true
    end
    if tr !== nothing
        for (column, ref) in tr.sets
            column[i] === _ref_value(ref, i) || reads(column) && return true
        end
    end
    return m.algebraic.switches[i] && !isempty(m.governed.states) && !isempty(m.governed.by_mode[ka][target])
end

"""
    _count_rounds!(rt, steps, t, budget)

Count the round at time `t` of each instance that takes a step in it (once for all its steps): one more
round at the instant it last switched where that was within `ZENO_ULPS` ulps of `t`, the
first round of a new instant otherwise. Throws when an instance needs more than `budget`
rounds at one instant (Zeno, or chattering about a threshold its own switching crosses):
a round still needed after the budget is spent, not the budget's being reached.
Names each such instance and the steps it takes.
"""
function _count_rounds!(rt::AutomatonRuntime, steps, t, budget::Int)
    over = false
    previous = nothing
    rt.rounds_taken[] += 1
    for (m, i, _, _, _) in steps
        (previous !== nothing && previous[1] === m && previous[2] == i) && continue
        previous = (m, i)
        last = m.last_switch[i]
        m.rounds[i] = abs(t - last) <= ZENO_ULPS * eps(max(abs(t), abs(last))) ? m.rounds[i] + 1 : 1
        m.last_switch[i] = t
        m.rounds[i] == 2 && push!(rt.multi_round, (Float64(t), m.name, i))
        over |= m.rounds[i] > budget
    end
    over || return nothing
    lines = String[]
    for (m, i, ka, tr, target) in steps
        m.rounds[i] > budget || continue
        a = m.automata[ka]
        from = a.automaton.modes[a.modes[i]].name
        to = a.automaton.modes[target].name
        how = tr === nothing ? "held clause" : "transition"
        push!(lines, "  $(m.name) ($(join(string.(m.keys[i]), ", "))): automaton $(a.automaton.name) " *
                     "$from -> $to ($how)")
    end
    throw(ErrorException(
        "Automata switch without end at t = $t: instances switched more than " *
        "$budget times at one instant, the budget of rounds of an instance (`max_rounds`; within $ZENO_ULPS ulps; Zeno, or " *
        "chattering about a threshold their own switching crosses). Still switching:\n" *
        join(lines, "\n")))
end

"""
    _held_at_zero!(rt::AutomatonRuntime, steps, stepped, delivers, rates, t, at_crossing)

After a round of transitions (`steps`, of the instances `stepped`): where a transition
fired on the zero of a value `surface` that is at its zero at this instant (located here,
held here, or decided as zero at rest), each value it holds at zero there by continuity
(`CompiledTransition` `continuity`) is decided by the sign of its rate in the discrete
state the transition entered, not by its value, which is zero to the rounding of the
reinitialization: it takes that sign as its side through the instant (`at_instant`, as a
located value does, `atom_holds`), and after it is measured from where the instant leaves
it (`rt.held`, `_instant_left!`). A transition's `crossing_continuity` values are decided the
same way, only where `at_crossing` (the round is the first of a located instant: the source mode was
integrated up to the crossing, so settled before the instant) and the surface is located (`at_instant`);
a surface at rest, a later round, a start or an event's reinitialization decides none of them. A rate that is zero decides it as the zero it is
(`zero`): at rest (a surface decided as zero), whose rates are all zero, without being
evaluated. `rates` (`AfterRates`) gives the rates of the values in the current discrete state
(`rates(cache, km, w)`: a matrix (value, instance) of model `km` with value `w` among them,
computed once per model and round in the `cache` of the round, `_rates_after`); it is
`nothing` where the integration has no rates, and a transition that needs them then throws.
"""
function _held_at_zero!(rt::AutomatonRuntime, steps, stepped, delivers, rates, t, at_crossing::Bool)
    computed = nothing
    for (s, ((km, i), (m, _, ka, tr, _))) in enumerate(zip(stepped, steps))
        tr === nothing && continue
        for (surface, held, kind2) in Iterators.flatten((
                ((s, h, false) for (s, h) in tr.continuity),
                at_crossing ? ((s, h, true) for (s, h) in tr.crossing_continuity) : ()))
            at_rest = !kind2 && m.zero[i, surface]
            (at_rest || m.at_instant[i, surface]) || continue
            for w in held
                _voided(steps, m, i, tr, w) && continue
                _moved_through_network(delivers, s, m, w) && continue
                rate = if at_rest
                    0.0
                else
                    rates === nothing && throw(ErrorException(
                        "Automaton settle at t = $t: $(m.name) ($(join(string.(m.keys[i]), ", "))), automaton " *
                        "$(m.automata[ka].automaton.name), holds $(keys(m.state.values)[w]) at zero by continuity, " *
                        "which is decided by its rate, and this integration evaluates no rates."))
                    computed === nothing && (computed = _rates_after(rates.rt, rates.before))
                    rates(computed, km, w)[w, i]
                end
                isfinite(rate) || throw(ErrorException(
                    "Automaton settle at t = $t: the rate of $(keys(m.state.values)[w]) of $(m.name) " *
                    "($(join(string.(m.keys[i]), ", "))), held at zero by continuity, is $rate."))
                # (a rate the integration does not resolve decides nothing: the value stays at its zero through the instant)
                unresolved = false
                if !at_rest
                    resolved = _resolved_rate(rates, computed, km, w, i, rate)
                    unresolved = iszero(resolved) && !iszero(rate)
                    rate = resolved
                end
                m.anchor[i, w] = 0.0
                m.pinned[i, w] = false
                m.returning[i, w] = false
                if iszero(rate)
                    _set_side!(m, i, w, Int8(0))
                    m.at_instant[i, w] = unresolved
                    _zero!(m, i, w, true)
                    unresolved && push!(rt.held, (km, i, w))
                else
                    _set_side!(m, i, w, rate > 0 ? Int8(1) : Int8(-1))
                    m.at_instant[i, w] = true
                    _zero!(m, i, w, false)
                    push!(rt.held, (km, i, w))
                end
                _measure!(m, i, w, true)
            end
        end
    end
    return nothing
end

# Whether the hold of the value `w` by the transition `tr` of instance `i` of `m` is voided: an automaton
# `tr.writers` lists for it fires in this round (`steps`), with `tr`: the fact the marker states is derived for the
# point as it was before such an automaton fired, and holds in the round in which none fires, the algebraic point having
# been solved again from the modes of the rounds before. The value is then decided by its own value, as one no marker names.
function _voided(steps, m::AutomatonModel, i::Int, tr::CompiledTransition, w::Int)
    for (v, writers) in tr.writers
        v == w || continue
        for (m2, i2, ka, _, _) in steps
            m2 === m && i2 == i && ka in writers && return true
        end
    end
    return false
end

# Whether the hold of the value `w` of `m` by the transition of step `s` of a round is not applied because the value
# reads a quantity the network solves (a variable without mass: a bus voltage, or what is computed from one) in a round in
# which another transition, of any instance, changes a quantity delivered to the network (`delivers`, `_deliveries!`): the
# fact a marker states reads the inputs of the model as unchanged across the transition, and what the other
# devices do reaches them through the network. The value is then decided by its own value, as one no marker names.
function _moved_through_network(delivers::Vector{Bool}, s::Int, m::AutomatonModel, w::Int)
    any(k -> k != s && delivers[k], eachindex(delivers)) || return false
    return reads_massless(m.reads, w, length(m.value_columns))
end

"""
    _release_zeros!(rt::AutomatonRuntime) -> Bool

At an instant that changes the equations (an event of the simulation, or a transition):
the values decided as zero (`AutomatonModel` `zero`) are decided by their values from here
on, without a side, as any value no event located. The rest that held them at zero is
over: the changed equations move what was at rest, and whether a value jumped at the
instant or moves from zero at its rate is not known here -- by its value it is right where
it jumped, and where it did not it is at zero still, its sign the rounding of the
instant's reinitialization, which decides at most whether a transition at that zero fires
now or as the value moves. Returns whether any was released.
"""
function _release_zeros!(rt::AutomatonRuntime)
    released = false
    for m in rt.models
        any(!iszero, m.n_zero) || continue
        for v in axes(m.zero, 2), i in axes(m.zero, 1)
            m.zero[i, v] || continue
            m.at_instant[i, v] && continue        # (a rate the integration did not resolve: at its zero until the instant ends)
            _zero!(m, i, v, false)
            _measure!(m, i, v, m.pinned[i, v] || !iszero(m.anchor[i, v]))
            released = true
        end
    end
    return released
end

# The index of the values column `column` among the values of `m` (`keys(state.values)`)
function _value_index(m::AutomatonModel, column::Vector{Float64})
    k = findfirst(c -> c === column, m.value_columns)
    isnothing(k) && error("A value read by the automata of $(m.name) is not one of its automaton values")
    return k
end

# Whether the transition `tr` of an automaton of `m` reads the value `column`: its guard, or
# an atom of its condition
_reads(m::AutomatonModel, tr::CompiledTransition, column::Vector{Float64}) =
    tr.guard === column || any(at -> m.atom_column[m.atom_index[at]] === column, atoms(tr.condition))

# The name of the values column `column` of model `m` (see `_value_index`)
_value_name(m::AutomatonModel, column::Vector{Float64}) = keys(m.state.values)[_value_index(m, column)]

# The discrete state of instance `i` of model `m`: the mode of each automaton, each
# discrete quantity and each memory
_discrete_state(m::AutomatonModel, i::Int) =
    [Float64[a.modes[i] for a in m.automata], Float64[c[i] for c in values(m.state.discrete)],
     Float64[c[i] for c in values(m.state.memory)]]

function _set_discrete_state!(m::AutomatonModel, i::Int, (modes, discrete, memory))
    for (a, v) in zip(m.automata, modes)
        a.modes[i] = Int(v)
    end
    for (c, v) in zip(values(m.state.discrete), discrete)
        c[i] = v
    end
    for (c, v) in zip(values(m.state.memory), memory)
        c[i] = v
    end
    return nothing
end

"""
    _instant_left!(rt::AutomatonRuntime, located, x, y, t)

Where an instant ends at `(x, y, t)` (the settles of a located event, of a fixed-time
event or of the start): each value `located` there (`LocatedValue`s), and each its
transitions held at zero by continuity with a side (`rt.held`, `_held_at_zero!`), is left
where the instant leaves it, its `anchor` its value there, with the side its crossing gave it
(`AutomatonModel` `side`; none where a transition changed what it reads, a jump,
`_quantity_changed!`), and so each value `pinned` since a refused crossing, which
is measured from there again (`_anchor`: from zero where it is on its side of it); the
atoms over them are decided from there on by their displacement from the anchor
(`atom_holds`), until a crossing is refused again. An
automaton whose mode the instant changed reads its values from the mode it entered, where
what an earlier instant left them at is not what its transitions read: the sides of the
values it reads that were not located here are cleared, and so is their zero (decided by
the value again). The
crossings recorded are forgotten.

Clears the automata switched (`switched_automata`).
"""
function _instant_left!(rt::AutomatonRuntime, located, x, y, t)
    for (km, i, ka) in rt.switched_automata
        m = rt.models[km]
        for v in m.automaton_values[ka]
            m.at_instant[i, v] && continue
            _set_side!(m, i, v, Int8(0))
            _zero!(m, i, v, false)
            _measure!(m, i, v, false)
            m.anchor[i, v] = 0.0
            m.pinned[i, v] = false
            m.returning[i, v] = false
        end
    end
    empty!(rt.switched_automata)
    # The values where the instant leaves them: of the models located and those with values
    # pinned
    fill!(rt.marked, false)
    for (km, _) in located
        rt.marked[km] = true
    end
    for (km, _, _) in rt.held
        rt.marked[km] = true
    end
    for (km, m) in enumerate(rt.models)
        any(m.pinned) && (rt.marked[km] = true)
    end
    any(rt.marked) || (rt.returning[] = 0; return nothing)
    _evaluate_marked!(rt, x, y, t)
    for (km, i, column) in located
        m = rt.models[km]
        v = _value_index(m, column)
        foreach(c -> m.crossed[c, i] = 0, axes(m.crossed, 1))
        # (a value whose discrete state the instant changed lost its side: by its value)
        m.side[i, v] == 0 && continue
        m.anchor[i, v] = _anchor(column[i], m.side[i, v])
        _measure!(m, i, v, !iszero(m.anchor[i, v]) || iszero(column[i]))
        m.at_instant[i, v] = false
    end
    # The values held at zero by continuity, as those located
    for (km, i, v) in rt.held
        m = rt.models[km]
        m.at_instant[i, v] || continue        # (decided otherwise since: a jump, a later hold)
        g = m.value_columns[v][i]
        m.anchor[i, v] = _anchor(g, m.side[i, v])
        _measure!(m, i, v, !iszero(m.anchor[i, v]) || iszero(g))
        m.at_instant[i, v] = false
        if m.zero[i, v]
            # (a value held at its zero through the instant: measured from where the instant leaves it, so that its
            # displacement is zero there whichever side of zero the rounding left it)
            _zero!(m, i, v, false)
            m.anchor[i, v] = g
            _measure!(m, i, v, true)
        end
    end
    empty!(rt.held)
    for m in rt.models, v in axes(m.pinned, 2), i in axes(m.pinned, 1)
        m.pinned[i, v] || continue
        m.pinned[i, v] = false
        m.returning[i, v] = false
        m.anchor[i, v] = _anchor(m.value_columns[v][i], m.side[i, v])
        _measure!(m, i, v, !iszero(m.anchor[i, v]) || iszero(m.value_columns[v][i]))
    end
    rt.returning[] = 0
    return nothing
end

"""
    _returned!(rt::AutomatonRuntime) -> Bool

At the start of a step, with the values of the automata at the state it starts from (the
integrated end of the step before): release each value pinned since it crossed against
its rate (`AutomatonModel` `returning`, `_moving_crossings!`) that is back on its side of
its anchor, measured from that anchor again. Returns whether any was released.
"""
function _returned!(rt::AutomatonRuntime)
    released = false
    for m in rt.models
        any(m.returning) || continue
        for v in axes(m.returning, 2), i in axes(m.returning, 1)
            m.returning[i, v] || continue
            d = m.value_columns[v][i] - m.anchor[i, v]
            (d > 0) == (m.side[i, v] > 0) && !iszero(d) || continue
            m.returning[i, v] = false
            m.pinned[i, v] = false
            rt.returning[] -= 1
            released = true
        end
    end
    return released
end

"""
    _anchor(g, side) -> Float64

Where a value with `side` (`+1` above its surface, `-1` below) is measured from after the
instant that left it at `g`: its surface, zero, where the instant left it on its side of
it, and `g` itself where the instant's settles and reinitializations (or the integration
of a value held still) left it across zero, so that its atoms hold as the crossing left
them there. A surface that followed the value to where every instant leaves it would make
the next step's dense output, which after an instant is of low order and may set off the
other way from the value's motion, cross it at once, and again after each such event.
"""
_anchor(g::Float64, side::Int8) = side > 0 ? min(g, 0.0) : max(g, 0.0)
