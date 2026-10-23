# ============================================================================
# Conditions That Change Inside One Step: Crossings Between the Sampled Points
# ============================================================================
#
# The automaton callback (automata.jl) locates a transition where a component, the signed
# truth of a condition, changes sign between two points DiffEq evaluates it at: the ends of a
# step and the `interp_points` points of `range(tprev, t, length = interp_points)` (the last is the
# end). A condition that holds or fails only between two of them goes unseen. That happens where
#
# - a value crosses its surface and comes back inside the step (an anti-windup limiter's
#   excursion past its limit that a coarse step spans, or one that a restart no longer cuts
#   short): the value has an extremum inside the step, beyond its surface, and
# - several atoms of a condition change their truths one after the other inside the step
#   (the state reaches its limit, then its rate turns outward: the condition holds between the
#   two).
#
# The atoms' signed truths at those points are what the callback computes anyway, and they are the
# prefilter (`_track!`, at each evaluation, per atom: the component is a least or greatest over its atoms, so
# the atom that binds it changes between samples and its own samples tell more). After the last sample, with
# no event at any (DiffEq stops sampling at the first), `_scan_step!` examines a component that can fire and has
# one sign at every sample (a different sign is an event DiffEq locates) only where
#
# - one of its atoms changed sign between samples while the component did not (another atom held the
#   condition: two atoms changing one after the other; exact, no rates needed), or
# - one of its atoms (or its guard) kept one sign at every sample and its sampled magnitude fell and then
#   rose (a turn toward zero; a turn away from it, or monotone samples, cannot cross),
#
# in the window between the samples around those changes, and the atoms that moved can give the condition the other
# truth (`_may_change`: some combination of their truths at the start and the other does; no rates needed) or the
# window is at an end of the step (the sampled values of an atom that did not move do not show a turn there).
#
# Only for such a component, in that window (`_examine_window!`):
#
# 1. each value its atoms (or its guard) read whose displacement from where the last instant
#    left it (`anchor`; the value itself without a side) is not decided as zero or pinned has
#    its crossings inside the window (`_crossing_count!`): one where its displacement has
#    different signs at the two ends, none or two where the signs agree. Those two come
#    from the extremum: where the rates at the ends point opposite ways, the first toward
#    zero (the rates are the derivative of the integrator's dense output carried through one
#    Dual evaluation of the values of each model, `_value_rates!`: no residual is evaluated, no
#    `F_ZZ` factorized), the extremum is located on the dense output (the Hermite of the ends
#    refined by a Newton step on the rate along the dense output, `_locate_peak`, shared with the
#    grazing of automaton_grazing.jl), and where the displacement there is beyond zero (strictly:
#    an extremum at the surface is a tangency, no crossing) the value crosses zero on each side of
#    it. Each crossing time is the zero of the displacement along the dense output (bisection);
# 2. the condition is evaluated between consecutive crossing times from the truths its atoms have
#    there, and where its truth differs from that at the start the component is folded: it is
#    reported as it is up to a time `t_m` inside that interval, and as its value there from `t_m`
#    on (`_fold!`), a continuous function that changes sign exactly where the component first
#    does. The last sample is reported folded, so that DiffEq sees the sign change; its own root
#    finding locates the time and moves the integrator back to it (`change_t_via_interpolation!`),
#    and the located event is that of any other (`_located_event!`): the atoms that changed their
#    truths there, the settle, the rates and the sliding check.
#
# What the prefilter cannot see: a turn inside the first or last tenth of the step, where the
# samples are monotone (their neighbour on one side is the end of the step), and two turns between
# adjacent samples that leave the samples monotone: an excursion shorter than a tenth of an accepted
# step that the sampled values do not show. Nothing else of the callback changes: a component with no
# turn costs a comparison of its sampled values.
# ============================================================================

using ForwardDiff: ForwardDiff

# ----------------------------------------------------------------------------
# The extremum of a function along a step (shared with the grazing of guards)
# ----------------------------------------------------------------------------

# Cubic Hermite on a step of length h (θ ∈ [0, 1]) with values w0, w1 and rates ẇ0, ẇ1
# (verified with SymPy): the value, its rate and its second derivative
@inline _hermite(w0, ẇ0, w1, ẇ1, h, θ) =
    (2θ^3 - 3θ^2 + 1) * w0 + (θ^3 - 2θ^2 + θ) * h * ẇ0 + (-2θ^3 + 3θ^2) * w1 + (θ^3 - θ^2) * h * ẇ1
@inline _hermite_rate(w0, ẇ0, w1, ẇ1, h, θ) =
    ((6θ^2 - 6θ) * w0 + (3θ^2 - 4θ + 1) * h * ẇ0 + (-6θ^2 + 6θ) * w1 + (3θ^2 - 2θ) * h * ẇ1) / h
@inline _hermite_second(w0, ẇ0, w1, ẇ1, h, θ) =
    ((12θ - 6) * w0 + (6θ - 4) * h * ẇ0 + (-12θ + 6) * w1 + (6θ - 2) * h * ẇ1) / h^2

# The zero of the rate of the Hermite from ẇ0 > 0 to ẇ1 ≤ 0 (bisection)
function _hermite_peak(w0, ẇ0, w1, ẇ1, h)
    a, b = 0.0, 1.0
    for _ in 1:60
        c = (a + b) / 2
        _hermite_rate(w0, ẇ0, w1, ẇ1, h, c) > 0 ? (a = c) : (b = c)
    end
    return (a + b) / 2
end

"""
    _locate_peak(rate_at, w0, ẇ0, w1, ẇ1, t0, t1) -> Union{Nothing, NamedTuple}

The maximum of a function `w` over the step from `t0` to `t1`, where its rate goes from
`ẇ0 > 0` to `ẇ1 ≤ 0` (`w0`, `w1` its values at the ends): the maximum of the Hermite of the
ends (`t_hermite`, with its curvature), refined by a Newton step on the rate `rate_at(t)` of
`w` itself along the integrator's dense output (`t`, the time of the maximum, within the step;
`shift` how far the step moved it). `nothing` where the Hermite has no maximum (its curvature
there is not negative).
"""
function _locate_peak(rate_at::F, w0, ẇ0, w1, ẇ1, t0, t1) where {F}
    h = t1 - t0
    θ = _hermite_peak(w0, ẇ0, w1, ẇ1, h)
    ẅ = _hermite_second(w0, ẇ0, w1, ẇ1, h, θ)
    ẅ < 0 || return nothing
    t_H = t0 + θ * h
    t = clamp(t_H - rate_at(t_H) / ẅ, t0, t1)
    return (; t, t_hermite = t_H, curvature = ẅ, shift = abs(t - t_H))
end

# Whether the displacement `d_x` at the extremum of a value with orientation `σ` (`σ d` below zero at the
# start of the step) is beyond its surface: strictly, so that an extremum at the surface (a tangency, which
# the value reaches and leaves on its own side) is no crossing
_beyond_surface(σ::Float64, d_x::Float64) = σ * d_x > 0

# ----------------------------------------------------------------------------
# The record of a step
# ----------------------------------------------------------------------------

"""
    ExtremumFold

A component of the automaton callback whose condition changes its truth inside the step being
checked (`StepExtrema`), though it has the same truth at the ends: the `component`, a time `t`
inside the interval the condition has the other truth in, and the component `s` there.
"""
struct ExtremumFold
    component::Int
    t::Float64
    s::Float64
end

# No crossings (shared, so that a value that does not cross allocates nothing)
const NO_CROSSINGS = Float64[]

# The points DiffEq samples a step at, and how many (`interp_points` of the callback)
const INTERP_POINTS = 10

"""
    StepExtrema(rt)

The record of the search inside the step the automaton callback is evaluating (see the header of
step_extremum.jl): the step (`tprev`, `t`), its sample points `ts` (those `filled`; `analyzed` once the
step is scanned), what the samples showed of each atom's signed truth and each guard (a cell of `first_neg`,
`best_mag`, `best_slot`, `flipped`, `flip_slot`: by `atom_cell` for the atoms the conditions read, `guard_cell` for the
guards; `cell_atom` the model, atom and instance of a cell (the component, for a guard), `atom_blocks` the
blocks that read an atom; the cells that `flips` and the `turns` of their magnitude toward zero, as
(cell, sample), and their signed values at the samples (`samples`); `queue` the components to examine) and the first sample of each component (`first_comp`), the
`folds` found, and the window
(`wa`, `wb`) being examined with the components at its start (`s`, scratch too). How often each value
crosses its surface inside the window is `crossings` (by model, value and instance, with the
displacements at the ends `start` and `finish` and the extremum of an excursion `peak`, `peak_d`) and
the times it does `times`. The rest is scratch: the values and rates of each model at its ends (`V0`, `R0`, `V1`, `R1`, computed for the
models in `ends`) and at a time within it (`V`, `R`), and the states and rates at those times.
"""
mutable struct StepExtrema
    tprev::Float64
    t::Float64
    ts::Vector{Float64}
    filled::Vector{Bool}
    analyzed::Bool
    atom_cell::Vector{Matrix{Int}}
    guard_cell::Vector{Int}
    first_neg::Vector{Bool}
    best_mag::Vector{Float64}
    best_slot::Vector{Int}
    flipped::Vector{Bool}
    turns::Vector{Tuple{Int, Int}}
    samples::Matrix{Float64}
    flips::Vector{Int}
    flip_slot::Vector{Int}
    cell_atom::Vector{Tuple{Int, Int, Int}}
    atom_blocks::Vector{Vector{Vector{Int}}}
    atom_base::Vector{Vector{Int}}
    queue::Vector{Int}
    queued::Vector{Bool}
    first_comp::Vector{Float64}
    folds::Vector{ExtremumFold}
    wa::Float64
    wb::Float64
    wlo::Int
    whi::Int
    kdeg::Int
    s::Vector{Float64}
    V0::Vector{Matrix{Float64}}
    R0::Vector{Matrix{Float64}}
    V1::Vector{Matrix{Float64}}
    R1::Vector{Matrix{Float64}}
    V::Vector{Matrix{Float64}}
    R::Vector{Matrix{Float64}}
    ends::Vector{Bool}
    crossings::Vector{Matrix{Int8}}
    start::Vector{Matrix{Float64}}
    finish::Vector{Matrix{Float64}}
    peak::Vector{Matrix{Float64}}
    peak_d::Vector{Matrix{Float64}}
    ua::Vector{Float64}
    ub::Vector{Float64}
    du0::Vector{Float64}
    du1::Vector{Float64}
    u::Vector{Float64}
    du::Vector{Float64}
    times::Dict{Tuple{Int, Int, Int}, Vector{Float64}}
    rounding::Vector{Matrix{Float64}}      # per model: the running error bound `e` of each value (value, instance) at the last sample
    rounding_known::Vector{Bool}           # per model: whether `rounding` holds a bound
    rounding_age::Vector{Int}              # per model: the steps since it was evaluated
    rounding_rounds::Vector{Int}           # per model: `rt.rounds_taken` when it was evaluated
    counts::ExtremumCounts
end

function StepExtrema(rt::AutomatonRuntime)
    n = rt.components.n
    matrices() = [zeros(length(m.state.values), n_instances(m)) for m in rt.models]
    # The cells the samples are followed in: each atom, per instance, and each guard
    comps = rt.components
    cells = 0
    atom_cell = [zeros(Int, n_instances(m), length(m.atoms)) for m in rt.models]
    # (every atom: which the conditions read, `live`, changes at settles)
    cell_atom = Tuple{Int, Int, Int}[]
    atom_base = [zeros(Int, length(m.atoms)) for m in rt.models]
    for (km, m) in enumerate(rt.models), k in eachindex(m.atoms)
        atom_base[km][k] = cells
        for i in 1:n_instances(m)
            atom_cell[km][i, k] = (cells += 1)
            push!(cell_atom, (km, k, i))
        end
    end
    guard_cell = zeros(Int, n)
    # (every block: which of them can fire, `conditions` and `crossings`, changes at settles)
    for b in eachindex(comps.model)
        comps.condition[b] === nothing || continue
        for i in 1:n_instances(rt.models[comps.model[b]])
            guard_cell[comps.offset[b] + i] = (cells += 1)
            push!(cell_atom, (0, 0, comps.offset[b] + i))
        end
    end
    atom_blocks = [[Int[] for _ in m.atoms] for m in rt.models]
    for b in eachindex(comps.model)
        c = comps.condition[b]
        c === nothing || foreach(k -> push!(atom_blocks[comps.model[b]][k], b), c.atoms)
    end
    return StepExtrema(NaN, NaN, zeros(INTERP_POINTS), fill(false, INTERP_POINTS), true, atom_cell, guard_cell,
                       fill(false, cells), fill(Inf, cells), zeros(Int, cells), fill(false, cells), Tuple{Int, Int}[],
                       zeros(cells, INTERP_POINTS), Int[], zeros(Int, cells),
                       cell_atom, atom_blocks, atom_base, Int[], fill(false, n), zeros(n), ExtremumFold[], NaN, NaN, 0, 0, NOT_POLYNOMIAL, zeros(n),
                       matrices(), matrices(), matrices(), matrices(),
                       matrices(), matrices(), fill(false, length(rt.models)),
                       [zeros(Int8, length(m.state.values), n_instances(m)) for m in rt.models],
                       matrices(), matrices(), matrices(), matrices(),
                       Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
                       Dict{Tuple{Int, Int, Int}, Vector{Float64}}(), matrices(), fill(false, length(rt.models)), zeros(Int, length(rt.models)), zeros(Int, length(rt.models)), rt.extrema)
end

_step_extrema!(::Nothing, out, rt, evaluators, host, jump, u, t, integrator) = nothing

"""
    _step_extrema!(ex, out, rt, evaluators, host, jump, u, t, integrator)

Called at every evaluation of the automaton callback's components `out` at `(u, t)` within the
integrator's current step, after the values and signed truths of the automata are computed there: it
follows the atoms and guards over the step's sample points (`_track!`), scans the step (`_scan_step!`) at
the last of them, when DiffEq has sampled all without finding an event, and folds `out` at that sample and at
every later evaluation where the scan found a condition that changes inside the step (`_fold!`).
"""
function _step_extrema!(ex::StepExtrema, out::AbstractVector, rt::AutomatonRuntime, evaluators::Tuple, host, jump,
                        u::AbstractVector, t::Real, integrator)
    t0, t1, t = Float64(integrator.tprev), Float64(integrator.t), Float64(t)
    t0 < t1 || return nothing
    if t == t0 && !(ex.tprev == t0 && ex.t == t1)
        # (the first evaluation of the step: DiffEq's root finding evaluates its start again)
        ex.tprev, ex.t = t0, t1
        # (the points DiffEq evaluates the step's components at, as it computes them)
        copyto!(ex.ts, range(t0, stop = t1, length = INTERP_POINTS))
        fill!(ex.filled, false)
        ex.analyzed = false
        empty!(ex.flips)
        empty!(ex.turns)
        ex.wa = ex.wb = NaN
        empty!(ex.folds)
        copyto!(ex.first_comp, out)
        _track!(ex, rt, out, 1)
        ex.filled[1] = true
    elseif ex.tprev == t0 && ex.t == t1
        k = t == t0 ? 0 : _sample_slot(ex.ts, t)
        if k > 1
            if k == INTERP_POINTS && ex.filled[k] && !ex.analyzed && all(@view ex.filled[1:(INTERP_POINTS - 1)])
                # (the last sample, after all the others and the end: nothing in the step was an event)
                ex.analyzed = true
                _track!(ex, rt, out, k)
                _derive_flags!(ex, rt)
                if !isempty(ex.flips) || !isempty(ex.turns)
                    started = time_ns()
                    _scan_step!(ex, rt, evaluators, host, jump, integrator)
                    ex.counts.nanoseconds += time_ns() - started
                end
            elseif !ex.filled[k]
                # (the end of the step is evaluated first, and again as the last sample, in order)
                k < INTERP_POINTS && _track!(ex, rt, out, k)
                ex.filled[k] = true
            end
        end
        isempty(ex.folds) || _fold!(out, ex, t)
    end
    return nothing
end

# A magnitude is less than the least before it, for a turn toward zero, where it is less by more than the rounding of
# the evaluation: the sampled values of a value that is constant (a limit, a parameter) differ by the rounding of the
# evaluations and are no motion toward its surface (every such step had an inner least, and was examined with a Dual evaluation
# of the model's values: two thirds of the time of a representative run). The rounding of a value `g` is bounded by `ROUNDING_ULPS u e`
# with `e` its running error bound (`RoundingBound`, rounding_bound.jl; `u = eps()/2`), and `e >= |g|` for every value: `ROUNDING_NOISE`
# times `|g|` is the part of it every value has, which the track applies at every sample with no evaluation of bounds;
# the rest (a value that cancels terms larger than itself: `e` from its terms) `_derive_flags!` applies, per cell with a turn
const ROUNDING_NOISE = ROUNDING_ULPS * (eps(Float64) / 2)

"""
    TURN_FILTER

Test hook: where `false` every decrease of a sampled magnitude is a turn (the scan before the rounding filter), to pair a run with the
same run filtered: the filter changes no decision, only how many windows are examined.
"""
const TURN_FILTER = Ref(true)

# What sample `slot` (1 at the start of the step, then in order) shows of the atoms the watched conditions read
# (their signed truths, `m.signed`) and the guards (the components `out`): their signed values (`samples`), whether
# each has changed sign since the first sample (`flipped`, and the sample it did at, `flip_slot`) and where its
# magnitude is least (`best_slot`: a least at an inner sample is a turn toward zero). Plain loops over the
# instances of an atom, whose cells are consecutive
function _track!(ex::StepExtrema, rt::AutomatonRuntime, out::AbstractVector, slot::Int)
    comps = rt.components
    noise = TURN_FILTER[] ? ROUNDING_NOISE : 0.0
    for (km, m) in enumerate(rt.models)
        signed, bases = m.signed, ex.atom_base[km]
        n = size(signed, 1)
        for k in comps.live[km]
            base = bases[k]
            if slot == 1
                @inbounds for i in 1:n
                    s = signed[i, k]
                    ex.samples[base + i, 1] = s
                    ex.first_neg[base + i] = signbit(s)
                    ex.best_mag[base + i] = abs(s)
                    ex.best_slot[base + i] = 1
                    ex.flipped[base + i] = false
                end
            else
                @inbounds for i in 1:n
                    s = signed[i, k]
                    c = base + i
                    ex.samples[c, slot] = s
                    newly = (signbit(s) != ex.first_neg[c]) & !ex.flipped[c]
                    ex.flip_slot[c] = ifelse(newly, slot, ex.flip_slot[c])
                    ex.flipped[c] |= newly
                    less = abs(s) + noise * ex.best_mag[c] < ex.best_mag[c]
                    ex.best_mag[c] = ifelse(less, abs(s), ex.best_mag[c])
                    ex.best_slot[c] = ifelse(less, slot, ex.best_slot[c])
                end
            end
        end
    end
    for b in comps.crossings, i in 1:n_instances(rt.models[comps.model[b]])
        j = comps.offset[b] + i
        c, s = ex.guard_cell[j], out[j]
        ex.samples[c, slot] = s
        if slot == 1
            ex.first_neg[c], ex.best_mag[c], ex.best_slot[c], ex.flipped[c] = signbit(s), abs(s), 1, false
        else
            newly = (signbit(s) != ex.first_neg[c]) & !ex.flipped[c]
            ex.flip_slot[c] = ifelse(newly, slot, ex.flip_slot[c])
            ex.flipped[c] |= newly
            less = abs(s) + noise * ex.best_mag[c] < ex.best_mag[c]
            ex.best_mag[c] = ifelse(less, abs(s), ex.best_mag[c])
            ex.best_slot[c] = ifelse(less, slot, ex.best_slot[c])
        end
    end
    return nothing
end

# After the last sample: the cells that changed sign (`flips`) and those that did not and whose magnitude is least at
# an inner sample (`turns`, as (cell, the sample after the least: the window is the three samples around it))
function _derive_flags!(ex::StepExtrema, rt::AutomatonRuntime)
    comps = rt.components
    ex.rounding_age .+= 1
    check!(c) = begin
        if ex.flipped[c]
            push!(ex.flips, c)
        elseif 1 < ex.best_slot[c] < INTERP_POINTS && _turn_beyond_rounding(ex, rt, c)
            push!(ex.turns, (c, ex.best_slot[c] + 1))
        end
        nothing
    end
    for (km, m) in enumerate(rt.models), k in comps.live[km]
        base = ex.atom_base[km][k]
        for i in 1:n_instances(m)
            check!(base + i)
        end
    end
    for b in comps.crossings, i in 1:n_instances(rt.models[comps.model[b]])
        check!(ex.guard_cell[comps.offset[b] + i])
    end
    return nothing
end

# The steps after which the running error bounds of a model's values are evaluated again (and at every round of an instant): the
# scale of a value's terms moves with the dynamics, over many steps, and a stale scale moves a noise floor of ~1e-16 by a factor
const ROUNDING_REFRESH_STEPS = 64

# Whether the least magnitude of cell `c` at an inner sample is a turn: the magnitude at each end of the step exceeds it by more than
# the rounding of the atom's value, `ROUNDING_ULPS u e` with `e` the running error bound of the value (`RoundingBound`) at the
# state of the last evaluation of its model (the last sample): its terms' magnitudes, which a value that cancels them
# has far above its own. The bounds of a model are evaluated where a cell of it needs them and the last evaluation is `ROUNDING_REFRESH_STEPS` steps old or
# before the last round of an instant. A guard
# (a component of the callback, with no value of its own here) is a turn where its magnitude is least at an inner sample.
function _turn_beyond_rounding(ex::StepExtrema, rt::AutomatonRuntime, c::Int)
    km, k, i = ex.cell_atom[c]
    (km == 0 || !TURN_FILTER[]) && return true
    m = rt.models[km]
    if !ex.rounding_known[km] || ex.rounding_age[km] >= ROUNDING_REFRESH_STEPS || ex.rounding_rounds[km] != rt.rounds_taken[]
        columns = values(rt.bounds[km](ex.t))
        for v in eachindex(columns), j in 1:n_instances(m)
            ex.rounding[km][v, j] = columns[v][j].e
        end
        ex.rounding_known[km] = true
        ex.rounding_age[km] = 0
        ex.rounding_rounds[km] = rt.rounds_taken[]
    end
    noise = ROUNDING_NOISE * ex.rounding[km][m.atom_value[k], i]
    best = ex.best_mag[c]
    return abs(ex.samples[c, 1]) > best + noise && abs(ex.samples[c, INTERP_POINTS]) > best + noise
end

# The index of the sample point `t` of the step, 0 where it is none
function _sample_slot(ts::Vector{Float64}, t::Float64)
    for k in eachindex(ts)
        ts[k] == t && return k
    end
    return 0
end

# The components of `ex.folds` folded: as they are up to their time `t`, their value there after it
function _fold!(out::AbstractVector, ex::StepExtrema, t::Float64)
    for f in ex.folds
        t > f.t && (out[f.component] = f.s)
    end
    return out
end

# The displacement of value `v` of instance `i` of `m` from where the last instant left it, which its
# atoms and guards read (`anchored_signed_truth`), for the value `g`
@inline _displacement(m::AutomatonModel, v::Int, i::Int, g::Float64) = g - m.anchor[i, v]

# Whether the displacement of value `v` of instance `i` of `m` is read at all: not where the value is
# decided as zero or pinned (`_signed_truths!`, `_guard_component`)
@inline _displaced(m::AutomatonModel, v::Int, i::Int) = !(m.zero[i, v] | m.pinned[i, v])

# The state at time `t` of the integrator's current step and its rate, from its dense output
function _step_state!(ex::StepExtrema, integrator, t::Float64)
    integrator(ex.u, t)
    integrator(ex.du, t, Val{1})
    return nothing
end

# The values and rates of model `km` at time `t` within the step (`ex.V[km]`, `ex.R[km]`)
function _step_rates!(ex::StepExtrema, rt::AutomatonRuntime, r::RateDAE, integrator, km::Int, t::Float64)
    _step_state!(ex, integrator, t)
    _value_rates!(ex.V[km], ex.R[km], r, rt.models[km], ex.u, ex.du, t)
    return nothing
end

# The displacement of value `v` of instance `i` of model `km` at time `t` within the step, from its dense output
function _displacement_at!(ex::StepExtrema, rt::AutomatonRuntime, integrator, km::Int, v::Int, i::Int, t::Float64)
    nd = rt.n_diff
    integrator(ex.u, t)
    rt.evaluators[km](view(ex.u, 1:nd), view(ex.u, (nd + 1):rt.n), t)
    m = rt.models[km]
    return _displacement(m, v, i, m.value_columns[v][i])
end

# The time (the first point beyond it) where the displacement of a value, `d_a` at `a`, changes sign
# between `a` and `b`, where it has the other sign: its zero along the dense output, by bisection to
# adjacent floats
function _crossing_time(ex::StepExtrema, rt::AutomatonRuntime, integrator, km::Int, v::Int, i::Int,
                        a::Float64, b::Float64, d_a::Float64)
    lo, hi = a, b
    for _ in 1:2000
        mid = lo + (hi - lo) / 2
        (mid <= lo || mid >= hi) && break
        d = _displacement_at!(ex, rt, integrator, km, v, i, mid)
        signbit(d) == signbit(d_a) ? (lo = mid) : (hi = mid)
    end
    return hi
end

# The values and rates of model `km` at the ends of the window (once per window)
function _end_rates!(ex::StepExtrema, rt::AutomatonRuntime, host, jump, integrator, km::Int)
    ex.ends[km] && return nothing
    ex.ends[km] = true
    ex.counts.dual_windows += 1
    r = _rate_dae(jump, rt, host, integrator)
    integrator(ex.ua, ex.wa)
    integrator(ex.ub, ex.wb)
    integrator(ex.du0, ex.wa, Val{1})
    integrator(ex.du1, ex.wb, Val{1})
    _value_rates!(ex.V0[km], ex.R0[km], r, rt.models[km], ex.ua, ex.du0, ex.wa)
    _value_rates!(ex.V1[km], ex.R1[km], r, rt.models[km], ex.ub, ex.du1, ex.wb)
    return nothing
end

# The degree in t of the dense output of the integrator's current step: for QNDF, p(Θ) = y₁ + Σ φ_j(Θ - 1) k[j] with
# φ_j of degree j and k[j] zero beyond the order of the step, so the index of the last k[j] that is not zero (the
# truncated step of a callback is the same polynomial rebased, `stiff_addsteps.jl`); `NOT_POLYNOMIAL` for another
# method, whose dense output this does not know
function _dense_degree(integrator)
    integrator.alg isa DiffEq.QNDF || return NOT_POLYNOMIAL
    k = integrator.k
    j = findlast(q -> any(!iszero, k[q]), eachindex(k))
    return j === nothing ? 0 : j
end

# The polynomial of degree INTERP_POINTS - 1 through the samples of a cell at the sample points of the step, in Newton
# form (its divided differences: nested multiplications, with no division per evaluation). An atom or guard that is a
# polynomial of degree d in the state and the time follows the dense output of the step, a polynomial of degree at
# most d k in t, k that of the dense output: this polynomial (the higher differences zero) where d k is below the
# number of samples
struct SamplePolynomial
    c::NTuple{INTERP_POINTS, Float64}
    x::NTuple{INTERP_POINTS, Float64}
end

function SamplePolynomial(ex::StepExtrema, cell::Int)
    work = [ex.samples[cell, j] for j in 1:INTERP_POINTS]
    coefficients = Vector{Float64}(undef, INTERP_POINTS)
    coefficients[1] = work[1]
    for order in 1:(INTERP_POINTS - 1)
        for j in 1:(INTERP_POINTS - order)
            work[j] = (work[j + 1] - work[j]) / (ex.ts[j + order] - ex.ts[j])
        end
        coefficients[order + 1] = work[1]
    end
    return SamplePolynomial(ntuple(j -> coefficients[j], INTERP_POINTS), ntuple(j -> ex.ts[j], INTERP_POINTS))
end

# The value of the polynomial at `x`
@inline function (p::SamplePolynomial)(x::Float64)
    r = p.c[INTERP_POINTS]
    @inbounds for j in (INTERP_POINTS - 1):-1:1
        r = p.c[j] + (x - p.x[j]) * r
    end
    return r
end

# The zero of the polynomial `p` between `a` and `b`, where it has different signs, to adjacent floats
function _sampled_root(p::SamplePolynomial, a::Float64, b::Float64)
    lo, hi = a, b
    s_lo = signbit(p(lo))
    for _ in 1:2000
        mid = lo + (hi - lo) / 2
        (mid <= lo || mid >= hi) && break
        signbit(p(mid)) == s_lo ? (lo = mid) : (hi = mid)
    end
    return hi
end

# The point of [a, b] where `σ` times the polynomial `p` is least (golden section: the polynomial of a sub-interval
# between samples has at most one turn)
function _sampled_least(p::SamplePolynomial, σ::Float64, a::Float64, b::Float64)
    ϕ = (sqrt(5.0) - 1) / 2
    x1, x2 = b - ϕ * (b - a), a + ϕ * (b - a)
    f1, f2 = σ * p(x1), σ * p(x2)
    for _ in 1:40
        if f1 < f2
            b, x2, f2 = x2, x1, f1
            x1 = b - ϕ * (b - a)
            f1 = σ * p(x1)
        else
            a, x1, f1 = x1, x2, f2
            x2 = a + ϕ * (b - a)
            f2 = σ * p(x2)
        end
    end
    return (a + b) / 2
end

# The times where the sampled polynomial of cell `c` crosses zero between the samples `lo` and `hi`: one in each
# interval between samples where it changes sign, and the two of an excursion beyond zero (its least point toward
# zero is beyond it, strictly) in each where it does not
function _sampled_crossings(ex::StepExtrema, c::Int, lo::Int, hi::Int)
    times = Float64[]
    p = SamplePolynomial(ex, c)
    for q in lo:(hi - 1)
        a, b = ex.ts[q], ex.ts[q + 1]
        fa, fb = ex.samples[c, q], ex.samples[c, q + 1]
        if signbit(fa) != signbit(fb)
            push!(times, _sampled_root(p, a, b))
        else
            # (a turn toward zero inside the interval shows as a sampled magnitude at one end that is not above its
            # neighbours; at an end of the step an end has no neighbour on one side)
            mag(r) = abs(ex.samples[c, r])
            turn_at_p = q == 1 || (mag(q) <= mag(q - 1) && mag(q) <= mag(q + 1))
            turn_at_q = q + 1 == INTERP_POINTS || (mag(q + 1) <= mag(q) && mag(q + 1) <= mag(q + 2))
            (turn_at_p || turn_at_q) || continue
            σ = fa > 0 ? 1.0 : -1.0
            x = _sampled_least(p, σ, a, b)
            if σ * p(x) < 0
                push!(times, _sampled_root(p, a, x), _sampled_root(p, x, b))
            end
        end
    end
    return times
end

"""
    _crossing_count!(ex, rt, host, jump, integrator, km, v, i, hint) -> Int

How many times value `v` of instance `i` of model `km` crosses its surface (the zero of its
displacement) inside the window, from what the samples showed of it (`hint`): one where it changed sign between
the samples around the window (`:flipped`; no rates needed); where it did not, none, or the two of an excursion
beyond the surface where it turned toward zero (`:turned`) or the window is at an end of the step (`:blind`), its
rates at the ends of the window point opposite ways, the first toward the surface, and its extremum, located on the
dense output, is beyond it (see the header of step_extremum.jl); none for a value that showed nothing (`:still`). A value that is a polynomial in the state and the time of degree d with d k below the number of samples (k the degree of the dense
output: `cell`, its samples) has its crossings from the samples alone
(`_sampled_crossings`). Kept for the window (`ex.crossings`, with the displacements at the ends and the extremum for `_crossing_times!`).
"""
function _crossing_count!(ex::StepExtrema, rt::AutomatonRuntime, host, jump, integrator, km::Int, v::Int, i::Int,
                          hint::Symbol, cell::Int)
    n = ex.crossings[km][v, i]
    n >= 0 && return Int(n)
    m = rt.models[km]
    n = 0
    if _displaced(m, v, i) && (cell > 0 || hint !== :still)
        if cell > 0
            # (a polynomial the samples determine, with no evaluation of the values)
            times = _sampled_crossings(ex, cell, ex.wlo, ex.whi)
            ex.times[(km, v, i)] = times
            n = length(times)
        elseif hint === :flipped
            d0 = _displacement_at!(ex, rt, integrator, km, v, i, ex.wa)
            d1 = _displacement_at!(ex, rt, integrator, km, v, i, ex.wb)
            ex.start[km][v, i], ex.finish[km][v, i] = d0, d1
            n = !(iszero(d0) || iszero(d1)) && signbit(d0) != signbit(d1) ? 1 : 0
        else
            _end_rates!(ex, rt, host, jump, integrator, km)
            d0, d1 = _displacement(m, v, i, ex.V0[km][v, i]), _displacement(m, v, i, ex.V1[km][v, i])
            ex.start[km][v, i], ex.finish[km][v, i] = d0, d1
            if iszero(d0) || iszero(d1)
                # (a value at its surface at an end: DiffEq reports no crossing from an exact zero)
            elseif signbit(d0) != signbit(d1)
                n = 1
            else
                n = _excursion!(ex, rt, host, jump, integrator, km, v, i, d0, d1) ? 2 : 0
            end
        end
    end
    ex.crossings[km][v, i] = n
    return n
end

# Whether value `v` of instance `i` of model `km` makes an excursion beyond its surface inside the window:
# its displacement `d0` at the start and `d1` at the end have one sign, its rates at the ends point
# opposite ways, the first toward zero, and its extremum (kept, with the displacement there, in
# `ex.peak` and `ex.peak_d`) is beyond it
function _excursion!(ex::StepExtrema, rt::AutomatonRuntime, host, jump, integrator, km::Int, v::Int, i::Int,
                     d0::Float64, d1::Float64)
    # w = σ d is below zero at the start, and its maximum is the extremum: heading up at the start, down by
    # the end
    σ = d0 > 0 ? -1.0 : 1.0
    ġ0, ġ1 = ex.R0[km][v, i], ex.R1[km][v, i]
    (σ * ġ0 > 0 && σ * ġ1 <= 0) || return false
    r = _rate_dae(jump, rt, host, integrator)
    peak = _locate_peak(σ * d0, σ * ġ0, σ * d1, σ * ġ1, ex.wa, ex.wb) do t_H
        _step_rates!(ex, rt, r, integrator, km, t_H)
        return σ * ex.R[km][v, i]
    end
    peak === nothing && return false
    d_x = _displacement_at!(ex, rt, integrator, km, v, i, peak.t)
    _beyond_surface(σ, d_x) || return false
    ex.peak[km][v, i], ex.peak_d[km][v, i] = peak.t, d_x
    return true
end

# The times inside the window at which value `v` of instance `i` of model `km` crosses its surface
# (`ex.crossings` of them, found once for the window, `ex.times`): the zeros of its displacement along
# the dense output on each side of its extremum, or between the ends
function _crossing_times!(ex::StepExtrema, rt::AutomatonRuntime, host, jump, integrator, km::Int, v::Int, i::Int)
    n = Int(ex.crossings[km][v, i])      # (counted: `_crossing_count!`)
    n <= 0 && return NO_CROSSINGS
    return get!(ex.times, (km, v, i)) do
        t0, t1, d0 = ex.wa, ex.wb, ex.start[km][v, i]
        n == 1 ? [_crossing_time(ex, rt, integrator, km, v, i, t0, t1, d0)] :
                 [_crossing_time(ex, rt, integrator, km, v, i, t0, ex.peak[km][v, i], d0),
                  _crossing_time(ex, rt, integrator, km, v, i, ex.peak[km][v, i], t1, ex.peak_d[km][v, i])]
    end
end

"""
    _scan_step!(ex, rt, evaluators, host, jump, integrator)

Scan the integrator's current step, whose atoms and guards `_track!` followed over its sample points
`ex.ts`, for the components that change their truths inside it though not at any sample (see the header of
step_extremum.jl), and keep their folds in `ex.folds`. A component that can fire, with one sign at every
sample, is examined (`_examine_window!`) where

- one of its atoms changed sign between samples (the condition did not: another atom held it), or
- one of its atoms (or its guard) kept one sign at every sample and its sampled magnitude turned toward zero,

in the window between the samples around those changes of its atoms (the two at a change of sign, the three at a
turn; the whole span they cover where there are several).
"""
function _scan_step!(ex::StepExtrema, rt::AutomatonRuntime, evaluators::Tuple, host, jump, integrator)
    comps = rt.components
    n = length(integrator.u)
    length(ex.u) == n || foreach(v -> resize!(v, n), (ex.u, ex.ua, ex.ub, ex.du, ex.du0, ex.du1))
    ex.counts.steps += 1
    ex.kdeg = _dense_degree(integrator)
    # The components that read an atom or guard that changed sign or turned
    empty!(ex.queue)
    for c in ex.flips
        _queue_components!(ex, rt, c)
    end
    for (c, _) in ex.turns
        _queue_components!(ex, rt, c)
    end
    examined = false
    for j in ex.queue
        ex.queued[j] = false
        iszero(ex.first_comp[j]) && continue
        b = searchsortedlast(comps.offset, j - 1)
        i = j - comps.offset[b]
        c = comps.condition[b]
        km = comps.model[b]
        # The samples around where its atoms changed sign (the two) and turned (the three), at the most
        lo, hi = INTERP_POINTS, 1
        for k in (c === nothing ? (0,) : c.atoms)
            cl = c === nothing ? ex.guard_cell[j] : ex.atom_cell[km][i, k]
            ex.flipped[cl] && (lo = min(lo, ex.flip_slot[cl] - 1); hi = max(hi, ex.flip_slot[cl]))
            for (tc, slot) in ex.turns
                tc == cl && !ex.flipped[cl] && (lo = min(lo, slot - 2); hi = max(hi, slot))
            end
        end
        # (an atom that did not move cannot have turned between two samples, but can inside the first or the
        # last interval of the step, where the samples before or after it are missing: its turn shows as no turn)
        blind = lo == 1 || hi == INTERP_POINTS
        (blind || _may_change(ex, c, km, i, j)) || continue
        examined = true
        ex.counts.candidates += 1
        _examine_window!(ex, rt, evaluators, host, jump, integrator, b, i, lo, hi, blind)
    end
    ex.counts.candidate_steps += examined
    return nothing
end

# Queue the components that can fire and read cell `c`
function _queue_components!(ex::StepExtrema, rt::AutomatonRuntime, c::Int)
    km, k, i = ex.cell_atom[c]
    comps = rt.components
    if km == 0
        queue!(ex, i)       # (a guard: its component)
    else
        for b in ex.atom_blocks[km][k]
            (comps.watch[b][i] & rt.models[km].active[i]) && queue!(ex, comps.offset[b] + i)
        end
    end
    return nothing
end

function queue!(ex::StepExtrema, j::Int)
    ex.queued[j] || (ex.queued[j] = true; push!(ex.queue, j))
    return nothing
end

# Whether component `j` (instance `i` of model `km`, condition `c`, `nothing` for a guard), of one sign at every
# sample, can have the other truth inside the step: with each of its atoms that changed sign or turned toward zero
# at its truth at the start or the other, and the rest at theirs at the start, some combination makes it so (a
# guard that turned changes its own truth)
function _may_change(ex::StepExtrema, c, km::Int, i::Int, j::Int)
    c === nothing && return true
    cells = ex.atom_cell[km]
    moving = 0                  # the atoms that moved, as bits of a mask over the atoms of the condition
    for (q, k) in enumerate(c.atoms)
        cl = cells[i, k]
        (ex.flipped[cl] || any(t -> t[1] == cl, ex.turns)) && (moving |= 1 << (q - 1))
    end
    return _flip_possible(ex, c, km, i, ex.first_comp[j] > 0, moving)
end

# Whether the condition `c` of instance `i` of model `km`, true at the start of the step as `start`, has the other
# truth with some of the atoms in `moving` (a mask over its atoms) at the other truth than at the start and the rest at
# it
function _flip_possible(ex::StepExtrema, c, km::Int, i::Int, start::Bool, moving::Int)
    cells = ex.atom_cell[km]
    mask = moving
    while true      # (the submasks of `moving`, each atom in it taking the other truth)
        sub = mask
        holds = condition_holds(c, k -> xor(!ex.first_neg[cells[i, k]], isodd(sub >> (findfirst(==(k), c.atoms) - 1))))
        holds != start && return true
        mask == 0 && return false
        mask = (mask - 1) & moving
    end
end

# Component `i` of block `b`, whose truth is the same at every sample, between `wa` and `wb`, around a turn of
# the sampled values of one of its atoms or where one of its atoms changed sign: whether its condition changes its truth inside the window, from
# where its atoms' values cross their surfaces (see the header of step_extremum.jl); its fold, kept in
# `ex.folds`
function _examine_window!(ex::StepExtrema, rt::AutomatonRuntime, evaluators::Tuple, host, jump, integrator, b::Int,
                          i::Int, lo::Int, hi::Int, blind::Bool)
    comps = rt.components
    km = comps.model[b]
    m = rt.models[km]
    c = comps.condition[b]
    j = comps.offset[b] + i
    # (what was found for a window stands for the next component examined over the same one)
    wa, wb = ex.ts[lo], ex.ts[hi]
    if (ex.wa, ex.wb) != (wa, wb)
        ex.wa, ex.wb, ex.wlo, ex.whi = wa, wb, lo, hi
        empty!(ex.times)
        fill!(ex.ends, false)
        foreach(c -> fill!(c, Int8(-1)), ex.crossings)
    end
    start = ex.first_comp[j] > 0
    # (the truths of the atoms at the start of the window are those at the first sample: an atom that changed sign
    # did so after it)
    atom_true(k) = !ex.first_neg[ex.atom_cell[km][i, k]]
    # The values the condition reads (its guard, or its atoms') and how often each crosses inside the window, from
    # what the samples showed of it
    guard = c === nothing ? (comps.transition[b]::CompiledTransition).guard_value : 0
    function hint(v)
        c === nothing && return any(t -> t[1] == ex.guard_cell[j], ex.turns) ? :turned : blind ? :blind : :still
        h = :still
        for k in c.atoms
            m.atom_value[k] == v || continue
            cl = ex.atom_cell[km][i, k]
            ex.flipped[cl] && return :flipped
            any(t -> t[1] == cl, ex.turns) ? (h = :turned) : (h === :turned || !blind || (h = :blind))
        end
        return h
    end
    function sampled_cell(v)       # the cell of an atom over the value `v` whose polynomial the samples determine, 0 if there is none
        deg = m.reads.degree[v]
        (deg < NOT_POLYNOMIAL && deg * ex.kdeg <= INTERP_POINTS - 1) || return 0
        c === nothing && return ex.guard_cell[j]
        for k in c.atoms
            m.atom_value[k] == v && return ex.atom_cell[km][i, k]
        end
        return 0
    end
    count(v) = _crossing_count!(ex, rt, host, jump, integrator, km, v, i, hint(v), sampled_cell(v))
    if c === nothing
        count(guard) == 0 && return false
    else
        # The condition can differ from its truth at the start only where, with each atom whose value
        # crosses at its truth at the start or the other, some combination makes it so
        moving = 0
        for (q, k) in enumerate(c.atoms)
            count(m.atom_value[k]) > 0 && (moving |= 1 << (q - 1))
        end
        (moving != 0 && _flip_possible(ex, c, km, i, start, moving)) || return false
    end
    # The truth of the condition after each crossing time, each atom's truth that at the start as often
    # flipped as its value has crossed
    times_of(v) = _crossing_times!(ex, rt, host, jump, integrator, km, v, i)
    times = sort!(unique!(c === nothing ? copy(times_of(guard)) :
                          reduce(vcat, [times_of(m.atom_value[k]) for k in c.atoms])))
    flipped(v, t) = isodd(Base.count(<=(t), times_of(v)))
    for (p, t_p) in enumerate(times)
        holds = c === nothing ? xor(start, flipped(guard, t_p)) :
                condition_holds(c, k -> xor(atom_true(k), flipped(m.atom_value[k], t_p)))
        holds == start && continue
        # It differs from the start from t_p to the next crossing (or the end of the window): a fold inside
        t_m = (t_p + (p < length(times) ? times[p + 1] : wb)) / 2
        _step_state!(ex, integrator, t_m)
        automaton_conditions!(ex.s, rt, evaluators, ex.u, t_m)
        # (a component at t_m with the start's sign is no change the integration can resolve: the interval is
        # at the rounding of the crossing times)
        signbit(ex.s[j]) == signbit(ex.first_comp[j]) && return false
        push!(ex.folds, ExtremumFold(j, t_m, ex.s[j]))
        ex.counts.located += 1
        return true
    end
    return false
end
