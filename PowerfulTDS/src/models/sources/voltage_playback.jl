"""VoltagePlayback: a bus driven by a recorded voltage (magnitude and angle as functions of time)

Hand-written: the recording is per-instance data (a time series), which the model
compiler does not describe.
"""

using InlineStrings
using OrderedCollections: OrderedDict
using ForwardDiff: ForwardDiff

export VoltagePlayback, PlaybackRecording, playback_voltage

# ============================================================================
# The recording and its interpolant
# ============================================================================

"""
    PlaybackRecording

The recorded voltage of a [`VoltagePlayback`](@ref): samples `t` (s), `v` (pu) and
`theta` (rad), cut at the `breaks` into pieces, each interpolated on its own by a
not-a-knot cubic spline (C²), held as the slopes `dv`, `dtheta` of its cubic Hermite
segments at the samples. Piece `k` holds the samples `first[k]:last[k]`; outside
its samples (up to the breaks around it) a piece continues its end segment. Read-only
once built.
"""
struct PlaybackRecording
    t::Vector{Float64}
    v::Vector{Float64}
    theta::Vector{Float64}
    dv::Vector{Float64}
    dtheta::Vector{Float64}
    first::Vector{Int}
    last::Vector{Int}
    breaks::Vector{Float64}
end

# The slopes at the samples of the not-a-knot cubic spline through y(x) (C²; the third
# derivative also continuous at the second and the last-but-one samples), in the Hermite
# slope form of SciPy's CubicSpline(bc_type = "not-a-knot"): three samples give the
# parabola through them, two the line
function _spline_slopes(x::AbstractVector{Float64}, y::AbstractVector{Float64})
    n = length(x)
    h = diff(x)
    δ = diff(y) ./ h
    n == 2 && return [δ[1], δ[1]]
    if n == 3
        A = [1.0 1.0 0.0; h[2] 2(h[1] + h[2]) h[1]; 0.0 1.0 1.0]
        return A \ [2δ[1], 3(h[2] * δ[1] + h[1] * δ[2]), 2δ[2]]
    end
    dl, dd, du, b = zeros(n - 1), zeros(n), zeros(n - 1), zeros(n)
    for k in 2:(n - 1)   # the second derivative continuous at sample k
        dl[k - 1] = h[k]
        dd[k] = 2(h[k - 1] + h[k])
        du[k] = h[k - 1]
        b[k] = 3(h[k] * δ[k - 1] + h[k - 1] * δ[k])
    end
    d0 = x[3] - x[1]
    dd[1], du[1] = h[2], d0
    b[1] = ((h[1] + 2d0) * h[2] * δ[1] + h[1]^2 * δ[2]) / d0
    d1 = x[n] - x[n - 2]
    dl[n - 1], dd[n] = d1, h[n - 2]
    b[n] = (h[n - 1]^2 * δ[n - 2] + (2d1 + h[n - 1]) * h[n - 2] * δ[n - 1]) / d1
    return Tridiagonal(dl, dd, du) \ b
end

"""
    PlaybackRecording(t, v, theta, breaks; name) -> PlaybackRecording

Check a recording and build its interpolant. `name` names the instance in the errors.

`t` is strictly increasing, except that a break time may appear twice in a row: its
first sample is the value just before the break (left limit), its second the value just
after (right limit). A break time that is not a sample time cuts between the samples
around it, and each side continues its end segment to the break. Errors (naming the
instance) for: series of unequal lengths, times that do not increase, a break time
sampled once (which side it belongs to is not known), breaks that do not increase or are
not inside the recording, a piece with fewer than two samples, and non-finite samples.
"""
function PlaybackRecording(t::AbstractVector{<:Real}, v::AbstractVector{<:Real},
                           theta::AbstractVector{<:Real}, breaks::AbstractVector{<:Real}; name::AbstractString)
    length(t) == length(v) == length(theta) || throw(ArgumentError(
        "$name: the recording has $(length(t)) times, $(length(v)) magnitudes v and $(length(theta)) " *
        "angles theta; they must have the same length"))
    length(t) >= 2 || throw(ArgumentError("$name: the recording has $(length(t)) samples; it needs at least two"))
    for (series, label) in ((t, "time"), (v, "magnitude v"), (theta, "angle theta"), (breaks, "break"))
        k = findfirst(!isfinite, series)
        isnothing(k) || throw(ArgumentError("$name: $label sample $k is $(series[k]); samples must be finite"))
    end
    tb = Float64.(breaks)
    for j in 2:length(tb)
        tb[j] > tb[j - 1] || throw(ArgumentError(
            "$name: breaks must be strictly increasing; break $j at $(tb[j]) s follows $(tb[j - 1]) s"))
    end
    for b in tb
        t[1] < b < t[end] || throw(ArgumentError(
            "$name: the break at $b s is not inside the recording ($(t[1]) s to $(t[end]) s)"))
    end

    # Pieces: the samples between consecutive breaks (a break sampled twice gives its first
    # sample to the piece before it and its second to the piece after it)
    n, nb = length(t), length(tb)
    first_k, last_k = Int[1], Int[]
    j = 1   # the next break
    for k in 2:n
        if t[k] == t[k - 1]
            (j <= nb && t[k] == tb[j]) || throw(ArgumentError(
                "$name: time sample $k repeats $(t[k]) s, which is not a break; times must be strictly " *
                "increasing (a break time may appear twice: the values before and after it)"))
            push!(last_k, k - 1); push!(first_k, k); j += 1
            continue
        end
        t[k] > t[k - 1] || throw(ArgumentError(
            "$name: times must be strictly increasing; sample $k at $(t[k]) s follows $(t[k - 1]) s"))
        (j <= nb && t[k] >= tb[j]) || continue
        if t[k] == tb[j]
            # a break at sample k: given twice, it is cut at the repeat
            (k < n && t[k + 1] == t[k]) || throw(ArgumentError(
                "$name: the break at $(tb[j]) s is sampled once (sample $k); give it twice, the values " *
                "before and after it, or no sample at it"))
            continue
        end
        # a break between samples k - 1 and k
        push!(last_k, k - 1); push!(first_k, k); j += 1
        j <= nb && tb[j] <= t[k] && throw(ArgumentError(
            "$name: no sample between the breaks at $(tb[j - 1]) s and $(tb[j]) s"))
    end
    push!(last_k, n)
    for p in eachindex(first_k)
        last_k[p] > first_k[p] || throw(ArgumentError(
            "$name: the piece of the recording from $(t[first_k[p]]) s has one sample; each piece " *
            "between breaks needs at least two"))
    end

    tt, vv, θθ = Float64.(t), Float64.(v), Float64.(theta)
    dv, dθ = zeros(n), zeros(n)
    for p in eachindex(first_k)
        r = first_k[p]:last_k[p]
        dv[r] = _spline_slopes(view(tt, r), view(vv, r))
        dθ[r] = _spline_slopes(view(tt, r), view(θθ, r))
    end
    return PlaybackRecording(tt, vv, θθ, dv, dθ, first_k, last_k, tb)
end

"""
    playback_piece(rec::PlaybackRecording, t::Real) -> Int

The piece of `rec` in which the playback is at time `t`: after every break at or before
`t` (at a break, the value after it).
"""
playback_piece(rec::PlaybackRecording, t::Real) = 1 + count(<=(t), rec.breaks)

# The float value of a time that may carry derivatives (nested Duals)
_primal(t::Real) = t
_primal(t::ForwardDiff.Dual) = _primal(ForwardDiff.value(t))

"""
    playback_segment(rec::PlaybackRecording, piece::Integer, s::Float64) -> Int

The segment `j` of piece `piece` (between samples `j` and `j + 1`) that holds the time
`s`, by binary search: before the first sample of the piece the first segment, after its
last sample the last segment. Pure (no state is kept between evaluations, which may go
back in time or run in parallel).
"""
@inline function playback_segment(rec::PlaybackRecording, piece::Integer, s::Float64)
    lo = rec.first[piece]
    hi = rec.last[piece] - 1
    t = rec.t
    # the last j in lo:hi with t[j] <= s (lo when s is before the piece)
    @inbounds while lo < hi
        mid = (lo + hi + 1) >>> 1
        if t[mid] <= s
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

# The cubic Hermite interpolant on the segment j of y (slopes d) at time t (Horner form of
# y0 h00 + h d0 h10 + y1 h01 + h d1 h11, s = (t - t_j)/h; SymPy-checked)
@inline function _hermite(tg::Vector{Float64}, y::Vector{Float64}, d::Vector{Float64}, j::Int, t::Real)
    @inbounds begin
        h = tg[j + 1] - tg[j]
        s = (t - tg[j]) / h
        Δ = y[j + 1] - y[j]
        c1 = h * d[j]
        c2 = 3Δ - h * (2d[j] + d[j + 1])
        c3 = h * (d[j] + d[j + 1]) - 2Δ
        return y[j] + s * (c1 + s * (c2 + s * c3))
    end
end

"""
    playback_voltage(rec::PlaybackRecording, piece::Integer, t::Real) -> (v, theta)

The recorded voltage magnitude `v` (pu) and angle `theta` (rad) at time `t`, within the
segment of the recording bounded by break `piece - 1` and break `piece` (`piece = 1`
before the first break, `playback_piece` finds the piece a given time falls in). Works
with a differentiable `t` (e.g. a Dual number), so it can be used inside sensitivity or
Jacobian computations.
"""
function playback_voltage(rec::PlaybackRecording, piece::Integer, t::Real)
    j = playback_segment(rec, piece, Float64(_primal(t)))
    return (_hermite(rec.t, rec.v, rec.dv, j, t), _hermite(rec.t, rec.theta, rec.dtheta, j, t))
end

# ============================================================================
# The model
# ============================================================================

struct VoltagePlayback{T} <: AbstractModel
    "Bus number the recorded voltage is played back at"
    BUS::Int32
    "Identifier of the power flow Slack unit it replaces"
    ID::String3
    "The recorded voltage and its interpolant"
    rec::PlaybackRecording
    "The piece of the recording the playback is on (set at the start, advanced by the break events)"
    piece::Int32
    "Model status: true=online, false=offline"
    status::Bool
end

"""
    VoltagePlayback(bus; id="1", t, v, theta, breaks=Float64[], status=true) -> VoltagePlayback{Float64}

A playback of the voltage recorded at bus `bus`: magnitude `v` (pu) and angle `theta`
(rad) at the times `t` (s), discontinuous at the times `breaks` (see
[`PlaybackRecording`](@ref) for the data and its checks). It replaces the power flow
Slack unit `id` at the bus.
"""
function VoltagePlayback(bus::Integer; id::AbstractString = "1", t::AbstractVector{<:Real},
                         v::AbstractVector{<:Real}, theta::AbstractVector{<:Real},
                         breaks::AbstractVector{<:Real} = Float64[], status::Bool = true)
    rec = PlaybackRecording(t, v, theta, breaks; name = _playback_name(bus, id))
    return VoltagePlayback{Float64}(Int32(bus), String3(id), rec, Int32(1), status)
end

_playback_name(bus, id) = "VoltagePlayback at bus $bus (id \"$id\")"

"""
    playback_voltage(pb::VoltagePlayback, t::Real) -> (v, theta)

The voltage `pb` plays back at time `t` (after a break at `t`): at the start of a
simulation, the voltage the power flow's slack must have at the bus.
"""
playback_voltage(pb::VoltagePlayback, t::Real) = playback_voltage(pb.rec, playback_piece(pb.rec, t), t)

@register_model VoltagePlayback

requires_initialization(::Type{<:VoltagePlayback}) = true

get_key_name(::Type{<:VoltagePlayback}) = [:BUS, :ID]

# It holds the bus in place of the power flow's slack unit
get_replacement_map(::Type{<:VoltagePlayback}) = OrderedDict(
    :Slack => (source_fields = (:BUS, :ID), target_fields = (:i, :id)),
)

const VOLTAGE_PLAYBACK_TERMINALS = (Terminal(:bus_group, :BUS),)
get_terminals(::Type{<:VoltagePlayback}) = VOLTAGE_PLAYBACK_TERMINALS

# The private residuals are paired with the model's own variables in declaration order
# (g_v with p, g_theta with q), which an out-of-service instance holds
function model_metadata(::Type{VoltagePlayback})
    ModelMetadata(
        name = :VoltagePlayback,
        vars = [
            ModelVar(:p, AlgebVar(); description = "Active power injected into the bus (pu, system base)", units = "pu"),
            ModelVar(:q, AlgebVar(); description = "Reactive power injected into the bus (pu, system base)", units = "pu"),
            ModelVar(:v, AlgebVar(), :bus_group, :v, :BUS; description = "Bus voltage magnitude (external)", units = "pu"),
            ModelVar(:theta, AlgebVar(), :bus_group, :theta, :BUS; description = "Bus voltage angle (external)", units = "rad"),
        ],
        residuals = [
            ModelResidual(:g_v, AlgebRes(), PrivateRes(); description = "v - V(t): the bus magnitude is the recorded one"),
            ModelResidual(:g_theta, AlgebRes(), PrivateRes(); description = "theta - Θ(t): the bus angle is the recorded one"),
            ModelResidual(:p_inj, AlgebRes(), SharedRes(:add), :bus_group, :p, :BUS;
                description = "Active power injection into the bus"),
            ModelResidual(:q_inj, AlgebRes(), SharedRes(:add), :bus_group, :q, :BUS;
                description = "Reactive power injection into the bus"),
        ],
        observable_vars = Symbol[],
        init_inputs = [
            InitInput(:p_gen, :from_observation, :Generator, :p_gen;
                description = "Active power of the replaced Slack unit in the power flow"),
            InitInput(:q_gen, :from_observation, :Generator, :q_gen;
                description = "Reactive power of the replaced Slack unit in the power flow"),
        ],
        status_field = :status,
    )
end

# ============================================================================
# Residuals: CPU kernel, GPU refusal, data adapter
# ============================================================================

@inline function voltage_playback_residual_kernel!(
    outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
    active::AbstractVector{Bool}, t::Real, ::Type{<:Vector}
)
    s = Float64(_primal(t))
    @inbounds for i in eachindex(active)
        mask = active[i] ? one(eltype(outputs.g_v)) : zero(eltype(outputs.g_v))
        rec = params.rec[i]
        j = playback_segment(rec, params.piece[i], s)
        V = _hermite(rec.t, rec.v, rec.dv, j, t)
        Θ = _hermite(rec.t, rec.theta, rec.dtheta, j, t)
        outputs.g_v[i] = mask * (vars.v[i] - V)
        outputs.g_theta[i] = mask * (vars.theta[i] - Θ)
        outputs.p_inj[i] = mask * vars.p[i]
        outputs.q_inj[i] = mask * vars.q[i]
    end
    nothing
end

# The recordings are per-instance tables of different lengths: they are not device data,
# and no TDS model runs on the GPU
PowerfulCore.gpu_support(::Type{<:VoltagePlayback}) = PowerfulCore.UnsupportedOnGPU(
    "VoltagePlayback has no GPU kernel: its recordings are interpolated on the CPU only")

function voltage_playback_residual_kernel!(outputs::NamedTuple, params::NamedTuple, vars::NamedTuple,
                                           active::AbstractVector{Bool}, t::Real, ::Type{AT}) where AT
    throw(ArgumentError("VoltagePlayback has no GPU kernel (array type $AT): its recordings are " *
                        "interpolated on the CPU only"))
end

@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT}, model::StructArray{<:VoltagePlayback},
    external::NamedTuple, ::PowerInjection, t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    voltage_playback_residual_kernel!(
        (g_v = ws.res_vals.g_v, g_theta = ws.res_vals.g_theta, p_inj = ws.res_vals.p_inj, q_inj = ws.res_vals.q_inj),
        (rec = model.rec, piece = model.piece),
        (p = ws.var_vals.p, q = ws.var_vals.q, v = ws.var_vals.v, theta = ws.var_vals.theta),
        active_mask, t, AT)
end

# ============================================================================
# Initialization: the power of the replaced Slack unit
# ============================================================================

function initialize_dynamic_states(model::StructArray{VoltagePlayback{T}}, ctx::InitContext) where T
    n = length(model)
    sources =ctx.links.source_idx[ctx.links.source_type .== :VoltagePlayback]
    for k in 1:n
        k in sources && continue
        throw(ArgumentError("$(_playback_name(model.BUS[k], model.ID[k])) replaces no power flow Slack " *
                            "unit: bus $(model.BUS[k]) is not the slack bus of the power flow, or it has no unit " *
                            "\"$(model.ID[k])\""))
    end
    sources == 1:n || error("VoltagePlayback: links out of instance order ($sources)")

    # Every slack unit at the bus is replaced (a unit left in place would hold the bus too)
    if haskey(ctx.sys.models, :Slack)
        slack = ctx.sys.models.Slack
        for k in 1:n, u in eachindex(slack)
            slack.i[u] == model.BUS[k] || continue
            throw(ArgumentError("$(_playback_name(model.BUS[k], model.ID[k])): the Slack unit \"$(slack.id[u])\" " *
                                "at the bus is not replaced by a dynamic model, and would hold the bus voltage too"))
        end
    end

    data = gather_init_data(ctx, VoltagePlayback)
    return (p = T.(data.p_gen), q = T.(data.q_gen))
end

# ============================================================================
# The start and the breaks of a problem
# ============================================================================

"""
    expand_voltage_playback!(sys::SystemModel, events::Vector, tspan, pf::NamedTuple) -> Vector

Start each `VoltagePlayback` of `sys` at `tspan[1]` and add its breaks to `events`:

- refuses a `tspan` outside an instance's recording;
- refuses a start where the power flow's voltage at the bus (`pf`, the power flow snapshot
  of the initialization) differs from the recording at `tspan[1]` by more than the
  power flow's residual (`pf.residual_norm`, which bounds the slack's own `v - vs` and
  `θ - va`), naming both values;
- sets the piece of the recording at `tspan[1]` (after a break at `tspan[1]`);
- adds a `FieldModification` of the piece at each break inside the span: the step before
  it ends on the value before the break, and the reinitialization after it starts from
  the value after it.

Called by the `TDSProblem` constructor.
"""
function expand_voltage_playback!(sys::SystemModel, events::Vector, tspan, pf::NamedTuple)
    haskey(sys.models, :VoltagePlayback) || return events
    model = sys.models.VoltagePlayback
    isempty(model) && return events
    t0, t1 = Float64(tspan[1]), Float64(tspan[2])
    expanded = Any[e for e in events]
    positions = pf_snapshot_positions(pf)
    for k in eachindex(model)
        rec, bus, id = model.rec[k], model.BUS[k], String(model.ID[k])
        name = _playback_name(bus, id)
        rec.t[1] <= t0 && t1 <= rec.t[end] || throw(ArgumentError(
            "$name: the time span ($t0 s, $t1 s) is not inside the recording ($(rec.t[1]) s to $(rec.t[end]) s)"))
        piece = playback_piece(rec, t0)
        model.piece[k] = piece
        V0, Θ0 = playback_voltage(rec, piece, t0)
        b = positions[bus]
        tol = pf.residual_norm
        abs(pf.v[b] - V0) <= tol && abs(pf.theta[b] - Θ0) <= tol || throw(ArgumentError(
            "$name: the power flow's voltage at the bus, v = $(pf.v[b]) pu, theta = $(pf.theta[b]) rad, is not the " *
            "recording's at the start t = $t0 s, v = $V0 pu, theta = $Θ0 rad (power flow residual $tol); set the " *
            "slack's voltage (Slack vs, Bus va) to the recording's start, `playback_voltage`"))
        for (j, tb) in enumerate(rec.breaks)
            t0 < tb < t1 || continue
            push!(expanded, FieldModification(tb, :VoltagePlayback, (Int(bus), id), :piece, :set, j + 1))
        end
    end
    return expanded
end
