# ============================================================================
# Rates of the Values of Automata Along the Trajectory
# ============================================================================
#
# The rates of the values of the automata of a model along the trajectory of a DAE
# `M u̇ = F(u, t)`, in the discrete state the runtime holds: `ġ = g_u u̇ + g_t`, from one
# Dual evaluation of the model's `calc_automaton_values!` along the rate `u̇` of the state
# (`_value_rates`), which is `F/M` on the rows with mass and, on the rows without mass,
# the derivative of their equations along the trajectory, `F_ZZ u̇_Z = -(F_ZW u̇_W +
# ∂F_Z/∂t)` (`_state_rate`). Shared by the DAEs whose rates the automata need (`RateDAE`):
#
# - the plain simulation's (`StateRates`, on the DAE of its host, `dae_residuals!`), for
#   the sliding check of located transitions (automata.jl, `_sliding_candidates`);
# - an augmented system such as a parameter-sensitivity integration, for the switching
#   times of transitions and the grazing of guards.
#
# The rate of a value that reads no variable without mass is that of the rows with mass
# alone, `F/M` (`_mass_rate`), which needs no solve; the rows without mass are completed
# (`F_ZZ`, factorized by KLU) only for a value that reads one (`reads_massless`). A
# factorization is kept (`RateFactor`) for the evaluations at the same state and time in
# an unchanged `F_ZZ`. An `F_ZZ` that is numerically singular has no rates: it throws
# (`SingularRatesError`) rather than let a decision read garbage.
# ============================================================================

using ForwardDiff: ForwardDiff
using KLU: KLU, klu
using LinearAlgebra: SingularException, svd
using PowerfulCore: get_typed_workspace!

"""
    RateDAE

A DAE `M u̇ = F(u, t)` along whose trajectory the rates of the values of automata are
evaluated (see the header of automaton_rates.jl): `StateRates`, or an augmented system
such as a parameter-sensitivity integration. Each has the fields

- `n`, `n_diff`: the lengths of `u = [x; y]` and of `x`;
- `mass`, `zero_mass`: the diagonal of the mass matrix, and the rows (and variables)
  without mass;
- `zz_pattern`, `zz_colors`, `zz_position`: the pattern of `F_ZZ`, the Jacobian of the
  rows without mass with respect to the variables without mass, a coloring of its
  columns, and the position in `zero_mass` of each variable (0 for one with mass);
- `u_d`, `du_d`: a Dual state and residual, whose element type's partials seed the
  colors of `F_ZZ` (`_zero_mass_factor`) and the rate (its first partial);

Extension point: every subtype implements the methods `_rate_host` (the host whose
`dae_residuals!` is `F`), `_dual_workspace` (the workspace of the element type of `u_d`,
in which the values of the automata are evaluated with Duals) and `_plain_values!` (makes
the quantities that move in that workspace hold their values without derivatives, before a
rate is evaluated: where the host sets them in the Float64 workspace and writes the Dual
workspaces too, as the plain simulation does, it only asserts that, under
`CHECK_MEMORY_SYNC`). `StateRates` is the plain simulation's; an augmentation that
integrates more than the DAE state (e.g. its sensitivities) defines its own subtype with
these three methods, and is handed back through `_rate_dae`.
"""
abstract type RateDAE end

# The Dual tag of the rates of a plain simulation whose host names no Dual workspace of its
# integration's Jacobian (`rate_dual_type`)
struct AutomatonRateTag end

"""
    rate_dual_type(host, width) -> Type

The Dual type in which the rates of the automata of the DAE of `host` are evaluated
(`StateRates`), in the host's workspace of that type (`dual_workspace`). A `TDSAdapter`'s
is that of its integration's Jacobian (`JacobianClosureFunctor`'s tag), whose workspace
the integration builds at its first Jacobian, before any step and so before any located
event: the rates build nothing. Sharing it is sound: the rates are evaluated at events,
never within an evaluation of the Jacobian, so no Dual of the one is a perturbation of the
other's (the confusion tags guard against), and every evaluation of the residuals gathers
what it reads into the workspace first. A host adds its method to name its Jacobian's
(PowerfulYbus does for the Y-bus DAE with the network simultaneous); one that names none
evaluates them with Duals of their own tag (`AutomatonRateTag`) and `width` partials, its
workspace of that type built at the first use.
"""
rate_dual_type(host, width::Int) = ForwardDiff.Dual{AutomatonRateTag, Float64, width}

function rate_dual_type(adapter::TDSAdapter, width::Int)
    for T in keys(adapter.adcache.dual_workspaces)
        T <: ForwardDiff.Dual{ForwardDiff.Tag{JacobianClosureFunctor, Float64}, Float64} && return T
    end
    throw(ErrorException(
        "The rates of the automata are evaluated in the Dual workspace of the integration's Jacobian, which " *
        "this integration has not built (no Jacobian of `JacobianClosureFunctor` was evaluated before the " *
        "located event that needs them)."))
end

"""
    StateRates(host, n_diff, mass, sparsity) <: RateDAE

The rates of the DAE of a simulation's `host` (`dae_residuals!`, a `TDSAdapter` or the
Y-bus DAE), `mass` the diagonal of its mass matrix and `sparsity` the pattern of its
Jacobian (the integrator's `jac_prototype`), evaluated with Duals of `rate_dual_type`
in the host's workspace of that type (`dual_workspace`), whose memory is the Float64
workspace's, without derivatives (everything that sets memory writes it there too,
`sync_automaton_memory!`; `CHECK_MEMORY_SYNC` asserts it).
"""
struct StateRates{H, D} <: RateDAE
    host::H
    n::Int
    n_diff::Int
    mass::Vector{Float64}
    zero_mass::Vector{Int}
    zz_pattern::SparseMatrixCSC{Float64, Int}
    zz_colors::Vector{Int}
    zz_position::Vector{Int}
    u_d::Vector{D}
    du_d::Vector{D}
end

function StateRates(host, n_diff::Int, mass::Vector{Float64}, sparsity::AbstractSparseMatrix)
    n = length(mass)
    size(sparsity) == (n, n) || throw(ArgumentError(
        "The rates of the automata need the pattern of the Jacobian of the DAE ($n × $n); got $(size(sparsity))"))
    zero_mass = findall(iszero, mass)
    zz_pattern, zz_colors, zz_position = _zero_mass_pattern(sparsity, zero_mass, n)
    # (of their own tag: as many partials as F_ZZ has colors, up to the width of a
    # ForwardDiff chunk)
    D = rate_dual_type(host, clamp(maximum(zz_colors; init = 1), 1, 8))
    return StateRates(host, n, n_diff, mass, zero_mass, zz_pattern, zz_colors, zz_position, zeros(D, n), zeros(D, n))
end

# The pattern of F_ZZ (rows and columns `zero_mass` of the Jacobian pattern `sparsity`), a
# coloring of its columns, and the position in `zero_mass` of each variable
function _zero_mass_pattern(sparsity::AbstractSparseMatrix, zero_mass::Vector{Int}, n::Int)
    zz_pattern = SparseMatrixCSC{Float64, Int}(sparsity[zero_mass, zero_mass])
    zz_colors = isempty(zero_mass) ? Int[] : matrix_colors(zz_pattern)
    zz_position = zeros(Int, n)
    zz_position[zero_mass] .= eachindex(zero_mass)
    return zz_pattern, zz_colors, zz_position
end

_rate_host(r::StateRates) = r.host
_dual_workspace(r::StateRates) = dual_workspace(r.host, eltype(r.u_d))

"""
    CHECK_MEMORY_SYNC

Test hook: where `true` (or the environment variable `POWERFUL_CHECK_MEMORY_SYNC` is set at load),
every rate evaluation asserts that the memory of the automata in the Dual workspace it evaluates
in equals the Float64 workspace's, without derivatives. Everything that sets memory writes the Dual
workspaces too (`sync_automaton_memory!`), so the rates read it as it is.
"""
const CHECK_MEMORY_SYNC = Ref(false)

function _plain_values!(r::StateRates)
    CHECK_MEMORY_SYNC[] && !automaton_memory_synced(_dual_workspace(r), dual_workspace(r.host, Float64)) &&
        throw(ErrorException("The memory of the automata in the Dual workspace of a rate evaluation differs from the " *
                             "Float64 workspace's: a path that sets memory did not `sync_automaton_memory!`."))
    return nothing
end

"""
    dual_workspace(host, ::Type{D})

The workspace of element type `D` of the DAE of a simulation's `host` (see
`dae_residuals!`), whose modes and discrete quantities are those of its Float64
workspace: a `SystemBuffers` (or the Y-bus DAE's) with `workspaces_by_name`,
`models_by_name` and `formulation`. A package that runs the automata of its models with
this runtime adds its method.
"""
dual_workspace(adapter::TDSAdapter, ::Type{D}) where D = get_typed_workspace!(adapter.adcache, D)

# ----------------------------------------------------------------------------
# Dual numbers
# ----------------------------------------------------------------------------

# A Dual number of type `D` with value `v` and partials `f(c)`, c = 1..width
@inline _dual_with(::Type{D}, v, f) where D =
    D(v, ForwardDiff.Partials(ntuple(c -> Float64(f(c)), Val(ForwardDiff.npartials(D)))))

# The Dual state: the values of u, and partial `c` of entry `i` from `partial(i, c)`
@inline function _load_dual!(u_d::Vector{D}, u::AbstractVector, partial) where D
    @inbounds for i in eachindex(u_d)
        u_d[i] = _dual_with(D, u[i], c -> partial(i, c))
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Rates
# ----------------------------------------------------------------------------

# The values of the automata of model `m` in the Dual workspace of `r` (`_dual_workspace`),
# at its Dual state `r.u_d` and the time `t` (a number or a Dual), with the moving values
# last written there; every variable gathered (those the values read among them)
function _dual_automaton_values(r::RateDAE, m::AutomatonModel, t)
    sys_ws = _dual_workspace(r)
    ws = sys_ws.workspaces_by_name[m.name]
    nd, n = r.n_diff, r.n
    evaluate_automaton_values!(ws, ws.var_gather_tuple, sys_ws.models_by_name[m.name], m.external,
                               view(r.u_d, 1:nd), view(r.u_d, (nd + 1):n), sys_ws.formulation, t,
                               sys_ws.metadata_by_name[m.name])
    return ws.automata.values
end

# The rate F/M of the rows with mass at the state u[1:n] and time t (zero elsewhere)
function _mass_rate(r::RateDAE, u::AbstractVector, t)
    n = r.n
    F = zeros(n)
    dae_residuals!(F, _rate_host(r), view(u, 1:n), t)
    return [iszero(r.mass[i]) ? 0.0 : F[i] / r.mass[i] for i in 1:n]
end

"""
    RateFactor

The KLU factorization of `F_ZZ` for the rates of a `RateDAE`: `klu`, kept from one factorization to
the next (the symbolic analysis, of the pattern of `F_ZZ`, is reused; the numeric factorization,
with pivoting, is made afresh each time: the modes of the automata change which entries are
nonzero), which holds the one built at the state `u[1:n]` and time `t` where `valid`. `reuse`
takes it for an evaluation at the same state and time where the equations of the rows without
mass are known to be the same (the caller's to decide: the discrete state that residuals
without mass read is unchanged). `built` counts the numeric factorizations.
"""
mutable struct RateFactor
    klu::Union{Nothing, KLU.KLUFactorization{Float64, Int}}
    valid::Bool
    u::Vector{Float64}
    t::Float64
    built::Int
end

RateFactor() = RateFactor(nothing, false, Float64[], NaN, 0)

# The factorization no longer stands for the current equations (the analysis is kept)
_forget!(f::RateFactor) = (f.valid = false; f)

"""
    SingularRatesError(rows, empty, variables, equations, rcond, t, cause, constraints)

`F_ZZ`, the Jacobian of the rows without mass with respect to the variables without mass,
is numerically singular at time `t`: the rates of the rows without mass, and with them the
rates of the values of the automata that read them, cannot be computed. `rows` (entries of `u`)
are the rows of `F_ZZ`; `empty` the rows and columns with no nonzero entry, named as the host names
them; `variables` and `equations` the largest entries of the null direction of its columns and of its
rows (inverse iteration on the factorization), the variables nothing determines and the equations that
do not: a combination of them is in the null space; `rcond` the reciprocal of KLU's estimate of its
condition number (`KLU.condest`; zero where the factorization failed), `cause` what was found, `constraints`
the constraints on the states alone that the host's reduction found in this discrete state and does not make index 1
(`UnsupportedConstraint`), the likely reason.
"""
struct SingularRatesError <: Exception
    rows::Vector{Int}
    empty::Vector{String}
    variables::Vector{String}
    equations::Vector{String}
    rcond::Float64
    t::Float64
    cause::String
    constraints::Vector{String}
end

function Base.showerror(io::IO, e::SingularRatesError)
    print(io, "SingularRatesError: the Jacobian F_ZZ of the $(length(e.rows)) rows without mass is numerically singular at t = ",
          e.t, " (reciprocal condition estimate ", e.rcond, ", ", e.cause, "), so the rates of the rows without mass, ",
          "and of the values of the automata that read them, cannot be computed.")
    isempty(e.empty) || print(io, " Rows or columns with no nonzero entry: ", join(first(e.empty, 12), ", "),
                              length(e.empty) > 12 ? ", ..." : "", ".")
    isempty(e.variables) || print(io, " The variables its null direction leaves undetermined: ", join(e.variables, ", "),
                                  "; the equations that do not determine them: ", join(e.equations, ", "), ".")
    isempty(e.constraints) || print(io, " The reduction found a constraint on the states alone in this discrete state that one ",
                                    "differentiation does not make index 1: ", join(e.constraints, "; "), ".")
    return nothing
end

# The reasons of the constraints the reduction of `host` found and does not reduce, in the modes it is in
_unreduced_constraints(host) = String[]
function _unreduced_constraints(adapter::TDSAdapter)
    red = adapter.reduction
    red === nothing && return String[]
    return [u.reason for u in red.unsupported if gate_open(u.gate)]
end

# How the host names the entry `i` of its state (`[x; y]`) as a variable and as an equation, for messages
_variable_name(host, i::Int, n_diff::Int) = "u[$i]"
_equation_name(host, i::Int, n_diff::Int) = "F[$i]"

function _variable_name(adapter::TDSAdapter, i::Int, n_diff::Int)
    am = get_address_manager(adapter.adcache)
    return string(i <= n_diff ? variable_name(am, StateVar(), UInt(i)) : variable_name(am, AlgebVar(), UInt(i - n_diff)))
end

function _equation_name(adapter::TDSAdapter, i::Int, n_diff::Int)
    am = get_address_manager(adapter.adcache)
    return string(i <= n_diff ? variable_name(am, StateRes(), UInt(i)) : variable_name(am, AlgebRes(), UInt(i - n_diff)))
end

"""
    _state_rate(r::RateDAE, u, t; factor=RateFactor(), reuse=false) -> Vector{Float64}

The rate of the whole state `u[1:n]` at time `t` in the current discrete state: `F/M`
on the rows with mass, and on the rows without mass the derivative of their equations
along the trajectory, `F_ZZ u̇_Z = -(F_ZW u̇_W + ∂F_Z/∂t)`. `factor` keeps the
factorization of `F_ZZ`, which `reuse` takes from it where it was built at this state and
time.
"""
_state_rate(r::RateDAE, u::AbstractVector, t; factor::RateFactor = RateFactor(), reuse::Bool = false) =
    _complete_direction!(_mass_rate(r, u, t), r, u, t; time = 1.0, factor, reuse)

"""
    _complete_direction!(d, r::RateDAE, u, t; time=0.0, factor=RateFactor(), reuse=false) -> d

Complete the direction `d` of the state `u[1:n]` (given on the rows with mass, zero on
the rows without) with the rows without mass that keep their equations: `F_ZZ d_Z =
-(F_ZW d_W + ∂F_Z/∂t time)`, the equations at time `t` moving with `time`. `F_ZZ` and the
derivative along `(d, time)` come from one evaluation of `F` with Duals where the colors of `F_ZZ`
leave a partial free for it (`_zero_mass_jacobian`); `factor` keeps the factorization, which
`reuse` takes at the same state and time (`_zero_mass_system!`).
"""
function _complete_direction!(d::Vector{Float64}, r::RateDAE, u::AbstractVector, t; time::Float64 = 0.0,
                              factor::RateFactor = RateFactor(), reuse::Bool = false)
    Z = r.zero_mass
    isempty(Z) && return d
    _plain_values!(r)
    K, along = _zero_mass_system!(factor, r, u, t, d, time; reuse)
    d[Z] .= .-(K \ along)
    return d
end

"""
    _zero_mass_system!(factor::RateFactor, r::RateDAE, u, t, d, time; reuse=false) -> (K, along)

The factorization `K` of `F_ZZ` at the state `u[1:n]` and time `t`, and the derivative `along` of
the rows without mass along the direction `d` of the state and `time` of the time: `K` that kept in
`factor` where `reuse` and it was built at exactly this state and time, otherwise built from the
Jacobian at this state (`_zero_mass_jacobian`) and kept there, its symbolic analysis reused
(`_factor_zero_mass`).
"""
function _zero_mass_system!(factor::RateFactor, r::RateDAE, u::AbstractVector, t, d::Vector{Float64}, time::Float64;
                            reuse::Bool = false)
    n = r.n
    if reuse && factor.valid && factor.t == t && length(factor.u) == n &&
       all(i -> @inbounds(factor.u[i] == u[i]), 1:n)
        return factor.klu, _derivative_along(r, u, t, d, time)
    end
    factor.valid = false
    J, along = _zero_mass_jacobian(r, u, t; direction = d, time)
    factor.klu = _factor_zero_mass(r, J, t, factor.klu)
    factor.u = collect(Float64, view(u, 1:n))
    factor.t = Float64(t)
    factor.built += 1
    factor.valid = true
    return factor.klu, along
end

"""
    _zero_mass_factor(r::RateDAE, u, t)

The KLU factorization of `F_ZZ`, the Jacobian of the rows without mass with respect to the
variables without mass, at the state `u[1:n]` and time `t` (`_zero_mass_jacobian`), checked
(`_factor_zero_mass`).
"""
_zero_mass_factor(r::RateDAE, u::AbstractVector, t) =
    _factor_zero_mass(r, first(_zero_mass_jacobian(r, u, t)), t, nothing)

# The derivative of the rows without mass of `F` at the state `u[1:n]` and time `t` along the
# direction `d` of the state and `time` of the time, by one evaluation with Duals
function _derivative_along(r::RateDAE, u::AbstractVector, t, d::AbstractVector, time::Float64)
    D = eltype(r.u_d)
    _load_dual!(r.u_d, u, (i, c) -> c == 1 ? d[i] : 0.0)
    dae_residuals!(r.du_d, _rate_host(r), r.u_d, _dual_with(D, t, c -> c == 1 ? time : 0.0))
    return [ForwardDiff.partials(r.du_d[i], 1) for i in r.zero_mass]
end

"""
    _zero_mass_jacobian(r::RateDAE, u, t; direction=nothing, time=0.0) -> (J, along)

`F_ZZ` at the state `u[1:n]` and time `t`, with the moving quantities at their values
(`_plain_values!` first: the caller's, which evaluates with them anyway), in the pattern
`r.zz_pattern`. Its columns are compressed by `r.zz_colors` and seeded, as many colors
per evaluation as the Dual type of `r.u_d` has partials, into evaluations of `F` in the
workspace of that type. With a `direction` of the state (and `time` of the time), `along` is the
derivative of the rows without mass along it (`_derivative_along`), taken from the partial of the
first evaluation that no color uses where there is one, from an evaluation of its own otherwise;
`nothing` without a `direction`.
"""
function _zero_mass_jacobian(r::RateDAE, u::AbstractVector, t; direction::Union{Nothing, AbstractVector} = nothing,
                             time::Float64 = 0.0)
    Z, P, colors, position = r.zero_mass, r.zz_pattern, r.zz_colors, r.zz_position
    D = eltype(r.u_d)
    J = copy(P)
    rows, vals = rowvals(J), nonzeros(J)
    width = ForwardDiff.npartials(D)
    n_colors = maximum(colors; init = 0)
    along = nothing
    for first_color in 1:width:n_colors
        n_pass = min(width, n_colors - first_color + 1)
        slot = (direction !== nothing && along === nothing && n_pass < width) ? n_pass + 1 : 0
        seed(i, c) = c == slot ? direction[i] : (k = position[i]; k != 0 && colors[k] == first_color - 1 + c ? 1.0 : 0.0)
        _load_dual!(r.u_d, u, seed)
        dae_residuals!(r.du_d, _rate_host(r), r.u_d, slot == 0 ? t : _dual_with(D, t, c -> c == slot ? time : 0.0))
        for j in axes(J, 2)
            c = colors[j] - first_color + 1
            1 <= c <= width || continue
            for q in nzrange(J, j)
                vals[q] = ForwardDiff.partials(r.du_d[Z[rows[q]]], c)
            end
        end
        slot == 0 || (along = [ForwardDiff.partials(r.du_d[i], slot) for i in Z])
    end
    direction !== nothing && along === nothing && (along = _derivative_along(r, u, t, direction, time))
    return J, along
end

"""
    _factor_zero_mass(r::RateDAE, J, t, previous) -> KLU.KLUFactorization

The KLU factorization (with pivoting, in full each time) of `F_ZZ` = `J` at time `t`: the numeric
factorization of `previous`, a factorization of a matrix of the same pattern, where there is one (its
symbolic analysis is the pattern's), a new one otherwise. Throws `SingularRatesError` where KLU finds
`J` singular, or estimates its condition number (`KLU.condest`) at or above `1 / eps()`: the rates
that would follow are noise. (`KLU.rcond`, the ratio of the extreme pivots, is not that estimate: for
an `F_ZZ` with an exact null direction it can stand at 1e-13 where the condition number is 1e17.)
"""
function _factor_zero_mass(r::RateDAE, J::SparseMatrixCSC{Float64, Int}, t, previous)
    K = try
        if previous === nothing
            klu(J)
        else
            copyto!(previous.nzval, nonzeros(J))
            KLU.klu_factor!(previous)
        end
    catch e
        e isa SingularException || rethrow()
        throw(_singular_rates(r, J, 0.0, t, "KLU found it singular", nothing))
    end
    κ = KLU.condest(K)
    κ * eps() >= 1 && throw(_singular_rates(r, J, 1 / κ, t, "KLU estimates its condition number at $κ", K))
    return K
end

"""
    SINGULAR_DIAGNOSIS_MAX_ROWS

The largest `F_ZZ` whose null direction a `SingularRatesError` finds by a dense SVD where KLU
found it exactly singular (the diagnosis only; above it the error says the direction was not
computed).
"""
const SINGULAR_DIAGNOSIS_MAX_ROWS = 3000

# The error for a singular `F_ZZ` `J` (rows `r.zero_mass`) at time `t`, factorized by `K` where it was: the
# rows and columns of `J` with no nonzero entry, and the largest entries of its null direction. With a
# factorization, by inverse iteration (`K \ b` is dominated by the right singular vector of the smallest
# singular value, `K' \ b` by the left one); where KLU found it exactly singular, by a dense SVD of `J` (its
# right and left singular vectors of the smallest singular value; the ratio of the extreme singular values
# is then the reciprocal condition reported)
function _singular_rates(r::RateDAE, J::SparseMatrixCSC, rcond::Float64, t, cause::String, K)
    Z, host, n_diff = r.zero_mass, _rate_host(r), r.n_diff
    empty_rows = setdiff(1:length(Z), rowvals(dropzeros(J)))
    empty_cols = [j for j in axes(J, 2) if all(iszero, view(nonzeros(J), nzrange(J, j)))]
    empty = vcat(["row " * _equation_name(host, Z[k], n_diff) for k in empty_rows],
                 ["column " * _variable_name(host, Z[k], n_diff) for k in empty_cols])
    variables, equations = String[], String[]
    largest(v) = (order = sortperm(abs.(v); rev = true); top = abs(v[first(order)]);
                  [k for k in first(order, 6) if abs(v[k]) >= 0.1 * top])
    if K !== nothing
        b = [sin(7.0 * k + 1.0) for k in eachindex(Z)]
        x, y = K \ b, K' \ b
        all(isfinite, x) && (variables = [_variable_name(host, Z[k], n_diff) for k in largest(x)])
        all(isfinite, y) && (equations = [_equation_name(host, Z[k], n_diff) for k in largest(y)])
    elseif length(Z) <= SINGULAR_DIAGNOSIS_MAX_ROWS
        F = svd(Matrix(J))
        rcond = F.S[end] / F.S[1]
        variables = [_variable_name(host, Z[k], n_diff) for k in largest(F.V[:, end])]
        equations = [_equation_name(host, Z[k], n_diff) for k in largest(F.U[:, end])]
        cause *= "; a dense SVD for the diagnosis gives the ratio of its extreme singular values, $rcond"
    else
        cause *= "; its null direction is not computed (a dense SVD of its $(length(Z)) rows exceeds " *
                 "SINGULAR_DIAGNOSIS_MAX_ROWS = $SINGULAR_DIAGNOSIS_MAX_ROWS)"
    end
    return SingularRatesError(copy(Z), empty, variables, equations, rcond, Float64(t), cause, _unreduced_constraints(host))
end

# The Dual columns of the values of the automata of model `m` at the state `u` moving at `rate`
# at time `t`, in the current discrete state
function _dual_value_columns(r::RateDAE, m::AutomatonModel, u::AbstractVector, rate::AbstractVector, t)
    D = eltype(r.u_d)
    _load_dual!(r.u_d, u, (a, c) -> c == 1 ? rate[a] : 0.0)
    _plain_values!(r)
    return values(_dual_automaton_values(r, m, _dual_with(D, t, c -> c == 1 ? 1.0 : 0.0)))
end

# The values of the automata of model `m` and their rates at the state `u` moving at
# `rate` at time `t`, in the current discrete state: matrices (value, instance)
_value_rates(r::RateDAE, m::AutomatonModel, u::AbstractVector, rate::AbstractVector, t) =
    _values_and_rates(_dual_value_columns(r, m, u, rate, t), n_instances(m))

# The values `V` and rates `R` alone, into matrices (value, instance)
function _value_rates!(V::Matrix{Float64}, R::Matrix{Float64}, r::RateDAE, m::AutomatonModel, u::AbstractVector,
                       rate::AbstractVector, t)
    _values_and_rates!(V, R, _dual_value_columns(r, m, u, rate, t))
    return nothing
end

# The values and first partials of the Dual columns (a function barrier: the workspace
# holding them is found by name, so their type is known at run time only)
_values_and_rates(columns::Tuple, n::Int) = _values_and_rates!(zeros(length(columns), n), zeros(length(columns), n), columns)

function _values_and_rates!(V::Matrix{Float64}, R::Matrix{Float64}, columns::Tuple)
    for k in eachindex(columns), i in axes(V, 2)
        V[k, i] = ForwardDiff.value(columns[k][i])
        R[k, i] = ForwardDiff.partials(columns[k][i], 1)
    end
    return V, R
end

# ----------------------------------------------------------------------------
# How well a rate is known: the bound of its error from the state's
# ----------------------------------------------------------------------------
#
# The rate `h(u, t) = ġ` of a value `g` along the DAE's trajectory (`_value_rates` along `_state_rate`) is a function of
# the state, which the integration knows to its tolerance: `|δu_j| <= tol_j = abstol_j + reltol_j |u_j|`. To first order
# `|δh| <= B = Σ_j |∂h/∂u_j| tol_j`, and the sign of `h` is the model's where `|h| > B`.
#
# With the weights `p = ∇_u g` (and `q = (p, ∂g/∂t)`), the rows with mass `W` (`u̇_W = F_W / M_W`) and the rows without
# mass `Z` (`J_ZZ u̇_Z = -(J_ZW u̇_W + F_Zt)`, so `d/dt F_Z = 0`), `h = q · (u̇, 1)`, and its gradient is, with the adjoint
# `λ = J_ZZ⁻ᵀ p_Z` (one solve with the factorization of `F_ZZ`),
#
#     ∇_u h = G + J_Wᵀ M_W⁻¹ (p_W - J_ZWᵀ λ) - D,
#     G = d/dε ∇_u g(u + ε u̇, t + ε),      D = d/dε J_Z(u + ε u̇, t + ε)ᵀ λ     (at ε = 0, λ fixed).
#
# `G` is the curvature of the value along the trajectory, `D` that of the rows without mass (the derivative of the
# rate of the variables without mass with respect to the state goes through the derivative of `J_Z u̇`, with the Hessian of
# `F_Z`), and the middle term the dependence of the rates of the states on the state. The derivation was checked with SymPy
# on a small DAE and the implementation against the difference quotients in every coordinate (test_automata.jl).

"""
    _rate_jacobian(r::RateDAE, u, t) -> SparseMatrixCSC

The Jacobian `∂F/∂u` of all the rows of the DAE at the state `u[1:n]` and time `t`, by Duals seeded `width` columns at a
time (the partials of the Dual type of `r.u_d`; `n / width` evaluations of `F`).
"""
function _rate_jacobian(r::RateDAE, u::AbstractVector, t)
    D = eltype(r.u_d)
    width, n = ForwardDiff.npartials(D), r.n
    rows, cols, vals = Int[], Int[], Float64[]
    for start in 1:width:n
        stop = min(start + width - 1, n)
        _load_dual!(r.u_d, u, (a, c) -> a == start + c - 1 ? 1.0 : 0.0)
        _plain_values!(r)
        dae_residuals!(r.du_d, _rate_host(r), r.u_d, t)
        for c in 1:(stop - start + 1), a in 1:n
            v = ForwardDiff.partials(r.du_d[a], c)
            iszero(v) || (push!(rows, a); push!(cols, start + c - 1); push!(vals, v))
        end
    end
    return sparse(rows, cols, vals, n, n)
end

"""
    _value_gradient(r::RateDAE, m::AutomatonModel, w, i, u, t; support = 1:n) -> Vector{Float64}

The gradient `∇_u g` of value `w` of instance `i` of model `m` at the state `u[1:n]` and time `t`, in the columns of `support`
(the others zero), by Duals seeded `width` columns at a time.
"""
function _value_gradient(r::RateDAE, m::AutomatonModel, w::Int, i::Int, u::AbstractVector, t;
                         support::AbstractVector{Int} = 1:r.n)
    D = eltype(r.u_d)
    width = ForwardDiff.npartials(D)
    p = zeros(r.n)
    for start in 1:width:length(support)
        chunk = view(support, start:min(start + width - 1, length(support)))
        seed = Dict(j => c for (c, j) in enumerate(chunk))
        _load_dual!(r.u_d, u, (a, c) -> get(seed, a, 0) == c ? 1.0 : 0.0)
        _plain_values!(r)
        columns = values(_dual_automaton_values(r, m, _dual_with(D, t, c -> 0.0)))
        for (c, j) in enumerate(chunk)
            p[j] = ForwardDiff.partials(columns[w][i], c)
        end
    end
    return p
end

"""
    _rate_gradient(r::RateDAE, m::AutomatonModel, w, i, u, t, rate) -> Vector{Float64}

`∇_u h`, the gradient with respect to the state `u[1:n]` of the rate `h = ġ` of value `w` of instance `i` of model `m` along
the state rate `rate` at time `t` (see the header of this section). The adjoint solve uses the factorization of `F_ZZ`
at `(u, t)` only where the value reads a variable without mass.
"""
function _rate_gradient(r::RateDAE, m::AutomatonModel, w::Int, i::Int, u::AbstractVector, t, rate::AbstractVector)
    n, Z = r.n, r.zero_mass
    p = _value_gradient(r, m, w, i, u, t)
    support = findall(!iszero, p)
    isempty(support) && return zeros(n)
    # the step along the trajectory: a state error of 1e-6 of the largest entry, and as far in time
    ε = 1e-6 * (1 + maximum(abs, view(u, 1:n))) / max(maximum(abs, view(rate, 1:n)), 1.0)
    shifted(sign) = (collect(Float64, view(u, 1:n)) .+ sign * ε .* view(rate, 1:n), t + sign * ε)
    (u⁺, t⁺), (u⁻, t⁻) = shifted(1.0), shifted(-1.0)
    G = (_value_gradient(r, m, w, i, u⁺, t⁺; support) .- _value_gradient(r, m, w, i, u⁻, t⁻; support)) ./ (2ε)
    J = _rate_jacobian(r, u, t)
    p_W = copy(p)
    λ = zeros(length(Z))
    D = zeros(n)
    if !isempty(Z) && any(!iszero, p[Z])
        K = _zero_mass_factor(r, u, t)
        λ = K' \ p[Z]
        p_W[Z] .= 0.0
        J_Zᵀλ = J[Z, :]' * λ
        p_W .-= J_Zᵀλ
        p_W[Z] .= 0.0
        D = (J_dual_rows(r, u⁺, t⁺, Z)' * λ .- J_dual_rows(r, u⁻, t⁻, Z)' * λ) ./ (2ε)
    end
    scaled = [iszero(r.mass[k]) ? 0.0 : p_W[k] / r.mass[k] for k in 1:n]
    return G .+ J' * scaled .- D
end

# The rows `Z` of the Jacobian at `(u, t)`
J_dual_rows(r::RateDAE, u, t, Z) = _rate_jacobian(r, u, t)[Z, :]

"""
    _rate_error_bound(r::RateDAE, m, w, i, u, t, rate, opts) -> Float64

`B = Σ_j |∂h/∂u_j| tol_j` for the rate `h` of value `w` of instance `i` of model `m` (`_rate_gradient`), `tol_j = abstol_j +
reltol_j |u_j|` the integration's tolerance of entry `j` (`opts` the integrator's options).
"""
function _rate_error_bound(r::RateDAE, m::AutomatonModel, w::Int, i::Int, u::AbstractVector, t, rate::AbstractVector, opts)
    a = _rate_gradient(r, m, w, i, u, t, rate)
    return sum(abs(a[j]) * _state_tolerance(opts, u, j) for j in 1:r.n)
end
