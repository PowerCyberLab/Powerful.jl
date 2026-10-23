# ============================================================================
# The Rounding of the Values of Automata: a Running Error Bound
# ============================================================================
#
# A value of an automaton is computed in floating point from terms its kernel adds and subtracts (an atom `a - b`, a
# guard `x - limit`, a gain times a difference), and where those terms cancel, the value is the rounding of them: a
# start computed onto a surface (`Qmax` set from a sum that `Qinv` reaches by another, `c - x` with x within an ulp of
# c) is off zero by an ulp of the largest term, with the sign of that rounding. `RoundingBound` carries, through the
# model's own `calc_automaton_values!`, a first-order bound on that rounding, the running error bound of Wilkinson
# (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §3.3): each number `v` is carried with `e` such
# that `|computed - exact| <= u e + O(u^2)`, `u = eps()/2` the unit roundoff, where
#
# - every number that enters (a variable, the time, a parameter, a setpoint, a discrete quantity, a memory, a literal of
#   the kernel) is its own rounding, `e = |v|` (`fl(x) = x (1 + d)`, `|d| <= u`);
# - `a ± b` has `e = e_a + e_b + |a ± b|`, `a b` has `e = |b| e_a + |a| e_b + |a b|`, `a / b` has
#   `e = (e_a + |a / b| e_b) / |b| + |a / b|` (the errors of the operands carried, plus the rounding of the result);
# - a smooth function `f` has `e = |f'(a)| e_a + |f(a)|` (`f2(a, b)` sums over its arguments; the derivatives by
#   ForwardDiff);
# - `abs`, `min`, `max` and `clamp` keep the bound of what they return (`sum`, `minimum` and `maximum` reduce with
#   `+`, `min` and `max`), `sign`, `floor`, `ceil` and `round` are their own rounding, and a comparison compares the
#   values.
#
# So `c - x` is bounded by `|c| + |x| + |c - x|`: the terms it subtracts, which is how the value's rounding is read
# from the kernel, operation by operation, without the kernel exposing its terms. The values carried are those of the
# Float64 evaluation's operations: a branch of the kernel takes the side the Float64 evaluation does.
#
# A kernel's operation without a method here throws a `MethodError` naming it: a value whose rounding is unknown is
# never taken as exact (and `Float64(::RoundingBound)` has no method, so a kernel cannot drop the bound by writing a
# value into a Float64 buffer).
# ============================================================================

using ForwardDiff: ForwardDiff

"""
    RoundingBound(v, e)

A number `v` computed in Float64 with the running error bound `e` of its rounding: `|v - exact| <= u e` to first
order, `u = eps(Float64) / 2` (see the header of rounding_bound.jl). `RoundingBound(x)` of a real number that enters
a computation is its own rounding, `e = |x|`.
"""
struct RoundingBound <: Real
    v::Float64
    e::Float64
end

RoundingBound(x::Real) = (v = Float64(x); RoundingBound(v, abs(v)))
RoundingBound(x::RoundingBound) = x

Base.promote_rule(::Type{RoundingBound}, ::Type{<:Real}) = RoundingBound
Base.convert(::Type{RoundingBound}, x::Real) = RoundingBound(x)
Base.convert(::Type{RoundingBound}, x::RoundingBound) = x
Base.zero(::Type{RoundingBound}) = RoundingBound(0.0, 0.0)
Base.one(::Type{RoundingBound}) = RoundingBound(1.0)
Base.zero(::RoundingBound) = zero(RoundingBound)
Base.one(::RoundingBound) = one(RoundingBound)
Base.float(x::RoundingBound) = x

# The rounding of a result `v` computed from operands whose bounds the operation carries as `carried`
_rounded(v::Float64, carried::Float64) = RoundingBound(v, carried + abs(v))

Base.:+(a::RoundingBound) = a
Base.:-(a::RoundingBound) = RoundingBound(-a.v, a.e)
Base.:+(a::RoundingBound, b::RoundingBound) = _rounded(a.v + b.v, a.e + b.e)
Base.:-(a::RoundingBound, b::RoundingBound) = _rounded(a.v - b.v, a.e + b.e)
Base.:*(a::RoundingBound, b::RoundingBound) = _rounded(a.v * b.v, abs(b.v) * a.e + abs(a.v) * b.e)
function Base.:/(a::RoundingBound, b::RoundingBound)
    q = a.v / b.v
    return _rounded(q, (a.e + abs(q) * b.e) / abs(b.v))
end
Base.inv(a::RoundingBound) = one(RoundingBound) / a
Base.muladd(a::RoundingBound, b::RoundingBound, c::RoundingBound) = a * b + c
Base.fma(a::RoundingBound, b::RoundingBound, c::RoundingBound) = a * b + c
# (a mask: multiplying by a Bool selects, as Base's `Bool * x` does)
Base.:*(a::Bool, b::RoundingBound) = a ? b : zero(RoundingBound)
Base.:*(a::RoundingBound, b::Bool) = b * a

Base.:^(a::RoundingBound, n::Integer) = _rounded(a.v^n, abs(n * a.v^(n - 1)) * a.e)
Base.literal_pow(::typeof(^), a::RoundingBound, ::Val{n}) where {n} = a^n
Base.:^(a::RoundingBound, b::RoundingBound) =
    (p = a.v^b.v; _rounded(p, abs(b.v * a.v^(b.v - 1)) * a.e + abs(p * log(abs(a.v))) * b.e))

# Smooth functions of one argument: the error of the argument through the derivative, and the rounding of the result
for f in (:sqrt, :cbrt, :exp, :exp2, :exp10, :expm1, :log, :log2, :log10, :log1p, :sin, :cos, :tan, :sec, :csc, :cot,
          :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh, :sinpi, :cospi, :deg2rad, :rad2deg, :abs2)
    @eval Base.$f(a::RoundingBound) = _rounded($f(a.v), abs(ForwardDiff.derivative($f, a.v)) * a.e)
end
# ... of two
for f in (:atan, :hypot)
    @eval begin
        function Base.$f(a::RoundingBound, b::RoundingBound)
            da = ForwardDiff.derivative(x -> $f(x, b.v), a.v)
            db = ForwardDiff.derivative(x -> $f(a.v, x), b.v)
            return _rounded($f(a.v, b.v), abs(da) * a.e + abs(db) * b.e)
        end
    end
end

# What selects one of its arguments keeps its bound; a sign, and an integer `floor`, `ceil`, `round`
# or `trunc` takes (a step of the value, a surface of its own, not its rounding), is its own rounding
Base.abs(a::RoundingBound) = RoundingBound(abs(a.v), a.e)
Base.sign(a::RoundingBound) = RoundingBound(sign(a.v))
Base.round(a::RoundingBound, r::RoundingMode) = RoundingBound(round(a.v, r))
Base.copysign(a::RoundingBound, b::RoundingBound) = RoundingBound(copysign(a.v, b.v), a.e)
Base.copysign(a::RoundingBound, b::Real) = RoundingBound(copysign(a.v, b), a.e)
for f in (:min, :max)
    @eval Base.$f(a::RoundingBound, b::RoundingBound) = $f(a.v, b.v) == a.v ? a : b
end
Base.clamp(a::RoundingBound, lo::Real, hi::Real) = max(min(a, hi), lo)

# (with another real number, these and the arithmetic above promote it, as Base does for any two reals)
for op in (:<, :<=, :(==), :isless, :isequal)
    @eval Base.$op(a::RoundingBound, b::RoundingBound) = $op(a.v, b.v)
end
for f in (:isfinite, :isnan, :isinf, :iszero, :signbit)
    @eval Base.$f(a::RoundingBound) = $f(a.v)
end

"""
    within_rounding(g::Float64, e::Float64) -> Bool

Whether the value `g`, whose rounding is bounded by `u e` (`RoundingBound`, `u = eps()/2`), is zero to its rounding:
`|g| <= ROUNDING_ULPS u e`. The bound is first order in the rounding of the kernel's own operations and of its
inputs taken as exact numbers rounded once; the inputs the kernel reads were computed themselves (the initialization's
states and limits, the power flow's voltages), each to a few ulps of its own terms, which the kernel cannot see:
`ROUNDING_ULPS` allows for them. An exact zero is within its rounding.
"""
within_rounding(g::Float64, e::Float64) = iszero(g) || abs(g) <= ROUNDING_ULPS * (eps(Float64) / 2) * e

"""The ulps of its running error bound within which a value is zero to its rounding (`within_rounding`)."""
const ROUNDING_ULPS = 4
