# ============================================================================
# The Polynomial Degree of the Values of Automata in the State
# ============================================================================
#
# Within a step of a multistep method the state is a polynomial in t of degree at most the method's order
# `k` (at most 5 for QNDF and FBDF), so a value that is a polynomial of degree `d` in the state and the time
# is one of degree at most `d k` in t: the samples of the automaton callback determine it exactly where
# `d k` is below their number (`step_extremum.jl`). `PolyDegree` carries the degree of such a polynomial through
# a kernel (`trace_value_degrees`): a variable and the time have degree 1, a parameter, a constant, a discrete
# quantity and a memory (constant within a step) degree 0; `+` takes the larger degree, `*` adds them, `/` by a
# constant keeps the degree; a power by a non-negative integer multiplies it; and everything else of a
# variable (a division by it, a root, an exponential, an absolute value, a comparison) is not a polynomial
# (`NOT_POLYNOMIAL`). A comparison of two constants cannot be decided without their values, and throws
# `PolyBranch`, as SparseConnectivityTracer's tracers do: the trace is then repeated with the parameters at
# their values, as `trace_value_reads` does.
# ============================================================================

"""The degree of a quantity that is not a polynomial in the state and the time."""
const NOT_POLYNOMIAL = 1_000

"""
    PolyDegree(d)

A number known by the degree `d` of the polynomial in the state and the time it is (`NOT_POLYNOMIAL` where it
is none).
"""
struct PolyDegree <: Real
    d::Int
end

PolyDegree(x::Real) = PolyDegree(0)
PolyDegree(x::PolyDegree) = x

"""A comparison or a branch on a quantity whose value the trace does not know."""
struct PolyBranch <: Exception end

_deg(x::PolyDegree) = x.d
_deg(::Real) = 0
_poly(d::Int) = PolyDegree(min(d, NOT_POLYNOMIAL))

Base.promote_rule(::Type{PolyDegree}, ::Type{<:Real}) = PolyDegree
Base.convert(::Type{PolyDegree}, x::Real) = PolyDegree(0)
Base.convert(::Type{PolyDegree}, x::PolyDegree) = x
Base.zero(::Type{PolyDegree}) = PolyDegree(0)
Base.one(::Type{PolyDegree}) = PolyDegree(0)
Base.zero(::PolyDegree) = PolyDegree(0)
Base.one(::PolyDegree) = PolyDegree(0)
Base.float(x::PolyDegree) = x
Base.eps(::Type{PolyDegree}) = PolyDegree(0)
Base.typemax(::Type{PolyDegree}) = PolyDegree(0)
Base.typemin(::Type{PolyDegree}) = PolyDegree(0)

Base.:+(a::PolyDegree, b::PolyDegree) = _poly(max(a.d, b.d))
Base.:-(a::PolyDegree, b::PolyDegree) = _poly(max(a.d, b.d))
Base.:-(a::PolyDegree) = a
Base.:+(a::PolyDegree) = a
Base.:*(a::PolyDegree, b::PolyDegree) = _poly(a.d + b.d)
Base.:/(a::PolyDegree, b::PolyDegree) = b.d == 0 ? a : PolyDegree(NOT_POLYNOMIAL)
Base.:^(a::PolyDegree, n::Integer) = n >= 0 ? _poly(a.d * n) : (a.d == 0 ? a : PolyDegree(NOT_POLYNOMIAL))
Base.:^(a::PolyDegree, b::PolyDegree) = (a.d == 0 && b.d == 0) ? a : PolyDegree(NOT_POLYNOMIAL)
Base.:^(a::PolyDegree, b::Real) = (b == round(b) && b >= 0) ? _poly(a.d * Int(b)) : (a.d == 0 ? a : PolyDegree(NOT_POLYNOMIAL))
Base.:^(a::Real, b::PolyDegree) = (b.d == 0 && _deg(a) == 0) ? PolyDegree(0) : PolyDegree(NOT_POLYNOMIAL)
Base.inv(a::PolyDegree) = a.d == 0 ? a : PolyDegree(NOT_POLYNOMIAL)
Base.literal_pow(::typeof(^), a::PolyDegree, ::Val{n}) where {n} = a ^ n

# Functions that are polynomials of nothing but constants
for f in (:sqrt, :cbrt, :exp, :exp2, :expm1, :log, :log2, :log10, :log1p, :sin, :cos, :tan, :asin, :acos, :atan, :sinh,
          :cosh, :tanh, :abs, :sign, :floor, :ceil, :round, :trunc, :deg2rad, :rad2deg)
    @eval Base.$f(a::PolyDegree) = a.d == 0 ? a : PolyDegree(NOT_POLYNOMIAL)
end
Base.abs2(a::PolyDegree) = _poly(2 * a.d)     # (a polynomial)
for f in (:atan, :hypot, :min, :max, :copysign, :rem, :mod, :div)
    @eval begin
        Base.$f(a::PolyDegree, b::PolyDegree) = (a.d == 0 && b.d == 0) ? PolyDegree(0) : PolyDegree(NOT_POLYNOMIAL)
        Base.$f(a::PolyDegree, b::Real) = $f(a, PolyDegree(b))
        Base.$f(a::Real, b::PolyDegree) = $f(PolyDegree(a), b)
    end
end
Base.clamp(a::PolyDegree, lo, hi) = (a.d == 0 && _deg(lo) == 0 && _deg(hi) == 0) ? a : PolyDegree(NOT_POLYNOMIAL)
Base.ifelse(c::Bool, a::PolyDegree, b::PolyDegree) = c ? a : b
Base.ifelse(c::Bool, a::PolyDegree, b::Real) = c ? a : PolyDegree(b)
Base.ifelse(c::Bool, a::Real, b::PolyDegree) = c ? PolyDegree(a) : b
Base.isfinite(a::PolyDegree) = throw(PolyBranch())
Base.isnan(a::PolyDegree) = throw(PolyBranch())
Base.iszero(a::PolyDegree) = throw(PolyBranch())
for op in (:<, :<=, :(==), :isless, :isequal)
    @eval Base.$op(a::PolyDegree, b::PolyDegree) = throw(PolyBranch())
end
