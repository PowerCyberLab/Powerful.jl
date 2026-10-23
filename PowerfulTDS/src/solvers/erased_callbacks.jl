"""
Callbacks with small types.

The functions of a callback built from the problem (the closures of the event callbacks, the
condition and affect of the automata, the `initializealg` of an event) capture the problem, its buffers and the automata, so their
types name all of them, and so does the type of every object built from them: the callback, the
`CallbackSet`, the keyword arguments of `solve`. Each of these is compiled for the whole type,
and a callback of 5 to 7 KB takes the compiler minutes (`CallbackSet(...)` of the two callbacks
of the ieee14 Y-bus case: 27 s in a process of its own; the first solve of each Y-bus network
treatment: 150 s; a callback that also carries a Dual workspace for derivatives through the
switches: over 300 s).

A callback built with `ErasedCall` instead holds one `Ref{Any}`: its type is a few hundred
characters whatever the problem. Each call goes through one dynamic dispatch into a function
(`_call`) that is compiled for the captured function, and from there everything is as concretely
typed as before. The dispatch costs about 100 ns per call, against the milliseconds of a step.
"""

using PowerfulCore.DiffEq: VectorContinuousCallback, NoInit

"""
    ErasedCall{N}

A function of `N` arguments, held behind a `Ref{Any}`: `(::ErasedCall)(args...)` calls the
function with one dynamic dispatch. `N` is the number of arguments only so that the type of the
callable does not depend on the function it holds.
"""
struct ErasedCall{N}
    f::Base.RefValue{Any}
end

# The function barrier: compiled for the held function's type, called through one dynamic dispatch
_call(f, args...) = f(args...)

(c::ErasedCall{1})(a) = _call(c.f[], a)
(c::ErasedCall{2})(a, b) = _call(c.f[], a, b)
(c::ErasedCall{3})(a, b, x) = _call(c.f[], a, b, x)
(c::ErasedCall{4})(a, b, x, y) = _call(c.f[], a, b, x, y)

"""
    erased(f, ::Val{N}) -> ErasedCall{N}

`f` behind an [`ErasedCall`](@ref) of `N` arguments.
"""
erased(f, ::Val{N}) where N = ErasedCall{N}(Ref{Any}(f))

"""
    erase_callback(cb::VectorContinuousCallback) -> VectorContinuousCallback

The callback `cb` with its condition, affect and initialization behind [`ErasedCall`](@ref)s, so
that its type does not name what they capture. Its other settings are those of `cb`.
"""
function erase_callback(cb::VectorContinuousCallback)
    return VectorContinuousCallback(
        erased(cb.condition, Val(4)), erased(cb.affect!, Val(2)), cb.len;
        initialize = erased(cb.initialize, Val(4)), rootfind = cb.rootfind,
        save_positions = cb.save_positions, initializealg = cb.initializealg)
end

"""
    ErasedInit <: DAEInitializationAlgorithm

The initialization algorithm `alg` of an event behind a `Ref{Any}`, for the `initializealg` of a
callback: it is a type parameter of the callback and `alg` names the adapter and the
automata. `_initialize_dae!` forwards to `alg` through one dynamic dispatch.
"""
struct ErasedInit <: DAEInitializationAlgorithm
    alg::Base.RefValue{Any}
end

_initialize_with(alg, integrator, prob, ::Val{iip}) where {iip} = _initialize_dae!(integrator, prob, alg, Val(iip))

_initialize_dae!(integrator, prob, alg::ErasedInit, ::Val{iip}) where {iip} =
    _initialize_with(alg.alg[], integrator, prob, Val(iip))

"""
    erase_init(alg)

`alg` behind an [`ErasedInit`](@ref), except `NoInit`, which has no type to hide.
"""
erase_init(alg::NoInit) = alg
erase_init(alg) = ErasedInit(Ref{Any}(alg))
