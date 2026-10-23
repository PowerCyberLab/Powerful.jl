"""
Precompilable Functors for TDS

This module defines callable struct types (functors) that replace anonymous closures
in the TDS solver. Named struct types can be precompiled, while anonymous closures
created at runtime get fresh type identities that invalidate precompilation.

# Why Functors?
Julia's precompilation caches method specializations by type. Anonymous closures
(like `x -> f(x, captured_var)`) get unique type identifiers each time they're
created, preventing precompilation reuse. Named structs with call methods have
stable types that persist across sessions.

# Type Barriers
We still use `Ref{Any}` to wrap complex types (like TDSAdapter) to prevent
type explosion in the ODEFunction type signature. This trades one runtime
dispatch per call for stable, precompilable types.

# Usage
Instead of:
```julia
ode_fn! = (du, u, p, t) -> evaluate!(du, adapter, u, t, ODEInterface())
```

Use:
```julia
ode_fn! = ODERHSFunctor(Ref{Any}(adapter))
```

The functor has stable type `ODERHSFunctor` regardless of what's inside the Ref.
"""

using DifferentiationInterface

# ============================================================================
# ODE Right-Hand-Side Functor
# ============================================================================

"""
    ODERHSFunctor

Callable struct for the ODE right-hand side: `du = f(u, p, t)`.

Replaces the anonymous closure in `create_diffeq_problem`. The struct type is
stable for precompilation while the adapter is type-erased via Ref{Any}.

# Fields
- `adapter_ref::Ref{Any}`: Type-erased TDSAdapter reference

# Call Signature
```julia
(f::ODERHSFunctor)(du, u, p, t)
```
"""
struct ODERHSFunctor
    adapter_ref::Ref{Any}
end

function (f::ODERHSFunctor)(du, u, p, t)
    adapter = f.adapter_ref[]::TDSAdapter
    evaluate!(du, adapter, u, t, ODEInterface())
    nothing
end

# ============================================================================
# Jacobian Closure Functor
# ============================================================================

"""
    JacobianClosureFunctor

Callable struct for Jacobian computation closure: `du = F(u)` at fixed time.

Used with DifferentiationInterface.jacobian! to compute sparse Jacobians.
The time reference is updated externally before each call.

# Fields
- `adapter_ref::Ref{Any}`: Type-erased TDSAdapter reference
- `t_ref::Ref{Float64}`: Mutable time reference

# Call Signature
```julia
(f::JacobianClosureFunctor)(du::AbstractVector{T}, u::AbstractVector{T}) where T
```
"""
struct JacobianClosureFunctor
    adapter_ref::Ref{Any}
    t_ref::Ref{Float64}
end

function (f::JacobianClosureFunctor)(du::AbstractVector{T}, u::AbstractVector{T}) where T
    _evaluate_packed!(du, f.adapter_ref[]::TDSAdapter, u, f.t_ref[], nothing)
end

# ============================================================================
# Tgrad Closure Functor
# ============================================================================

"""
    TgradClosureFunctor

Callable struct for time gradient computation: `du = F(u, t)` as function of t.

Used with DifferentiationInterface.derivative! to compute ∂F/∂t.
The state reference is updated externally before each call.

# Fields
- `adapter_ref::Ref{Any}`: Type-erased TDSAdapter reference
- `u_ref::Ref{Vector{Float64}}`: Mutable state reference

# Call Signature
```julia
(f::TgradClosureFunctor)(du::AbstractVector{T}, t::T) where T
```
"""
struct TgradClosureFunctor
    adapter_ref::Ref{Any}
    u_ref::Ref{Vector{Float64}}
end

function (f::TgradClosureFunctor)(du::AbstractVector{T}, t::T) where T
    # u is Float64; unpacking into T buffers gives Duals with zero partials,
    # so only explicit time dependence contributes to ∂F/∂t
    _evaluate_packed!(du, f.adapter_ref[]::TDSAdapter, f.u_ref[], t, nothing)
end

# ============================================================================
# Sparse Jacobian Functor
# ============================================================================

"""
    SparseJacobianFunctor

Callable struct for computing sparse Jacobians with DifferentiationInterface.

Wraps the Jacobian closure, AD preparation, and backend in a single callable
that updates time before computing.

# Fields
- `jac_closure_ref::Ref{Any}`: Type-erased JacobianClosureFunctor
- `jac_prep_ref::Ref{Any}`: Type-erased jacobian preparation from DI
- `jac_backend_ref::Ref{Any}`: Type-erased AD backend
- `t_ref::Ref{Float64}`: Time reference (shared with closure)
- `du_temp::Vector{Float64}`: Temporary output buffer

# Call Signature
```julia
(f::SparseJacobianFunctor)(J, u, p, t)
```
"""
struct SparseJacobianFunctor
    jac_closure_ref::Ref{Any}
    jac_prep_ref::Ref{Any}
    jac_backend_ref::Ref{Any}
    t_ref::Ref{Float64}
    du_temp::Vector{Float64}
end

function (f::SparseJacobianFunctor)(J, u, p, t)
    f.t_ref[] = t
    DifferentiationInterface.jacobian!(f.jac_closure_ref[], f.du_temp, J, f.jac_prep_ref[], f.jac_backend_ref[], u)
    nothing
end

# ============================================================================
# Sparse Tgrad Functor
# ============================================================================

"""
    SparseTgradFunctor

Callable struct for computing time gradients with DifferentiationInterface.

Wraps the tgrad closure, AD preparation, and backend in a single callable
that updates state before computing.

# Fields
- `tgrad_closure_ref::Ref{Any}`: Type-erased TgradClosureFunctor
- `tgrad_prep_ref::Ref{Any}`: Type-erased derivative preparation from DI
- `tgrad_backend_ref::Ref{Any}`: Type-erased AD backend
- `u_ref::Ref{Vector{Float64}}`: State reference (shared with closure)
- `du_temp::Vector{Float64}`: Temporary output buffer

# Call Signature
```julia
(f::SparseTgradFunctor)(dT, u, p, t)
```
"""
struct SparseTgradFunctor
    tgrad_closure_ref::Ref{Any}
    tgrad_prep_ref::Ref{Any}
    tgrad_backend_ref::Ref{Any}
    u_ref::Ref{Vector{Float64}}
    du_temp::Vector{Float64}
end

function (f::SparseTgradFunctor)(dT, u, p, t)
    # Update state reference (non-allocating copy)
    u_storage = f.u_ref[]
    @inbounds for i in eachindex(u)
        u_storage[i] = u[i]
    end
    # Compute ∂f/∂t via AD
    DifferentiationInterface.derivative!(f.tgrad_closure_ref[], f.du_temp, dT, f.tgrad_prep_ref[], f.tgrad_backend_ref[], t)
    nothing
end

# ============================================================================
# Exports
# ============================================================================

export ODERHSFunctor, JacobianClosureFunctor, TgradClosureFunctor
export SparseJacobianFunctor, SparseTgradFunctor
