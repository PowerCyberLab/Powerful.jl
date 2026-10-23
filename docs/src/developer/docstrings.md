# Docstring Guidelines

A three-tier framework for consistent, maintainable docstrings.

## Tier 1: Public API

**When:** Exported from `Powerful`, used in guides and examples.

**Include:**
- One-line summary
- Fields (for types)
- Constructor variants
- Usage example

**Example:**
````julia
"""
    ACPowerFlowProblem <: PowerFlowProblem

AC power flow problem with full nonlinear equations.

# Fields
- `sys`: Power system (single `SystemModel` or `Vector{SystemModel}` for multiple islands)
- `formulation`: `PowerInjection`, `CurrentInjection` or `DCPowerFlow`

# Example
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())
```
"""
````

## Tier 2: Developer API

**When:** Used by model and solver implementers, exported from sub-packages.

**Include:**
- One-line summary
- Brief arguments (only if non-obvious)
- Return type (only if non-obvious)

**Example:**
```julia
"""
    setup_powerflow_newton(prob; use_gpu=false, ordering=:AMD, status, verbose=false) -> NewtonWorkspace

Set up the workspace for Newton power flow iterations.

Prepares the sparse Jacobian (sparsity pattern and coloring) and allocates the
iteration buffers.
"""
```

**Anti-pattern (too verbose):**
```julia
"""
    setup_powerflow_newton(prob; use_gpu=false, ...)

Set up workspace for Newton-Raphson power flow iterations.

# Arguments
- `prob::ACPowerFlowProblem`: The power flow problem to solve
- `use_gpu::Bool=false`: Whether to use the GPU linear solver (default: false)

# Returns
- `NewtonWorkspace`: Workspace containing...

# Implementation Notes
This function performs the following steps:
1. Creates the autodiff cache...
2. Computes the sparsity pattern...
...
"""
```

## Tier 3: Internal Helpers

**When:** Private functions (`_` prefix), implementation details.

**Include:**
- One-line summary with the key behavior
- Critical side effects or assumptions (about 15 words at most)

**Pattern:** `Summary. [Behavior note.]`

**Examples:**
```julia
"""Gather each variable from x (states) or y (algebraic) into ws.var_vals."""

"""Scatter one residual entry by mode (:private/:add/:set). Skips inactive instances."""

"""Build (addr_pairs, val_pairs) for metadata items. Allocates on the workspace array type."""

"""Fault reactance (pu) below which a fault warns unless loads are constant impedance."""
```

**Key behavior notes to include:**
- **Allocation** - allocates, or works on pre-allocated buffers
- **Skips inactive** - filtering behavior
- **Handles GPU** - platform behavior
- **Idempotent** - safe for repeated calls
- **Throws if** - error conditions
- **Mutates** - side effects

## Quick Reference

| Tier | Scope | Length | Includes |
|------|-------|--------|----------|
| 1 | Public API | 10-30 lines | Summary, fields, example |
| 2 | Developer API | 3-8 lines | Summary, args/returns if needed |
| 3 | Internal | 1-2 lines | Summary + key behavior |

## Deciding the Tier

```
Is it exported from Powerful?
  → Yes: Tier 1 (full docs)
  → No: Is it exported from a sub-package?
    → Yes: Tier 2 (moderate docs)
    → No: Tier 3 (minimal docs)
```

## DocStringExtensions

Use `$TYPEDEF` and `$TYPEDFIELDS` for struct field documentation:

```julia
"""
    $TYPEDEF

Static PQ load data records.

# Fields
$TYPEDFIELDS
"""
@kwdef mutable struct PQ{Tv} <: AbstractModel
    "Bus number or extended bus name enclosed in single quotes"
    i::Int32
    ...
end
```

## Physics Functions

Physics and math functions (like `calc_zip_load`) warrant fuller documentation
regardless of tier, since the formulas are essential context:

```julia
"""
    calc_zip_load(v, pl, ql, ip, iq, yp, yq, vmin) -> (p=..., q=...)

ZIP load model: P = pl·f(v) + ip·v + yp·v², Q = ql·f(v) + iq·v + yq·v²,
where f = constant_power_factor(v, vmin) is 1 at or above vmin and falls off
smoothly to zero below it.
"""
```
