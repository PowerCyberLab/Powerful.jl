# Docstring Guidelines

Docstrings in Powerful.jl follow one of three tiers. The tier is selected by the visibility of the documented symbol, and it determines the length and the sections of the docstring.

| Tier | Scope | Length | Contents |
|------|-------|--------|----------|
| 1 | Public API | 10-30 lines | Summary, fields, constructors, example |
| 2 | Developer API | 3-8 lines | Summary, arguments and return value where needed |
| 3 | Internal helpers | 1-2 lines | Summary and key behavior |

A symbol exported from `Powerful` receives a Tier 1 docstring. A symbol exported only from a sub-package receives a Tier 2 docstring. All remaining symbols receive a Tier 3 docstring.

## Tier 1: Public API

Tier 1 applies to symbols that are exported from `Powerful` and used in the guides and examples. The docstring contains a one-line summary, the fields of a type, the constructor variants, and a usage example.

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

Tier 2 applies to symbols used by model and solver implementers and exported from sub-packages. The docstring contains a one-line summary, a brief description of the arguments when they are not self-explanatory, and the return type when it is not obvious. Sections such as `# Arguments`, `# Returns` and `# Implementation Notes` are omitted for signatures that the summary already explains.

```julia
"""
    setup_powerflow_newton(prob; use_gpu=false, ordering=:AMD, status, verbose=false) -> NewtonWorkspace

Set up the workspace for Newton power flow iterations.

Prepares the sparse Jacobian (sparsity pattern and coloring) and allocates the
iteration buffers.
"""
```

## Tier 3: Internal Helpers

Tier 3 applies to private functions (names with a `_` prefix) and other implementation details. The docstring is a one-line summary, followed by at most one short sentence of about 15 words for a critical side effect or assumption, in the pattern `Summary. [Behavior note.]`.

```julia
"""Gather each variable from x (states) or y (algebraic) into ws.var_vals."""

"""Scatter one residual entry by mode (:private/:add/:set). Skips inactive instances."""

"""Build (addr_pairs, val_pairs) for metadata items. Allocates on the workspace array type."""

"""Fault reactance (pu) below which a fault warns unless loads are constant impedance."""
```

The behavior notes worth stating are the following:

- Allocation: whether the function allocates or works on pre-allocated buffers.
- Skipped instances: which instances are filtered out, for example inactive ones.
- Platform: whether the function handles GPU arrays.
- Idempotence: whether repeated calls are safe.
- Errors: the conditions under which the function throws.
- Mutation: the arguments or workspace fields that are modified.

## DocStringExtensions

Struct fields are documented with `$TYPEDEF` and `$TYPEDFIELDS`, with a string preceding each field:

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

Physics and math functions (such as `calc_zip_load`) receive fuller documentation at any tier, because the formulas are essential context for the reader:

```julia
"""
    calc_zip_load(v, pl, ql, ip, iq, yp, yq, vmin) -> (p=..., q=...)

ZIP load model: P = pl·f(v) + ip·v + yp·v², Q = ql·f(v) + iq·v + yq·v²,
where f = constant_power_factor(v, vmin) is 1 at or above vmin and falls off
smoothly to zero below it.
"""
```
