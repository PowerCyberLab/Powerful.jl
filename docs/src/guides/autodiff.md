# Automatic Differentiation

How Powerful.jl computes derivatives, and how to choose the AD backend.

## Where Derivatives Are Used

Every model's equations are written once, as residual functions that accept any number
type. Derivatives of these residuals are computed by automatic differentiation (AD)
wherever a method needs them:

- **Power flow**: the Newton Jacobian ``\partial g / \partial y``.
- **Time-domain simulation**: the Jacobian used by the implicit integrator (QNDF, FBDF),
  in every simulation mode.
- **Small-signal analysis**: the state matrix for [eigenvalues](eigenvalue.md) and the
  input matrix for [frequency response](bode.md) (`linearize`).

Nothing has to be set up for this: models provide no hand-written Jacobians, and the
Jacobians are exact (no finite differences).

## How Jacobians Are Computed

Power system Jacobians are large and sparse: each equation involves only a few
variables. Three steps keep the cost proportional to the number of nonzeros:

1. **Sparsity pattern from the model structure.** Each model declares which variables
   its residuals read (its own variables and those of the buses it connects to). The
   pattern of the Jacobian follows from these declarations and the variable addresses;
   it is built once, without evaluating the equations.
2. **Coloring.** Columns that never share a nonzero row are grouped into one color
   (greedy coloring from SparseMatrixColorings.jl). All columns of a color are
   differentiated in one pass, so a Jacobian costs one residual evaluation per color
   rather than one per variable.
3. **Forward-mode AD.** Each pass evaluates the residuals with ForwardDiff.jl dual
   numbers and writes the derivatives into the sparse matrix.
   DifferentiationInterface.jl orchestrates the passes.

The verbose output of the Newton solver shows the size and density of the Jacobian:

```@example autodiff
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton(); verbose = true)
nothing # hide
```

## Choosing the AD Backend

`ACPowerFlowProblem` and `TDSProblem` take a `backend` keyword:

| Backend | Description |
|---------|-------------|
| `ForwardDiffBackend()` | ForwardDiff.jl dual numbers (default) |
| `EnzymeBackend()` | Enzyme.jl, differentiating the compiled code |

Both produce the same Jacobian:

```@example autodiff
sol_fd = solve(ACPowerFlowProblem(sys; backend = ForwardDiffBackend()), Newton())
sol_ez = solve(ACPowerFlowProblem(sys; backend = EnzymeBackend()), Newton())
(sol_ez.iterations, maximum(abs, sol_ez.dae.y .- sol_fd.dae.y))
```

For time-domain simulation:

```julia
prob = TDSProblem(init; tspan = (0.0, 10.0), backend = EnzymeBackend())
```

## Limitations

- Differentiating through `solve` itself (for example `ForwardDiff.gradient` of a
  function that calls `solve`) is not supported.
- For dynamic systems, `linearize` provides the state-space matrices around an
  equilibrium; see [Frequency Response](bode.md) and [Eigenvalue Analysis](eigenvalue.md).

## See Also

- [Power Flow Analysis](powerflow.md): Newton solver and formulations
- [Time-Domain Simulation](tds.md): integrators and solver options
