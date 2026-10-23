# PowerfulPF

Power flow models and solvers for Powerful.jl. Its public API is available through
`using Powerful`.

## Models (`src/models/`)

- `Bus` - Network buses
- `PQ` - Loads (constant power, current, and admittance components)
- `PV` - Voltage-controlled generators
- `Slack` - Slack (swing) generators
- `Branch` - Transmission lines
- `Transformer2W` - Two-winding transformers
- `Shunt` - Shunt compensators
- Area, zone, and owner records: `AreaInterchange`, `Zone`, `Owner`

## Usage

```julia
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))

# Create the problem
prob = ACPowerFlowProblem(sys)                           # PowerInjection (default)
prob_ci = ACPowerFlowProblem(sys, CurrentInjection())    # Alternative formulation

# Solve with different algorithms
sol = solve(prob, Newton())            # Native Newton-Raphson
sol = solve(prob, RobustNR())          # NonlinearSolve.jl backend with globalization
sol = solve(prob, FastDecoupledXB())   # Fast decoupled (PowerInjection only)
sol = solve(prob, FastDecoupledBX())

# Results
get_var(sol, :Bus, :v, (1,))
print_model_table(sol, :Branch)
```

## Formulations

- `PowerInjection()` - Polar voltages (θ, V) with power balance equations (default)
- `CurrentInjection()` - Rectangular voltages (V_r, V_i) with current balance equations
- `DCPowerFlow()` - Linearized DC power flow

## Features

- Sparse Jacobians by automatic differentiation (ForwardDiff.jl), factorized with KLU
- Multi-island detection and solving
- PSS/E RAW file parsing (via PowerFlowData.jl)
- Result tables in physical units

See the [Power Flow guide](../docs/src/guides/powerflow.md).

## Dependencies

- PowerfulCore
- ForwardDiff.jl, DifferentiationInterface.jl, SparseConnectivityTracer.jl
- KLU.jl
- NonlinearSolve.jl (optional, for `RobustNR`)
