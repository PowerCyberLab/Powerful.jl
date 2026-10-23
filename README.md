# Powerful.jl

[![CI](https://github.com/PowerCyberLab/Powerful.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PowerCyberLab/Powerful.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://powercyberlab.github.io/Powerful.jl/dev/)

Fast, composable power system simulation in Julia: power flow, time-domain
simulation, contingency analysis and Y-bus time-domain simulation, built on one
sparse differential-algebraic equation (DAE) system whose Jacobians come from
automatic differentiation.

- **Power flow**: Newton-Raphson, fast decoupled (XB/BX) and a NonlinearSolve.jl
  backend, in the power-injection or current-injection formulation, plus DC power
  flow. Multi-island systems are detected and solved automatically.
- **Time-domain simulation**: generators, exciters, governors and renewable models,
  solved with DifferentialEquations.jl, with eigenvalue and frequency-response
  analysis.
- **Contingency analysis**: N-1 and N-k AC screening, and DC screening with PTDF/LODF
  (`PowerfulCA`).
- **Y-bus simulation**: an alternative time-domain formulation with a sparse Y-bus
  network solve (`PowerfulYbus`).
- **Automatic differentiation**: every residual is AD-ready, with sparsity detection
  and matrix coloring from the model structure.

See the [documentation](https://powercyberlab.github.io/Powerful.jl/dev/) for guides,
a developer's architecture overview and the full API reference.

## Status

Research software: the API may change without notice between releases.

## Installation

Powerful.jl requires Julia 1.11 or later. Clone the repository and instantiate it;
`[sources]` in `Project.toml` resolves the package forks it depends on:

```bash
git clone https://github.com/PowerCyberLab/Powerful.jl
cd Powerful.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick Example

```julia
using Powerful
using PowerfulCases

# Power flow on the IEEE 14-bus case
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
sol = solve(ACPowerFlowProblem(sys), Newton())
println("Converged: ", sol.converged, " in ", sol.iterations, " iterations")

# A three-phase fault, cleared after 0.1 s
config = SimConfig(solver = (TDS = TDSConfig(
    convert_pq_to_zip = true, p_zip = (1.0, 0.0, 0.0), q_zip = (1.0, 0.0, 0.0)),))
dyr_data = parse_dynamics(case.dyr, sys)
init = setup_tds(sol, dyr_data; config = config)
prob = TDSProblem(init; tspan = (0.0, 2.0),
    events = [BusFault(0.5, 0.6, 7)])
tds_sol = solve(prob)
```

## Citation

If you use Powerful.jl in your research, please cite:

```bibtex
@software{powerful_jl,
  title   = {Powerful.jl},
  author  = {Cui, Hantao},
  url     = {https://github.com/PowerCyberLab/Powerful.jl},
  year    = {2026}
}
```

## License

GNU Affero General Public License v3.0 - see [LICENSE.md](LICENSE.md).

## Contributing

Contributions are welcome! Please open an issue or pull request on
[GitHub](https://github.com/PowerCyberLab/Powerful.jl).
