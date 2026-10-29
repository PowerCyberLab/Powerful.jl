```@meta
CurrentModule = Powerful
```

# Powerful.jl

**Fast, composable power system simulation in Julia**

Powerful.jl is a power systems analysis package for power flow, time-domain
simulation, and related studies. It builds every study on one sparse
differential-algebraic equation (DAE) system whose Jacobians come from automatic
differentiation.

## Key Features

- **Problem/algorithm API**: DifferentialEquations.jl-style interface: build a problem, then `solve` it with an algorithm
- **Sparse solvers**: Sparse Jacobians from automatic differentiation, factorized with KLU
- **Multiple formulations**: Power injection (polar θ, V) and current injection (rectangular V_r, V_i) AC power flow, plus DC power flow
- **Power flow algorithms**: Native Newton-Raphson, fast decoupled (XB and BX), and a NonlinearSolve.jl backend
- **Automatic differentiation**: ForwardDiff.jl throughout, usable for sensitivities and gradients
- **Island detection**: Disconnected networks are detected and each island is solved
- **Time-domain simulation**: Dynamic models (generators, exciters, governors, renewables) solved with DifferentialEquations.jl
- **More studies**: Eigenvalue and frequency-response analysis<!--GENERATED:more-studies-->

## Quick Start

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/cuihantao/Powerful.jl")
```

### Basic Power Flow

```@example index
using Powerful
using PowerfulCases

# Load a test case and build the system
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# Create the problem and solve it
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())
```

Check convergence and read results by the external bus number. Keys are tuples:

```@example index
println("Converged: ", sol.converged, " in ", sol.iterations, " iterations")

v = get_var(sol, :Bus, :v, (1,))       # Bus 1 voltage magnitude (pu)
θ = get_var(sol, :Bus, :theta, (1,))   # Bus 1 voltage angle (rad)
(v, θ)
```

### Formulations and Algorithms

The same system can be solved with different formulations and algorithms:

```@example index
# Current injection formulation (rectangular voltages)
sol_ci = solve(ACPowerFlowProblem(sys, CurrentInjection()), Newton())

# Fast decoupled power flow (requires the default PowerInjection formulation)
sol_fd = solve(prob, FastDecoupledXB())

# NonlinearSolve.jl backend with automatic globalization
sol_robust = solve(prob, RobustNR())

(sol_ci.converged, sol_fd.converged, sol_robust.converged)
```

<!--GENERATED:more-examples-->
### Generating Reports

Result tables are converted to physical units (MW, MVAR, degrees):

```@example index
using DataFrames

# Branch flows as a table
DataFrame(collect_model_table(sol, :Branch))
```

A full text report:

```@example index
report = build_powerflow_report(sol, Default)
print(format_text(report))
```

## Architecture Overview

Powerful.jl uses a layered architecture that separates concerns:

1. **Foundation**: Type system, metadata, traits
2. **Indexing**: External IDs ↔ internal indices ↔ DAE addresses
3. **Models**: Component implementations (buses, loads, generators, branches)
4. **System**: Coordinates models into a unified DAE system
5. **Solvers**: Numerical methods using the DAE interface

This design lets different formulations and studies reuse the same model
implementations, and lets new models and solvers be added independently. See the
[Developer Documentation](developer/overview.md) for details.

## Packages

| Package | Description |
|---------|-------------|
| `Powerful` | Umbrella package with the public API (power flow, TDS, eigenvalues) |
| `PowerfulCore` | Core infrastructure: types, DAE, address management |
| `PowerfulPF` | Power flow models and solvers |
| `PowerfulTDS` | Time-domain simulation |
<!--GENERATED:packages-table-->
| `PowerfulCases` | Test case data (IEEE, synthetic, MATPOWER cases) |

## Documentation Structure

- [User Guides](guides/index.md) - Task-oriented guides, starting with [Getting Started](guides/getting_started.md)
- [API Reference](api.md) - Complete API documentation
- [Index](index_page.md) - Function and type index
- [Developer Documentation](developer/overview.md) - Architecture and internals

## Command-Line Interface

A development checkout includes a CLI that keeps a Julia daemon running for fast
repeated commands:

```bash
# Power flow report for a PSS/E RAW file (daemon auto-starts if needed)
./bin/powerful pf path/to/case.raw --output report.txt

# Run tests
./bin/powerful test                                       # Run all tests
./bin/powerful test Powerful/test/tds/test_dyr_trait.jl   # Run one test file

# Run Julia code or scripts
./bin/powerful run examples/tds/tds_flat_run.jl
```

See [bin/README.md](https://github.com/cuihantao/Powerful.jl/blob/main/bin/README.md) for complete CLI documentation.

## Citation

If you use Powerful.jl in your research, please cite:

```bibtex
@software{powerful_jl,
  author = {Cui, Hantao},
  title = {Powerful.jl: Fast Power System Simulation in Julia},
  year = {2025},
  url = {https://github.com/cuihantao/Powerful.jl}
}
```

## License

GNU Affero General Public License v3.0 - see LICENSE file for details.

## Contributing

Contributions are welcome! Please open an issue or pull request on [GitHub](https://github.com/cuihantao/Powerful.jl).

## Acknowledgements

This material is based upon work supported by the National Science Foundation under
Award Nos. ECCS-2514755 and ECCS-2525884, and by the National Academies of Sciences,
Engineering, and Medicine under Award No. SCON-10001546. Any opinions, findings, and
conclusions or recommendations expressed in this material are those of the author(s)
and do not necessarily reflect the views of the National Science Foundation or the
National Academies of Sciences, Engineering, and Medicine.
