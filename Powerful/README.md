# Powerful

Umbrella package with the public API of Powerful.jl: power flow, time-domain
simulation, eigenvalue analysis, and reporting.

## Usage

```julia
using Powerful
using PowerfulCases

# Load system
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# Power flow
sol = solve(ACPowerFlowProblem(sys), Newton())

# Access results (keys are tuples)
v = get_var(sol, :Bus, :v, (1,))
print_model_table(sol, :Bus)

# Time-domain simulation
init = setup_tds(sol, parse_dynamics(case.dyr, sys))
tds_sol = solve(TDSProblem(init; tspan = (0.0, 5.0)))
t, v1 = get_series(tds_sol, :Bus, :v, (1,))
```

## Exported API

- `SystemModel`, `parse_network` - Build a system from a PSS/E RAW file
- `ACPowerFlowProblem`, `PowerInjection`, `CurrentInjection`, `DCPowerFlow` - Power flow problems and formulations
- `solve` - Unified solve interface
- `Newton`, `RobustNR`, `FastDecoupledXB`, `FastDecoupledBX` - Power flow algorithms
- `get_var`, `get_all_vars`, `get_vars_dict` - Solution variable access
- `collect_model_table`, `print_model_table`, `ReportConfig` - Result tables in physical units
- `analyze_topology`, `split_islands`, `MultiIslandSolution` - Island handling
- `parse_dynamics`, `setup_tds`, `TDSProblem`, `TDSSolution`, `get_series` - Time-domain simulation
- `BusFault`, `BranchTrip`, `BranchClose`, `FieldModification` - Simulation events
- `compute_eigenvalues`, `linearize` - Small-signal analysis
- `SimConfig`, `TDSConfig`, `load_config` - Configuration

## Model Types

Model types are in the `Powerful.Models` submodule:

```julia
using Powerful.Models: Bus, PQ, PV, Branch, GENCLS, TGOV1
```

## Related Packages

`PowerfulYbus` (Y-bus simulation) is loaded separately. Test cases come from
`PowerfulCases`.

## Dependencies

- PowerfulCore
- PowerfulPF
- PowerfulTDS
- PowerfulCases
