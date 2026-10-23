# PowerfulTDS

Phasor-domain time-domain simulation (TDS) for Powerful.jl. Its public API is
available through `using Powerful`.

## Models

### Generators (`src/models/generators/`)
- `GENCLS` - Classical generator
- `GENROU` - Round-rotor generator
- `GENSAL` - Salient-pole generator
- `GENTPJ1` - Round-rotor generator with saturation on both axes

### Governors (`src/models/governors/`)
- `TGOV1` - Steam turbine governor

### Exciters (`src/models/exciters/`)
- `SEXS` - Simplified excitation system

### Renewables
- `REGCA1` - Renewable generator/converter
- `REECA1` - Renewable electrical control

## Usage

```julia
using Powerful
using PowerfulCases

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# Power flow provides the initial operating point
pf_sol = solve(ACPowerFlowProblem(sys), Newton())

# Dynamic models from the DYR file, initialized from the power flow
dyr = parse_dynamics(case.dyr, sys)
init = setup_tds(pf_sol, dyr)

# Simulate a 0.1 s fault at bus 5
events = [BusFault(1.0, 1.1, 5, "F1"; xf = 0.01)]
prob = TDSProblem(init; tspan = (0.0, 5.0), events = events)
sol = solve(prob)

t, v5 = get_series(sol, :Bus, :v, (5,))
```

## Features

- Initialization from a power flow solution
- DYR file parsing (PSS/E dynamic data)
- Integration with DifferentialEquations.jl
- Events: bus faults, branch trips and closures, parameter changes
- ZIP load conversion
- Eigenvalue analysis and linearization

See the [Time-Domain Simulation guide](../docs/src/guides/tds.md).

## Dependencies

- PowerfulCore
- PowerfulPF
- OrdinaryDiffEq components
- PowerDynData.jl (DYR parsing)
- ControlSystemsBase.jl (frequency response)
