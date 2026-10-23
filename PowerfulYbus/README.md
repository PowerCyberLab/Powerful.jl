# PowerfulYbus

Time-domain simulation with a sparse Y-bus network model for Powerful.jl. Machines
are Norton sources in the network admittance matrix, and the network can be solved
together with the machine equations or eliminated by a network solve at each step.
Load it separately with `using PowerfulYbus`.

## Usage

A `YBusTDSProblem` wraps a phasor `TDSProblem`:

```julia
using Powerful
using PowerfulCases
using PowerfulYbus

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
pf_sol = solve(ACPowerFlowProblem(sys), Newton())
init = setup_tds(pf_sol, parse_dynamics(case.dyr, sys))

events = [BusFault(1.0, 1.1, 5, "F1"; xf = 0.01)]
prob = TDSProblem(init; tspan = (0.0, 5.0), events = events)

sol = solve(YBusTDSProblem(prob; events = events))
t, v5 = get_series(sol, :Bus, :v, (5,))
```

## Features

- Sparse Y-bus network with Norton machine sources
- DifferentialEquations.jl solvers with the network equations simultaneous or eliminated
- Fixed-step partitioned trapezoidal integration
- Topology changes from events

See the [Y-Bus Simulation guide](../docs/src/guides/ybus.md).

## Dependencies

- PowerfulCore
- PowerfulPF
- PowerfulTDS
- KLU.jl
