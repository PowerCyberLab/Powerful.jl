# Getting Started

This guide walks through the basic power flow workflow step by step and explains what
each piece does. All examples use the IEEE 14-bus case.

## The Basic Workflow

```@example getting_started
using Powerful
using PowerfulCases

# 1. Load system data
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# 2. Create the problem and solve it
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# 3. Check convergence
println("Converged: ", sol.converged)
println("Iterations: ", sol.iterations)
println("Final residual: ", sol.residual_norm)

# 4. Extract results
v = get_var(sol, :Bus, :v, (1,))          # Bus 1 voltage magnitude (pu)
θ = get_var(sol, :Bus, :theta, (1,))      # Bus 1 voltage angle (rad)
nothing # hide
```

The rest of this guide breaks down each step.

## Step 1: Loading System Data

Test cases are distributed through `PowerfulCases`. `load` returns a case bundle
that knows where the case files are:

```@example getting_started
case = load("ieee14")
case.raw     # Path to the PSS/E RAW file
case.dyr     # Path to the default dynamic data (DYR) file
basename(case.raw), basename(case.dyr)
```

List all available cases with `cases()`. Some files come in variants, such as the
two-island version of the IEEE 14-bus case:

```@example getting_started
basename(file(case, :raw; variant = "island"))
```

`parse_network` reads a PSS/E RAW file (from `PowerfulCases` or any path on disk)
into a `Network` with one table per data category:

```@example getting_started
network = parse_network(case.raw)
```

## Step 2: Creating a SystemModel

```@example getting_started
sys = SystemModel(network)
```

`SystemModel` converts the parsed data into a computational system:

1. **Organizes models** - Splits the data into model types (`Bus`, `Branch`,
   `Transformer2W`, `PQ`, `PV`, `Slack`, ...)
2. **Maps external IDs** - Builds a lookup from external keys (bus numbers, load and
   generator IDs) to internal array indices
3. **Converts to per-unit** - Stores per-unit values for computation and keeps the
   conversion coefficients

### Models as StructArrays

Each model is a StructArray (struct-of-arrays), so each field is a vector:

```@example getting_started
buses = sys.models.Bus
propertynames(buses)
```

```@example getting_started
using DataFrames

# i: bus number, ide: bus type (1 = PQ, 2 = PV, 3 = slack),
# vm: voltage magnitude from the case file (initial value, pu)
DataFrame(i = buses.i, ide = buses.ide, vm = buses.vm)
```

### External IDs and Internal Indices

Models refer to each other by **external IDs** from the data file. The system maps
external keys to **internal indices** (positions in the model arrays). Keys are
tuples: `(bus,)` for buses, `(bus, id)` for loads and generators, and
`(from_bus, to_bus, circuit)` for branches.

```@example getting_started
sys.order.key_indices[:Bus][(14,)]          # Internal index of bus 14
```

```@example getting_started
sys.order.key_indices[:PQ][(4, "1")]        # Internal index of load "1" at bus 4
```

A third layer maps internal indices to positions in the DAE vectors. It is created
with each problem, because different formulations have different variables. You do
not need to handle these layers directly: `get_var` takes the external key and
translates it. See [Three-Layer Indexing](../developer/indexing.md) for details.

### Per-Unit Values

`sys.models` holds per-unit values. The original physical values are kept in
`sys.parsed_data`, and the conversion coefficients in `sys.pu_coeffs`, with

```
per-unit = physical × coefficient
physical = per-unit / coefficient
```

For example, the active power of the first load, in MW and in per-unit on the
100 MVA system base:

```@example getting_started
p_mw = sys.parsed_data.models.PQ.pl[1]
coeff = sys.pu_coeffs[:PQ][:pl][1]
p_pu = sys.models.PQ.pl[1]
(p_mw, coeff, p_pu, p_mw * coeff ≈ p_pu)
```

Reports use these coefficients to convert results back to physical units.

## Step 3: Creating a Problem

```@example getting_started
prob = ACPowerFlowProblem(sys)
```

The problem specifies *what* to solve: the system, the formulation (default
`PowerInjection()`, with voltage angle and magnitude at each bus as variables), and
the variable layout. Creating a problem does not solve anything yet. It also detects
disconnected islands (see [Multi-Island Systems](multi_island.md)).

## Step 4: Solving

```@example getting_started
sol = solve(prob, Newton())
```

The algorithm specifies *how* to solve. At each Newton iteration the solver:

1. Evaluates the residuals g(y) (power balance equations)
2. Computes the sparse Jacobian J = ∂g/∂y by automatic differentiation
3. Solves the linear system J·Δy = -g
4. Updates y ← y + Δy

until the largest residual is below the tolerance (default `1e-8`). Tolerance and
iteration limit are keywords of the algorithm:

```@example getting_started
sol = solve(prob, Newton(abstol = 1e-10, maxiters = 20))
sol.residual_norm
```

### What's in the Solution?

```julia
sol.converged        # true/false
sol.iterations       # Number of iterations
sol.residual_norm    # Final residual norm
sol.retcode          # Return code
sol.dae.y            # Final state vector
sol.prob             # The problem (includes sol.prob.sys)
```

## Step 5: Extracting Results

`get_var` reads a variable by model, variable name, and external key:

```@example getting_started
v1 = get_var(sol, :Bus, :v, (1,))          # Voltage magnitude (pu)
θ14 = get_var(sol, :Bus, :theta, (14,))    # Angle (rad)
(v1, rad2deg(θ14))
```

`get_var` reads solution *variables*. Power injections, load consumption and branch
flows are computed from the solution and are available as tables:

```@example getting_started
DataFrame(collect_model_table(sol, :Branch))
```

See [Accessing Results](accessing_variables.md) and
[Report Customization](reports.md) for more.

## Design: Problems and Algorithms

Powerful.jl separates *what* to solve from *how* to solve it:

```@example getting_started
# WHAT to solve
prob = ACPowerFlowProblem(sys)

# HOW to solve
sol_nr = solve(prob, Newton())             # Native Newton-Raphson
sol_rb = solve(prob, RobustNR())           # NonlinearSolve.jl backend
sol_fd = solve(prob, FastDecoupledXB())    # Fast decoupled

[s.iterations for s in (sol_nr, sol_rb, sol_fd)]
```

The same problem is solved with different algorithms without rebuilding `sys`.

## Common Patterns

### Modify and Re-solve

Model parameters can be changed in place; create a new problem to pick up the
change. The third load (internal index 3) is load `(4, "1")` at bus 4, which has
no generator, so its voltage responds to the load change:

```@example getting_started
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# Increase load 3 by 20% (values are per-unit)
sys.models.PQ.pl[3] *= 1.2
sys.models.PQ.ql[3] *= 1.2

sol_new = solve(ACPowerFlowProblem(sys), Newton())

v_base = get_var(sol, :Bus, :v, (4,))
v_new = get_var(sol_new, :Bus, :v, (4,))
println("Voltage change at bus 4: ", v_new - v_base, " pu")
```

Buses with a voltage-controlling generator (PV buses) keep their voltage magnitude,
so load changes there do not change the voltage.

### Warm Starting

Pass a previous solution's state vector as the initial guess:

```@example getting_started
prob_warm = ACPowerFlowProblem(sys; initial_guess = sol.dae.y)
sol_warm = solve(prob_warm, Newton())
(sol_new.iterations, sol_warm.iterations)
```

The initial guess must come from a problem with the same formulation and variable
layout.

### Multiple Formulations

```@example getting_started
# Polar formulation: θ and V at each bus
sol_pi = solve(ACPowerFlowProblem(sys, PowerInjection()), Newton())

# Rectangular formulation: V_r and V_i at each bus
sol_ci = solve(ACPowerFlowProblem(sys, CurrentInjection()), Newton())

v_pi = get_var(sol_pi, :Bus, :v, (14,))
v_r = get_var(sol_ci, :Bus, :v_r, (14,))
v_i = get_var(sol_ci, :Bus, :v_i, (14,))
(v_pi, hypot(v_r, v_i))
```

Both formulations give the same operating point. Their variables differ, so the
variable names passed to `get_var` depend on the formulation.

## Next Steps

- **[Power Flow Guide](powerflow.md)** - Formulations, algorithms, advanced options
- **[Accessing Results](accessing_variables.md)** - Working with solutions
- **[Report Customization](reports.md)** - Tables in physical units
- **[Time-Domain Simulation](tds.md)** - Dynamic simulation from a power flow solution

## Summary

| Concept | What It Is | Why It Matters |
|---------|-----------|----------------|
| **SystemModel** | Computational representation of the network | Organizes data, maps IDs, stores per-unit values |
| **Problem** | What to solve | Separates the problem from the algorithm |
| **Algorithm** | How to solve | Affects speed and robustness |
| **Formulation** | Choice of variables and equations | Affects numerical behavior, not the physical result |
| **Solution** | Converged result | Holds the final state; read with `get_var` and report tables |
