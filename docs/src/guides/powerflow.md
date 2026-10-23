# Power Flow Analysis

Guide to AC power flow in Powerful.jl: formulations, algorithms, and options.

## Basic Usage

```@example powerflow
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())
```

See [Getting Started](getting_started.md) for an explanation of each step.

## Formulations

The power flow equations can be written in different forms. The formulation is the
second argument of `ACPowerFlowProblem`.

### PowerInjection (Default)

Polar coordinates with power balance equations.

**Variables** at each bus: voltage angle `theta` and magnitude `v`.

**Equations** at each bus i, with θᵢⱼ = θᵢ - θⱼ and Yᵢⱼ = Gᵢⱼ + jBᵢⱼ the
admittance matrix entries:

```
Pᵢ - Vᵢ ∑ⱼ Vⱼ (Gᵢⱼ cos θᵢⱼ + Bᵢⱼ sin θᵢⱼ) = 0
Qᵢ - Vᵢ ∑ⱼ Vⱼ (Gᵢⱼ sin θᵢⱼ - Bᵢⱼ cos θᵢⱼ) = 0
```

where Pᵢ and Qᵢ are the net active and reactive power injected at bus i by
generators and loads.

```@example powerflow
prob = ACPowerFlowProblem(sys)                     # Default
prob = ACPowerFlowProblem(sys, PowerInjection())   # Explicit
nothing # hide
```

**When to use:** The default. Residuals are power mismatches, and the fast decoupled
algorithms require this formulation.

### CurrentInjection

Rectangular coordinates with current balance equations.

**Variables** at each bus: the real and imaginary voltage parts `v_r` and `v_i`.

**Equations** at each bus: the real and imaginary parts of the current balance. The
network current I = Y·V is linear in the rectangular voltages, so the network
equations contain no trigonometric functions.

```@example powerflow
prob_ci = ACPowerFlowProblem(sys, CurrentInjection())
sol_ci = solve(prob_ci, Newton())
```

The bus variables of this formulation are `v_r` and `v_i`; compute the magnitude
from them:

```@example powerflow
v_r = get_var(sol_ci, :Bus, :v_r, (14,))
v_i = get_var(sol_ci, :Bus, :v_i, (14,))
hypot(v_r, v_i), get_var(sol, :Bus, :v, (14,))
```

**When to use:** An alternative when PowerInjection has convergence difficulty.

### Comparison

| Aspect | PowerInjection | CurrentInjection |
|--------|----------------|------------------|
| **Coordinates** | Polar (θ, V) | Rectangular (V_r, V_i) |
| **Residuals** | Power mismatches | Current mismatches |
| **Network equations** | Trigonometric | Linear in voltage |
| **Fast decoupled algorithms** | Supported | Not supported |

Both give the same physical result; they differ in numerical behavior during the
iterations.

### Formulation Fallback Example

```@example powerflow
sol_pi = solve(ACPowerFlowProblem(sys, PowerInjection()), Newton())

if !sol_pi.converged
    @warn "PowerInjection did not converge, trying CurrentInjection"
    sol_pi = solve(ACPowerFlowProblem(sys, CurrentInjection()), Newton())
end
sol_pi.converged
```

### DC Power Flow

`DCPowerFlow()` is the linearized power flow (angles only, active power only) for
fast approximations and screening:

```@example powerflow
sol_dc = solve(ACPowerFlowProblem(sys, DCPowerFlow()))
get_var(sol_dc, :Bus, :theta, (14,))
```

See [DC Analysis](dc_analysis.md) for DC power flow, PTDF/LODF, and DC contingency
screening.

## Algorithms

The algorithm is the second argument of `solve`.

### Newton (Default)

Native Newton-Raphson with a sparse Jacobian from automatic differentiation,
factorized with KLU.

```@example powerflow
sol = solve(prob, Newton())
```

At each iteration:
1. Evaluate the residuals g(y)
2. Build the sparse Jacobian J = ∂g/∂y
3. Solve J·Δy = -g with KLU
4. Update y ← y + Δy

until the largest residual is below `abstol`.

Keyword arguments:

```julia
Newton(abstol = 1e-8,        # Convergence tolerance on the residual norm
       maxiters = 50,        # Maximum iterations
       ordering = :AMD)      # KLU fill-reducing ordering: :AMD, :COLAMD, :natural
```

**When to use:** The default. Fast and reliable for most systems.

### RobustNR

Newton through the NonlinearSolve.jl backend, with automatic globalization (trust
region, line search).

```@example powerflow
sol = solve(prob, RobustNR())
```

Keywords: `abstol`, `reltol`, `maxiters`.

**When to use:** When `Newton()` fails to converge, for example from a poor initial
guess on a heavily loaded system.

### FastDecoupledXB and FastDecoupledBX

Fast decoupled power flow. The P-θ and Q-V equations are decoupled and solved with
constant susceptance matrices B' and B'', which are factorized once and reused:

```
ΔP/V ≈ B' · Δθ
ΔQ/V ≈ B'' · ΔV
```

The XB and BX variants differ in how branch resistance enters B' and B''.

```@example powerflow
sol_xb = solve(prob, FastDecoupledXB())
sol_bx = solve(prob, FastDecoupledBX())
(sol_xb.iterations, sol_bx.iterations)
```

Keywords: `abstol` (default `1e-6`), `maxiters` (default 100), `ordering`.

Each iteration is cheaper than a Newton iteration, but more iterations are needed
and convergence is linear. The approximations assume small angle differences and
low resistance-to-reactance ratios, so these algorithms suit transmission networks
and are less reliable for distribution networks or heavily loaded systems.
Whether they are faster than `Newton()` depends on the system; time both on your
case.

Fast decoupled algorithms require the `PowerInjection` formulation.

**When to use:** Large transmission systems where many solves are needed.

### Algorithm Comparison

```@example powerflow
algorithms = ["Newton" => Newton(), "RobustNR" => RobustNR(),
              "FastDecoupledXB" => FastDecoupledXB(), "FastDecoupledBX" => FastDecoupledBX()]

for (name, alg) in algorithms
    s = solve(prob, alg)
    println(rpad(name, 16),
            " converged = ", s.converged,
            ", iterations = ", s.iterations,
            ", residual = ", s.residual_norm)
end
```

## Advanced Options

### Variable Ordering

The `addressing` keyword controls the layout of the variables in the DAE vectors.
`AMDAddressing()` reorders them with approximate minimum degree ordering:

```@example powerflow
prob_amd = ACPowerFlowProblem(sys, PowerInjection(); addressing = AMDAddressing())
sol_amd = solve(prob_amd, Newton())
```

The default is `NaturalAddressing()`. The result is the same; only the layout of the
state vector and Jacobian changes. KLU applies its own fill-reducing ordering
(the `ordering` keyword of `Newton`) in either case.

### Warm Starting

Use a previous solution's state vector as the initial guess:

```@example powerflow
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# Increase load (4, "1") at bus 4 by 20%
sys.models.PQ.pl[3] *= 1.2
sys.models.PQ.ql[3] *= 1.2

sol_cold = solve(ACPowerFlowProblem(sys), Newton())
sol_warm = solve(ACPowerFlowProblem(sys; initial_guess = sol.dae.y), Newton())
(sol_cold.iterations, sol_warm.iterations)
```

The initial guess must have the layout of the new problem (same formulation and
addressing). For systems with several islands, create one problem per island (see
[Multi-Island Systems](multi_island.md)).

**When to use:** Parameter sweeps and time series.

```@example powerflow
# Restore the load for the examples below
sys.models.PQ.pl[3] /= 1.2
sys.models.PQ.ql[3] /= 1.2
nothing # hide
```

## Convergence Diagnostics

### Check Convergence

```@example powerflow
if sol.converged
    println("Converged in $(sol.iterations) iterations, residual $(sol.residual_norm)")
else
    @warn "Did not converge" sol.iterations sol.residual_norm
end
```

Pass `verbose = true` to `solve` to print the setup steps and the residual at each
iteration.

### Common Convergence Issues

**Divergence (residual increasing)**

Causes:
- Infeasible operating point (too much load, not enough generation)
- Operating point near voltage collapse
- Poor initial guess

Remedies:
1. Check the power balance: total generation ≈ total load + losses
2. Try the `CurrentInjection` formulation
3. Try the `RobustNR()` algorithm
4. Increase the loading in steps, warm starting each solve

**Slow convergence (many iterations)**

Causes:
- Poorly conditioned Jacobian
- Unusual data (very small impedances, unrealistic setpoints)

Remedies:
1. Check for very small branch impedances
2. Check generator voltage setpoints and limits

## Multi-Island Systems

Disconnected networks are detected when the problem is created, and each island is
solved separately. The solution is a `MultiIslandSolution` that is read with the
same functions:

```@example powerflow
sys_island = SystemModel(parse_network(file(load("ieee14"), :raw; variant = "island")))
sol_island = solve(ACPowerFlowProblem(sys_island), Newton())
get_var(sol_island, :Bus, :v, (14,))
```

See [Multi-Island Systems](multi_island.md) for details.

## Accessing Results

```@example powerflow
v = get_var(sol, :Bus, :v, (1,))          # Voltage magnitude at bus 1
θ = get_var(sol, :Bus, :theta, (1,))      # Angle at bus 1 (rad)
all_v = get_all_vars(sol, :Bus, :v)       # All bus voltages (internal order)
v_dict = get_vars_dict(sol, :Bus, :v)     # Dict: (bus,) => voltage
v_dict[(14,)]
```

Power injections and branch flows are available as tables in physical units:

```@example powerflow
using DataFrames
DataFrame(collect_model_table(sol, :PV))
```

See [Accessing Results](accessing_variables.md) and
[Report Customization](reports.md).

## Examples

### Load Scaling Sweep

Scale all loads and record the lowest bus voltage, warm starting each solve from
the previous one:

```@example powerflow
function load_sweep(sys, factors)
    pl0 = copy(sys.models.PQ.pl)
    ql0 = copy(sys.models.PQ.ql)
    results = NamedTuple[]
    y_prev = nothing

    for factor in factors
        sys.models.PQ.pl .= factor .* pl0
        sys.models.PQ.ql .= factor .* ql0

        sol = solve(ACPowerFlowProblem(sys; initial_guess = y_prev), Newton())
        if !sol.converged
            @warn "Did not converge at factor = $factor"
            break
        end

        push!(results, (factor = factor, v_min = minimum(get_all_vars(sol, :Bus, :v)),
                        iterations = sol.iterations))
        y_prev = sol.dae.y
    end

    # Restore the original loads
    sys.models.PQ.pl .= pl0
    sys.models.PQ.ql .= ql0
    return results
end

DataFrame(load_sweep(sys, 1.0:0.2:2.0))
```

### Contingency Analysis

For contingency analysis, use the standalone `PowerfulCA` package:

```julia
using Powerful
using PowerfulCases
using PowerfulCA

sys = SystemModel(parse_network(load("ieee14").raw))

contingencies = [
    Contingency("LINE_1_5", [BranchOutage(1, 5, "1")]),
    Contingency("LINE_2_3", [BranchOutage(2, 3, "1")]),
    Contingency("GEN_3", [GeneratorOutage(3, "1")]),
]

result = solve_accc(sys, contingencies)
print_accc_report(result; vmin = 0.95, vmax = 1.05)
```

See [Contingency Analysis](accc.md).

## Further Reading

- [Getting Started](getting_started.md) - Understanding SystemModel
- [Accessing Results](accessing_variables.md) - Working with solutions
- [Multi-Island Systems](multi_island.md) - Disconnected networks
- [API Reference: Power Flow](../api/powerflow.md) - Complete API documentation
