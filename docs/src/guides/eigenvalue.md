# Eigenvalue Analysis

Small-signal stability analysis through the eigenvalues of the linearized DAE system.

## Quick Start

```@example eig
using Powerful
using PowerfulCases

# Load system and solve power flow
case = load("ieee14")
sys_pf = SystemModel(parse_network(case.raw))
sol_pf = solve(ACPowerFlowProblem(sys_pf), Newton())

# Dynamic models: GENCLS generators, SEXS exciters, TGOV1 governors
dyr_data = parse_dynamics(case.dyr, sys_pf)

# Create TDS problem at the power flow operating point
init = setup_tds(sol_pf, dyr_data)
prob = TDSProblem(init; tspan = (0.0, 10.0))

# Compute eigenvalues
result = compute_eigenvalues(prob)
result.n_positive, result.n_zeros, result.n_negative
```

```@example eig
# Print the report (the participation factor table is omitted above 10 states)
print_eigenvalue_report(result; max_modes = 10)
```

## Understanding the Results

The `EigenvalueResult` struct contains:

| Field | Type | Description |
|-------|------|-------------|
| `eigenvalues` | `Vector` | Complex eigenvalues (modes) |
| `right_eigenvectors` | `Matrix` | Mode shapes (columns) |
| `left_eigenvectors` | `Matrix` | Left eigenvectors |
| `participation_factors` | `Matrix{Float64}` | State-mode associations (states × modes) |
| `state_names` | `Vector{Symbol}` | Name of each state, e.g. `:GENCLS_omega_2` |
| `n_positive` | `Int` | Unstable modes (Re > tol) |
| `n_zeros` | `Int` | Marginal modes (\|Re\| <= tol) |
| `n_negative` | `Int` | Stable modes (Re < -tol) |
| `zero_mass_indices` | `Vector{Int}` | States with a zero mass coefficient (reduced) |
| `n_states` | `Int` | Number of states in the state matrix |

State names are `{Model}_{state}_{n}`, where `n` is the position of the instance in
the model (the order of `prob.sys.models.GENCLS`, for example), not its bus number:

```@example eig
prob.sys.models.GENCLS.BUS   # bus of GENCLS instance 1, 2, ...
```

### Stability Interpretation

| n_positive | Meaning |
|------------|---------|
| 0 | System is stable or marginally stable |
| > 0 | System has unstable modes |

When every generator is modeled dynamically and none is an infinite bus, the rotor
angles are only defined relative to each other, which gives one eigenvalue at zero (the
marginal mode in the Quick Start).

Oscillatory modes appear as complex conjugate pairs. The imaginary part gives the
oscillation frequency and the real part the damping:

```@example eig
λ = result.eigenvalues
osc = findall(x -> imag(x) > 0, λ)           # one of each conjugate pair
freq_hz = imag.(λ[osc]) ./ (2π)               # oscillation frequency
damping = -real.(λ[osc]) ./ abs.(λ[osc])      # damping ratio
[freq_hz damping]
```

### Participation Factors

Participation factors indicate which states contribute most to each mode. Each
column (mode) sums to 1:

```@example eig
# Most associated state for each oscillatory mode
[result.state_names[argmax(result.participation_factors[:, j])] for j in osc]
```

### Eigenvalue Plot

```@example eig
using Plots
scatter(real.(λ), imag.(λ); xlabel = "Real (1/s)", ylabel = "Imaginary (rad/s)",
        label = "Eigenvalues")
```

## Reporting

### Console Output

```julia
print_eigenvalue_report(result)                    # stdout; participation factors up to 50 states
print_eigenvalue_report(result; max_modes = 100)   # participation factors up to 100 states
print_eigenvalue_report(result; io = IOBuffer())   # Custom IO
```

The report lists every mode; `max_modes` is the largest number of modes for which the
participation factor table is printed.

### File Output

```julia
write_eigenvalue_report(result, "stability.txt")
```

## Theory

The DAE system is linearized at the steady state:

```
T · dx/dt = f(x, y)    (differential equations)
        0 = g(x, y)    (algebraic equations)
```

Eliminating the algebraic variables gives the state matrix:

```
As = T⁻¹ · (fx - fy · gy⁻¹ · gx)
```

where `fx`, `fy`, `gx`, `gy` are the Jacobian blocks (computed with automatic
differentiation) and `T` is the diagonal mass matrix of the differential states (for
example `2H` for rotor speeds).

### Zero Mass Coefficient Handling

A differential state whose mass coefficient is zero (for example a lag block with a
zero time constant) is an algebraic equation in disguise. Such states are:
1. Identified from the mass matrix diagonal
2. Moved to the end of the state vector
3. Eliminated together with the algebraic variables

Their indices are listed in `result.zero_mass_indices`.

## Complete Example

Synchronous machines with converter-interfaced renewables, loaded from a DYR variant of
the IEEE 14-bus case (GENROU generators with SEXS and TGOV1, REGCA1 converters with
REECA1 electrical controls):

```@example eig
sys2 = SystemModel(parse_network(case.raw))
sol2 = solve(ACPowerFlowProblem(sys2), Newton())
dyr2 = parse_dynamics(file(case, :dyr; variant = "reeca1"), sys2)
keys(dyr2)
```

```@example eig
prob2 = TDSProblem(setup_tds(sol2, dyr2); tspan = (0.0, 10.0))
result2 = compute_eigenvalues(prob2)

if result2.n_positive > 0
    println("System has $(result2.n_positive) unstable modes")
else
    println("System is stable ($(result2.n_negative) stable, $(result2.n_zeros) marginal)")
end
```

`variants(case, :psse_dyr)` lists the dynamic data variants of a case.

## Sensitivity to Model Parameters

`eigenvalue_sensitivity` gives the derivatives of the eigenvalues with respect to
parameters of the dynamic models (controller gains, time constants, inertias, ...), for
tuning controllers to a required damping with a gradient-based optimizer. It takes the
power flow and the dynamic models (as `setup_tds`) and the parameters as
`(model, field, instances)` tuples:

```@example eig
params = [(:SEXS, :K, 1:length(dyr_data.SEXS))]   # every exciter gain
s = eigenvalue_sensitivity(sol_pf, dyr_data, params)
size(s.dlambda)   # eigenvalues × parameters
```

A parameter moves the eigenvalues through the Jacobian and through the operating
point: the power flow is held, and the states and setpoints of the dynamic models
follow the parameter through the initialization (for example the voltage reference of
an exciter whose gain changes). Both are differentiated exactly with forward-mode AD.

`damping_ratio_sensitivity` converts the derivatives to those of the damping ratios
`ζ = -σ/|λ|`. The gradient of the least-damped oscillatory mode:

```@example eig
λ = s.eigen.eigenvalues
osc = findall(>(0.1), imag.(λ))
i = osc[argmin(-real.(λ[osc]) ./ abs.(λ[osc]))]
dζ_dK = damping_ratio_sensitivity(s)[i, :]
```

The derivative of an eigenvalue is defined while it is distinct; for an objective over
several modes (such as the minimum damping ratio), use a smooth form such as a soft
minimum, since the least-damped mode can change as the parameters move.

## Building Dynamic Models in Code

Dynamic models can also be constructed directly instead of parsed from a DYR file.
Pass them to `setup_tds` as a NamedTuple of StructArrays:

```julia
using Powerful.Models: GENCLS
using StructArrays

slack_bus = sys_pf.models.Slack.i[1]
gencls = StructArray([
    GENCLS(BUS = slack_bus, ID = "1", H = 5.0, D = 0.0, xd_p = 0.3, r_a = 0.0)
])

init = setup_tds(sol_pf, (GENCLS = gencls,))
result = compute_eigenvalues(TDSProblem(init))
```

Generators without a dynamic model stay static (PV or Slack), so this system has a
single machine against constant-voltage buses.

## API Reference

### `compute_eigenvalues`

```julia
compute_eigenvalues(prob::TDSProblem; tol = 1e-6, verbose = false) -> EigenvalueResult
```

Computes eigenvalues and participation factors at the problem's initial (steady-state)
operating point.

**Arguments:**
- `prob`: TDSProblem with an initialized system
- `tol`: Tolerance for classifying eigenvalues as zero (default: 1e-6)
- `verbose`: Print progress information (default: false)

### `print_eigenvalue_report`

```julia
print_eigenvalue_report(result::EigenvalueResult; io = stdout, max_modes = 50)
```

Prints the stability statistics, all modes, and the participation factors for up to
`max_modes` states.

### `write_eigenvalue_report`

```julia
write_eigenvalue_report(result::EigenvalueResult, filename::String; max_modes = 50)
```

Writes the same report to a file.

### `eigenvalue_sensitivity`

```julia
eigenvalue_sensitivity(sol_pf, dyn_models, params; config = SimConfig(),
                       formulation = PowerInjection(), addressing = NaturalAddressing(),
                       tol = 1e-6, verbose = false) -> EigenvalueSensitivity
```

Eigenvalues of the system initialized from `sol_pf` with `dyn_models`, and their
derivatives `s.dlambda[i, j]` with respect to the parameters `params`
(`(model, field, instances)` tuples). `damping_ratio_sensitivity(s)` gives the
derivatives of the damping ratios.

## Further Reading

- [Time-Domain Simulation](tds.md): TDS problem setup
- [Frequency Response](bode.md): Linearized input-output models
