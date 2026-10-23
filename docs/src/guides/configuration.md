# Configuration

Powerful.jl separates what the system *is* from how it is simulated:

- **System settings** (frequency base, MVA base, constant-power load breakpoint) belong
  to the `SystemModel` (`sys.settings`, a `SystemConfig`), so every analysis of the
  system uses the same values.
- **`SimConfig`** holds the solver settings. Models take no runtime options: their data
  set what they do.

## Quick Start

```@example config
using Powerful
using PowerfulCases

# Simulation settings, programmatically (or from a TOML file with load_config)
config = SimConfig(solver = (TDS = TDSConfig(abstol = 1e-8),))

# System settings are set when the system is built
case = load("ieee14")
sys = SystemModel(parse_network(case.raw); frequency_base = 50.0)
pf_sol = solve(ACPowerFlowProblem(sys), Newton())

# The config enters at TDS setup and travels with the problem
init = setup_tds(pf_sol, parse_dynamics(case.dyr, sys); config)
prob = TDSProblem(init; tspan = (0.0, 10.0))
```

## TOML Format

A configuration file has `[solver.<Solver>]` sections. Arrays become tuples, e.g. for
the ZIP load fractions:

```toml
[solver.TDS]
convert_pq_to_zip = true
p_zip = [0.3, 0.3, 0.4]  # 30% Z, 30% I, 40% P
q_zip = [0.3, 0.3, 0.4]
abstol = 1e-8
```

```julia
config = load_config("simulation.toml")
```

`load_config_from_string` parses the same format from a string:

```@example config
config = load_config_from_string("""
[solver.TDS]
convert_pq_to_zip = true
p_zip = [0.3, 0.3, 0.4]
q_zip = [0.3, 0.3, 0.4]
abstol = 1e-8
""")
```

## System Settings

System-wide settings are owned by the `SystemModel`. `SystemModel(case)` takes the
frequency base from the case header and uses the per-unit base (`system_base`) as the
MVA base; both can be set when the system is built:

```@example config
sys = SystemModel(parse_network(case.raw); frequency_base = 50.0)
sys.settings.frequency_base, sys.settings.omega_base   # 50 Hz, 2π · 50 rad/s
```

### Constant-Power Load Breakpoint

`constant_power_vmin` (default 0.7 pu) is the voltage below which the constant-power
part of PQ loads (`pl`, `ql`) is reduced. At or above it the load draws constant
power; below it the drawn power is `pl·f(V)` with the raised cosine
`f(V) = sin(πV/(2·vmin))²`, which is zero at zero voltage and has a continuous
derivative at the breakpoint. Every solver uses it:

```@example config
sys = SystemModel(parse_network(case.raw); constant_power_vmin = 0.8)
sys.settings.constant_power_vmin
```

Systems derived from it (TDS, islands, multi-island copies) keep the same settings. A
`[system]` section in a TOML config is an error: system settings are keywords of
`SystemModel`.

## Solver Settings

### TDSConfig

Time-domain simulation settings:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `convert_pq_to_zip` | Bool | false | Convert PQ loads to ZIP loads at TDS setup |
| `p_zip` | NTuple{3,Float64} | (0, 0, 1) | Active power (Z, I, P) fractions |
| `q_zip` | NTuple{3,Float64} | (0, 0, 1) | Reactive power (Z, I, P) fractions |
| `flat_start` | Bool | true | Reinitialize after events from V = 1, θ = 0 (false: pre-fault voltages) |
| `abstol` | Float64 | 1e-6 | Absolute tolerance of the DAE solver |
| `reltol` | Float64 | 1e-3 | Relative tolerance of the DAE solver |
| `maxiters` | Int | 100 | Maximum Newton iterations of the post-event reinitialization |
| `adjust_limits` | Bool | true | Move a limit to an initial value outside it (false: error) |

`solve(::TDSProblem)` uses these as defaults; explicit `abstol`/`reltol` keywords and
an explicit `initializealg` take precedence.

**Initial values outside their limits:**
The dynamic models initialize from the power flow, which knows nothing of their limits:
a governor can start with its valve beyond `Vmax` when the power flow dispatches the
unit above it. `setup_tds` and any simulation mode built on it check every limited
variable at its initial value: those a
model declares as a `VariableLimit` (the anti-windup limits of SEXS and TGOV1 and the
REECA1 power order, which their automata enforce). With `adjust_limits` (the
default) the violated limit moves to the initial value, so the simulation starts in
steady state with the device at its limit, and a warning lists each adjustment; with
`adjust_limits = false` the violations are an error. Clamping the value instead would
start the simulation away from the power flow solution. Limits written inside a
model's equations rather than declared as `VariableLimit`s (the REECA1 limits other than its
power order) are not checked; REGCA1 refuses a start outside its LVPL and its
high-voltage management instead, naming what the start requires.

**ZIP model:**
The tuple `(Z, I, P)` splits each load's consumption at its power flow voltage into:
- **Z** (constant impedance): power varies with V²
- **I** (constant current): power varies with V
- **P** (constant power): power independent of voltage

Each fraction lies in [0, 1] and the three sum to 1, so the conversion keeps the
load's consumption at the operating voltage. The conversion modifies the loads of the
power flow system passed to `setup_tds`.

**Example: 100% constant-impedance loads**
```julia
config = SimConfig(
    solver = (TDS = TDSConfig(
        convert_pq_to_zip = true,
        p_zip = (1.0, 0.0, 0.0),  # 100% Z
        q_zip = (1.0, 0.0, 0.0)
    ),)
)
```

### Accessing Solver Settings

```@example config
config = SimConfig(solver = (TDS = TDSConfig(convert_pq_to_zip = true),))
get_solver_config(config, :TDS, TDSConfig)
```

## Validation

Settings are validated when they are loaded or constructed.

**Typo detection:**
```@example config
try
    load_config_from_string("""
    [solver.TDS]
    convert_pq_to_zi = true
    """)
catch err
    showerror(stdout, err)
end
```

**ZIP validation:**
```@example config
try
    TDSConfig(p_zip = (0.5, 0.5, 0.5))   # sums to 1.5
catch err
    showerror(stdout, err)
end
```

## Complete Example

```@example config
# Build the configuration (or load it with load_config("simulation.toml"))
config = SimConfig(
    solver = (TDS = TDSConfig(
        convert_pq_to_zip = true,
        p_zip = (1.0, 0.0, 0.0),
        q_zip = (1.0, 0.0, 0.0)
    ),)
)

# System settings, e.g. a 50 Hz grid, are set on the system
sys = SystemModel(parse_network(case.raw); frequency_base = 50.0)
pf_sol = solve(ACPowerFlowProblem(sys), Newton())
dyr_data = parse_dynamics(case.dyr, sys)

init = setup_tds(pf_sol, dyr_data; config)
prob = TDSProblem(init; tspan = (0.0, 10.0))
sol = solve(prob)
```

See [Developer: Models](../developer/models.md) for the model implementation guide.
