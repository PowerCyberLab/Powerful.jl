# Supported Models

Reference list of the power flow and dynamic models, their input data and their
configuration types.

## Finding Models

Power flow and dynamic model types are collected in the `Powerful.Models` submodule,
together with their configuration types:

```@example models
using Powerful
names(Powerful.Models)
```

Import the ones you need by name:

```julia
using Powerful.Models: Bus, PQ, GENROU, TGOV1
```

Each model's docstring lists its fields, e.g. `?GENROU` in the REPL.

## Input Formats

- **Network data**: PSS/E RAW files, version 33, read with `parse_network`.
- **Dynamic data**: PSS/E DYR files, read with `parse_dynamics(dyr_file, sys)`. Records
  are matched to the models below by their DYR model name; records of other models
  are skipped with a warning listing them.

```@example models
using PowerfulCases

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
dyr = parse_dynamics(case.dyr, sys)
map(length, dyr)
```

The test cases shipped with PowerfulCases include DYR variants that exercise the
different dynamic models:

```@example models
variants(case, :psse_dyr)
```

```@example models
dyr_genrou = parse_dynamics(file(case, :dyr; variant = "genrou"), sys)
map(length, dyr_genrou)
```

## Power Flow Models

Created from the RAW file by `SystemModel(parse_network(raw_file))`:

| Model | RAW data | Description |
|-------|----------|-------------|
| `Bus` | Bus | Network bus. Holds the bus voltage variables (magnitude `v` and angle `theta` in the default formulation) |
| `PQ` | Load | Load as a ZIP combination of constant power, constant current and constant admittance. Constant-power parts are reduced smoothly below `constant_power_vmin` (a `SystemModel` keyword, default 0.7 pu) |
| `PV` | Generator (not at the swing bus) | Generator holding its active power and voltage setpoint |
| `Slack` | Generator (at the swing bus, type 3) | Swing generator holding the voltage magnitude and angle reference |
| `Branch` | Non-transformer branch | Pi-model line with charging and line shunts at both ends |
| `Transformer2W` | Transformer | Two-winding transformer with off-nominal tap and phase shift. Three-winding transformers are converted to three two-winding transformers joined at an added star-point bus |
| `Shunt` | Fixed shunt | Fixed shunt admittance (capacitor or reactor) |
| `SwitchedShunt` | Switched shunt | Switched shunt held at its initial admittance `binit` (a fixed susceptance, the physics of `Shunt`); its control (`modsw`) is not modeled, and the RAW reader warns of those in service with a control mode. Its control data and blocks are kept |
| `AreaInterchange` | Area | Area records (data only, no equations) |
| `Zone` | Zone | Zone records (data only) |
| `Owner` | Owner | Owner records (data only) |

The remaining RAW sections (two-terminal and multi-terminal DC lines, VSC DC lines, FACTS
devices, impedance correction tables, multi-section line groups, inter-area transfers, ...)
are not read. Each table is classed as power-carrying or data only (`PowerfulPF.RAW_TABLE_POWER`).
A case with records that carry power and that no model reads is refused, naming the table,
the counts and the records: a two-terminal DC line not blocked and at a nonzero setpoint, an
in-service VSC or multi-terminal DC line or FACTS device, or an impedance correction table
that an in-service transformer references would change its power flow. The records that carry none (out of service, or at a zero setpoint) and those of
data-only tables are warned of with their counts.

Bus and device keys used by `get_var` and the reports are tuples: `(bus,)` for `Bus` and `SwitchedShunt`
(one switched shunt per bus in RAW revision 33), `(bus, id)` for `PQ`, `PV`, `Slack` and `Shunt`, and `(from, to, circuit)` for
branches and transformers.

## Dynamic Models

Read from the DYR file by `parse_dynamics` and used by [time-domain simulation](tds.md)
and [Y-bus simulation](ybus.md), among other studies. Generator and converter models
replace the static `PV` or `Slack` generator with the same bus and ID; controllers
attach to the generator with the same bus and ID.

| Model | DYR name | Kind | Description |
|-------|----------|------|-------------|
| `GENCLS` | `GENCLS` | Synchronous machine | Classical model: swing equation with a constant voltage behind transient reactance |
| `GENROU` | `GENROU` | Synchronous machine | Round-rotor machine with quadratic saturation |
| `GENSAL` | `GENSAL` | Synchronous machine | Salient-pole machine with quadratic saturation |
| `GENTPJ1` | `GENTPJ1` | Synchronous machine | Round-rotor machine with quadratic saturation and the `Kis` term |
| `SEXS` | `SEXS` | Exciter | Simplified excitation system; controls a machine |
| `TGOV1` | `TGOV1` | Governor | Steam turbine governor with lead-lag; controls a machine |
| `REGCA1` | `REGCA1` | Converter | Renewable generator/converter (current source, WECC REGC_A): low voltage power logic, low-voltage active current and high-voltage reactive current management; publishes the power it delivers (`Pe`, `Qe`) |
| `REGFM_A1` | `REGFM_A1` | Converter | Droop-controlled grid-forming converter (WECC REGFM_A1, PNNL-35110): a voltage source behind a coupling reactance with P-f and Q-V droop, voltage control, overload mitigation and a reactive power limit controller; publishes the power it delivers (`Pe`, `Qe`) |
| `REECA1` | `REECA1` | Converter control | Renewable electrical control (reactive power and voltage control, current limits); controls a `REGCA1`, measuring its `Pe` and `Qe` |

A machine or converter has no rating of its own in the DYR record: it is per unit on the
MVA base `MBASE` of the generator it replaces, and Powerful converts to the system base
only where it meets the network. REGCA1 and REGFM_A1 take their base the same way as the
machines (the [converter table](../developer/miscellaneous.md) lists who converts what):
their states, `Pe` and `Qe` are on the converter's base, and their injections into the bus
are `Mbase/Sbase` times those. A rating is read when the system is built, so a
`FieldModification` of it (the generator's `mbase`) is refused.

One more dynamic model is not read from files:

| Model | Created by | Description |
|-------|------------|-------------|
| `Fault` | `BusFault` event | Balanced three-phase fault at a bus, as a shunt impedance switched on and off by the event |

## Configuration Types

Solver options are set through `SimConfig(solver = (...))`; models take no runtime
options. See [Configuration](configuration.md).

| Type | Applies to | Settings |
|------|------------|----------|
| `TDSConfig` | Time-domain simulation (`solver.TDS`) | Load conversion (`convert_pq_to_zip`, `p_zip`, `q_zip`), `flat_start`, default `abstol`, `reltol`, `maxiters`, `adjust_limits` |

System-wide settings (frequency base, `constant_power_vmin`) are keywords of
`SystemModel`.

```@example models
config = SimConfig(solver = (TDS = TDSConfig(convert_pq_to_zip = true),))
pf_sol = solve(ACPowerFlowProblem(sys), Newton())
init = setup_tds(pf_sol, dyr; config)
get_solver_config(config, :TDS, TDSConfig)
```
