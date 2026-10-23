# Report Customization

This guide shows how to generate result tables and text reports from power flow
solutions, with control over units, formatting, and column selection.

## The Basic Workflow

```@example reports
using Powerful
using PowerfulCases
using DataFrames

sys = SystemModel(parse_network(load("ieee14").raw))
sol = solve(ACPowerFlowProblem(sys), Newton())

bus_data = collect_model_table(sol, :Bus)
DataFrame(bus_data)
```

`collect_model_table(sol, :Model)` returns a `NamedTuple` of columns: the instance
key columns, then the model's outputs, then its parameters. It works with any
Tables.jl consumer (DataFrames, CSV, ...). The default configuration gives:

- Power in MW/MVAR
- Voltage magnitudes in per-unit
- Angles in degrees
- Column names with units

`print_model_table` prints the same table, rounded for display:

```@example reports
print_model_table(sol, :Branch)
```

## Understanding Report Configuration

Reports are controlled by `ReportConfig`, which has three parts:

```@example reports
config = ReportConfig(
    units = UnitPreferences(power = :physical, voltage = :pu),
    formatting = FormattingPreferences(float_precision = 4),
    columns = ColumnPreferences(include_id = true, use_descriptions = true),
)
nothing # hide
```

Pass it with the `config` keyword:

```julia
collect_model_table(sol, :Bus; config = config)
print_model_table(sol, :Bus; config = config)
```

### 1. Unit Preferences

```julia
UnitPreferences(
    power = :physical,      # :physical (MW/MVAR) or :pu
    voltage = :pu,          # :pu or :kv
    angle = :deg,           # :deg or :rad
    impedance = :pu,        # :pu or :ohm
    base_mva = 100.0,       # System base power
)
```

The values shown are the defaults.

```@example reports
# Voltages in kV
config = ReportConfig(units = UnitPreferences(voltage = :kv))
DataFrame(collect_model_table(sol, :Bus; config = config))
```

```@example reports
# Power in per-unit, angles in radians
config = ReportConfig(units = UnitPreferences(power = :pu, angle = :rad))
DataFrame(collect_model_table(sol, :Branch; config = config))
```

### 2. Formatting Preferences

```julia
FormattingPreferences(
    float_precision = 4,    # Decimal places
)
```

Formatting applies to printed tables (`print_model_table`, text reports).
`collect_model_table` always returns full precision for data export.

### 3. Column Preferences

```julia
ColumnPreferences(
    use_descriptions = true,      # Descriptive names vs. field names
    include_units = true,         # Append units, e.g. "(MW)"
    column_style = :short,        # :compact, :short, :verbose
    include_id = true,            # Include the key columns
    include_vars = nothing,       # Vector{Symbol}: select specific outputs
    exclude_vars = Symbol[],      # Exclude these outputs
    include_params = true,        # Include model parameters
    exclude_params = Symbol[],    # Exclude these parameters
)
```

**Column styles:**

```@example reports
for style in (:compact, :short, :verbose)
    style_config = ReportConfig(columns = ColumnPreferences(column_style = style))
    println(style, ": ", keys(collect_model_table(sol, :Branch; config = style_config)))
end
```

## Preset Configurations

Five presets cover common cases:

| Preset | Units | Columns |
|--------|-------|---------|
| `Default` | MW/MVAR, pu voltage, degrees | Excludes the input voltages `vm`, `va` |
| `AllPerUnit` | Everything in per-unit, radians | All |
| `AllPhysical` | MW/MVAR, kV, degrees | All |
| `Minimal` | Default units | Outputs only, no parameters |
| `Compact` | Default units | All, 2 decimal places when printed |

`collect_model_table` and `print_model_table` use `Default` unless another
configuration is given:

```@example reports
DataFrame(collect_model_table(sol, :Bus; config = AllPhysical))
```

```@example reports
DataFrame(collect_model_table(sol, :PQ; config = Minimal))
```

```@example reports
print_model_table(sol, :PV; config = Compact)
```

## Selecting Outputs

### Discovering Available Outputs

```@example reports
using Powerful.Models: Bus, Branch, PQ, PV

get_available_outputs(Branch)      # All outputs a model can report
```

```@example reports
get_default_outputs(Branch)        # Outputs shown by default
```

```@example reports
get_available_outputs(PQ), get_available_outputs(PV)
```

Branch losses are available but not shown by default.

### Include Only Specific Outputs

```@example reports
config = ReportConfig(columns = ColumnPreferences(include_vars = [:p_from, :p_to, :p_loss]))
DataFrame(collect_model_table(sol, :Branch; config = config))
```

### Include Everything Available

```@example reports
config = ReportConfig(columns = ColumnPreferences(include_vars = get_available_outputs(Branch)))
keys(collect_model_table(sol, :Branch; config = config))
```

### Exclude Specific Outputs

```@example reports
config = ReportConfig(columns = ColumnPreferences(exclude_vars = [:theta]))
keys(collect_model_table(sol, :Bus; config = config))
```

## Working with Parameters

Model parameters with unit metadata are included after the outputs. For loads, these
are the constant-power, constant-current, and constant-admittance components:

```@example reports
keys(collect_model_table(sol, :PQ))
```

### Exclude All Parameters

```@example reports
config = ReportConfig(columns = ColumnPreferences(include_params = false))
keys(collect_model_table(sol, :PQ; config = config))
```

### Exclude Specific Parameters

The bus parameters `vm` and `va` are the voltages from the case file (the initial
guess). `Default` excludes them; other configurations include them unless excluded:

```@example reports
keys(collect_model_table(sol, :Bus; config = AllPhysical))
```

```@example reports
config = ReportConfig(AllPhysical; columns = ColumnPreferences(exclude_params = [:vm, :va]))
keys(collect_model_table(sol, :Bus; config = config))
```

## Modifying a Configuration

`ReportConfig(base; kwargs...)` copies a configuration and replaces the given parts:

```@example reports
config = ReportConfig(AllPhysical;
                      formatting = FormattingPreferences(float_precision = 2))
print_model_table(sol, :Bus; config = config)
```

## Exporting Reports

### To DataFrame

```@example reports
df = DataFrame(collect_model_table(sol, :Bus))

# Buses below 1.0 pu, sorted by voltage
sort(filter(row -> row["V (pu)"] < 1.0, df), "V (pu)")
```

### To CSV

```julia
using CSV

CSV.write("bus_results.csv", collect_model_table(sol, :Bus; config = AllPhysical))
CSV.write("branch_flows.csv", collect_model_table(sol, :Branch))
```

### Full Text Report

`build_powerflow_report` collects the summary, convergence information, and the bus,
generator, load, and branch tables; `format_text` renders it:

```@example reports
report = build_powerflow_report(sol, Compact)
print(format_text(report; config = Compact))
```

The command-line tool `./bin/powerful pf case.raw` writes this report (see the
[CLI documentation](https://github.com/cuihantao/Powerful.jl/blob/main/bin/README.md)).

## Multi-Island Systems

For a `MultiIslandSolution`, `collect_model_table`, `print_model_table` and
`build_powerflow_report` merge the tables of all islands:

```@example reports
sys_island = SystemModel(parse_network(file(load("ieee14"), :raw; variant = "island")))
sol_island = solve(ACPowerFlowProblem(sys_island), Newton())

DataFrame(collect_model_table(sol_island, :Bus))
```

## Formulations

Reports read bus voltages through `bus_voltages(sol)`, which gives the magnitude and
angle of every bus in any formulation, so a `CurrentInjection` solution reports the
same tables as a `PowerInjection` one. A `DCPowerFlow` solution reports what the DC
model solves: bus angles and active power; the voltage magnitude and reactive power
columns are left out. Devices that are out of service report no flow or consumption.

## Unit Conversion

Reports convert values in two ways:

- **Model parameters** use the per-unit coefficients stored when the system was
  built (`sys.pu_coeffs`): physical = per-unit / coefficient. This is the exact
  inverse of the conversion at load time.
- **Computed outputs** (flows, injections, angles) are converted by their unit type:
  power × `base_mva` for MW/MVAR, radians to degrees, and per-unit voltage × the bus
  base voltage (`Bus.basekv`) for kV.

## Validation and Error Messages

The report system validates the configuration and lists the valid names:

```@example reports
config = ReportConfig(columns = ColumnPreferences(include_vars = [:v, :theta, :nonexistent_var]))
try
    collect_model_table(sol, :Bus; config = config)
catch e
    println(sprint(showerror, e))
end
```

```@example reports
# Branch outputs requested from the Bus model
config = ReportConfig(columns = ColumnPreferences(include_vars = [:p_from, :q_from]))
try
    collect_model_table(sol, :Bus; config = config)
catch e
    println(sprint(showerror, e))
end
```

## Summary

| Feature | Configuration | Example |
|---------|--------------|---------|
| **Units** | `UnitPreferences` | `power = :physical, voltage = :kv` |
| **Formatting** | `FormattingPreferences` | `float_precision = 2` |
| **Column selection** | `ColumnPreferences` | `include_vars = [:p_from, :p_to]` |
| **Presets** | Ready-made configurations | `Default`, `AllPerUnit`, `AllPhysical`, `Minimal`, `Compact` |
| **Discovery** | Available outputs | `get_available_outputs(Branch)` |
| **Export** | Tables.jl interface | `DataFrame(table)`, `CSV.write(file, table)` |

## Next Steps

- **[Accessing Results](accessing_variables.md)** - Reading solution variables by key
- **[API Reference: Power Flow](../api/powerflow.md)** - Reporting function reference
