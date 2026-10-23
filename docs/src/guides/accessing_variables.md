# Accessing Variables and Results

After solving a power flow, results are read in two ways:

- **`get_var` and friends** read solution *variables* (the unknowns of the power
  flow, such as bus voltages) by external key, in per-unit.
- **Report tables** (`collect_model_table`) give computed quantities such as power
  injections, load consumption, and branch flows, converted to physical units.

The examples use the IEEE 14-bus case:

```@example access
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))
sol = solve(ACPowerFlowProblem(sys), Newton())
nothing # hide
```

## The Generic Access Pattern

```julia
value = get_var(solution, :ModelName, :variable, external_key)
```

- **solution**: A `PowerFlowSolution` or `MultiIslandSolution`
- **ModelName**: Model type as a symbol (`:Bus`, `:PQ`, `:PV`, ...)
- **variable**: Variable name as a symbol (`:v`, `:theta`, ...)
- **external_key**: Key from the data file, as a tuple

Keys are tuples:

| Model | Key | Example |
|-------|-----|---------|
| `Bus` | `(bus,)` | `(14,)` |
| `PQ`, `PV`, `Slack` | `(bus, id)` | `(4, "1")` |
| `Branch`, `Transformer2W` | `(from_bus, to_bus, circuit)` | `(1, 2, "1")` |

## Bus Variables

With the default `PowerInjection` formulation, the bus variables are the voltage
magnitude `:v` (pu) and angle `:theta` (rad):

```@example access
v1 = get_var(sol, :Bus, :v, (1,))
θ14 = get_var(sol, :Bus, :theta, (14,))
println("Bus 1 voltage: $v1 pu")
println("Bus 14 angle: $(rad2deg(θ14))°")
```

Variables depend on the formulation. With `CurrentInjection`, the bus variables are
the real and imaginary voltage parts `:v_r` and `:v_i`:

```@example access
sol_ci = solve(ACPowerFlowProblem(sys, CurrentInjection()), Newton())
v_r = get_var(sol_ci, :Bus, :v_r, (14,))
v_i = get_var(sol_ci, :Bus, :v_i, (14,))
hypot(v_r, v_i)
```

Asking for a name that is not a variable of the model raises an error that lists
the available variables:

```@example access
try
    get_var(sol, :Bus, :p, (2,))
catch e
    println(sprint(showerror, e))
end
```

## Variables of Other Models

Loads, generators, and branches hold the voltages of the buses they connect to as
variables, so their terminal voltages can be read by their own keys:

```@example access
get_var(sol, :PQ, :v, (4, "1"))          # Voltage at the bus of load (4, "1")
```

```@example access
get_var(sol, :PV, :v, (2, "1"))          # Voltage at the bus of generator (2, "1")
```

```@example access
# Branch (1, 2, "1"): voltage at the from-bus and angle at the to-bus
get_var(sol, :Branch, :v_i, (1, 2, "1")), get_var(sol, :Branch, :theta_j, (1, 2, "1"))
```

## Powers and Flows: Report Tables

Power injections, load consumption, and branch flows are computed from the solution.
`collect_model_table` returns them, with each instance's key, in physical units
(MW, MVAR, degrees). The result is a `NamedTuple` of columns, which works with any
Tables.jl consumer such as DataFrames:

```@example access
using DataFrames

branches = DataFrame(collect_model_table(sol, :Branch))
```

```@example access
DataFrame(collect_model_table(sol, :PQ))
```

```@example access
DataFrame(collect_model_table(sol, :PV))
```

The slack generator is in `:Slack`. Key columns (`I`, `J`, `Ckt`, `Id`) are strings.
Select rows with DataFrame operations, for example the flow on branch 1-2:

```@example access
row = only(filter(r -> r.I == "1" && r.J == "2" && r.Ckt == "1", branches))
row["P_from (MW)"], row["P_to (MW)"]
```

The losses on a branch are the sum of the two ends:

```@example access
loss_mw = row["P_from (MW)"] + row["P_to (MW)"]
```

Use `get_available_outputs` to list the quantities a model can report. See
[Report Customization](reports.md) for units, column selection, and presets.

## Getting All Values of a Variable

### Vectors (Internal Order)

```@example access
all_voltages = get_all_vars(sol, :Bus, :v)
```

Values are in *internal order* (the order of the model arrays), not sorted by
external key. The external keys in the same order are the model's key fields, e.g.
`sys.models.Bus.i`.

### Dictionaries (External Keys)

```@example access
v_dict = get_vars_dict(sol, :Bus, :v)
v_dict[(14,)]
```

Iterate in external key order:

```@example access
for (key, voltage) in sort(collect(v_dict))
    println("Bus $(key[1]): $(round(voltage, digits = 4)) pu")
end
```

## Multi-Island Solutions

For a system with several islands, `solve` returns a `MultiIslandSolution`. The same
functions work on it and find the island that holds each key:

```@example access
sys_island = SystemModel(parse_network(file(load("ieee14"), :raw; variant = "island")))
sol_island = solve(ACPowerFlowProblem(sys_island), Newton())

get_var(sol_island, :Bus, :v, (1,)), get_var(sol_island, :Bus, :v, (14,))
```

See [Multi-Island Systems](multi_island.md).

## Direct Index Access

When the internal index is known, `get_var_by_index` skips the key lookup:

```@example access
bus_idx = sys.order.key_indices[:Bus][(14,)]    # External key → internal index
get_var_by_index(sol, :Bus, :v, bus_idx)
```

**Use case**: Loops that already work with internal indices.

## Common Patterns

### Voltage Profile Analysis

```@example access
v_dict = get_vars_dict(sol, :Bus, :v)

# Buses outside 0.95-1.05 pu
violations = filter(p -> !(0.95 <= p.second <= 1.05), v_dict)

# Lowest and highest voltage (value, key)
min_v, min_bus = findmin(v_dict)
max_v, max_bus = findmax(v_dict)
println("Lowest voltage: $min_v pu at bus $(min_bus[1])")
println("Highest voltage: $max_v pu at bus $(max_bus[1])")
```

### Angle Profile Plot

```julia
using Plots

θ_dict = get_vars_dict(sol, :Bus, :theta)
keys_sorted = sort(collect(keys(θ_dict)))
bus_nums = [k[1] for k in keys_sorted]
angles_deg = [rad2deg(θ_dict[k]) for k in keys_sorted]

plot(bus_nums, angles_deg,
     xlabel = "Bus Number", ylabel = "Angle (degrees)",
     marker = :circle, legend = false)
```

### Branch Loading Analysis

Compare branch flows with the ratings in the model data. Branch ratings (`rate_a`)
are in MVA:

```@example access
branch_model = sys.models.Branch
flows = collect_model_table(sol, :Branch)

s_from = hypot.(flows[Symbol("P_from (MW)")], flows[Symbol("Q_from (MVAR)")])
loading = 100 .* s_from ./ branch_model.rate_a

DataFrame(from = branch_model.i, to = branch_model.j, ckt = branch_model.ckt,
          s_mva = s_from, rating = branch_model.rate_a, loading_pct = loading)
```

The report table and `sys.models.Branch` list branches in the same (internal) order.

## Error Handling

A key that does not exist raises an error:

```@example access
try
    get_var(sol, :Bus, :v, (9999,))
catch e
    println(sprint(showerror, e))
end
```

Check whether a model is present with `haskey`:

```@example access
haskey(sys.models, :Transformer2W), haskey(sys.models, :Shunt)
```

## Further Reading

- [Report Customization](reports.md) - Units, columns, and presets for result tables
- [API Reference: Power Flow](../api/powerflow.md) - Variable accessors and reporting functions
