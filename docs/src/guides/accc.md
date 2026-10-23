# AC Contingency Analysis (ACCC)

Run N-1 and N-k contingency analysis with full AC power flows to assess system
security.

## Basic Usage

```@example accc
using Powerful
using PowerfulCA
using PowerfulCases

# Load system
sys = SystemModel(parse_network(load("ieee14").raw))

# Define contingencies
contingencies = [
    Contingency("LINE_1_5", [BranchOutage(1, 5, "1")]),
    Contingency("LINE_2_3", [BranchOutage(2, 3, "1")]),
    Contingency("GEN_3", [GeneratorOutage(3, "1")]),
]

# Run ACCC
result = solve_accc(sys, contingencies)

# Print violation report
report = print_accc_report(result);
nothing # hide
```

## Defining Contingencies

A `Contingency` has a name and a list of outages that are applied together.

### Branch Outages

Open a transmission line (`Branch`):

```julia
# BranchOutage(from_bus, to_bus, circuit_id)
BranchOutage(1, 2, "1")   # Line from bus 1 to bus 2, circuit "1"
BranchOutage(2, 4, "1")   # Line from bus 2 to bus 4
```

### Transformer Outages

Open a two-winding transformer (`Transformer2W`):

```julia
# Transformer2WOutage(from_bus, to_bus, circuit_id)
Transformer2WOutage(4, 9, "1")
```

### Generator Outages

Trip a generator:

```julia
# GeneratorOutage(bus, unit_id)
GeneratorOutage(3, "1")   # Generator at bus 3, unit "1"
```

The keys must match the data exactly: an outage of an element that does not exist is
an error.

### N-1 Contingencies

Single element outage:

```julia
contingencies = [
    Contingency("LINE_1_5", [BranchOutage(1, 5, "1")]),
    Contingency("GEN_3", [GeneratorOutage(3, "1")]),
]
```

### N-2 Contingencies

Multiple simultaneous outages:

```julia
# Combined branch + generator outage
Contingency("N2_LINE_GEN", [
    BranchOutage(1, 5, "1"),
    GeneratorOutage(3, "1")
])

# Double line outage
Contingency("DOUBLE_LINE", [
    BranchOutage(1, 5, "1"),
    BranchOutage(2, 3, "1")
])
```

## Understanding Results

`solve_accc` returns a NamedTuple with the base case solution (`base`) and one result
per contingency (`contingencies`), each with `name`, `converged`, `iterations`,
`residual_norm`, the solution vector `y`, `out_of_service` and `singular`.

`out_of_service` holds the instances the contingency switched off, by model: the
outaged devices and, when the outages cut buses off from every slack bus, those
buses (`:Bus`) and every device on them. A de-energized bus keeps its base-case
voltage in `y`. A contingency whose Jacobian is singular has `singular = true`,
`converged = false` and `y = nothing`; the other contingencies still run.

### Convergence Check

```@example accc
println("Base case: ", result.base.converged ? "OK" : "FAILED")

for r in result.contingencies
    status = r.converged ? "OK" : "FAILED"
    println("$(r.name): $status in $(r.iterations) iterations")
end
```

### Violation Report

`print_accc_report` reports:

- **Voltage violations**: energized buses outside the voltage limits
- **Branch overloads**: in-service lines above their thermal rating
- **Isolated buses**: buses a contingency de-energized
- **Non-converged cases**: contingencies that failed to solve

```julia
report = print_accc_report(result;
    vmin = 0.95,       # Low voltage limit (pu)
    vmax = 1.05,       # High voltage limit (pu)
    rating = :rate_a   # Rating set (:rate_a, :rate_b, :rate_c)
)
```

The report returns the violations as data:

- `report.voltage_violations`: `contingency`, `bus`, `basekv`, `v_cont`, `v_base`,
  `v_min`, `v_max`, `violation_type` (`"LOW"` or `"HIGH"`)
- `report.branch_violations`: `contingency`, `from_bus`, `to_bus`, `circuit`,
  `flow_mva`, `rating_mva`, `pct_loading`, `base_flow_mva`
- `report.isolated_buses`: `contingency`, `bus`, `basekv`
- `report.nonconverged`: `name`, `iterations`, `residual`, `singular`

## Example: Full N-1 Branch Screening

```@example accc
# Build an N-1 contingency for every line
branches = sys.models.Branch
contingencies = Contingency[]
for k in eachindex(branches.i)
    i, j, ckt = Int(branches.i[k]), Int(branches.j[k]), String(branches.ckt[k])
    push!(contingencies, Contingency("BR_$(i)_$(j)_$(ckt)", [BranchOutage(i, j, ckt)]))
end

result = solve_accc(sys, contingencies)
report = print_accc_report(result; vmin = 0.95, vmax = 1.05)
nothing # hide
```

```@example accc
# Worst overloads
for b in sort(report.branch_violations; by = b -> b.pct_loading, rev = true)
    println("$(b.contingency): $(b.from_bus)-$(b.to_bus) at $(round(b.pct_loading, digits = 1))%")
end
```

```@example accc
# Number of non-converged contingencies
count(r -> !r.converged, result.contingencies)
```

## How It Works

1. **Solve the base case**: full AC power flow, used as the warm start
2. **For each contingency:**
   - Mark the outaged elements out of service (the system itself is not modified)
   - Solve with Newton's method from the base case solution
   - Store the solution vector
3. **Report**: check voltages and branch flows of the stored solutions

## Parallel Execution

For large contingency lists, use multi-threaded solving:

```julia
using Powerful
using PowerfulCA
using PowerfulCases

sys = SystemModel(parse_network(load("synthetic/ACTIVSg2000").raw))
contingencies = [...]  # Your contingency list

# Parallel solve using all available threads
result = solve_accc_parallel(sys, contingencies; verbose = true)
```

Start Julia with several threads (`julia --threads=auto` or `JULIA_NUM_THREADS`).
Each thread owns its own power flow workspace. For a handful of contingencies the
sequential `solve_accc` can be faster.

Both solvers take the keywords `abstol` (default 1e-8), `maxiters` (default 20),
`damping` (default 1.0) and `verbose`.

## Limitations

- **Islands**: an island left without a slack bus is de-energized, not solved on its
  own (its generation is lost with it).
- **Branch overloads** are checked for lines (`Branch`); transformer loading is not
  reported.

For fast linear screening with PTDF/LODF sensitivities before running the AC
analysis, see [DC Analysis](dc_analysis.md).

## Further Reading

- [DC Analysis](dc_analysis.md): PTDF, LODF and DC contingency screening
- [Power Flow Analysis](powerflow.md): Formulations and algorithms
- [Accessing Results](accessing_variables.md): Working with solution data
