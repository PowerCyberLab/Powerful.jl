# PowerfulCA

Contingency analysis for Powerful.jl: AC contingency analysis (ACCC) with warm-started
Newton solves, and DC screening with PTDF/LODF sensitivity factors. Load it
separately with `using PowerfulCA`.

## Usage

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

## Features

- AC contingency analysis: `solve_accc`, and `solve_accc_parallel` for multithreaded runs
- Outage types: `BranchOutage`, `GeneratorOutage`, `Transformer2WOutage`, combined into
  N-1 or N-k `Contingency` definitions
- DC screening: `compute_ptdf`, `compute_lodf`, `screen_contingencies_dc`
- Voltage and branch-loading violation reports

See the [Contingency Analysis guide](../docs/src/guides/accc.md) and the
[DC Analysis guide](../docs/src/guides/dc_analysis.md).

## Dependencies

- PowerfulCore
- PowerfulPF
- KLU.jl
