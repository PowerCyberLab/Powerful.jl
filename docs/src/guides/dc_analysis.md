# DC Analysis

DC power flow and DC contingency screening with sensitivity factors (PTDF and LODF).

The DC approximation linearizes the active power flow equations:

- all bus voltage magnitudes are 1.0 pu,
- angle differences are small, so ``\sin(\theta_i - \theta_j) \approx \theta_i - \theta_j``,
- branch resistance is neglected, so the flow on a branch with reactance ``x`` is
  ``P_{ij} = (\theta_i - \theta_j) / x``,
- reactive power and voltage magnitudes are not modeled.

The result is a linear system ``B \theta = P`` in the bus angles. It is solved in a
single linear solve, which makes it useful for screening many outages quickly. It says
nothing about voltages or reactive power; use [AC power flow](powerflow.md) and
[AC contingency analysis](accc.md) for those.

## DC Power Flow

`DCPowerFlow()` is a formulation of `ACPowerFlowProblem`, next to `PowerInjection()`
and `CurrentInjection()`. The only variables are the bus angles:

```@example dc
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))

prob_dc = ACPowerFlowProblem(sys, DCPowerFlow())
sol_dc = solve(prob_dc, Newton())
```

The equations are linear, so Newton needs no more than a couple of iterations.

Read angles (in radians) with `get_var`, keyed by the bus number tuple:

```@example dc
get_var(sol_dc, :Bus, :theta, (14,))
```

Compare with the AC solution of the same system:

```@example dc
sol_ac = solve(ACPowerFlowProblem(sys), Newton())

for bus in (2, 5, 9, 14)
    θ_dc = rad2deg(get_var(sol_dc, :Bus, :theta, (bus,)))
    θ_ac = rad2deg(get_var(sol_ac, :Bus, :theta, (bus,)))
    println("Bus $bus: DC $(round(θ_dc, digits = 2))°, AC $(round(θ_ac, digits = 2))°")
end
```

The DC solution has no voltage magnitude variable: `get_var(sol_dc, :Bus, :v, ...)`
raises an error.

!!! note
    The report functions (`collect_model_table`, `print_model_table`) currently
    support AC solutions only. For a DC solution, read the angles with `get_var`.

## Contingency Screening with PowerfulCA

DC sensitivity factors and contingency screening live in the standalone
`PowerfulCA` package, which also provides [AC contingency analysis](accc.md):

```@example dc
using PowerfulCA
```

### Power Transfer Distribution Factors (PTDF)

`PTDF[l, k]` is the change of active power flow on branch `l` for a 1 pu injection at
bus `k`, withdrawn at the slack bus. `compute_ptdf` builds the full dense matrix
(branches by buses):

```@example dc
ptdf = compute_ptdf(sys)
size(ptdf.matrix)
```

- **Rows** are the in-service branches and two-winding transformers.
  `ptdf.branch_indices[l]` gives `(model, instance)` for row `l`, e.g.
  `(:Branch, 3)` is the third entry of `sys.models.Branch`.
- **Columns** are the buses in the order of `sys.models.Bus`. `ptdf.slack_bus` is the
  position of the slack bus; its column is zero.

For example, the flow change on every branch when 0.1 pu (10 MW on a 100 MVA base) is
injected at bus 14 and withdrawn at bus 2:

```@example dc
k_from = findfirst(==(14), sys.models.Bus.i)
k_to = findfirst(==(2), sys.models.Bus.i)
Δflow = 0.1 .* (ptdf.matrix[:, k_from] .- ptdf.matrix[:, k_to])

l = argmax(abs.(Δflow))
model, idx = ptdf.branch_indices[l]
m = getproperty(sys.models, model)
println("Largest change on $model $(m.i[idx])-$(m.j[idx]): $(round(Δflow[l], digits = 4)) pu")
```

For large systems, where the dense L×N matrix does not fit in memory, `LazyPTDF` keeps
only the sparse factorization and computes rows, columns or single entries on demand:

```@example dc
lptdf = LazyPTDF(sys)

row = get_ptdf_row(lptdf, 1)        # sensitivities of branch 1 to every bus
col = get_ptdf_column(lptdf, 14)    # sensitivities of every branch to bus 14
get_ptdf_element(lptdf, 1, 14) ≈ ptdf.matrix[1, 14]
```

### Line Outage Distribution Factors (LODF)

`LODF[m, l]` is the fraction of the pre-outage flow on branch `l` that moves to branch
`m` when `l` is taken out. The diagonal is -1 (the outaged branch loses its flow):

```@example dc
lodf = compute_lodf(sys)
size(lodf.matrix)
```

Rows and columns use the same branch ordering as the PTDF (`lodf.branch_indices`).
Outaging a branch whose removal splits the network (a radial branch) has no LODF: its
column is `NaN`.

Post-outage flows follow from the base-case flows:
``P^{\text{post}}_m = P_m + \text{LODF}[m, l] \, P_l``. The screening function below
applies this to whole contingency lists; [Estimating a Single Outage](#Estimating-a-Single-Outage)
shows the direct use.

### Screening a Contingency List

`screen_contingencies_dc` solves the base-case DC power flow once and then estimates
the post-contingency flows of every contingency from the sensitivity factors, without
solving a power flow per contingency. Contingencies are defined as for
[AC contingency analysis](accc.md):

```@example dc
br = sys.models.Branch
contingencies = [
    Contingency("BR_$(br.i[k])_$(br.j[k])", [BranchOutage(br.i[k], br.j[k], br.ckt[k])])
    for k in eachindex(br)
]
push!(contingencies, Contingency("GEN_2", [GeneratorOutage(2, "1")]))

result = screen_contingencies_dc(sys, contingencies; loading_threshold = 70.0)
print_dc_screening_report(result; max_violations = 3)
```

How each outage type is handled:

- **Branch and transformer outages** (`BranchOutage`, `Transformer2WOutage`) remove
  the branch with the LODF. Several branch outages in one contingency are removed
  together (multiple-outage LODF).
- **Generator outages** (`GeneratorOutage`) remove the unit's injection; the slack bus
  picks up the lost generation.
- An outage that splits the network is flagged as radial (`r.is_radial`); its flows
  are not computed and it is listed as critical.

In the IEEE 14-bus system, bus 8 is connected only through the 8-7 transformer:

```@example dc
radial = screen_contingencies_dc(sys, [Contingency("XF_8_7", [Transformer2WOutage(8, 7, "1")])])
radial.results[1].is_radial
```

Keyword arguments:

| Keyword | Default | Meaning |
|---------|---------|---------|
| `loading_threshold` | `100.0` | Report branches loaded above this percentage of their rating |
| `use_lazy_ptdf` | `false` | Compute PTDF columns on demand instead of the full matrix (large systems) |
| `verbose` | `false` | Print progress and timing |

Loading is `|flow| / rating × 100`, with the rating from RATE A of the RAW file
(branches) or RATA1 (transformers). Branches with a zero rating are not checked.

### Working with the Result

`screen_contingencies_dc` returns a `DCContingencyScreenResult`:

| Field | Content |
|-------|---------|
| `base_flows` | Base-case DC flow of every branch (pu on the system base) |
| `branch_indices` | `(model, instance)` of every branch row |
| `branch_limits` | Rating of every branch row (pu) |
| `results` | One `DCScreeningResult` per contingency |
| `critical_contingencies` | Indices (into `results`) of contingencies with violations or a network split |

Each `DCScreeningResult` has `name`, `max_loading` (%), `is_radial`,
`outaged_branches`, `outaged_generators` and `violations`. A `BranchViolation` has
`branch_idx` (row in the branch ordering), `model_name`, `instance_idx`, `flow` (pu),
`limit` (pu) and `loading` (%).

Listing the violations with bus numbers:

```@example dc
for idx in result.critical_contingencies
    r = result.results[idx]
    for v in r.violations
        m = getproperty(sys.models, v.model_name)
        k = v.instance_idx
        println(rpad(r.name, 10), " $(v.model_name) $(m.i[k])-$(m.j[k]): ",
                round(v.loading, digits = 1), "%")
    end
end
```

### Estimating a Single Outage

`estimate_post_outage_flows` gives the post-outage flows of all branches for the outage
of one branch row, from the base-case flows (`result.base_flows`, pu) and the LODF:

```@example dc
l = 1                                    # branch row to take out
model, idx = result.branch_indices[l]    # (:Branch, 1)
post = estimate_post_outage_flows(result.base_flows, lodf, l)

post[l]   # the outaged branch carries no flow
```

`estimate_post_outage_flows(result.base_flows, lptdf, l)` does the same with a
`LazyPTDF`, computing only the LODF column it needs.

## From DC Screening to AC Contingency Analysis

DC screening is fast but approximate: it neglects losses and reactive power, and it
cannot detect voltage violations or non-convergence. A common workflow screens a long
contingency list with DC factors and runs the full AC analysis only on the
contingencies it flags:

```@example dc
critical = contingencies[result.critical_contingencies]
ac_result = solve_accc(sys, critical)

n_converged = count(r -> r.converged, ac_result.contingencies)
println("AC analysis: $n_converged of $(length(critical)) contingencies converged")
```

The AC result is checked for voltage and thermal violations as described in
[Contingency Analysis](accc.md). Contingencies that DC screening passes can still have
voltage problems in AC, so choose the loading threshold of the screen with margin.

## Limitations

- The DC model neglects voltage magnitudes, reactive power and losses.
- Generator outages are balanced by the slack bus only.
- Thermal checks use RATE A (branches) and RATA1 (transformers).
- Report functions (`collect_model_table`) do not support DC power flow solutions yet.

## See Also

- [Power Flow Analysis](powerflow.md): AC formulations and algorithms
- [Contingency Analysis](accc.md): AC contingency analysis and violation reports
