# API Reference

Reference documentation for the public API of every Powerful.jl package, generated
from the docstrings in the source code.

`Powerful` is the umbrella package: `using Powerful` brings in the public API of
PowerfulCore, PowerfulPF and PowerfulTDS.
<!--GENERATED:standalone-intro-->
Model types are imported from `Powerful.Models`.

## Packages

| Page | Contents |
|:-----|:---------|
| [Powerful](api/powerful.md) | Umbrella module and the `Powerful.Models` submodule |
| [Core: Types and Variables](api/foundation.md) | `AbstractModel`, variable and residual declarations, address spaces |
| [Core: Metadata and Traits](api/metadata.md) | Model metadata, traits, groups, base quantities |
| [Core: Formulations, Units and Registry](api/formulations.md) | Formulations, data formats, per-unit kinds, bus type codes, model registry |
| [Core: DAE and Addressing](api/dae.md) | DAE vectors, address manager, indexing, residual buffers |
| [Core: Jacobians and Newton](api/ad.md) | Sparsity detection, automatic differentiation, Newton engine |
| [Core: Time Series and Injection Recording](api/recording.md) | `TDSSolution`, time series access, injection recording |
| [Core: System](api/system.md) | `SystemModel`, per-unit conversion, configuration, topology |
| [Core: Islands and Outages](api/outages.md) | Island splitting, contingency and outage types |
| [Network Models](api/models.md) | Bus, branch, transformer, load, generator and shunt models (PowerfulPF) |
| [Power Flow](api/powerflow.md) | Problems, algorithms and solutions (PowerfulPF) |
| [Power Flow Results and Reports](api/results.md) | Variable accessors, parameter system, reports (PowerfulPF) |
| [Dynamic Models](api/dynamic_models.md) | Generators, exciters, governors, renewable models, faults, DYR parsing (PowerfulTDS) |
| [Time-Domain Simulation](api/tds.md) | Setup, initialization, solving, events (PowerfulTDS) |
| [Validation and Small-Signal Analysis](api/small_signal.md) | Initialization validation, eigenvalues, linearization (PowerfulTDS) |
<!--GENERATED:api-table-->
| [Test Cases](api/cases.md) | Case loading and file access (PowerfulCases) |

## Quick Start

```julia
using Powerful
using PowerfulCases

# Load a system and solve the power flow
sys = SystemModel(parse_network(load("ieee14").raw))
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# Access results by external key
v = get_var(sol, :Bus, :v, (1,))
```

The [Index](index_page.md) lists every documented name alphabetically.
