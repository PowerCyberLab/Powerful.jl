# PowerfulCore

Core infrastructure for Powerful.jl: the type system, DAE containers, address
management, and the `SystemModel`. Most users work through `Powerful`; this package
is for developers of models, formulations, and solvers.

## Components

### Foundation (`src/foundation/`)
- Type system and model metadata
- Model traits and formulation traits
- Solution interfaces

### DAE System (`src/dae/`)
- `DAE` - Differential-algebraic equation state and residual vectors
- `AddressManager` - Maps model variables to DAE vector positions
- `ResidualBuffers` - Workspace for residual computation
- Layout strategies: `NaturalAddressing`, `AMDAddressing`, `CustomPermutation`
- Sparse Jacobians via automatic differentiation

### System (`src/system/`)
- `SystemModel` - Container for models, external IDs, and per-unit data
- Per-unit conversion with stored coefficients
- Topology analysis and island splitting
- Configuration (`SystemConfig`, `SimConfig`)
- Contingency types (`Contingency`, `BranchOutage`, `GeneratorOutage`, `Transformer2WOutage`)

## Key Types

```julia
using Powerful
using PowerfulCases
using PowerfulCore: AddressManager, NaturalAddressing, setup_with_formulation!, get_addresses, DAE

sys = SystemModel(parse_network(load("ieee14").raw))

# Assign DAE addresses for a formulation
am = AddressManager(strategy = NaturalAddressing())
setup_with_formulation!(am, sys, PowerInjection())
v_addrs = get_addresses(am, :Bus, :v)     # Positions of the bus voltages in dae.y

# DAE vectors sized from the address manager
dae = DAE(am)
length(dae.y)
```

See the [developer documentation](../docs/src/developer/overview.md) for the
architecture.

## Dependencies

- AMD.jl, KLU.jl (sparse ordering and factorization)
- StructArrays.jl (struct-of-arrays layout)
- ForwardDiff.jl, DifferentiationInterface.jl, SparseMatrixColorings.jl (automatic differentiation)
- PowerFlowData.jl (PSS/E RAW parsing)
- OrdinaryDiffEq components (time integration)
