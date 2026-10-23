# System Architecture Overview

A guide to the architecture of Powerful.jl for contributors.

## The Big Picture

**The workflow:**

1. **Parse** - PSS/E RAW file → `SystemModel` (data only)
2. **Setup** - Create `ACPowerFlowProblem` (picks formulation and variable layout)
3. **Solve** - Newton iterations until convergence
4. **Extract** - Results addressed by device IDs from the data file

**Key principle:** Data structures don't know about solvers. Solvers don't know about topology.

## Package Structure (Monorepo)

Powerful.jl is organized as a monorepo with sub-packages:

| Package | Purpose | Key Contents |
|---------|---------|--------------|
| `PowerfulCore/` | Foundation infrastructure | Types, traits, metadata, DAE, `AddressManager`, `SystemModel` |
| `PowerfulPF/` | Power flow | Bus, PQ, PV, Slack, Branch, Transformer models; Newton and fast-decoupled solvers |
| `PowerfulTDS/` | Time-domain simulation | GENCLS, GENROU, GENSAL, GENTPJ1, SEXS, TGOV1, REGCA1, REECA1 models; TDS solver, events |
| `PowerfulYbus/` | Y-bus time-domain simulation | Sparse network solve with trapezoidal machine integration |
| `PowerfulCA/` | Contingency analysis | ACCC solvers, DC screening (PTDF/LODF), violation reporting |
| `PowerfulCases/` | Test cases | IEEE 14, 39, 118 bus and other cases |
| `Powerful/` | Umbrella package | Curated public API |

**Dependencies:**

```
PowerfulCore → PowerfulPF → PowerfulTDS → Powerful (also depends on PowerfulCases)
                    │            └→ PowerfulYbus  (standalone)
                    └→ PowerfulCA   (standalone)
```

`Powerful` depends on `PowerfulCore`, `PowerfulPF`, `PowerfulTDS` and `PowerfulCases`.
`PowerfulCA` and `PowerfulYbus` are standalone packages that `Powerful` does not
re-export; load them explicitly (e.g. `using PowerfulCA`). Other standalone packages
extend the same layers without changing them.

**Import patterns:**
```julia
# User code - use the umbrella package
using Powerful
using PowerfulCases
using PowerfulCA  # Contingency analysis (standalone)

# Access model types via the Models submodule
using Powerful.Models: Bus, PQ, PV, TGOV1

# Developer code - import from the defining package
using PowerfulCore: model_metadata, AddressManager
using PowerfulPF: setup_powerflow_newton
```

## Architectural Layers

Layers build bottom-up. Each layer knows nothing about layers above it.

### Foundation Layer (`PowerfulCore/src/foundation/`)

**Defines:**
- `AbstractModel` - Base type for all device models
- `ModelMetadata` - The variables, residuals and parameters of each model
- `FormulationTrait` - `PowerInjection`, `CurrentInjection`, `DCPowerFlow`
- Address types - `AlgebVar`, `StateVar`, `AlgebRes`, `StateRes`, `ObservableVar`

**Changes rarely** - These are fundamental abstractions.

### Indexing Layer (`PowerfulCore/src/dae/`)

**Maps between three coordinate systems:**

1. **External IDs** → Bus 4001, 4002 (user-facing, from file)
2. **Internal indices** → 1, 2, 3... (contiguous array positions)
3. **DAE addresses** → y\[42\], g\[58\] (equation vector positions)

**Two mappers:**
- `IndexOrder` - External → Internal (built by `ext2int`)
- `AddressManager` - Internal → DAE

**Why separate:** Different formulations have different variables, so the
`AddressManager` is created per problem rather than stored in system data. See
[Three-Layer Indexing](indexing.md).

### Model Layer (`PowerfulPF/src/models/`, `PowerfulTDS/src/models/`)

**Each device type** (PQ load, PV generator, Branch):
- Struct definition (pure data)
- Metadata (variables, residuals, parameters)
- Residual kernels (compute power injections or currents)

**No cross-references in struct fields.** Models store external IDs of the devices they
connect to and declare dependencies via metadata. The system resolves them during setup.

### System Layer (`PowerfulCore/src/system/`)

**`SystemModel` contains:**
- Model collections (StructArrays, per-unit)
- External parameters retrieved from other models
- External-to-internal index mappings (`IndexOrder`)
- Per-unit conversion coefficients and the original parsed data
- System-wide settings (`SystemConfig`: frequency, MVA base, constant-power breakpoint)

**No computational state** - No DAE vectors, no Jacobians. Those live in workspaces.

### Problem Layer (`PowerfulPF/src/solvers/`)

**`ACPowerFlowProblem` specifies:**
- Which `SystemModel` to solve
- Which formulation (`PowerInjection`, `CurrentInjection` or `DCPowerFlow`)
- Variable ordering strategy (`addressing=` keyword: `NaturalAddressing()`, `AMDAddressing()`, `CustomPermutation(perm)`)
- Autodiff backend (`backend=` keyword)

**Creates the AddressManager** - Assigns every variable and residual to a DAE vector position.

### Solver Layer (`PowerfulPF/src/solvers/`, `PowerfulTDS/src/solvers/`)

**Sets up workspaces and runs iterations:**

1. Create the autodiff cache (`create_adcache`), which builds the `SystemBuffers`
2. Prepare the sparse Jacobian (sparsity pattern from model structure, coloring)
3. **Iteration loop:**
   - Evaluate residuals
   - Compute Jacobian
   - Solve linear system
   - Update state
   - Check convergence

**Delegates physics to workspaces** - The solver does numerical methods only.

## Core Design Principles

### Separation of Data from Computation

- `SystemModel` = data only, no computation
- Workspaces = computation, read from system
- **Benefit:** Create multiple problems (different formulations) from the same system

### Metadata-Driven Modeling

- Models declare their interface via `ModelMetadata`
- The system introspects metadata to allocate variables
- **Benefit:** Add new models without modifying core code

### Type Stability

- Concrete types and type parameters throughout
- NamedTuples for hot-path data in `ResidualBuffers` (`var_vals`, `res_vals`)
- `Vector{Function}` + `Dict` in `SystemBuffers` (avoids type explosion with many models)
- Model structs are parametric (`T<:Number` fields) so they work with `ForwardDiff.Dual`

### Pre-allocated Buffers

- Workspaces are allocated once during setup
- Residual evaluation fills existing buffers
- Jacobian evaluation reuses the prepared sparsity pattern and coloring

### Automatic Differentiation Backends

Powerful.jl uses [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
to compute sparse Jacobians with pluggable AD backends. The backend types
`ForwardDiffBackend()` (default) and `EnzymeBackend()` are defined in
`PowerfulCore/src/dae/ad_backend.jl` and passed with the `backend` keyword:

```julia
using PowerfulCore: EnzymeBackend, ForwardDiffBackend

prob = ACPowerFlowProblem(sys; backend = EnzymeBackend())
tds_prob = TDSProblem(init; tspan = (0.0, 10.0), backend = ForwardDiffBackend())
```

| Aspect | `ForwardDiffBackend()` | `EnzymeBackend()` |
|--------|------------------------|-------------------|
| **Mechanism** | Dual number propagation | LLVM-level differentiation with shadow arrays |
| **Workspace numeric type** | `Dual` workspaces cached per Dual type | `Float64` throughout |
| **Cache type** | `ForwardDiffCache` | `EnzymeCache` |

**ForwardDiff** propagates derivatives by replacing `Float64` with `Dual{Tag,Float64,N}`:
- Each `Dual` carries a value and N partial derivatives
- Requires type-specialized workspaces for each `Dual` configuration
- The first evaluation with a given `Dual` type creates its workspace; later calls reuse it

**Enzyme** computes derivatives at the LLVM IR level:
- All computation stays in `Float64` - no special numeric types
- More restrictive about which code patterns it can differentiate; for example, the
  target function must not allocate arrays, so `EnzymeCache` pre-allocates its empty
  vectors

#### The Closure Type Problem

Julia closures get unique compiler-generated types. Enzyme's batch mode
(`BatchDuplicated`) requires the function passed to `prepare_jacobian` and to
`jacobian!` to have the same type:

```julia
f1 = (g, y) -> adcache(g, y, status, 0.0)  # Type: var"#1#2"{...}
f2 = (g, y) -> adcache(g, y, status, 0.0)  # Type: var"#3#4"{...}  <- different type

prep = DifferentiationInterface.prepare_jacobian(f1, g, backend, y)
DifferentiationInterface.jacobian!(f2, g, J, prep, backend, y)  # fails with Enzyme
```

#### The ResidualFunctor Pattern

Power flow wraps the residual in a callable struct, whose type is stable:

```julia
# PowerfulPF/src/solvers/powerflow.jl
struct ResidualFunctor{C, S}
    adcache::C
    status::S
end

function (f::ResidualFunctor)(g::AbstractVector, y::AbstractVector)
    f.adcache(g, y, f.status, 0.0)   # a power flow's state: the operating point, at t = 0
    return nothing
end
```

The functor is created once in setup, stored in the workspace and used for both
preparation and every `jacobian!` call:

| Module | Storage Location | Notes |
|--------|------------------|-------|
| Power flow | `NewtonWorkspace.residual_functor` | Used in Newton iterations |
| ACCC | `NewtonWorkspace.residual_functor` | One workspace per worker for parallel contingencies |
| TDS reinitialization | `ReinitCache.closure` | Closure created once, stored in the cache |

The power flow setup prepares the Jacobian with `strict=Val(false)` because parallel
ACCC evaluates the same functor type with different status dicts; the sparsity
pattern comes from the model structure and is a superset of the non-zeros under any
status.

TDS uses `JacobianClosureFunctor` (`PowerfulTDS/src/solvers/functors.jl`), a callable
struct created once per problem. Its adapter is held behind a `Ref{Any}` so that the
functor type does not embed the full system type:

```julia
struct JacobianClosureFunctor
    adapter_ref::Ref{Any}
    t_ref::Ref{Float64}
end
```

**Rules for Jacobian code:**

1. Create the functor once during setup and prepare the Jacobian with it.
2. Reuse the same functor and preparation in every iteration.
3. For parallel code, create one functor and preparation per thread or worker.

#### Implementation Files

| File | Purpose |
|------|---------|
| `PowerfulCore/src/dae/ad_backend.jl` | Backend types, `ADCacheCore`, `create_di_backend` |
| `PowerfulCore/src/dae/forwarddiff_cache.jl` | `ForwardDiffCache` with Dual workspace caching |
| `PowerfulCore/src/dae/enzyme_cache.jl` | `EnzymeCache` with the Float64-only path |
| `PowerfulCore/src/dae/autodiff_cache.jl` | `prepare_sparse_jacobian` |
| `PowerfulPF/src/solvers/powerflow.jl` | `ResidualFunctor`, Newton workspace |
| `PowerfulTDS/src/solvers/functors.jl` | `JacobianClosureFunctor` for TDS |

#### Testing Backend Equivalence

Both backends must produce the same residuals and Jacobians. See
`Powerful/test/performance/test_backend_equivalence.jl`:

```julia
using PowerfulCore: AddressManager, setup_with_formulation!, create_adcache
using PowerfulCore: ForwardDiffBackend, EnzymeBackend

cache_fd = create_adcache(am, sys, PowerInjection(), ForwardDiffBackend())
cache_ez = create_adcache(am, sys, PowerInjection(), EnzymeBackend())

cache_fd(g_fd, y, status, 0.0)
cache_ez(g_ez, y, status, 0.0)
@test g_fd ≈ g_ez rtol=1e-14
```

### Formulation as a First-Class Concept

- Different formulations = different types (`PowerInjection`, `CurrentInjection`, `DCPowerFlow`)
- Not config options or special cases
- The `AddressManager` is formulation-specific
- **Benefit:** Add new formulations without touching existing code

### Composition Over Inheritance

- Multi-island systems = `Vector{SystemModel}`
- Other array types = the system's `array_type`, not subclasses
- **Benefit:** Natural extensibility

## Data Flow

### 1. Parsing

```julia
using Powerful, PowerfulCases

case = parse_network(load("ieee14").raw)
# Returns a PowerFlowData.Network (buses, branches, loads, generators)
```

### 2. System Construction

```julia
sys = SystemModel(case)
```

**Steps:**
- `ext2int` - Map external keys (bus numbers, device IDs) to internal indices
- Retrieve external parameters needed before conversion (base voltages, bus types)
- Convert to per-unit (MW, kV → pu) and store the coefficients
- Retrieve external parameters that must be in per-unit

### 3. Problem Creation

```julia
prob = ACPowerFlowProblem(sys, PowerInjection())
```

**Steps:**
- Detect islands (each island is solved separately)
- Create the `AddressManager`
- Assign addresses based on formulation
  - PowerInjection: θ, V per bus
  - CurrentInjection: V_r, V_i per bus
- Apply the addressing strategy (e.g. AMD reordering)

### 4. Workspace Setup

```julia
using PowerfulPF: setup_powerflow_newton

ws = setup_powerflow_newton(prob)
```

**Steps:**
- Build the autodiff cache with a `ResidualBuffers` for each model
- Compute the Jacobian sparsity pattern and coloring
- Allocate the DAE vectors `y` (variables) and `g` (residuals)
- Fill the initial guess through the per-model `set_algebraic_guess!` hooks

### 5. Iteration

**Each Newton step:**

1. **Gather** - Copy `y` → local buffers using addresses
2. **Compute** - Each model's kernel computes residuals
3. **Scatter** - Copy local buffers → `g` using addresses
4. **Jacobian** - AD computes ∂g/∂y
5. **Solve** - J·Δy = -g (sparse linear solve)
6. **Update** - y += Δy
7. **Check** - Converged when ‖g‖ is small enough

### 6. Result Extraction

```julia
v = get_var(sol, :Bus, :v, (4001,))  # Voltage magnitude (pu) at bus 4001
bus_table = collect_model_table(sol, :Bus)  # Tables.jl-compatible NamedTuple
```

`get_var` accepts variable names only; asking for a residual name raises an error
that lists the model's variables. `collect_model_table` converts values to the units
selected by its `ReportConfig`.

## Extension Points

### New Device Models

**Implement:**
1. Struct subtyping `AbstractModel`
2. `model_metadata` (variables, residuals, parameters)
3. Residual kernels and a `calc_residuals!` method for each supported formulation
4. Register with `@register_model`

See [Implementing Models](models.md).

### New Formulations

**Implement:**
1. A new `FormulationTrait` subtype
2. Formulation-specific variables and residuals in each model's metadata
3. `calc_residuals!` methods for each model

The `AddressManager` and workspaces handle the rest. See [Multi-Physics Formulations](formulations.md).

### New Addressing Strategies

Addressing strategies are subtypes of `AddressingStrategy`
(`PowerfulCore/src/foundation/core_types.jl`), selected with the `addressing=` keyword of
`ACPowerFlowProblem` and `TDSProblem`. After natural addresses are assigned,
`apply_addressing_strategy!` (`PowerfulPF/src/solvers/power_flow_api.jl`) computes the
permutation for the strategy and remaps the addresses with `remap_addresses!`. Any
precomputed ordering can be supplied as `CustomPermutation(perm)`; a new strategy that
computes its own ordering needs a new type and a branch in `apply_addressing_strategy!`.

### New Solver Algorithms

**Implement:**
1. Create workspaces
2. Run your iteration method

**Reuse `ResidualBuffers`, `SystemBuffers` and the AD infrastructure.**

### New Result Access Patterns

Add convenience functions in `PowerfulPF/src/access/variable_accessors.jl`, built on
`get_var` and the solution's problem and system references.

## Reading Paths by Role

### Model Developers

1. This overview
2. [Implementing Models](models.md) - Implementation checklist
3. [Per-Unit & Multi-Island](miscellaneous.md) - Per-unit details
4. Reference: `PowerfulPF/src/models/pq.jl`

### Solver Developers

1. This overview
2. [Computation Architecture](workspaces.md) - Computational structure
3. [Multi-Physics Formulations](formulations.md)
4. Reference: `PowerfulPF/src/solvers/powerflow.jl`

### Architecture Explorers

1. This overview
2. [Three-Layer Indexing](indexing.md)
3. [Multi-Physics Formulations](formulations.md)

### Performance Optimizers

1. [Computation Architecture](workspaces.md) - Pre-allocated workspaces
2. `PowerfulCore/src/dae/residual_workspace.jl` - gather/scatter
3. [Implementing Models](models.md) - Residual kernels specialized per array type

## Logging and Verbose Output

A two-tier logging pattern: `@vprint` for user-facing output, `@debug` for developer diagnostics.

### `@vprint` - User-Facing Verbose Output

For solver progress, iteration details and timing information:

```julia
using PowerfulCore: @vprint

function solve(prob; verbose::Bool = false)
    @vprint "Setting up problem..."           # Only prints if verbose=true
    @vprint "Size: $(length(data)) elements"  # Interpolated only when printing

    verbose && print("Step 1/5: Computing... ")  # Mid-line (no newline)
    # ... work ...
    @vprint "done ($(round(elapsed, digits=3))s)"  # Completes the line
end
```

**Requirements:**
- A `verbose::Bool` variable must exist in the calling scope (the macro reads `verbose` from the caller)
- Import via `using PowerfulCore: @vprint`
- When `verbose=false` the string is never built

### `@debug` - Developer Diagnostics

For internal diagnostics:

```julia
@debug "PowerFlowInit: converged in $iter iterations"
@debug "Event at t=$time: applying fault"
```

**Enable with:** `ENV["JULIA_DEBUG"] = "all"` or `ENV["JULIA_DEBUG"] = "PowerfulTDS"`

### Comparison

| Aspect | `@vprint` | `@debug` |
|--------|-----------|----------|
| Purpose | User-facing progress | Developer diagnostics |
| Default | Off (`verbose=false`) | Off (logging level) |
| Enable | Pass `verbose=true` | `ENV["JULIA_DEBUG"]` |
| Location shown | No | Yes (file:line) |

### CLI `--verbose` Pass-Through

The `bin/powerful` CLI extracts `--verbose` wherever it appears and passes it on to
the Julia script:

```bash
./bin/powerful --verbose pf case.raw
./bin/powerful pf case.raw --verbose
```

CLI scripts must forward the flag, as `bin/powerflow.jl` does:

```julia
using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "--verbose", "-v"
        action = :store_true
end
args = parse_args(s)

sol = solve(prob, solver; verbose=args["verbose"])
```

Scripts run with `./bin/powerful run` receive `--verbose` in `ARGS`:
```julia
verbose_mode = "--verbose" in ARGS
```

### When to Use Each

| Pattern | Use For |
|---------|---------|
| `@vprint` | Solver iterations, setup steps, timing, convergence status |
| `@debug` | Internal state dumps, event triggers, diagnostic info |
| `@info`/`@warn` | Sparingly - important notices only, avoid in hot paths |

### Implementation Reference

- `PowerfulPF/src/solvers/powerflow.jl` - Newton solver verbose output
- `PowerfulPF/src/solvers/fastdecoupled.jl` - Fast decoupled verbose output
- `PowerfulTDS/src/solvers/solve.jl` - TDS solver verbose output
