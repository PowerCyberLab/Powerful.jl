# System Architecture Overview

This page describes the architecture of Powerful.jl for contributors. It introduces the package structure, the layers of the design, the data flow of a power flow solution, and the points at which the package is extended.

## Workflow

A power flow study proceeds in four steps:

1. **Parse**: a PSS/E RAW file is read into a `SystemModel` that holds data only.
2. **Setup**: an `ACPowerFlowProblem` selects the formulation and the variable layout.
3. **Solve**: Newton iterations run until convergence.
4. **Extract**: results are addressed by the device IDs of the data file.

The data structures hold no solver logic, and the solvers hold no topology logic. This separation allows several problems with different formulations to be created from one system.

## Package Structure

Powerful.jl is a monorepo of sub-packages:

| Package | Purpose | Key Contents |
|---------|---------|--------------|
| `PowerfulCore/` | Foundation infrastructure | Types, traits, metadata, DAE, `AddressManager`, `SystemModel` |
| `PowerfulPF/` | Power flow | Bus, PQ, PV, Slack, Branch, Transformer models; Newton and fast-decoupled solvers |
| `PowerfulTDS/` | Time-domain simulation | GENCLS, GENROU, GENSAL, GENTPJ1, SEXS, TGOV1, REGCA1, REECA1 models; TDS solver, events |
| `PowerfulYbus/` | Y-bus time-domain simulation | Sparse network solve with trapezoidal machine integration |
| `PowerfulCA/` | Contingency analysis | ACCC solvers, DC screening (PTDF/LODF), violation reporting |
| `PowerfulCases/` | Test cases | IEEE 14, 39, 118 bus and other cases |
| `Powerful/` | Umbrella package | Curated public API |

The dependencies between the packages are:

```
PowerfulCore → PowerfulPF → PowerfulTDS → Powerful (also depends on PowerfulCases)
                    │            └→ PowerfulYbus  (standalone)
                    └→ PowerfulCA   (standalone)
```

`Powerful` depends on `PowerfulCore`, `PowerfulPF`, `PowerfulTDS` and `PowerfulCases`. `PowerfulCA` and `PowerfulYbus` are standalone packages that `Powerful` does not re-export; they are loaded explicitly (for example `using PowerfulCA`). Further standalone packages extend the same layers without changing them.

User code imports the umbrella package and reaches model types through the `Models` submodule. Developer code imports each symbol from the package that defines it:

```julia
# User code
using Powerful
using PowerfulCases
using PowerfulCA  # Contingency analysis (standalone)

using Powerful.Models: Bus, PQ, PV, TGOV1

# Developer code
using PowerfulCore: model_metadata, AddressManager
using PowerfulPF: setup_powerflow_newton
```

## Architectural Layers

The layers build bottom-up, and each layer is independent of the layers above it.

### Foundation Layer

The foundation layer (`PowerfulCore/src/foundation/`) defines the fundamental abstractions, which change rarely:

- `AbstractModel`, the base type of all device models.
- `ModelMetadata`, which lists the variables, residuals and parameters of each model.
- `FormulationTrait`, with the subtypes `PowerInjection`, `CurrentInjection` and `DCPowerFlow`.
- The address types `AlgebVar`, `StateVar`, `AlgebRes`, `StateRes` and `ObservableVar`.

### Indexing Layer

The indexing layer (`PowerfulCore/src/dae/`) maps between three coordinate systems:

1. External IDs, such as bus 4001 or 4002, as they appear in the data file.
2. Internal indices (1, 2, 3, ...), the contiguous positions in the model arrays.
3. DAE addresses, such as `y[42]` or `g[58]`, the positions in the equation vectors.

Two mappers connect them. `IndexOrder`, built by `ext2int`, maps external IDs to internal indices. `AddressManager` maps internal indices to DAE addresses. Different formulations have different variables, so the `AddressManager` is created per problem and is not stored in the system data. See [Three-Layer Indexing](indexing.md).

### Model Layer

The model layer (`PowerfulPF/src/models/`, `PowerfulTDS/src/models/`) contains one entry per device type (PQ load, PV generator, Branch). Each entry consists of a struct definition that holds pure data, the metadata that declares variables, residuals and parameters, and the residual kernels that compute power injections or currents.

Struct fields hold the external IDs of the devices a model connects to, and the metadata declares the resulting dependencies. The system resolves these dependencies during setup.

### System Layer

The system layer (`PowerfulCore/src/system/`) defines `SystemModel`, which contains:

- The model collections (StructArrays, in per-unit).
- The external parameters retrieved from other models.
- The mappings from external to internal indices (`IndexOrder`).
- The per-unit conversion coefficients and the original parsed data.
- The system-wide settings (`SystemConfig`: frequency, MVA base, constant-power breakpoint).

A `SystemModel` holds no computational state. The DAE vectors and the Jacobians live in workspaces.

### Problem Layer

The problem layer (`PowerfulPF/src/solvers/`) defines `ACPowerFlowProblem`, which specifies:

- The `SystemModel` to solve.
- The formulation (`PowerInjection`, `CurrentInjection` or `DCPowerFlow`).
- The variable ordering, selected with the `addressing=` keyword (`NaturalAddressing()`, `AMDAddressing()` or `CustomPermutation(perm)`).
- The autodiff backend, selected with the `backend=` keyword.

Creating the problem also creates the `AddressManager`, which assigns every variable and residual to a position in the DAE vectors.

### Solver Layer

The solver layer (`PowerfulPF/src/solvers/`, `PowerfulTDS/src/solvers/`) sets up the workspaces and runs the iterations. The setup creates the autodiff cache (`create_adcache`), which builds the `SystemBuffers`, and prepares the sparse Jacobian from the model structure and its coloring. Each iteration evaluates the residuals and the Jacobian, solves the linear system, updates the state, and checks convergence. The physics is evaluated by the workspaces, and the solver implements the numerical methods only.

## Design Principles

Four principles shape the code base.

**Data and computation are separate.** A `SystemModel` stores data, and workspaces read from it to perform computation. Several problems with different formulations can therefore share one system.

**Models are described by metadata.** A model declares its interface through `ModelMetadata`, and the system introspects the metadata to allocate variables. A new model is added without modifying the core.

**Formulations are types.** `PowerInjection`, `CurrentInjection` and `DCPowerFlow` are distinct Julia types, so they participate in method dispatch, and the `AddressManager` is formulation-specific. A new formulation is added without changing existing code.

**Composition provides extensibility.** A multi-island system is a `Vector{SystemModel}`, and other array types (such as GPU arrays) are selected through the `array_type` of the system.

Performance is supported by three implementation conventions. Concrete types and type parameters are used throughout; model structs are parametric (`T<:Number` fields) so that they work with `ForwardDiff.Dual`. The hot-path data of each model are held in NamedTuples in `ResidualBuffers` (`var_vals`, `res_vals`), while `SystemBuffers` holds one closure per model in a `Vector{Function}` to bound the compilation cost of many model types. Finally, workspaces are allocated once during setup, residual evaluation fills the existing buffers, and Jacobian evaluation reuses the prepared sparsity pattern and coloring.

## Automatic Differentiation Backends

Powerful.jl uses [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl) to compute sparse Jacobians with pluggable AD backends. The backend types `ForwardDiffBackend()` (the default) and `EnzymeBackend()` are defined in `PowerfulCore/src/dae/ad_backend.jl` and passed with the `backend` keyword:

```julia
using PowerfulCore: EnzymeBackend, ForwardDiffBackend

prob = ACPowerFlowProblem(sys; backend = EnzymeBackend())
tds_prob = TDSProblem(init; tspan = (0.0, 10.0), backend = ForwardDiffBackend())
```

| Aspect | `ForwardDiffBackend()` | `EnzymeBackend()` |
|--------|------------------------|-------------------|
| Mechanism | Dual number propagation | LLVM-level differentiation with shadow arrays |
| Workspace numeric type | `Dual` workspaces cached per Dual type | `Float64` throughout |
| Cache type | `ForwardDiffCache` | `EnzymeCache` |

ForwardDiff replaces `Float64` with `Dual{Tag,Float64,N}`, where each `Dual` carries a value and N partial derivatives. It needs a type-specialized workspace for each `Dual` configuration; the first evaluation with a given `Dual` type creates the workspace, and later calls reuse it.

Enzyme differentiates at the LLVM IR level, so all computation stays in `Float64`. It restricts the code patterns that can be differentiated. For example, the target function must not allocate arrays, so `EnzymeCache` pre-allocates its empty vectors.

### Closures and Callable Structs

Julia closures receive unique compiler-generated types. Enzyme's batch mode (`BatchDuplicated`) requires the function passed to `prepare_jacobian` and the function passed to `jacobian!` to have the same type. Two separately created closures differ in type, and the second call fails:

```julia
f1 = (g, y) -> adcache(g, y, status, 0.0)  # Type: var"#1#2"{...}
f2 = (g, y) -> adcache(g, y, status, 0.0)  # Type: var"#3#4"{...}  <- different type

prep = DifferentiationInterface.prepare_jacobian(f1, g, backend, y)
DifferentiationInterface.jacobian!(f2, g, J, prep, backend, y)  # fails with Enzyme
```

Power flow therefore wraps the residual in a callable struct with a stable type:

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

The functor is created once in setup, stored in the workspace, and used for the preparation and for every `jacobian!` call:

| Module | Storage Location | Notes |
|--------|------------------|-------|
| Power flow | `NewtonWorkspace.residual_functor` | Used in Newton iterations |
| ACCC | `NewtonWorkspace.residual_functor` | One workspace per worker for parallel contingencies |
| TDS reinitialization | `ReinitCache.closure` | Closure created once, stored in the cache |

The power flow setup prepares the Jacobian with `strict=Val(false)`, because parallel ACCC evaluates the same functor type with different status dicts. The sparsity pattern comes from the model structure and is a superset of the non-zeros under any status.

TDS uses `JacobianClosureFunctor` (`PowerfulTDS/src/solvers/functors.jl`), a callable struct created once per problem. Its adapter is held behind a `Ref{Any}`, which keeps the full system type out of the functor type:

```julia
struct JacobianClosureFunctor
    adapter_ref::Ref{Any}
    t_ref::Ref{Float64}
end
```

Jacobian code follows three rules:

1. The functor is created once during setup, and the Jacobian is prepared with it.
2. The same functor and preparation are reused in every iteration.
3. Parallel code creates one functor and one preparation per thread or worker.

### Implementation Files

| File | Purpose |
|------|---------|
| `PowerfulCore/src/dae/ad_backend.jl` | Backend types, `ADCacheCore`, `create_di_backend` |
| `PowerfulCore/src/dae/forwarddiff_cache.jl` | `ForwardDiffCache` with Dual workspace caching |
| `PowerfulCore/src/dae/enzyme_cache.jl` | `EnzymeCache` with the Float64-only path |
| `PowerfulCore/src/dae/autodiff_cache.jl` | `prepare_sparse_jacobian` |
| `PowerfulPF/src/solvers/powerflow.jl` | `ResidualFunctor`, Newton workspace |
| `PowerfulTDS/src/solvers/functors.jl` | `JacobianClosureFunctor` for TDS |

### Testing Backend Equivalence

Both backends must produce the same residuals and Jacobians. `Powerful/test/performance/test_backend_equivalence.jl` checks this:

```julia
using PowerfulCore: AddressManager, setup_with_formulation!, create_adcache
using PowerfulCore: ForwardDiffBackend, EnzymeBackend

cache_fd = create_adcache(am, sys, PowerInjection(), ForwardDiffBackend())
cache_ez = create_adcache(am, sys, PowerInjection(), EnzymeBackend())

cache_fd(g_fd, y, status, 0.0)
cache_ez(g_ez, y, status, 0.0)
@test g_fd ≈ g_ez rtol=1e-14
```

## Data Flow

### Parsing

```julia
using Powerful, PowerfulCases

case = parse_network(load("ieee14").raw)
# Returns a PowerFlowData.Network (buses, branches, loads, generators)
```

### System Construction

```julia
sys = SystemModel(case)
```

The constructor performs the following steps:

- `ext2int` maps external keys (bus numbers, device IDs) to internal indices.
- External parameters needed before conversion (base voltages, bus types) are retrieved.
- The data are converted to per-unit (MW, kV → pu) and the coefficients are stored.
- External parameters that must be in per-unit are retrieved.

### Problem Creation

```julia
prob = ACPowerFlowProblem(sys, PowerInjection())
```

Problem creation performs the following steps:

- Islands are detected, and each island is solved separately.
- The `AddressManager` is created.
- Addresses are assigned according to the formulation: θ and V per bus for `PowerInjection`, V_r and V_i per bus for `CurrentInjection`.
- The addressing strategy (for example AMD reordering) is applied.

### Workspace Setup

```julia
using PowerfulPF: setup_powerflow_newton

ws = setup_powerflow_newton(prob)
```

The setup performs the following steps:

- The autodiff cache is built with a `ResidualBuffers` for each model.
- The Jacobian sparsity pattern and coloring are computed.
- The DAE vectors `y` (variables) and `g` (residuals) are allocated.
- The initial guess is filled through the per-model `set_algebraic_guess!` hooks.

### Iteration

Each Newton step consists of seven operations:

1. Gather: `y` is copied into local buffers using the addresses.
2. Compute: each model's kernel computes its residuals.
3. Scatter: the local buffers are copied into `g` using the addresses.
4. Jacobian: AD computes ∂g/∂y.
5. Solve: the sparse linear system J·Δy = -g is solved.
6. Update: `y += Δy`.
7. Check: the iteration has converged when ‖g‖ is sufficiently small.

### Result Extraction

```julia
v = get_var(sol, :Bus, :v, (4001,))  # Voltage magnitude (pu) at bus 4001
bus_table = collect_model_table(sol, :Bus)  # Tables.jl-compatible NamedTuple
```

`get_var` accepts variable names; a residual name raises an error that lists the variables of the model. `collect_model_table` converts values to the units selected by its `ReportConfig`.

## Extension Points

**New device models.** A model is a struct that subtypes `AbstractModel`, a `model_metadata` method (variables, residuals, parameters), residual kernels with a `calc_residuals!` method for each supported formulation, and a registration with `@register_model`. See [Implementing Models](models.md).

**New formulations.** A formulation is a new `FormulationTrait` subtype together with formulation-specific variables and residuals in the metadata of each model and a `calc_residuals!` method for each model. The `AddressManager` and the workspaces handle the remaining work. See [Multi-Physics Formulations](formulations.md).

**New addressing strategies.** Addressing strategies are subtypes of `AddressingStrategy` (`PowerfulCore/src/foundation/core_types.jl`), selected with the `addressing=` keyword of `ACPowerFlowProblem` and `TDSProblem`. After natural addresses are assigned, `apply_addressing_strategy!` (`PowerfulPF/src/solvers/power_flow_api.jl`) computes the permutation of the strategy and remaps the addresses with `remap_addresses!`. Any precomputed ordering can be supplied as `CustomPermutation(perm)`. A strategy that computes its own ordering needs a new type and a branch in `apply_addressing_strategy!`.

**New solver algorithms.** A solver creates its workspaces and runs its iteration method, reusing `ResidualBuffers`, `SystemBuffers` and the AD infrastructure.

**New result accessors.** Convenience functions are added in `PowerfulPF/src/access/variable_accessors.jl`, built on `get_var` and on the problem and system references of the solution.

## Further Reading

- Model developers: [Implementing Models](models.md), [Per-Unit & Multi-Island](miscellaneous.md), and the reference implementation `PowerfulPF/src/models/pq.jl`.
- Solver developers: [Computation Architecture](workspaces.md), [Multi-Physics Formulations](formulations.md), and the reference implementation `PowerfulPF/src/solvers/powerflow.jl`.
- Architecture: [Three-Layer Indexing](indexing.md) and [Multi-Physics Formulations](formulations.md).
- Performance: [Computation Architecture](workspaces.md), the gather and scatter code in `PowerfulCore/src/dae/residual_workspace.jl`, and the array-type specialization of residual kernels in [Implementing Models](models.md).

## Logging and Verbose Output

Powerful.jl provides two logging tiers. `@vprint` (from `PowerfulCore`) prints user-facing progress such as solver iterations, setup steps and timing. It requires a `verbose::Bool` variable in the calling scope, prints only when `verbose` is `true`, and builds its string only in that case. `@debug` carries developer diagnostics and is enabled with `ENV["JULIA_DEBUG"] = "all"` or `ENV["JULIA_DEBUG"] = "PowerfulTDS"`. `@info` and `@warn` are reserved for important notices and are avoided in hot paths.

```julia
using PowerfulCore: @vprint

function solve(prob; verbose::Bool = false)
    @vprint "Setting up problem..."
    @vprint "Size: $(length(data)) elements"

    @debug "Event at t=$time: applying fault"
end
```

The `bin/powerful` CLI extracts `--verbose` wherever it appears and passes it on to the Julia script. A CLI script forwards the flag to the solver, as `bin/powerflow.jl` does with `solve(prob, solver; verbose=args["verbose"])`. Scripts started with `./bin/powerful run` find the flag in `ARGS` (`"--verbose" in ARGS`).
