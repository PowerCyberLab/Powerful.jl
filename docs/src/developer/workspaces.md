# Computation Architecture

How residual evaluation is organized so that Newton and time-domain iterations work on
pre-allocated buffers.

## The Workspace Hierarchy

**Three nested layers with clear responsibilities:**

```
SystemBuffers                    # Coordinates all models
├── ResidualBuffers (PQ)         # Data movement for PQ loads
│   └── Residual kernels         # Pure math
├── ResidualBuffers (PV)         # Data movement for PV generators
│   └── Residual kernels         # Pure math
└── ResidualBuffers (Branch)     # Data movement for branches
    └── Residual kernels         # Pure math
```

### SystemBuffers (Coordination Layer)

Defined in `PowerfulCore/src/system/system_buffers.jl`:

```julia
struct SystemBuffers{T<:Number, F<:FormulationTrait}
    formulation::F                              # Formulation trait
    evaluators::Vector{Function}                # Hot path: (f, g, x, y, t, status) -> nothing
    evaluators_with_recorder::Vector{Function}  # Same, with an InjectionRecorder argument
    evaluators_accumulate::Vector{Function}     # Cold path: evaluate_generator_contributions!
    workspaces_by_name::Dict{Symbol, Any}       # Model name → ResidualBuffers
    metadata_by_name::Dict{Symbol, Any}         # Model name → formulation-specific ModelMetadata
    own_models::Vector{Symbol}                  # Models with own algebraic or zero-mass rows
    own_rows::Vector{Union{Nothing, OwnRows}}   # Their held rows (built on first use)
    zero_mass::Vector{Vector{Pair{Symbol, BitVector}}}  # Own states with zero-mass instances
    out_of_service::Vector{Bool}                # Per own model: out-of-service instances?
    hold::Vector{Float64}                       # Held values of out-of-service own variables
    hold_x::Vector{Float64}                     # Held values of their own zero-mass states
end
```

**Design:** `SystemBuffers` stores one closure per model in a `Vector{Function}` instead of
a NamedTuple of typed workspaces. With a typed NamedTuple, every model type would be a
type parameter of the system workspace, and compilation time grows quickly with the
number of models. Each closure captures its model's concrete types internally, so the
code inside it is fully specialized; calling it through the `Vector{Function}` costs one
dynamic dispatch per model per evaluation.

**Job:** Loop through models in evaluation order and call each model's evaluation

**Knows:** System structure, evaluation order (`evaluation_pass` from metadata)

**Doesn't know:** Model-specific calculations

### ResidualBuffers (Data Movement Layer)

Defined in `PowerfulCore/src/dae/residual_workspace.jl`. For each model it holds:

- `var_addrs`, `var_vals` - NamedTuples of DAE addresses and gathered values per variable (e.g. `theta`, `v`)
- `res_addrs`, `res_vals` - NamedTuples of DAE addresses and computed values per residual (e.g. `p`, `q`)
- `var_gather_tuple`, `res_scatter_tuple` - Precomputed tuples used by gather/scatter
- `n_instances` - Instance count
- `setpoints` - Initialization setpoints of dynamic models (empty NamedTuple otherwise)
- `system` - The `SystemConfig` of the system (`sys.settings`: frequency, MVA base, constant-power breakpoint)

The in-service mask is **not** stored in the workspace. It is passed to
`calc_residuals!` and `scatter_residuals!` as `active_mask`, taken from a status dict
(`Dict{Symbol, Vector{Bool}}`, built by `build_status_dict(sys)`). Parallel contingency
analysis relies on this: each thread can evaluate the same workspace structure with its
own status.

**Out-of-service instances.** An out-of-service instance contributes nothing, but its
own algebraic variables (variables without a source model, such as a bus's `theta`
and `v` or a machine's `P_m`) are still entries of `y`, and nothing else writes their
residual rows. Each model's own algebraic residuals share addresses with its own
algebraic variables (checked the first time an instance of the model is out of
service), so after all models are
evaluated `hold_out_of_service!` sets each such row `a` to `y[a] - hold[a]`: an identity
Jacobian row that keeps the variable at its held value. The solver sets `hold` with
`set_hold!` from its starting point (the initial guess of a power flow, the initialized
values of a TDS). This is what keeps a bus de-energized by a contingency, or a machine
tripped during a simulation, from making the Jacobian singular. Own states whose mass
is zero for a particular instance (a state whose mass coefficient reduces to zero there) are algebraic too:
their rows `f[a]` are held the same way, at `hold_x`; states with mass keep
`M·ẋ = 0`, their values. Models therefore just mask their equations with the in-service
mask. Each model's evaluator
records whether it has out-of-service instances, so the pass costs nothing when every
instance is in service.

**Job:** Gather variables, compute residuals, scatter results

**Doesn't know:** Mathematical details

### Residual Kernels (Computation Layer)

**Pure functions:**
```julia
kernel!(outputs::NamedTuple, params::NamedTuple, vars::NamedTuple, active, AT)
```

**Job:** Loop over arrays, compute residuals from parameters and variables

**Knows:** Physics equations only

**Doesn't know:** Workspaces, addresses, system structure

## The Gather-Compute-Scatter Pattern

Each model's evaluator closure runs:

```julia
gather_variables!(ws, x, y)              # y (and x for states) → ws.var_vals
gather_residuals!(ws, f, g)              # accumulated values of SharedRes(:set) residuals
calc_residuals!(ws, model, external, formulation, t; active_mask)
scatter_residuals!(ws, f, g, active_mask)
```

`x`/`f` are the differential states and residuals (empty in power flow), `y`/`g` the
algebraic variables and residuals.

### Phase 1: Gather

**Copy from the global DAE vectors to the workspace:**
```julia
# For each variable of the model, e.g. the bus voltage seen by a PQ load:
ws.var_vals.v[i] = y[ws.var_addrs.v[i]]
```

State variables are gathered from `x`, algebraic variables from `y`.

**Result:** Contiguous arrays ready for SIMD loops or GPU kernels

### Phase 2: Compute

**The kernel operates on workspace arrays** (the PQ load kernel from `PowerfulPF/src/models/pq.jl`):
```julia
@inbounds @simd for i in eachindex(outputs.p)
    mask = active[i]
    power = calc_zip_load(vars.v[i], params.pl[i], params.ql[i],
                          params.ip[i], params.iq[i], params.yp[i], params.yq[i], params.vmin)
    # Load consumption = negative injection into the bus
    outputs.p[i] = mask * (-power.p)
    outputs.q[i] = mask * (-power.q)
end
```

**Result:** Residual values in workspace arrays

### Phase 3: Scatter

**Copy from the workspace to the global residual vectors:**
```julia
# SharedRes(:add) and PrivateRes - accumulate / write
g[ws.res_addrs.p[i]] += ws.res_vals.p[i]

# SharedRes(:set) - overwrite (second pass, after all :add contributions of this model)
g[ws.res_addrs.v_ctrl[i]] = ws.res_vals.v_ctrl[i]
```

Only active instances are scattered.

**Result:** Global residual vector updated

## Why NamedTuples?

**Problem with `Dict{Symbol, Any}`:**
- Runtime type checks on every access
- Dictionary lookups on the hot path
- Boxed values

**Benefits of NamedTuples:**

### 1. Type Stability
```julia
ws.var_vals.v  # The compiler knows the array type
```

### 2. Direct Field Access
Field access on a NamedTuple compiles to a direct load, like struct field access.

### 3. Compile-Time Specialization
```julia
# Each model gets a specialized code path
calc_residuals!(ws::ResidualBuffers{Float64, Vector{Float64}, ...}, pq::StructArray{PQ{Float64}}, ...)
```

**Trade-off:** Fields can't be added at runtime. Workspaces are created once at setup
with field names known from the formulation-specific metadata.

## CPU and GPU Dispatch

**Same kernel interface, different implementations.** The last argument `AT` is the
workspace array type.

### CPU Kernel

In the model file (e.g. `PowerfulPF/src/models/pq.jl`):

```julia
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    # @inbounds @simd loop shown above
    nothing
end
```

**Key:** `::Type{<:Vector}` matches CPU arrays

### GPU Kernel

A device function and a launcher method of the same kernel, for the GPU extension's
array type:

```julia
# Device function: one thread per load
function pq_residual_cuda_kernel!(
    outputs_p, outputs_q,
    params_pl, params_ql, params_ip, params_iq, params_yp, params_yq, vmin,
    vars_v, active, n
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= n
        mask = active[i] ? one(eltype(outputs_p)) : zero(eltype(outputs_p))
        power = calc_zip_load(vars_v[i], params_pl[i], params_ql[i],
                              params_ip[i], params_iq[i], params_yp[i], params_yq[i], vmin)
        @inbounds outputs_p[i] = mask * (-power.p)
        @inbounds outputs_q[i] = mask * (-power.q)
    end
    return nothing
end

# Launcher
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:CUDA.CuArray}  # GPU dispatch
)
    n = length(outputs.p)
    (; threads, blocks) = cuda_kernel_config(n)
    CUDA.@cuda threads=threads blocks=blocks pq_residual_cuda_kernel!(
        outputs.p, outputs.q,
        params.pl, params.ql, params.ip, params.iq, params.yp, params.yq, params.vmin,
        vars.v, active, n
    )
    nothing  # kernels run in order on the default stream
end
```

**Key:** `::Type{<:CuArray}` matches GPU arrays. The physics function `calc_zip_load` is
shared by both versions.

### Automatic Dispatch

The array type comes from the system (`sys.array_type`). Moving the system or problem
to a different array type selects the matching kernels; no workspace option is needed,
and no manual branching: the type system routes to the right hardware.

## How Workspaces Are Created

**During solver setup** (`setup_powerflow_newton`, TDS `solve`), the autodiff cache
(`create_adcache`) builds the `SystemBuffers`.

### 1. SystemBuffers Construction
```julia
using PowerfulCore: SystemBuffers

sys_ws = SystemBuffers(Float64, am, sys, formulation)
```

**Steps:**
- Order the models by `evaluation_pass`, skipping models without instances, variables or residuals
- Get the formulation-specific metadata of each model
- Create a `ResidualBuffers` for each model
- Build the evaluator closures (hot path, recorder path, accumulate path)

### 2. ResidualBuffers Construction
```julia
using PowerfulCore: ResidualBuffers

ws = ResidualBuffers(Float64, metadata, am, model_vec, sys.array_type;
                     system = sys.settings, setpoints = NamedTuple(), config = SimConfig())
```

**Steps:**
- Allocate arrays for each variable and residual
- Extract address vectors from the AddressManager
- Allocate arrays for observable variables and setpoints (or reference the given setpoints)
- Precompute the gather and scatter tuples

### 3. Type Parameter Propagation
```julia
ResidualBuffers{Float64, Vector{Float64}, ...}    # CPU
ResidualBuffers{Float64, CuVector{Float64}, ...}  # GPU
```

**The array type flows through the hierarchy via type parameters.**

### 4. Cold-Path Access
```julia
sys_ws.metadata_by_name[:PQ]    # ModelMetadata for PQ (formulation-specific)
sys_ws.workspaces_by_name[:PQ]  # ResidualBuffers for PQ
```

### 5. One-Time Allocation

**After setup:**
- Array storage is not reallocated
- Iterations fill existing buffers
- All indexing uses fixed address vectors

## How Model References Flow Through Workspaces

**Reference chain:**

```
ACPowerFlowProblem
  ├── sys::SystemModel         # Reference (not copy)
  └── address_manager

SystemBuffers
  ├── evaluators::Vector{Function}  # Closures capture model_vec references
  └── workspaces_by_name::Dict      # Cold path access

Each evaluator closure captures:
  ├── ws::ResidualBuffers           # Pre-allocated workspace
  ├── model_vec::StructArray        # Reference to e.g. sys.models.PQ (not a copy)
  └── external::NamedTuple          # sys.external for that model
```

**Design:** Each evaluator closure captures a reference to the model's StructArray.
Changing a model field in place (`sys.models.PQ.pl[3] = 0.5`) is visible at the next
residual evaluation.

### Modify Parameters Without Rebuilding the Workspace
```julia
using PowerfulPF: setup_powerflow_newton, solve_powerflow_newton!

ws = setup_powerflow_newton(prob)
y0 = copy(ws.iter_dae.y)
base_pl = copy(sys.models.PQ.pl)

for scale in [0.9, 1.0, 1.1]
    sys.models.PQ.pl .= base_pl .* scale   # in place
    sol = solve_powerflow_newton!(ws, copy(y0))
end
```

Only in-place changes to the arrays the closures hold are seen. Replacing a model
collection, adding instances or changing topology requires a new problem and workspace.
External parameters (`sys.external`) are retrieved at system construction and are not
recomputed when model fields change.

## Performance Characteristics

### Pre-allocated Evaluation

- Gather-compute-scatter uses pre-allocated buffers
- Array storage is allocated at setup, not per iteration
- Allocation regressions are guarded by tests in `Powerful/test/performance/` and
  `Powerful/test/perf/`

### Cache-Friendly Access Patterns

**Gather creates contiguous arrays from scattered DAE vectors:**
- Kernels then access memory sequentially
- SIMD processes several elements per instruction
- This is also why models are stored as StructArrays (struct of arrays)

### Type-Stable Kernels

**Inside each evaluator closure:**
- Every call is type-stable
- Code is specialized for the concrete model, workspace and formulation types
- The only dynamic dispatch is the call to the closure itself

### GPU Parallelism

**When using GPU arrays:**
- Gather, compute and scatter run as device kernels, one thread per element
- Scatter of `SharedRes(:add)` residuals uses atomic adds, since several devices
  contribute to the same bus equation

## Extension Points

### Custom Gather/Scatter

The default gather/scatter loops over address arrays. Specialized versions for other
array types can be added by dispatching on the array type in the same way.

### Alternative Dispatch Strategies

**Current:** Dispatch on the array type (`Vector` vs `CuArray`)

**Possible:** Dispatch on element type or memory layout:
```julia
kernel!(outputs, params, vars, active, ::Type{Vector{Float32}})
kernel!(outputs, params, vars, active, ::Type{Vector{Float64}})
```

### Workspace Reuse Across Solves

```julia
ws = setup_powerflow_newton(prob)

# Solve several times with different initial guesses
sol1 = solve_powerflow_newton!(ws, y0_flat)
sol2 = solve_powerflow_newton!(ws, y0_warm)
```

**Amortizes setup cost for repeated similar problems.**

### Parallel Model Evaluation

**Current:** Models are evaluated sequentially, in the order given by `evaluation_pass`.

**Possible:** Parallelize within a pass; `SharedRes(:add)` scatter then needs atomic
operations, as on the GPU.

## Architecture Summary

**Key principles:**

1. **Three-layer separation** - Coordination, data movement, computation
2. **NamedTuples inside each model workspace** - Type stability and specialization
3. **Type-erased closures across models** - Bounded compilation time
4. **Gather-compute-scatter pattern** - Enables SIMD and GPU parallelism
5. **Reference-based system access** - Change parameters without recreating workspaces
6. **Type parameters for dispatch** - Automatic CPU/GPU selection

## Quick Reference

### Workspace Creation
```julia
using PowerfulPF: setup_powerflow_newton

prob = ACPowerFlowProblem(sys, PowerInjection())
ws = setup_powerflow_newton(prob)

# GPU arrays; use_gpu selects the GPU linear solver
ws_gpu = setup_powerflow_newton(to_gpu(prob); use_gpu = true)
```

### Residual Evaluation
```julia
using PowerfulCore: SystemBuffers, build_status_dict, evaluate_system_residuals!

sys_ws = SystemBuffers(Float64, prob.address_manager, prob.sys, prob.formulation)
status = build_status_dict(prob.sys)
evaluate_system_residuals!(sys_ws, dae, status, 0.0)      # power flow (DAE container), at t = 0
evaluate_system_residuals!(sys_ws, f, g, x, y, status, t)  # differential-algebraic
```

### Parameter Modification
```julia
sys.models.PQ.pl[3] = 0.5   # Visible at the next evaluation
sol = solve_powerflow_newton!(ws, y0)
```
