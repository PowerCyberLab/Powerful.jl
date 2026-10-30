# Computation Architecture

Residual evaluation in Newton and time-domain iterations runs on pre-allocated buffers. This page describes how those buffers are organized, created and reused.

## The Workspace Hierarchy

Residual evaluation is organized in three nested layers.

```
SystemBuffers                    # Coordinates all models
├── ResidualBuffers (PQ)         # Data movement for PQ loads
│   └── Residual kernels         # Pure math
├── ResidualBuffers (PV)         # Data movement for PV generators
│   └── Residual kernels         # Pure math
└── ResidualBuffers (Branch)     # Data movement for branches
    └── Residual kernels         # Pure math
```

### SystemBuffers

`SystemBuffers` is defined in `PowerfulCore/src/system/system_buffers.jl`. It loops over the models in evaluation order (`evaluation_pass` from the metadata) and calls the evaluation of each. It has no knowledge of model-specific calculations.

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

`SystemBuffers` stores one closure per model in a `Vector{Function}`. A NamedTuple of typed workspaces would make every model type a type parameter of the system workspace, and compilation time would grow quickly with the number of models. Each closure captures the concrete types of its model, so the code inside it is fully specialized. Calling a closure through the `Vector{Function}` costs one dynamic dispatch per model per evaluation.

### ResidualBuffers

`ResidualBuffers` is defined in `PowerfulCore/src/dae/residual_workspace.jl`. It gathers the variables of one model, passes them to the kernel and scatters the results; it contains no model mathematics. For each model it holds:

- `var_addrs`, `var_vals`: NamedTuples of DAE addresses and gathered values per variable (e.g. `theta`, `v`)
- `res_addrs`, `res_vals`: NamedTuples of DAE addresses and computed values per residual (e.g. `p`, `q`)
- `var_gather_tuple`, `res_scatter_tuple`: precomputed tuples used by gather and scatter
- `n_instances`: the instance count
- `setpoints`: initialization setpoints of dynamic models (an empty NamedTuple for other models)
- `system`: the `SystemConfig` of the system (`sys.settings`: frequency, MVA base, constant-power breakpoint)

The in-service mask is passed to `calc_residuals!` and `scatter_residuals!` as `active_mask`. It is taken from a status dict (`Dict{Symbol, Vector{Bool}}`, built by `build_status_dict(sys)`). Parallel contingency analysis relies on this arrangement, since each thread can evaluate the same workspace structure with its own status.

#### Out-of-service instances

An out-of-service instance contributes nothing to the residuals. Its own algebraic variables remain entries of `y`; these are the variables without a source model, such as the `theta` and `v` of a bus or the `P_m` of a machine. No other model writes their residual rows, so the following mechanism keeps the Jacobian nonsingular:

1. The own algebraic residuals of a model share addresses with its own algebraic variables. This is checked the first time an instance of the model is out of service.
2. After all models are evaluated, `hold_out_of_service!` sets each such row `a` to `y[a] - hold[a]`. The row is an identity row of the Jacobian that keeps the variable at its held value.
3. The solver sets `hold` with `set_hold!` from its starting point: the initial guess of a power flow or the initialized values of a time-domain simulation.

This mechanism keeps a bus de-energized by a contingency, or a machine tripped during a simulation, from making the Jacobian singular. Own states whose mass coefficient reduces to zero for a particular instance are algebraic for that instance, and their rows `f[a]` are held in the same way at `hold_x`. States with mass keep `M·ẋ = 0`, which holds their values.

Models therefore only mask their equations with the in-service mask. Each evaluator records whether it has out-of-service instances, so the pass costs nothing when every instance is in service.

### Residual Kernels

A residual kernel is a pure function that loops over arrays and computes residuals from parameters and variables. It contains the physics equations of one model and has no access to workspaces, addresses or the system structure.

```julia
kernel!(outputs::NamedTuple, params::NamedTuple, vars::NamedTuple, active, AT)
```

## The Gather-Compute-Scatter Pattern

The evaluator closure of each model runs four calls in sequence:

```julia
gather_variables!(ws, x, y)              # y (and x for states) → ws.var_vals
gather_residuals!(ws, f, g)              # accumulated values of SharedRes(:set) residuals
calc_residuals!(ws, model, external, formulation, t; active_mask)
scatter_residuals!(ws, f, g, active_mask)
```

Here `x` and `f` are the differential states and residuals (empty in power flow), and `y` and `g` are the algebraic variables and residuals. The calls form three phases.

**Gather.** Each variable of the model is copied from the global DAE vectors to the workspace, which yields contiguous arrays for SIMD loops or GPU kernels. State variables are gathered from `x` and algebraic variables from `y`.

```julia
# For each variable of the model, e.g. the bus voltage seen by a PQ load:
ws.var_vals.v[i] = y[ws.var_addrs.v[i]]
```

**Compute.** The kernel operates on the workspace arrays and fills the residual values. The PQ load kernel from `PowerfulPF/src/models/pq.jl` is:

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

**Scatter.** The residual values are copied from the workspace to the global residual vectors, for active instances only. `SharedRes(:add)` and `PrivateRes` residuals are accumulated or written first. `SharedRes(:set)` residuals overwrite in a second pass, after all `:add` contributions of the model.

```julia
# SharedRes(:add) and PrivateRes: accumulate / write
g[ws.res_addrs.p[i]] += ws.res_vals.p[i]

# SharedRes(:set): overwrite
g[ws.res_addrs.v_ctrl[i]] = ws.res_vals.v_ctrl[i]
```

## NamedTuples in the Workspace

The variable and residual arrays of a workspace are stored in NamedTuples. This gives type stability (`ws.var_vals.v` has a known array type), compiles field access to a direct load as for a struct field, and lets the compiler generate a specialized code path for each model:

```julia
calc_residuals!(ws::ResidualBuffers{Float64, Vector{Float64}, ...}, pq::StructArray{PQ{Float64}}, ...)
```

The field names are fixed when the workspace is created, at setup, from the formulation-specific metadata.

## CPU and GPU Dispatch

CPU and GPU implementations of a kernel share one interface. The last argument `AT` is the workspace array type.

The CPU kernel is a method on `::Type{<:Vector}`. In the model file (e.g. `PowerfulPF/src/models/pq.jl`):

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

The GPU version consists of a device function and a launcher method of the same kernel, defined for the array type of the GPU extension (`::Type{<:CUDA.CuArray}`):

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

The physics function `calc_zip_load` is shared by both versions.

The array type comes from the system (`sys.array_type`). Moving the system or problem to another array type selects the matching kernels through ordinary method dispatch, and the workspace needs no option for it. A GPU power flow is set up with `use_gpu` selecting the GPU linear solver:

```julia
ws_gpu = setup_powerflow_newton(to_gpu(prob); use_gpu = true)
```

## How Workspaces Are Created

During solver setup (`setup_powerflow_newton`, or `solve` for a time-domain problem), the autodiff cache (`create_adcache`) builds the `SystemBuffers`. Construction proceeds in two levels.

`SystemBuffers` is constructed from the address manager, the system and the formulation:

```julia
using PowerfulCore: SystemBuffers

sys_ws = SystemBuffers(Float64, am, sys, formulation)
```

1. The models are ordered by `evaluation_pass`; models without instances, variables or residuals are skipped.
2. The formulation-specific metadata of each model is retrieved.
3. A `ResidualBuffers` is created for each model.
4. The evaluator closures are built for the hot path, the recorder path and the accumulate path.

`ResidualBuffers` is constructed for each model:

```julia
using PowerfulCore: ResidualBuffers

ws = ResidualBuffers(Float64, metadata, am, model_vec, sys.array_type;
                     system = sys.settings, setpoints = NamedTuple(), config = SimConfig())
```

1. Arrays are allocated for each variable and residual.
2. The address vectors are extracted from the `AddressManager`.
3. Arrays are allocated for observable variables and setpoints, or the given setpoints are referenced.
4. The gather and scatter tuples are precomputed.

The array type flows through the hierarchy as a type parameter:

```julia
ResidualBuffers{Float64, Vector{Float64}, ...}    # CPU
ResidualBuffers{Float64, CuVector{Float64}, ...}  # GPU
```

The workspaces and metadata of a model can be reached by name on the cold path:

```julia
sys_ws.metadata_by_name[:PQ]    # ModelMetadata for PQ (formulation-specific)
sys_ws.workspaces_by_name[:PQ]  # ResidualBuffers for PQ
```

After setup the array storage is allocated once. Iterations fill the existing buffers, and all indexing uses fixed address vectors. A residual evaluation outside a solver, for example in a test, is made with `evaluate_system_residuals!`:

```julia
using PowerfulCore: SystemBuffers, build_status_dict, evaluate_system_residuals!

sys_ws = SystemBuffers(Float64, prob.address_manager, prob.sys, prob.formulation)
status = build_status_dict(prob.sys)
evaluate_system_residuals!(sys_ws, dae, status, 0.0)      # power flow (DAE container), at t = 0
evaluate_system_residuals!(sys_ws, f, g, x, y, status, t)  # differential-algebraic
```

## How Model References Flow Through Workspaces

The problem holds a reference to the system and the address manager. The `SystemBuffers` holds the evaluator closures and, for cold-path access, the workspaces by name.

```
ACPowerFlowProblem
  ├── sys::SystemModel         # Reference to the system
  └── address_manager

SystemBuffers
  ├── evaluators::Vector{Function}  # Closures capture model_vec references
  └── workspaces_by_name::Dict      # Cold path access

Each evaluator closure captures:
  ├── ws::ResidualBuffers           # Pre-allocated workspace
  ├── model_vec::StructArray        # Reference to e.g. sys.models.PQ
  └── external::NamedTuple          # sys.external for that model
```

Because each closure captures a reference to the StructArray of its model, a model field changed in place (`sys.models.PQ.pl[3] = 0.5`) is visible at the next residual evaluation. A parameter sweep therefore reuses one workspace:

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

Only in-place changes to the arrays held by the closures are seen. Replacing a model collection, adding instances or changing the topology requires a new problem and workspace. External parameters (`sys.external`) are retrieved at system construction and are not recomputed when model fields change.

## Performance Characteristics

Gather, compute and scatter use buffers that are allocated at setup, and no array storage is allocated per iteration. Tests in `Powerful/test/performance/` and `Powerful/test/perf/` guard against allocation regressions.

Gather copies the scattered entries of the DAE vectors into contiguous arrays, so kernels access memory sequentially and SIMD processes several elements per instruction. The same property motivates the storage of models as StructArrays.

Inside each evaluator closure every call is type-stable, and the code is specialized for the concrete model, workspace and formulation types. The call to the closure itself is the only dynamic dispatch.

With GPU arrays, gather, compute and scatter run as device kernels with one thread per element. The scatter of `SharedRes(:add)` residuals uses atomic adds, since several devices contribute to the same bus equation.

## Custom Gather and Scatter

The default gather and scatter loop over address arrays. Specialized versions for other array types are added by dispatching on the array type in the same way.
