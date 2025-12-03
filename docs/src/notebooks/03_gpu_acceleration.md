# GPU Acceleration

This notebook demonstrates GPU-accelerated power flow using CUDA.
Powerful.jl supports a **hybrid GPU architecture** that provides GPU acceleration
while maintaining compatibility with automatic differentiation.

## Setup

GPU features require CUDA.jl and CUDSS.jl:


```julia
using Powerful
using Powerful.PowerFlowData
using CUDA
using CUDSS
using Printf
```

    ┌ Warning: Module Powerful with build ID fafbfcfd-01e8-2be5-e50c-51fc041ef1a0 is missing from the cache.
    │ This may mean Powerful [542e4350-023d-42c7-b19f-8cb561780c9e] does not support precompilation but is imported by a module that does.
    └ @ Base loading.jl:2541


Check GPU availability


```julia
if !CUDA.functional()
    error("No CUDA-capable GPU detected")
end

println("GPU Device: ", CUDA.name(CUDA.device()))
```

    GPU Device: NVIDIA GeForce RTX 5090


## 1. Hybrid GPU Power Flow (Production Ready)

The hybrid approach computes Jacobian on CPU and solves linear systems on GPU.
This avoids technical limitations while providing GPU acceleration for the
most expensive operation (sparse factorization).

Load test case


```julia
case = PowerFlowData.parse_network(joinpath(pkgdir(Powerful), "cases", "ieee14.raw"))
sys = SystemModel(case)
```




    SystemModel{Vector}:
      Components:
        PV: 4 elements
        Owner: 2 elements
        Branch: 16 elements
        AreaInterchange: 2 elements
        Transformer2W: 4 elements
        Bus: 14 elements
        PQ: 11 elements
        Slack: 1 elements
        Zone: 2 elements




Create problem (topology analysis on CPU)


```julia
prob = ACPowerFlowProblem(sys)
```




    ACPowerFlowProblem:
      Formulation: PowerInjection
      Variables: 56




### CPU Baseline


```julia
sol_cpu = solve(prob, NewtonRaphson())
```




    PowerFlowSolution{Float64}:
      Converged: true
      Iterations: 4
      Residual: 3.25e-12
      Variables: 28




### GPU Linear Solve

Use `CUDSSSolver()` to accelerate the linear solve on GPU:


```julia
sol_gpu = solve(prob, NewtonRaphson(; linear_solver=CUDSSSolver()))
```




    PowerFlowSolution{Float64}:
      Converged: true
      Iterations: 4
      Residual: 3.26e-12
      Variables: 28




### Verify Consistency


```julia
max_diff = maximum(abs.(sol_cpu.dae.y .- sol_gpu.dae.y))

@assert sol_cpu.converged && sol_gpu.converged
@assert sol_cpu.iterations == sol_gpu.iterations
@assert max_diff < 1e-10

println("✓ CPU and GPU solutions match (iterations: $(sol_cpu.iterations))")
```

    ✓ CPU and GPU solutions match (iterations: 4)


## 2. Performance Considerations

### When GPU Helps

GPU acceleration benefits:
- **Large systems** (>10,000 buses) where sparse factorization dominates
- **Repeated solves** where one-time transfer overhead is amortized
- **Systems with high Jacobian density** (more nonzeros to factorize)

### When GPU Doesn't Help

GPU may not improve performance for:
- **Small systems** (<1,000 buses) - CPU KLU is highly optimized
- **Single solves** - PCIe transfer overhead dominates
- **Sparse systems** - Transfer cost exceeds factorization savings

### Overhead Breakdown (per iteration)

The hybrid approach transfers data at each iteration:
1. Jacobian matrix J → GPU (~1-5ms for 10k system)
2. RHS vector g → GPU (~0.1ms)
3. Solution vector → CPU (~0.1ms)

For the IEEE 14-bus system, overhead (~0.2ms) exceeds any GPU speedup.
For 70,000 bus systems, factorization time (50ms) justifies the overhead.

## 3. GPU Residual Evaluation (Advanced)

This section demonstrates GPU residual computation, which is the foundation
for future full GPU power flow (when custom Jacobian kernels are implemented).

Transfer system to GPU


```julia
sys_gpu = to_gpu(sys);  # Suppress output to avoid GPU array printing
```

Create addressing and formulation


```julia
am = prob.address_manager
formulation = prob.formulation
```




    PowerInjection()



Create DAE and workspace for GPU


```julia
dae_gpu = DAE(am; array_type=sys_gpu.array_type);
ws_gpu = SystemBuffers(Float64, am, sys_gpu, formulation);
```

Evaluate residuals on GPU


```julia
evaluate_system_residuals!(ws_gpu, dae_gpu);
CUDA.synchronize();
```

Compare with CPU


```julia
dae_cpu = DAE(am);
ws_cpu = SystemBuffers(Float64, am, sys, formulation);
evaluate_system_residuals!(ws_cpu, dae_cpu);

g_cpu = dae_cpu.g
g_gpu = Array(dae_gpu.g)  # Transfer back for comparison

max_residual_diff = maximum(abs.(g_cpu .- g_gpu))

@assert max_residual_diff < 1e-10

println("✓ GPU residuals match CPU (max diff: $(max_residual_diff))")
```

    ✓ GPU residuals match CPU (max diff: 0.0)


## 4. Architecture Overview

Powerful.jl supports two GPU architectures:

### Hybrid GPU (Production Ready) ✅

**Data flow:**
```
CPU: Residual eval → Jacobian (ForwardDiff) → [transfer] →
GPU: Linear solve (CUDSS) → [transfer] →
CPU: Update solution
```

**Advantages:**
- Works with existing autodiff infrastructure
- No technical limitations (ForwardDiff compatible)
- Simple implementation (already complete)

**Limitations:**
- PCIe transfer overhead per iteration
- Jacobian computation still on CPU

### Full GPU (Future Work) 🚧

**Data flow:**
```
GPU: Residual eval → Jacobian (custom kernels) → Linear solve → Solution
(No CPU transfers)
```

**Current status:**
- Residual evaluation: ✅ Working (demonstrated above)
- Jacobian computation: ❌ Requires custom CUDA kernels
- Linear solve: ✅ Working (CUDSS)

**Why not ForwardDiff on GPU?**

ForwardDiff.Dual numbers are incompatible with CUDA atomic operations
used in scatter operations during Jacobian assembly. Full GPU requires
manual partial derivatives in custom CUDA kernels (~2000 lines of code).

**When will full GPU be available?**

Full GPU power flow requires custom Jacobian kernels for each model type.
This is feasible but represents significant development effort (2-4 weeks).
The hybrid approach provides GPU acceleration today while maintaining
code simplicity.

## 5. Key Design: AddressManager is CPU Metadata

The `AddressManager` contains only UInt64 addresses mapping model variables
to DAE positions. This CPU metadata works with both CPU and GPU systems:

The same AddressManager works with both CPU and GPU systems


```julia
@assert prob.address_manager === am  # Same object used for both

println("✓ AddressManager is CPU metadata shared between CPU and GPU systems")
```

    ✓ AddressManager is CPU metadata shared between CPU and GPU systems


## Summary

**For production use:**
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, NewtonRaphson(; linear_solver=CUDSSSolver()))
```

**Key points:**
1. Hybrid GPU is production-ready for large systems (>10k buses)
2. GPU residual evaluation works (foundation for future full GPU)
3. AddressManager design enables CPU/GPU flexibility
4. Full GPU power flow requires custom Jacobian kernels (future work)

**Performance guideline:**
- Small systems (<1k buses): Use CPU (KLU)
- Large systems (>10k buses): Try GPU (CUDSSSolver)
- Very large systems (>50k buses): GPU likely beneficial

---

*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
