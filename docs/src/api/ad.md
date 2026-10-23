# Core: Jacobians and Newton (PowerfulCore)

Jacobian sparsity detection, automatic differentiation backends and caches, and the
Newton engine with its linear solver interface. These names are defined in
`PowerfulCore` and `PowerfulCore.DAEInfra`.

## Sparsity Detection

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/sparsity_pattern.jl"]
Private = false
```

## Automatic Differentiation

```@autodocs
Modules = [PowerfulCore]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/ad_backend.jl", "dae/autodiff_cache.jl", "dae/forwarddiff_cache.jl", "dae/enzyme_cache.jl"]
Private = false
```

## Newton Solver

```@autodocs
Modules = [PowerfulCore]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/newton.jl"]
Private = false
```
