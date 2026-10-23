# Core: DAE and Addressing (PowerfulCore)

DAE vectors, address management, indexing and residual evaluation. These names are
defined in `PowerfulCore.DAEInfra` and exported by `PowerfulCore`. Jacobians and the
Newton engine are on [Jacobians and Newton](ad.md); time series solutions and
injection recording on [Time Series and Injection Recording](recording.md).

## Module

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/DAEInfra.jl"]
Private = false
```

## DAE Vectors

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/vectors.jl"]
Private = false
```

## Address Management

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/address_manager.jl"]
Private = false
```

## External-to-Internal Indexing

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/ext2int.jl"]
Private = false
```

## Gather and Scatter

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/gather_scatter.jl"]
Private = false
```

## Residual Buffers

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/residual_workspace.jl"]
Private = false
```

## Automaton State

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/automaton_state.jl"]
Private = false
```
