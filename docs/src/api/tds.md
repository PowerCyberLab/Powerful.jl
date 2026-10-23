# Time-Domain Simulation (PowerfulTDS)

Setup, initialization, solution and events for phasor time-domain simulation. DYR
parsing is on [Dynamic Models](dynamic_models.md).
Initialization validation, eigenvalue and frequency response analysis are on
[Validation and Small-Signal Analysis](small_signal.md).

## Module

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulTDS/src/PowerfulTDS.jl"]
Private = false
```

## Problem Setup

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/tds.jl", "solvers/linking.jl"]
Private = false
```

## Solving

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/solve.jl"]
Private = true
Filter = t -> t === solve
```

## Initialization

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/initialization.jl", "solvers/reinit.jl"]
Private = false
```

## Events

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/events.jl", "solvers/initial_limits.jl"]
Private = false
```

## Solver Interface

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/adapter.jl", "solvers/functors.jl", "solvers/index_reduction.jl"]
Private = false
```

## Parameter Overrides

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/parameter_override.jl"]
Private = false
```
