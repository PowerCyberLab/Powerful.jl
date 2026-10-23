# Validation and Small-Signal Analysis (PowerfulTDS)

Validation of the initialized state, eigenvalue analysis and linearization for
frequency response, computed from an initialized time-domain simulation problem.

## Initialization Validation

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/validation.jl"]
Private = false
```

## Eigenvalue Analysis

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/eigenvalue.jl"]
Private = false
```

## Eigenvalue Sensitivities

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/eigen_sensitivity.jl"]
Private = false
```

## Linearization and Frequency Response

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["solvers/bode.jl", "solvers/linearization_utils.jl"]
Private = false
```
