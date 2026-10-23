# Contingency Analysis (PowerfulCA)

AC contingency analysis and DC screening with PTDF and LODF. PowerfulCA is a
standalone package: load it with `using PowerfulCA`. Contingency and outage types
(`Contingency`, `BranchOutage`, ...) are documented under
[Outages and Contingencies](outages.md#Outages-and-Contingencies).

## Module

```@autodocs
Modules = [PowerfulCA]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulCA/src/PowerfulCA.jl"]
Private = false
```

## AC Contingency Analysis

```@autodocs
Modules = [PowerfulCA]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/accc.jl", "src/report.jl"]
Private = false
```

## DC Contingency Screening

```@autodocs
Modules = [PowerfulCA]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/dccc.jl"]
Private = false
```

## Power Transfer Distribution Factors

```@autodocs
Modules = [PowerfulCA]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/ptdf.jl"]
Private = false
```

## Line Outage Distribution Factors

```@autodocs
Modules = [PowerfulCA]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/lodf.jl"]
Private = false
```
