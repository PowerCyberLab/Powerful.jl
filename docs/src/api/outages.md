# Core: Islands and Outages (PowerfulCore)

Splitting a system into islands, and the contingency and outage types used by
contingency analysis and time-domain events. These names are defined in
`PowerfulCore.System` and exported by `PowerfulCore`.

## Multi-Island Systems

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/multi_island.jl"]
Private = false
```

## Outages and Contingencies

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/outages.jl"]
Private = false
```
