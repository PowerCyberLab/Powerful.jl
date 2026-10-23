# Network Models (PowerfulPF)

Static network models used by the power flow and, through replacement, by the
time-domain simulation. Import model types from `Powerful.Models`.

## Buses

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/bus1ph.jl", "models/bus_types.jl"]
Private = false
```

## Branches

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/branch.jl"]
Private = false
```

## Transformers

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/transformer.jl", "models/transformer_conversion.jl"]
Private = false
```

## Loads

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/pq.jl"]
Private = false
```

## Generators

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/pv.jl", "models/slack.jl"]
Private = false
```

## Shunts

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/shunt.jl", "models/switched_shunt.jl"]
Private = false
```

## Areas, Zones and Owners

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["topological/area_interchange.jl", "topological/zone.jl", "topological/owner.jl", "topological/topological.jl"]
Private = false
```

## Model Groups and Post-Processing

```@autodocs
Modules = [PowerfulPF]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/group_definitions.jl", "models/post_processing.jl", "models/dc_physics.jl"]
Private = false
```
