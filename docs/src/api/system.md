# Core: System (PowerfulCore)

System construction, per-unit conversion, configuration and topology. Islands and
outages are on [Islands and Outages](outages.md). These names are defined in `PowerfulCore.System` and exported by
`PowerfulCore`; most are re-exported by `Powerful`.

## Module

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/System.jl"]
Private = false
```

## System Construction

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/system_builder.jl"]
Private = true
Filter = t -> t === SystemModel
```

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/assembly.jl", "system/model_utils.jl"]
Private = false
```

## System Buffers

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/system_buffers.jl"]
Private = false
```

## Per-Unit Conversion

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/per_unit.jl"]
Private = false
```

## Configuration

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/config_types.jl"]
Private = false
```

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/config.jl"]
Private = false
```

## Topology

```@autodocs
Modules = [PowerfulCore.System]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["system/topology.jl"]
Private = false
```
