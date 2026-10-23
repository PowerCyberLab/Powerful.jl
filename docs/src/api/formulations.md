# Core: Formulations, Units and Registry (PowerfulCore)

Formulation and data format traits, per-unit kinds, bus type codes, the solution
interface and the model registry. These names are defined in
`PowerfulCore.Foundation` (or `PowerfulCore` itself) and exported by `PowerfulCore`.

## Formulations

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/formulation_traits.jl"]
Private = false
```

## Data Formats

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/format_traits.jl"]
Private = false
```

## Per-Unit Kinds

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/core_types.jl"]
Private = false
Filter = t -> t isa Type && (t <: PUType || t <: RatingType || t <: ReferenceType || t <: RetrievalPhase)
```

## Bus Type Codes

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/constants.jl"]
Private = false
```

## System Settings

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/system_config.jl"]
Private = false
```

## Solution Interface

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/solution_interface.jl"]
Private = false
```

## Model Registry

```@autodocs
Modules = [PowerfulCore]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulCore/src/registry.jl"]
Private = false
```

## Verbose Output

```@autodocs
Modules = [PowerfulCore]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulCore/src/PowerfulCore.jl"]
Private = false
Filter = t -> !(t isa Module) && t !== PowerfulCore.InjectionRecorder
```
