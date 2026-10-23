# Core: Types and Variables (PowerfulCore)

Model supertype, variable and residual declarations and address spaces shared by every model and solver. Metadata and traits are on
[Metadata and Traits](metadata.md); per-unit kinds are on
[Formulations, Units and Registry](formulations.md). These names are defined in
`PowerfulCore.Foundation` and exported by `PowerfulCore`.

## Modules

```@autodocs
Modules = [PowerfulCore, PowerfulCore.Foundation, PowerfulCore.DiffEq]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulCore/src/PowerfulCore.jl", "foundation/Foundation.jl"]
Private = false
Filter = t -> t isa Module
```

## Types and Variables

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/core_types.jl", "foundation/model_var.jl"]
Private = false
Filter = t -> !(t isa Type && (t <: PUType || t <: RatingType || t <: ReferenceType || t <: RetrievalPhase))
```
