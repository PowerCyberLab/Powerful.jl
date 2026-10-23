# Core: Metadata and Traits (PowerfulCore)

Model metadata, traits, groups and base quantities. Formulations, data formats and
the model registry are on [Formulations, Units and Registry](formulations.md). These
names are defined in `PowerfulCore.Foundation` (or `PowerfulCore` itself) and
exported by `PowerfulCore`.

## Model Metadata

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/model_metadata.jl", "foundation/metadata_queries.jl"]
Private = false
```

## Automata

The discrete behaviour of a model (modes, transitions, memory), declared in
`ModelMetadata.automata` and `ModelMetadata.memory`. These names are public in
`PowerfulCore` and not exported: import them by name (`using PowerfulCore: Automaton, Mode, ...`).

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/automata.jl"]
Private = true
```

## Model Traits

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/model_traits.jl"]
Private = false
```

## Model Groups

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/groups.jl"]
Private = false
```

## Base Quantities

```@autodocs
Modules = [PowerfulCore.Foundation]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["foundation/base_quantity_traits.jl"]
Private = false
```
