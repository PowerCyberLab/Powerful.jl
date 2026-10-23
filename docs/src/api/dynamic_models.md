# Dynamic Models (PowerfulTDS)

Dynamic models for time-domain simulation, parsed from PSS/E DYR files with
`parse_dynamics`. Generators and converter models replace the power flow generator
with the same bus number and ID; exciters, governors and electrical controllers attach
to the device they control. Import model types from `Powerful.Models`.

## Synchronous Generators

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/docstrings.jl"]
Private = false
Filter = t -> t in (GENCLS, GENROU, GENSAL, GENTPJ1)
```

## Exciters

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/docstrings.jl"]
Private = false
Filter = t -> t === SEXS
```

## Governors

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/docstrings.jl"]
Private = false
Filter = t -> t === TGOV1
```

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["governors/tgov1.jl"]
Private = false
```

## Renewable Energy Models

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/docstrings.jl"]
Private = false
Filter = t -> t in (AbstractREGCA1, REGCA1, REGFM_A1, REECA1)
```

## Sources

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/docstrings.jl"]
Private = false
Filter = t -> t === VoltagePlayback
```

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/sources/voltage_playback.jl"]
Private = false
```

## Faults

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["models/fault.jl"]
Private = false
```

## DYR Parsing

```@autodocs
Modules = [PowerfulTDS]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["formats/psse_dyr.jl"]
Private = false
```
