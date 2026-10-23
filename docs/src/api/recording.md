# Core: Time Series and Injection Recording (PowerfulCore)

Time series solutions of time-domain simulations and the recording of per-model
contributions to shared residuals. These names are defined in `PowerfulCore.DAEInfra`
and exported by `PowerfulCore`.

## Time Series Solutions

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/tds_solution.jl"]
Private = false
```

## Injection Recording

```@autodocs
Modules = [PowerfulCore.DAEInfra]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["dae/injection_recorder.jl"]
Private = false
```

```@autodocs
Modules = [PowerfulCore]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulCore/src/PowerfulCore.jl"]
Private = false
Filter = t -> t === PowerfulCore.InjectionRecorder
```
