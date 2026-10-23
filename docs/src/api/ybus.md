# Y-Bus Time-Domain Simulation (PowerfulYbus)

Time-domain simulation with the network represented by its sparse admittance matrix.
PowerfulYbus is a standalone package: load it with `using PowerfulYbus`. The network
solution, machine sources and system buffers are on
[Y-Bus Network Solution](ybus_network.md), the model traits on
[Y-Bus Model Interface](ybus_traits.md).

## Module

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["PowerfulYbus/src/PowerfulYbus.jl"]
Private = false
```

## Problem and Solvers

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["integration/ybus_tds_problem.jl", "integration/ybus_residual.jl"]
Private = false
```

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["integration/ybus_solve.jl"]
Private = true
Filter = t -> t === solve
```

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["integration/trapezoidal_solver.jl"]
Private = false
```

## Events and Topology Changes

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["integration/ybus_event_processor.jl", "integration/ybus_event_schedule.jl", "ybus/topology_cache.jl"]
Private = false
```

## Y Matrix

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["ybus/ybus_builder.jl"]
Private = false
```
