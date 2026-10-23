# Y-Bus Model Interface (PowerfulYbus)

Traits and functions through which models take part in the Y-bus formulation: their
network partition role, nonconforming load type, admittance contributions and
residuals.

## Network Partition Roles

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["traits/network_partition_traits.jl"]
Private = false
```

## Nonconforming Loads

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["traits/nc_load_traits.jl"]
Private = false
```

## Admittances and Residuals

```@autodocs
Modules = [PowerfulYbus]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["traits/ybus_interface.jl", "traits/ybus_residuals.jl"]
Private = false
```
