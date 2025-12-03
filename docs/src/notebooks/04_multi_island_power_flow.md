# Multi-Island Power Flow

Interactive walkthrough of solving power flow for disconnected islands. The API is transparent - you don't need to manually track islands!

## Setup


```julia
using Powerful
using Powerful.PowerFlowData
using DataFrames
```

## 1. Load System with Islands


```julia
# Load IEEE 14-bus system modified to have disconnected islands
case = PowerFlowData.parse_network(joinpath(pkgdir(Powerful), "cases", "ieee14_island.raw"))
sys = SystemModel(case)
```




    SystemModel{Vector}:
      Components:
        PV: 4 elements
        Owner: 2 elements
        Branch: 14 elements
        AreaInterchange: 2 elements
        Transformer2W: 2 elements
        Bus: 14 elements
        PQ: 10 elements
        Slack: 2 elements
        Zone: 2 elements




## 2. Check System Topology (Optional)


```julia
# The solver automatically detects islands, but we can inspect them manually
result = split_islands(sys)

if result isa SystemModel
    println("System is fully connected (1 island)")
    islands = [result]
else
    islands = result
    println("System has $(length(islands)) disconnected islands")
end
```

    [ Info: Splitting system into 2 islands
    System has 2 disconnected islands


## 3. Inspect Individual Islands


```julia
for (idx, island) in enumerate(islands)
    buses = island.models[:Bus]
    slack_idx = findfirst(buses.ide .== 3)

    println("\nIsland $idx:")
    println("  Buses: $(length(buses))")
    println("  Bus numbers: $(collect(buses.i))")
    if slack_idx !== nothing
        println("  Slack bus: $(buses.i[slack_idx])")
    end
end
```

    
    Island 1:
      Buses: 6
      Bus numbers: Int32[1, 2, 3, 4, 5, 6]
      Slack bus: 1
    
    Island 2:
      Buses: 8
      Bus numbers: Int32[7, 8, 9, 10, 11, 12, 13, 14]
      Slack bus: 8


## 4. Solve Power Flow - Transparent API

**Key insight:** The API is identical whether you have one island or multiple islands!
Just call `solve()` - the solver handles island detection automatically.


```julia
prob = ACPowerFlowProblem(sys)  ## auto_split=true by default
sol = solve(prob, NewtonRaphson())
```

    [ Info: Splitting system into 2 islands





    MultiIslandSolution{Float64}:
      Islands: 2
      Converged: 2/2
      Max iterations: 4
      Max residual: 1.46e-9
      Total variables: 28




## 5. Check Solution Status


```julia
if sol isa MultiIslandSolution
    println("✓ Multi-island solution")
    println("  Islands: $(length(sol.solutions))")
    println("  All converged: $(sol.converged)")
    println("  Max iterations: $(sol.iterations)")
    println("  Max residual: $(round(sol.residual_norm, digits=8))")

    # Show per-island details
    println("\n  Per-island results:")
    for (idx, island_sol) in enumerate(sol.solutions)
        println("    Island $idx: converged=$(island_sol.converged), iters=$(island_sol.iterations)")
    end
else
    println("✓ Single-island solution")
    println("  Converged: $(sol.converged)")
    println("  Iterations: $(sol.iterations)")
end
```

    ✓ Multi-island solution
      Islands: 2
      All converged: true
      Max iterations: 4
      Max residual: 0.0
    
      Per-island results:
        Island 1: converged=true, iters=4
        Island 2: converged=true, iters=4


## 6. Access Variables - Transparent API

**No need to track which island contains which bus!** Just use the bus number.


```julia
# Access bus voltages from different islands
v1 = get_var(sol, :Bus, :v, (1,))  ## Bus 1 is in island 1
println("Voltage at bus 1: $(round(v1, digits=4)) pu")

v8 = get_var(sol, :Bus, :v, (8,))  ## Bus 8 is in island 2
println("Voltage at bus 8: $(round(v8, digits=4)) pu")
```

    Voltage at bus 1: 1.03 pu
    Voltage at bus 8: 1.03 pu


## 7. Get All Voltages (Merged from All Islands)


```julia
# get_vars_dict automatically merges results from all islands
voltages = get_vars_dict(sol, :Bus, :v)
println("Total buses: $(length(voltages))")

println("\nSample voltages (first 5 buses):")
for bus_key in sort(collect(keys(voltages)))[1:5]
    # Unwrap single-element tuple for display
    display_key = bus_key isa Tuple && length(bus_key) == 1 ? bus_key[1] : bus_key
    v = voltages[bus_key]
    println("  Bus $display_key: $(round(v, digits=4)) pu")
end
```

    Total buses: 14
    
    Sample voltages (first 5 buses):
      Bus 1: 1.03 pu
      Bus 2: 1.03 pu
      Bus 3: 1.01 pu
      Bus 4: 1.0149 pu
      Bus 5: 1.0223 pu


## 8. Generate Reports

Reports work identically for single or multi-island systems.


```julia
# Collect bus results from all islands
bus_data = collect_model_table(sol, :Bus)
df = DataFrame(bus_data)

# Show first few buses
println("\nBus results (first 8 buses):")
println(first(df, 8))
```

    
    Bus results (first 8 buses):
    8×3 DataFrame
     Row │ I       Theta (deg)  V (pu)
         │ String  Float64      Float64
    ─────┼──────────────────────────────
       1 │ 1          0.0       1.03
       2 │ 2         -0.699549  1.03
       3 │ 3         -1.83086   1.01
       4 │ 4         -2.17056   1.01495
       5 │ 5         -1.54592   1.02231
       6 │ 6          0.504863  1.03
       7 │ 7         -3.67755   1.0123
       8 │ 8         -1.3945    1.03


## 9. Key Takeaways

1. **Automatic island detection**: Just call `ACPowerFlowProblem(sys)`
2. **Transparent API**: No manual island tracking needed
3. **Unified solution**: `get_var()`, `get_vars_dict()` work across all islands
4. **Same code**: Works identically for single or multi-island systems

**Recommended pattern:**

```julia
## This works for any system (single or multi-island)
prob = ACPowerFlowProblem(sys)
sol = solve(prob, NewtonRaphson())

## Access variables transparently
v = get_var(sol, :Bus, :v, (bus_number,))
voltages = get_vars_dict(sol, :Bus, :v)
```

## Next Steps

- **[Getting Started Guide](../guides/getting_started.md)** - Basic workflow
- **[Power Flow Analysis](../guides/powerflow.md)** - Formulations and algorithms
- **[Accessing Results](../guides/accessing_variables.md)** - Complete API reference

---

*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
