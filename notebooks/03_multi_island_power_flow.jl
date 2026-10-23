# # Multi-Island Power Flow
#
# Interactive walkthrough of solving power flow for disconnected islands. The API is transparent - you don't need to manually track islands!

# ## Setup

using Powerful
using PowerfulCases
using DataFrames

# ## 1. Load System with Islands

## The IEEE 14-bus case ships an "island" variant of its RAW file in which
## the network is split into disconnected parts
case_bundle = load("ieee14")
case = parse_network(file(case_bundle, :raw; variant="island"))
sys = SystemModel(case)

# ## 2. Check System Topology (Optional)

## The solver automatically detects islands, but we can inspect them manually
result = split_islands(sys)

if result isa SystemModel
    println("System is fully connected (1 island)")
    islands = [result]
else
    islands = result
    println("System has $(length(islands)) disconnected islands")
end

# ## 3. Inspect Individual Islands

for (idx, island) in enumerate(islands)
    buses = island.models[:Bus]
    slack_idx = findfirst(==(BUS_SLACK), buses.ide)  # bus type code of the slack bus

    println("\nIsland $idx:")
    println("  Buses: $(length(buses))")
    println("  Bus numbers: $(Int.(buses.i))")
    if slack_idx !== nothing
        println("  Slack bus: $(buses.i[slack_idx])")
    end
end

# ## 4. Solve Power Flow - Transparent API
#
# **Key insight:** The API is identical whether you have one island or multiple islands!
# Just call `solve()` - the solver handles island detection automatically.

prob = ACPowerFlowProblem(sys)  # islands are detected automatically (auto_split=true)
sol = solve(prob, Newton())

# ## 5. Check Solution Status

if sol isa MultiIslandSolution
    println("✓ Multi-island solution")
    println("  Islands: $(length(sol.solutions))")
    println("  All converged: $(sol.converged)")
    println("  Max iterations: $(sol.iterations)")
    println("  Max residual: $(sol.residual_norm)")

    ## Show per-island details
    println("\n  Per-island results:")
    for (idx, island_sol) in enumerate(sol.solutions)
        println("    Island $idx: converged=$(island_sol.converged), iterations=$(island_sol.iterations), residual=$(island_sol.residual_norm)")
    end
else
    println("✓ Single-island solution")
    println("  Converged: $(sol.converged)")
    println("  Iterations: $(sol.iterations)")
end

# ## 6. Access Variables - Transparent API
#
# **No need to track which island contains which bus!** Just use the bus number.

## Access bus voltages from different islands
v1 = get_var(sol, :Bus, :v, (1,))  # Bus 1 is in island 1
println("Voltage at bus 1: $(round(v1, digits=4)) pu")

v8 = get_var(sol, :Bus, :v, (8,))  # Bus 8 is in island 2
println("Voltage at bus 8: $(round(v8, digits=4)) pu")

# ## 7. Get All Voltages (Merged from All Islands)

## get_vars_dict automatically merges results from all islands
voltages = get_vars_dict(sol, :Bus, :v)
println("Total buses: $(length(voltages))")

println("\nSample voltages (first 5 buses):")
for bus_key in sort(collect(keys(voltages)))[1:5]
    ## Keys are tuples: (bus_number,)
    println("  Bus $(bus_key[1]): $(round(voltages[bus_key], digits=4)) pu")
end

# ## 8. Generate Reports
#
# Reports work identically for single or multi-island systems.

## Collect bus results from all islands
bus_data = collect_model_table(sol, :Bus)
df = DataFrame(bus_data)

## Show the first 8 buses (both islands are included)
first(df, 8)

# ## 9. Key Takeaways
#
# 1. **Automatic island detection**: Just call `ACPowerFlowProblem(sys)`
# 2. **Transparent API**: No manual island tracking needed
# 3. **Unified solution**: `get_var()`, `get_vars_dict()` work across all islands
# 4. **Same code**: Works identically for single or multi-island systems
#
# **Recommended pattern:**
#
# ```julia
# # This works for any system (single or multi-island)
# prob = ACPowerFlowProblem(sys)
# sol = solve(prob, Newton())
#
# # Access variables transparently
# v = get_var(sol, :Bus, :v, (bus_number,))
# voltages = get_vars_dict(sol, :Bus, :v)
# ```

# ## Next Steps
#
# - **[Getting Started Guide](../guides/getting_started.md)** - Basic workflow
# - **[Power Flow Analysis](../guides/powerflow.md)** - Formulations and algorithms
# - **[Accessing Results](../guides/accessing_variables.md)** - Complete API reference
