# Island Power Flow Demo
#
# Demonstrates solving power flow for disconnected islands using transparent API
# Users don't need to manually handle islands - the solver does it automatically!

using Powerful
using PowerfulCases

println("=" ^ 70)
println("Island Power Flow Demo - Transparent Multi-Island API")
println("=" ^ 70)

# Load system with multiple islands
println("\n1. Loading system with multiple islands...")
case = parse_network(file(load("ieee14"), :raw; variant = "island"))
sys = SystemModel(case)
println("   ✓ Loaded ieee14_island.raw")

# Check topology (optional - solver does this automatically)
println("\n2. Checking system topology...")
result = split_islands(sys)
if result isa SystemModel
    println("   System is fully connected (1 island)")
else
    islands = result
    println("   System has ", length(islands), " disconnected islands")

    # Inspect individual islands
    println("\n3. Inspecting individual islands...")
    for (idx, island) in enumerate(islands)
        buses = island.models[:Bus]
        slack_bus_idx = findfirst(buses.ide .== 3)
        println("\n   Island $idx:")
        println("     Buses: ", length(buses))
        println("     Bus numbers: ", collect(buses.i))
        println("     Slack bus (ide=3): ", slack_bus_idx === nothing ? "None" : "bus $(buses.i[slack_bus_idx])")
    end
end

# Solve power flow - API is identical regardless of islands!
# The solver automatically detects islands and solves each independently
println("\n4. Solving power flow (auto-handles islands)...")
prob = ACPowerFlowProblem(sys)  # auto_split=true by default
sol = solve(prob, RobustNR())

# Check if we got multi-island solution
if sol isa MultiIslandSolution
    println("   ✓ Solved multi-island system")
    println("   Total islands: ", length(sol.islands))
    println("   All converged: ", sol.converged)
    println("   Max iterations: ", sol.iterations)
    println("   Max residual: ", sol.residual_norm)

    # Show per-island details
    println("\n   Per-island results:")
    for (idx, island_sol) in enumerate(sol.solutions)
        println("     Island $idx: converged=", island_sol.converged,
                ", iters=", island_sol.iterations,
                ", |g|=", round(island_sol.residual_norm, digits=8))
    end
else
    println("   ✓ Solved single-island system")
    println("   Converged: ", sol.converged)
    println("   Iterations: ", sol.iterations)
end

# Access variables - API is identical for single or multi-island!
println("\n5. Accessing variables (transparent API)...")
v1 = get_var(sol, :Bus, :v, (1,))
println("   Voltage at bus 1: ", round(v1, digits=4), " pu")

v8 = get_var(sol, :Bus, :v, (8,))
println("   Voltage at bus 8: ", round(v8, digits=4), " pu")

# Get all voltages across all islands
println("\n6. Getting all voltages (merged from all islands)...")
voltages = get_vars_dict(sol, :Bus, :v)
println("   Total buses: ", length(voltages))
println("\n   Sample voltages (first 5 buses):")
for bus_key in sort(collect(keys(voltages)))[1:5]
    # Unwrap single-element tuple for display
    display_key = bus_key isa Tuple && length(bus_key) == 1 ? bus_key[1] : bus_key
    println("     Bus $display_key: ", round(voltages[bus_key], digits=4), " pu")
end

println("\n" * "=" ^ 70)
println("✓ Key Takeaway: The API is identical for single or multi-island!")
println("  - Just call solve(ACPowerFlowProblem(sys), algorithm)")
println("  - get_var(), get_vars_dict() work transparently")
println("  - No need to manually track islands or solutions!")
println("=" ^ 70)
