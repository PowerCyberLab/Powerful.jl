"""
Show the power flow solution for IEEE 14-bus system
"""

using Powerful
using PowerfulCases
using Printf
using DataFrames

println("="^80)
println("IEEE 14-BUS POWER FLOW SOLUTION")
println("="^80)

# Load case
println("\nLoading IEEE 14-bus case...")
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

# Solve power flow
println("Solving power flow...")
prob = ACPowerFlowProblem(sys)
result = solve(prob, RobustNR())

if result.converged
    println("Converged in ", result.iterations, " iterations")
    println(@sprintf("  Final residual norm: %.6e", result.residual_norm))

    # Show bus results using reporting infrastructure
    println("\n" * "="^80)
    println("BUS RESULTS")
    println("="^80)

    bus_data = collect_model_table(result, :Bus)
    df_bus = DataFrame(bus_data)
    println(df_bus)

    # Show branch flows
    println("\n" * "="^80)
    println("BRANCH FLOWS")
    println("="^80)

    branch_data = collect_model_table(result, :Branch)
    df_branch = DataFrame(branch_data)
    println(df_branch)

    # Show PQ loads
    println("\n" * "="^80)
    println("PQ LOADS")
    println("="^80)

    pq_data = collect_model_table(result, :PQ)
    df_pq = DataFrame(pq_data)
    println(df_pq)

    # Show PV generators
    println("\n" * "="^80)
    println("PV GENERATORS")
    println("="^80)

    pv_data = collect_model_table(result, :PV)
    df_pv = DataFrame(pv_data)
    println(df_pv)

else
    println("Did NOT converge")
    println(@sprintf("  Final residual norm: %.6e", result.residual_norm))
    println("  Iterations: ", result.iterations)
end

println("\n" * "="^80)
