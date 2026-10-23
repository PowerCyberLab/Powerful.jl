#!/usr/bin/env julia
"""
Minimal IEEE 14-bus power flow example with default (contiguous) addressing.
"""

using Powerful
using PowerfulCases
using DataFrames

# Load case and solve
println("Starting IEEE 14-bus power flow example...")

println("Parsing case...")
case = parse_network(load("ieee14").raw)

println("Loading case...")
sys = SystemModel(case)

println("Solving power flow...")
prob = ACPowerFlowProblem(sys)
result = solve(prob, Newton())

# Display results
if result.converged
    println("✓ Converged in $(result.iterations) iterations (residual: $(result.residual_norm))\n")

    # Collect and display bus results using reporting infrastructure
    bus_data = collect_model_table(result, :Bus)
    df = DataFrame(bus_data)
    println(df)

    branch_data = collect_model_table(result, :Branch)
    df = DataFrame(branch_data)
    println(df)
else
    println("✗ Did not converge ($(result.iterations) iterations)")
end
