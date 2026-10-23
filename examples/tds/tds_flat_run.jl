"""
TDS flat run example - demonstrates state trajectory extraction API
"""

using Powerful
using PowerfulCases
using PowerfulCore: get_series
using CSV, DataFrames

# Load and solve power flow
case = parse_network(load("ieee14").raw)
sys_pf = SystemModel(case)
pf_prob = ACPowerFlowProblem(sys_pf)
pf_sol = solve(pf_prob, Newton())

# Parse dynamics and create TDS problem
dyr_data = parse_dynamics(load("ieee14").dyr, sys_pf)
init = setup_tds(pf_sol, dyr_data)
prob = TDSProblem(init; tspan=(0.0, 1.0))

# Solve TDS
sol = solve(prob)

println("TDS simulation completed")
println("  Return code: $(sol.retcode)")
println("  Time steps: $(length(sol.t))")
println("  Final time: $(sol.t[end])")

# Extract generator state variables using get_series
# Get all GENCLS omega values (returns NamedTuple with t and all instances)
data_omega = get_series(sol, :GENCLS, :omega)
data_delta = get_series(sol, :GENCLS, :delta)

# Merge into single DataFrame
results = DataFrame(data_omega)
for (k, v) in pairs(data_delta)
    k == :t && continue  # Skip duplicate time column
    results[!, k] = v
end

# Write to CSV
csv_file = joinpath(@__DIR__, "tds_flat_run.csv")
CSV.write(csv_file, results)
println("\nResults written to: $csv_file")
println("Columns: ", join(names(results), ", "))
