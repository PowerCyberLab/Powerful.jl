"""
TDS run example for ACTIVSg2000 (545 generators, 2000 buses)

Usage:
    julia --project examples/tds/tds_activsg2000_run.jl
"""

using Powerful
using PowerfulCases

# Load and solve power flow
println("Loading ACTIVSg2000 case...")
case = load("synthetic/ACTIVSg2000")
sys_pf = SystemModel(parse_network(case.raw))

println("Solving power flow...")
pf_prob = ACPowerFlowProblem(sys_pf)
pf_sol = solve(pf_prob, Newton())
println("  PF converged: $(pf_sol.converged)")

# Parse dynamics and create TDS problem
println("Parsing DYR file...")
dyr_data = parse_dynamics(case.dyr, sys_pf)
println("  Models: ", join(["$k: $(length(v))" for (k,v) in pairs(dyr_data)], ", "))

println("Creating TDS problem...")
init = setup_tds(pf_sol, dyr_data)
prob = TDSProblem(init; tspan=(0.0, 0.1))

# Validate initialization
println("\n=== Validating TDS Initialization ===")
result = validate_tds_initialization(prob)
print_validation_report(result)
