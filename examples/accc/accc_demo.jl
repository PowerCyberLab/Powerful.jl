# AC Contingency Calculation (ACCC) Demo
#
# This example demonstrates how to run AC contingency analysis using Powerful.jl.
# ACCC solves power flow for a set of contingencies (N-1 or N-k outages) starting
# from a converged base case solution.

using Powerful
using PowerfulCases
using PowerfulCA

# Load IEEE 14-bus test case
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

# Define contingencies
# Each contingency has a name and a list of outage actions
contingencies = [
    # N-1 branch outages
    Contingency("LINE_1_5", [BranchOutage(1, 5, "1")]),
    Contingency("LINE_2_3", [BranchOutage(2, 3, "1")]),
    Contingency("LINE_4_5", [BranchOutage(4, 5, "1")]),

    # N-1 generator outage
    Contingency("GEN_3", [GeneratorOutage(3, "1")]),

    # N-2 combined outage (branch + generator)
    Contingency("N2_LINE_GEN", [BranchOutage(1, 5, "1"), GeneratorOutage(3, "1")]),
]

# Run ACCC
# This:
# 1. Saves the base case status
# 2. Solves base case power flow
# 3. For each contingency:
#    - Applies outages (sets status to false)
#    - Solves power flow with warm start from base case
#    - Stores the full solution state vector (y)
#    - Restores base case status
result = solve_accc(sys, contingencies)

# Print violation report
# - vmin/vmax: voltage limits for violation checking (default: 0.95-1.05)
# - rating: which rating set to use for branch limits (:rate_a, :rate_b, :rate_c)
report = print_accc_report(result; vmin=0.95, vmax=1.05, rating=:rate_a)

# The report function returns the violation data for further processing
println("\n--- Programmatic Access ---")
println("Number of voltage violations: ", length(report.voltage_violations))
println("Number of branch overloads: ", length(report.branch_violations))
println("Number of non-converged: ", length(report.nonconverged))

# Access individual contingency results
println("\n--- Individual Results ---")
for r in result.contingencies
    status = r.converged ? "OK" : "FAILED"
    println("$(rpad(r.name, 15)) $(rpad(status, 8)) iters=$(r.iterations)")
end
