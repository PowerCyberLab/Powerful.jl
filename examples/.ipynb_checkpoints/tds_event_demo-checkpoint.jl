# TDS Event Demo: Branch Trip and Reclose
#
# This script demonstrates TDS with branch trip/reclose events.
# The IEEE 14-bus system has 5 GENCLS generators (buses 1, 2, 3, 6, 8).
#
# Run with: julia --project examples/tds_event_demo.jl

using Powerful
using PowerFlowData
using DifferentialEquations
using CSV
using DataFrames
using Plots

# =============================================================================
# Step 1: Load system and solve power flow
# =============================================================================
println("Loading IEEE 14-bus system...")
case = PowerFlowData.parse_network(joinpath(pkgdir(Powerful), "cases", "ieee14.raw"))
sys = SystemModel(case)

println("Solving power flow...")
prob_pf = ACPowerFlowProblem(sys)
sol_pf = solve(prob_pf, NewtonRaphson())

if !sol_pf.converged
    error("Power flow did not converge!")
end
println("  Power flow converged ✓")

# =============================================================================
# Step 2: Parse dynamic models
# =============================================================================
println("\nLoading dynamic models...")
dyr_data = parse_dyr(joinpath(pkgdir(Powerful), "cases", "ieee14.dyr"), sys)
println("  GENCLS generators: ", length(dyr_data.GENCLS))
println("  TGOV1 governors: ", length(dyr_data.TGOV1))

# =============================================================================
# Step 3: Define events - trip line 1-5 at t=1.0 (from slack bus - bigger impact)
# =============================================================================
events = [
    BranchTrip(1.0, 1, 5, "1"),
]

println("\nEvents defined:")
for e in events
    action = e.value == 0.0 ? "TRIP" : "CLOSE"
    println("  t=$(e.time)s: $action Branch $(e.key)")
end

# =============================================================================
# Step 4: Create and solve TDS problem
# =============================================================================
println("\nCreating TDS problem...")

# Convert all PQ loads to 100% constant impedance (Z) for TDS
# ZIP tuple format: (Z%, I%, P%)
config = SimConfig(
    solver = Dict{Symbol, Any}(:TDS => TDSConfig(
        convert_pq_to_zip = true,
        p_zip = (1.0, 0.0, 0.0),  # 100% constant Z for active power
        q_zip = (1.0, 0.0, 0.0)   # 100% constant Z for reactive power
    ))
)
println("  Load conversion: 100% constant impedance (Z)")

prob_tds = TDSProblem(
    sys, sol_pf,
    (GENCLS = dyr_data.GENCLS, TGOV1 = dyr_data.TGOV1);
    tspan = (0.0, 5.0),
    events = events,
    config = config
)

println("Solving TDS with FBDF (NoInit for event handling)...")
sol = solve(prob_tds; alg=FBDF(), initializealg=NoInit())

println("  Simulation completed: $(sol.retcode)")
println("  Time steps: $(length(sol.t))")
println("  Final time: $(sol.t[end])s")

# =============================================================================
# Step 5: Extract and save results
# =============================================================================
println("\nExtracting results...")
data = collect_state_trajectory(sol, prob_tds)

# Save to CSV
outfile = joinpath(pkgdir(Powerful), "examples", "tds_event_demo_results.csv")
CSV.write(outfile, DataFrame(data))
println("  Results saved to: $outfile")

# Plot generator speeds
# Note: GENCLS instances are ordered by dyr file, NOT by bus number
# Instance 1→Bus 2, Instance 2→Bus 3, Instance 3→Bus 6, Instance 4→Bus 8, Instance 5→Bus 1
println("  Generating plot...")
p = plot(data.t, data.GENCLS_omega_5, label="Gen 5 (bus 1, slack)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_1, label="Gen 1 (bus 2)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_2, label="Gen 2 (bus 3)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_3, label="Gen 3 (bus 6)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_4, label="Gen 4 (bus 8)", linewidth=1.5)
vline!([1.0], linestyle=:dash, color=:red, label="Trip at t=1.0", linewidth=1)
xlabel!("Time (s)")
ylabel!("Speed (pu)")
title!("Generator Speeds: Line 1-5 Trip at t=1.0")
plotfile = joinpath(pkgdir(Powerful), "examples", "tds_event_demo_speeds.png")
savefig(p, plotfile)
println("  Plot saved to: $plotfile")

# =============================================================================
# Step 6: Summary
# =============================================================================
println("\n" * "="^60)
println("DEMO COMPLETE")
println("="^60)

# Show delta and omega ranges
if length(data.t) > 0
    delta_range = extrema(data.GENCLS_delta_1)
    omega_range = extrema(data.GENCLS_omega_1)
    println("\nGenerator 1 state ranges:")
    println("  Delta: ", delta_range, " rad")
    println("  Omega: ", omega_range, " pu")
end

println("""

Results saved to: examples/tds_event_demo_results.csv
Plot saved to: examples/tds_event_demo_speeds.png
""")
