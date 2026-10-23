# TDS Event Demo: Branch Trip and Reclose
#
# This script demonstrates TDS with branch trip/reclose events.
# The IEEE 14-bus system has 5 GENCLS generators (buses 1, 2, 3, 6, 8).
#
# Run with: ./bin/powerful run examples/tds/tds_event_demo.jl

using Powerful
using PowerfulCases
using PowerfulTDS
using PowerfulCore: get_series
using PowerfulCore.DiffEq: FBDF, KLUFactorization
using CSV
using DataFrames
using Plots


# =============================================================================
# Step 1: Load system and solve power flow
# =============================================================================
println("Loading IEEE 14-bus system...")
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

println("Solving power flow...")
prob_pf = ACPowerFlowProblem(sys)
sol_pf = solve(prob_pf, Newton())

if !sol_pf.converged
    error("Power flow did not converge!")
end
println("  Power flow converged ✓")

# =============================================================================
# Step 2: Parse dynamic models
# =================================================================e===========
println("\nLoading dynamic models...")
dyr_data = parse_dynamics(load("ieee14").dyr, sys)
println("  Loaded models: ", keys(dyr_data))

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
    solver = (TDS = TDSConfig(
        convert_pq_to_zip = true,
        p_zip = (1.0, 0.0, 0.0),  # 100% constant Z for active power
        q_zip = (1.0, 0.0, 0.0)   # 100% constant Z for reactive power
    ),)
)
println("  Load conversion: 100% constant impedance (Z)")

init = setup_tds(sol_pf, dyr_data; config = config)
prob_tds = TDSProblem(init; tspan = (0.0, 5.0), events = events)

# Validate initialization (residuals should be ~0 at t=0)
result = validate_tds_initialization(prob_tds)
print_validation_report(result)

# The algebraic variables are reinitialized after the trip (PowerFlowInit, the default)
println("Solving TDS with FBDF...")
@time sol = solve(prob_tds, FBDF(linsolve=KLUFactorization()); dtmax=0.03);

println("  Simulation completed: $(sol.retcode)")
println("  Time steps: $(length(sol.t))")
println("  Final time: $(sol.t[end])s")

# =============================================================================
# Step 5: Extract and save results
# =============================================================================
println("\nExtracting results...")
data = get_series(sol, :GENCLS, :omega)

# Save to CSV
outfile = joinpath(@__DIR__, "tds_event_demo_results.csv")
CSV.write(outfile, DataFrame(data))
println("  Results saved to: $outfile")

# Plot generator speeds
# Note: External key format is bus_id (e.g., 1_1 for bus 1, id "1")
println("  Generating plot...")
p = plot(data.t, data.GENCLS_omega_1_1, label="Gen (bus 1, slack)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_2_1, label="Gen (bus 2)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_3_1, label="Gen (bus 3)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_6_1, label="Gen (bus 6)", linewidth=1.5)
plot!(data.t, data.GENCLS_omega_8_1, label="Gen (bus 8)", linewidth=1.5)
vline!([1.0], linestyle=:dash, color=:red, label="Trip at t=1.0", linewidth=1)
xlabel!("Time (s)")
ylabel!("Speed (pu)")
title!("Generator Speeds: Line 1-5 Trip at t=1.0")
plotfile = joinpath(@__DIR__, "tds_event_demo_speeds.png")
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
    _, delta_vals = get_series(sol, :GENCLS, :delta, (1, "1"))
    _, omega_vals = get_series(sol, :GENCLS, :omega, (1, "1"))
    delta_range = extrema(delta_vals)
    omega_range = extrema(omega_vals)
    println("\nGenerator 1 state ranges:")
    println("  Delta: ", delta_range, " rad")
    println("  Omega: ", omega_range, " pu")
end

println("""

Results saved to: $outfile
Plot saved to: $plotfile
""")
