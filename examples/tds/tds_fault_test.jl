#!/usr/bin/env julia
#= Phasor TDS fault test =#

using Powerful
using PowerfulCases
using PowerfulCore.DiffEq: FBDF, KLUFactorization
using PowerfulCore: get_series
using Plots

println("=" ^ 60)
println("Phasor TDS Fault Demo - Bus 10 Fault")
println("=" ^ 60)

# Load and solve power flow
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

prob = ACPowerFlowProblem(sys)
pf_sol = solve(prob, Newton())
println("\n[1] Power flow: ", pf_sol.converged ? "converged" : "FAILED")

# Parse dynamic models
dyr_data = parse_dynamics(load("ieee14").dyr, sys)

# BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
events = [
    BusFault(0.5, 0.6, 10, "F1"; xf=0.003, rf=0.0),
]

println("\n[2] Fault Events:")
println("    t=0.5s: Apply 3-phase fault at bus 10")
println("    t=0.6s: Clear fault")

# 100% constant impedance loads for numerical stability during fault
config = SimConfig(solver = (TDS = TDSConfig(
    convert_pq_to_zip=true, p_zip=(1.0,0.0,0.0), q_zip=(1.0,0.0,0.0)),))

# Create TDS problem (phasor mode)
init = setup_tds(pf_sol, dyr_data; config=config)
tds_prob = TDSProblem(init; tspan=(0.0, 2.0), events=events)

println("\n[3] TDSProblem: tspan=$(tds_prob.tspan)")

# Solve with PowerFlowInit for DAE reinitialization after events
println("\n[4] Solving Phasor TDS...")

# Track allocations during solve
println("\n=== Starting allocation tracking ===")
GC.gc()
result = @timed solve(tds_prob, FBDF(linsolve=KLUFactorization()); initializealg=PowerFlowInit(; initial_guess=:pre_fault),
                  dtmax=0.01)
sol = result.value

# Report allocation stats
println("\n=== Allocation Report ===")
println("  Total time: $(round(result.time, digits=2)) s")
println("  Total allocations: $(round(result.bytes / 1e9, digits=3)) GB")
println("  Timesteps: $(length(sol.t))")
println("  Bytes per timestep: $(round(result.bytes / length(sol.t) / 1e3, digits=1)) KB")

if sol.retcode == :Success
    println("\n[5] Result: SUCCESS")
    println("    Accepted steps: $(sol.stats.naccept)")
    println("    Rejected steps: $(sol.stats.nreject)")
    println("    Function evals: $(sol.stats.nf)")

    # Collect trajectory data - all GENCLS omega values using get_series
    data = get_series(sol, :GENCLS, :omega)

    # Plot generator speeds
    println("\n[6] Plotting...")

    p = plot(title="Phasor TDS - Generator Speeds (Bus 10 Fault, xf=0.3)",
             xlabel="Time (s)", ylabel="Speed (pu)",
             legend=:topright, size=(900, 500))

    for (var, vals) in pairs(data)
        var == :t && continue  # Skip time column
        var_str = string(var)
        # External key format: GENCLS_omega_bus_id
        parts = split(var_str, "_")
        gen_key = join(parts[3:end], "_")  # e.g., "1_1" for bus 1, id "1"
        plot!(p, data.t, vals, label="Gen $gen_key", linewidth=1.5)
    end
    vline!(p, [0.5, 0.6], linestyle=:dash, color=:red, label="Fault", linewidth=1)

    display(p)
else
    println("\n[5] Result: $(sol.retcode) - TDS failed at t=$(sol.t[end])")
end

println("\n" * "=" ^ 60)
println("Phasor TDS Demo Complete!")
println("=" ^ 60)
