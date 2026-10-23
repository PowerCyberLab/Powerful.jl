#!/usr/bin/env julia
"""
TDS fault benchmark for ACTIVSg70k (singlegen variant).
Compares QNDF (default) vs Trapezoid method performance.

Expect roughly 20-25 minutes and 7 GiB of peak memory: two warmup solves and two
benchmarked solves (QNDF, then Trapezoid) on a 70k-bus system, where Trapezoid's
per-step cost runs over ten times QNDF's.
"""

using Powerful
using PowerfulCore: get_series
using PowerfulCore.DiffEq: Trapezoid, QNDF, KLUFactorization
using PowerfulCases
using Plots

println("="^60)
println("TDS Fault Benchmark: ACTIVSg70k (singlegen variant)")
println("="^60)

# Load case using collection/case format to avoid ambiguity
case = load("synthetic/ACTIVSg70k")
raw_path = file(case, :raw, variant="singlegen")
dyr_path = file(case, :dyr, variant="singlegen")
sys = SystemModel(parse_network(raw_path))

pf_sol = solve(ACPowerFlowProblem(sys), Newton(), verbose = true)
println("Power flow: ", pf_sol.converged ? "converged" : "FAILED")

# Parse dynamics
dyr_data = parse_dynamics(dyr_path, sys)

# BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
events = [BusFault(0.1, 0.2, 845, "F1"; xf=0.0001, rf=0.0)]

# 100% constant impedance loads
config = SimConfig(
    solver = (
        TDS = TDSConfig(
            convert_pq_to_zip = true,
            p_zip = (1.0, 0.0, 0.0),
            q_zip = (1.0, 0.0, 0.0),
        ),
    ),
)

# Create TDS problem
init = setup_tds(pf_sol, dyr_data; config=config)
tds_prob = TDSProblem(init; tspan=(0.0, 0.5), events=events)

# Common solve options
solve_opts = (
    initializealg = PowerFlowInit(; initial_guess = :pre_fault),
    dtmax = 0.01,
    verbose = false,
)

println("\n" * "="^60)
println("Benchmarking QNDF (default) vs Trapezoid")
println("="^60)

# Warmup runs
println("\nWarmup runs...")
_ = solve(tds_prob, QNDF(linsolve=KLUFactorization()); solve_opts...)
_ = solve(tds_prob, Trapezoid(linsolve=KLUFactorization()); solve_opts...)
println("Warmup complete.")

# Benchmark QNDF
println("\nBenchmarking QNDF...")
t_qndf = @elapsed sol_qndf = solve(tds_prob, QNDF(linsolve=KLUFactorization()); solve_opts...)
steps_qndf = sol_qndf.stats.naccept
println("  QNDF time: $(round(t_qndf * 1000, digits=2)) ms")
println("  Steps: $steps_qndf")
println("  Time/step: $(round(t_qndf * 1000 / steps_qndf, digits=3)) ms")
println("  Status: $(sol_qndf.retcode)")

# Benchmark Trapezoid
println("\nBenchmarking Trapezoid...")
t_trap = @elapsed sol_trap = solve(tds_prob, Trapezoid(linsolve=KLUFactorization()); solve_opts...)
steps_trap = sol_trap.stats.naccept
println("  Trapezoid time: $(round(t_trap * 1000, digits=2)) ms")
println("  Steps: $steps_trap")
println("  Time/step: $(round(t_trap * 1000 / steps_trap, digits=3)) ms")
println("  Status: $(sol_trap.retcode)")

# Summary
println("\n" * "="^60)
println("Summary")
println("="^60)
println("  QNDF:      $(round(t_qndf * 1000, digits=2)) ms, $steps_qndf steps, $(round(t_qndf * 1000 / steps_qndf, digits=3)) ms/step")
println("  Trapezoid: $(round(t_trap * 1000, digits=2)) ms, $steps_trap steps, $(round(t_trap * 1000 / steps_trap, digits=3)) ms/step")
println("  Speedup:   $(round(t_qndf / t_trap, digits=2))x (Trapezoid vs QNDF)")

# Use QNDF solution for plotting
sol1 = sol_qndf
println("\nTDS solution status: ", sol1.retcode)

# Plot fault bus voltage
if sol1.retcode == :Success
    println("\nPlotting fault bus voltage...")

    fault_bus  = 845
    plot_buses = [844, 845, 846]

    # Create voltage plot
    p = plot(
        title = "Fault Bus Voltage (Bus $fault_bus) – QNDF: $(round(t_qndf, digits=2)) s, $steps_qndf steps",
        xlabel = "Time (s)",
        ylabel = "Voltage (pu)",
        legend = :bottomright,
        size = (800, 400),
    )

    for bnum in plot_buses
        t_b, v_b = get_series(sol1, :Bus, :v, bnum)
        lw    = bnum == fault_bus ? 2.5 : 1.5
        label = bnum == fault_bus ? "Bus $bnum (fault)" : "Bus $bnum"
        plot!(p, t_b, v_b; label, linewidth = lw)
    end

    # Add fault timing markers
    vline!(p, [0.1, 0.2], linestyle = :dash, color = :red, label = "Fault on/off", linewidth = 1)

    display(p)
    savefig(p, joinpath(@__DIR__, "activsg70k_fault_voltage.png"))
    println("Plot saved to activsg70k_fault_voltage.png")
else
    println("\nTDS FAILED: $(sol1.retcode)")
end
