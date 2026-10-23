#!/usr/bin/env julia
"""
TDS fault benchmark for ACTIVSg2000 (singlegen variant).
"""

using Powerful
using PowerfulCases
using PowerfulCore: get_series
using Plots
using Printf

println("="^60)
println("TDS Fault Benchmark: ACTIVSg2000 (singlegen variant)")
println("="^60)

# Load case using collection/case format to avoid ambiguity
case = load("synthetic/ACTIVSg2000")
raw_path = file(case, :raw, variant="singlegen")
dyr_path = file(case, :dyr, variant="singlegen")
sys = SystemModel(parse_network(raw_path))

pf_sol = solve(ACPowerFlowProblem(sys), Newton(), verbose = true)
println("Power flow: ", pf_sol.converged ? "converged" : "FAILED")

# Parse dynamics
dyr_data = parse_dynamics(dyr_path, sys)

# BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
events = [BusFault(0.1, 0.2, 5001, "F1"; xf=0.0001, rf=0.0)]

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

timed_result = @timed solve(
    tds_prob;
    initializealg = PowerFlowInit(; initial_guess = :pre_fault),
    dtmax = 0.01,
    verbose = true,
)
sol1 = timed_result.value
solve_time = timed_result.time

println("TDS solution status: ", sol1.retcode)
println(@sprintf("Solve time: %.3f s | Timesteps: %d", solve_time, length(sol1.t)))

# ── Voltage plot ──────────────────────────────────────────────────────────────
const FAULT_BUS   = 5001
const FAULT_START = 0.1
const FAULT_END   = 0.2
const N_NEARBY    = 4   # most-affected buses to add alongside fault bus

# Get all bus voltages
data_v = get_series(sol1, :Bus, :v)

# Find N buses with the biggest voltage dip (excluding fault bus)
bus_vars = filter(x -> x != :t, collect(propertynames(data_v)))
bus_dips = [(parse(Int, join(split(string(v), "_")[3:end], "_")),
             minimum(getproperty(data_v, v)))
            for v in bus_vars]
sort!(bus_dips, by = x -> x[2])   # ascending: most-affected first

plot_buses = [FAULT_BUS]
for (bnum, _) in bus_dips
    bnum == FAULT_BUS && continue
    push!(plot_buses, bnum)
    length(plot_buses) >= 1 + N_NEARBY && break
end

p = plot(
    title  = @sprintf("ACTIVSg2000 Bus Voltages – Fault at Bus %d (t=%.1f–%.1fs)\nSolve time: %.3f s  |  Timesteps: %d",
                      FAULT_BUS, FAULT_START, FAULT_END, solve_time, length(sol1.t)),
    xlabel = "Time (s)",
    ylabel = "Voltage Magnitude (pu)",
    legend = :bottomright,
    size   = (900, 500),
)

for bnum in plot_buses
    t_bus, v_bus = get_series(sol1, :Bus, :v, bnum)
    lw    = bnum == FAULT_BUS ? 2.5 : 1.5
    label = bnum == FAULT_BUS ? "Bus $bnum (fault)" : "Bus $bnum"
    plot!(p, t_bus, v_bus; label, linewidth = lw)
end

vline!(p, [FAULT_START, FAULT_END];
       linestyle = :dash, color = :red, label = "Fault window", linewidth = 1.5)

output_file = joinpath(@__DIR__, "activsg2000_fault_voltages.png")
savefig(p, output_file)
println("Plot saved: ", abspath(output_file))
display(p)