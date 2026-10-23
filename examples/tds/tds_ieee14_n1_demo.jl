#!/usr/bin/env julia
"""
TDS Demo: Brute-Force N-1 Line Trip Screening on IEEE 14-Bus
=============================================================

Runs a full time-domain simulation for EVERY line outage (N-1) and
overlays bus voltages and generator speeds across all scenarios.
Designed to showcase fast repeated transient simulation.

Dynamic models: GENCLS + SEXS + TGOV1 (from ieee14.dyr).
Each scenario: line trips at t=1.0s, simulate to t=5.0s.

# Usage
    ./bin/powerful run examples/tds/tds_ieee14_n1_demo.jl
"""

# Pop GUI windows when a display is available, else fall back to headless
if get(ENV, "DISPLAY", "") == "" && get(ENV, "WAYLAND_DISPLAY", "") == ""
    ENV["GKSwstype"] = "100"  # headless plotting
else
    delete!(ENV, "GKSwstype")  # use GR's Qt window (gksqt)
end

using Powerful
using PowerfulCore: get_addresses
using PowerfulCore.DiffEq: QNDF, KLUFactorization
using PowerfulCases
using Plots
using Printf

# =============================================================================
# Setup (once)
# =============================================================================
println("=" ^ 64)
println("IEEE 14-Bus Brute-Force N-1 Transient Screening")
println("=" ^ 64)

t_load = @elapsed begin
    case = load("ieee14")
    sys = SystemModel(parse_network(case.raw))
end
println("\n1. Network loaded in $(round(t_load, digits=2)) s")

t_pf = @elapsed pf_sol = solve(ACPowerFlowProblem(sys), Newton())
pf_sol.converged || error("Base case power flow failed")
println("2. Base power flow converged in $(round(t_pf, digits=3)) s")

t_setup = @elapsed begin
    dyr_data = parse_dynamics(case.dyr, sys)
    init = setup_tds(pf_sol, dyr_data)
end
println("3. TDS setup in $(round(t_setup, digits=2)) s")

# =============================================================================
# Enumerate N-1 line-trip scenarios
# =============================================================================
branches = init.sys.models.Branch
orig_status = copy(branches.status)
scenarios = [(Int(branches.i[k]), Int(branches.j[k]), String(branches.ckt[k]))
             for k in eachindex(branches.status) if branches.status[k]]

tspan = (0.0, 5.0)
t_trip = 1.0
println("\nScenarios: $(length(scenarios)) line trips (t=$(t_trip)s), " *
        "each simulated for $(tspan[2]) s")

# Buses/generators to record
plot_buses = [1, 3, 4, 8, 9, 14]
gen_buses = Int.(init.sys.models[:GENCLS].BUS)

# =============================================================================
# Run all scenarios
# =============================================================================
t_plot = range(tspan..., length=501)
v_traces = Dict(b => Vector{Vector{Float64}}() for b in plot_buses)
w_traces = Dict(b => Vector{Vector{Float64}}() for b in gen_buses)
labels = String[]
times = Float64[]
failed = String[]

for (i, j, ckt) in scenarios
    name = "$(i)-$(j)"
    events = [BranchTrip(t_trip, i, j, ckt)]
    tds_prob = TDSProblem(init; tspan=tspan, events=events)

    t_run = @elapsed sol = solve(tds_prob, QNDF(linsolve=KLUFactorization());
        initializealg = PowerFlowInit(; initial_guess = :pre_fault),
        dtmax = 0.05,
        verbose = false,
    )
    branches.status .= orig_status  # restore tripped line for next scenario

    if sol.retcode != :Success
        push!(failed, name)
        @printf("   trip %-7s  FAILED (%s)\n", name, sol.retcode)
        continue
    end

    # Record trajectories on common time grid
    am = sol.address_manager
    n_diff = sol.n_diff
    v_addrs = get_addresses(am, :Bus, :v)
    w_addrs = get_addresses(am, :GENCLS, :omega)
    bus_idx = init.sys.order.key_indices[:Bus]

    for b in plot_buses
        addr = v_addrs[bus_idx[(b,)]] + n_diff
        push!(v_traces[b], [sol.sol(t)[addr] for t in t_plot])
    end
    for (g, b) in enumerate(gen_buses)
        push!(w_traces[b], [sol.sol(t)[w_addrs[g]] for t in t_plot])
    end
    push!(labels, name)
    push!(times, t_run)
    @printf("   trip %-7s  %6.3f s  (%d steps)\n", name, t_run, length(sol.t))
end

# =============================================================================
# Performance summary
# =============================================================================
println("\n" * "=" ^ 64)
println("PERFORMANCE SUMMARY")
println("=" ^ 64)
n_ok = length(times)
@printf("Scenarios solved:        %d / %d\n", n_ok, length(scenarios))
@printf("Total solve time:        %.2f s (incl. JIT on first solve)\n", sum(times))
if n_ok > 1
    warm = times[2:end]
    @printf("Warm avg per scenario:   %.3f s  (%.1f s simulated -> %.0fx real time)\n",
            sum(warm)/length(warm), tspan[2], tspan[2] / (sum(warm)/length(warm)))
end
isempty(failed) || println("Failed scenarios: ", join(failed, ", "))

# =============================================================================
# Plots: overlay all scenarios
# =============================================================================
println("\nGenerating plots...")

const PLOT_SIZE = (1200, 680)
const WINDOW_POS = "+30+80"  # window position on screen (X11 geometry offsets)

"""Position the live window: set Qt's -geometry flag on the command GR uses to
spawn its gksqt viewer, replacing any stale flag left in the daemon env by other
demos. A running viewer would be reused at its old position/size, so close it
(GR launches gksqt through a `sh` wrapper, so match the full command line)."""
function position_window!()
    haskey(ENV, "GKS_QT") || return  # GR only sets it when a display exists
    geometry = " -geometry $(PLOT_SIZE[1])x$(PLOT_SIZE[2])$(WINDOW_POS)"
    ENV["GKS_QT"] = replace(ENV["GKS_QT"], r" -geometry \S+" => "") * geometry
    success(`pkill -f bin/gksqt`)  # ok if none running
    return nothing
end

warm_avg = n_ok > 1 ? sum(times[2:end]) / (n_ok - 1) : times[1]
timing_str = @sprintf("%d x %.0f s scenarios solved in %.2f s wall (%.0f ms each, %.0fx real time) — QNDF + KLU",
                      n_ok, tspan[2], sum(times), 1000 * warm_avg, tspan[2] / warm_avg)

vplots = [begin
    p = plot(title="Bus $b", titlefontsize=10, ylabel="V (pu)",
             xlabel="Time (s)", legend=false)
    for tr in v_traces[b]
        plot!(p, t_plot, tr, linewidth=1.0, alpha=0.7)
    end
    vline!(p, [t_trip], linestyle=:dash, color=:black, linewidth=0.8)
    p
end for b in plot_buses]
pv = plot(vplots..., layout=(2, 3), size=PLOT_SIZE,
          left_margin=6Plots.mm, bottom_margin=5Plots.mm,
          plot_title="Bus Voltages: all $(n_ok) N-1 line trips overlaid\n$(timing_str)",
          plot_titlefontsize=13, plot_titlevspan=0.11)
savefig(pv, joinpath(@__DIR__, "ieee14_n1_voltages.png"))
println("   Saved ieee14_n1_voltages.png")

wplots = [begin
    p = plot(title="Gen at Bus $b", titlefontsize=10, ylabel="ω (pu)",
             xlabel="Time (s)", legend=false)
    for tr in w_traces[b]
        plot!(p, t_plot, tr, linewidth=1.0, alpha=0.7)
    end
    vline!(p, [t_trip], linestyle=:dash, color=:black, linewidth=0.8)
    p
end for b in gen_buses]
pw = plot(wplots..., layout=(2, 3), size=PLOT_SIZE,
          left_margin=6Plots.mm, bottom_margin=5Plots.mm,
          plot_title="Generator Speeds: all $(n_ok) N-1 line trips overlaid\n$(timing_str)",
          plot_titlefontsize=13, plot_titlevspan=0.11)
savefig(pw, joinpath(@__DIR__, "ieee14_n1_speeds.png"))
println("   Saved ieee14_n1_speeds.png")

if !haskey(ENV, "GKSwstype")  # display available: pop the plot window
    position_window!()
    gui(pv)
    println("   Displayed voltage plot in GR window")
end

println("\nDone.")
