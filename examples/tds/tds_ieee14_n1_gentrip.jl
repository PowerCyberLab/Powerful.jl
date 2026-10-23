#!/usr/bin/env julia
"""
Trial: N-1 generator trips on IEEE 14-bus (GENCLS + SEXS + TGOV1).

Trips each GENCLS at t=1.0s via FieldModification(:GENCLS, key, :status).
Checks convergence and plots frequency response to see if instability shows.
"""

ENV["GKSwstype"] = "100"

using Powerful
using PowerfulCore: get_addresses
using PowerfulCore.DiffEq: QNDF, KLUFactorization
using PowerfulTDS: FieldModification
using PowerfulCases
using Plots
using Printf

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
pf_sol = solve(ACPowerFlowProblem(sys), Newton())
pf_sol.converged || error("PF failed")
dyr_data = parse_dynamics(case.dyr, sys)
init = setup_tds(pf_sol, dyr_data)

gens = init.sys.models[:GENCLS]
orig_status = copy(gens.status)
gen_keys = [(Int(gens.BUS[k]), String(gens.ID[k])) for k in eachindex(gens.status)]
println("Generators: ", gen_keys)

tspan = (0.0, 5.0)
t_trip = 1.0
t_plot = range(tspan..., length=501)
gen_buses = Int.(gens.BUS)

w_traces = Dict{Int, Vector{Vector{Float64}}}(b => [] for b in gen_buses)
trace_labels = String[]

for (bus, id) in gen_keys
    events = [FieldModification(t_trip, :GENCLS, (Int32(bus), id), :status, :set, 0.0)]
    tds_prob = TDSProblem(init; tspan=tspan, events=events)
    t_run = @elapsed sol = solve(tds_prob, QNDF(linsolve=KLUFactorization());
        initializealg = PowerFlowInit(; initial_guess = :pre_fault),
        dtmax = 0.05, verbose = false,
    )
    gens.status .= orig_status

    if sol.retcode != :Success
        @printf("trip G%-3d  FAILED  (%s, reached t=%.3f)\n", bus, sol.retcode, sol.t[end])
        continue
    end

    am = sol.address_manager
    w_addrs = get_addresses(am, :GENCLS, :omega)
    for (g, b) in enumerate(gen_buses)
        push!(w_traces[b], [sol.sol(t)[w_addrs[g]] for t in t_plot])
    end
    push!(trace_labels, "trip G$bus")
    wmin, wmax = extrema(vcat([tr for (b, trs) in w_traces for tr in trs[end:end]]...))
    @printf("trip G%-3d  %6.3f s  (%d steps)  omega range [%.4f, %.4f]\n",
            bus, t_run, length(sol.t), wmin, wmax)
end

wplots = [begin
    p = plot(title="Gen at Bus $b", titlefontsize=10, ylabel="ω (pu)",
             xlabel="Time (s)", legend=(b == gen_buses[1]), legendfontsize=6)
    for (s, tr) in enumerate(w_traces[b])
        plot!(p, t_plot, tr, linewidth=1.2, label=trace_labels[s])
    end
    vline!(p, [t_trip], linestyle=:dash, color=:black, linewidth=0.8, label="")
    p
end for b in gen_buses]
pw = plot(wplots..., layout=(2, 3), size=(1200, 660),
          left_margin=6Plots.mm, bottom_margin=5Plots.mm,
          plot_title="Generator Speeds: N-1 generator trips")
savefig(pw, joinpath(@__DIR__, "ieee14_n1_gentrip_speeds.png"))
println("Saved ieee14_n1_gentrip_speeds.png")
