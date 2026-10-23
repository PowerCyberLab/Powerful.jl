#!/usr/bin/env julia
#= REGCA1 Fault Test - Mixed synchronous and renewable generators =#

using Powerful
using PowerfulCases
using PowerfulCore.DiffEq: FBDF, KLUFactorization
using PowerfulCore: get_series
using Plots

println("REGCA1 Fault Test - Bus 4 Fault")
println("=" ^ 50)

# Load system
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

# Solve power flow
prob = ACPowerFlowProblem(sys)
pf_sol = solve(prob, Newton())
println("Power flow: ", pf_sol.converged ? "OK" : "FAILED")

# Parse dynamic models with REGCA1
dyr_file = joinpath(dirname(load("ieee14").raw), "ieee14_regca1.dyr")
dyr_data = parse_dynamics(dyr_file, sys)
println("Models: GENROU=$(length(dyr_data.GENROU)), REGCA1=$(length(dyr_data.REGCA1))")

# BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
events = [BusFault(0.5, 0.6, 4, "F1"; xf=0.01, rf=0.0)]

# 100% constant impedance loads
config = SimConfig(solver = (TDS = TDSConfig(
    convert_pq_to_zip=true, p_zip=(1.0,0.0,0.0), q_zip=(1.0,0.0,0.0)),))

# Create and solve TDS
init = setup_tds(pf_sol, dyr_data; config=config)
tds_prob = TDSProblem(init; tspan=(0.0, 3.0), events=events)
println("Solving TDS...")

sol = solve(tds_prob, FBDF(linsolve=KLUFactorization());
            initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.01)

if sol.retcode == :Success
    println("Result: SUCCESS ($(sol.stats.naccept) steps)")

    # Plot generator speeds (GENROU only - REGCA1 has no omega)
    data_gen = get_series(sol, :GENROU, :omega)
    data_ren = get_series(sol, :REGCA1, :Ip)

    p1 = plot(title="Generator Speeds", xlabel="Time (s)", ylabel="ω (pu)", legend=:topright)
    for (var, vals) in pairs(data_gen)
        var == :t && continue
        plot!(p1, data_gen.t, vals, label="Gen $(split(string(var), "_")[3])", linewidth=1.5)
    end
    vline!(p1, [0.5, 0.6], linestyle=:dash, color=:red, label="Fault", linewidth=1)

    p2 = plot(title="REGCA1 Active Current", xlabel="Time (s)", ylabel="Ip (pu)", legend=:topright)
    for (var, vals) in pairs(data_ren)
        var == :t && continue
        plot!(p2, data_ren.t, vals, label="REGCA1 $(split(string(var), "_")[3])", linewidth=1.5)
    end
    vline!(p2, [0.5, 0.6], linestyle=:dash, color=:red, label="Fault", linewidth=1)

    display(plot(p1, p2, layout=(2,1), size=(900, 600)))
else
    println("Result: $(sol.retcode) - failed at t=$(sol.t[end])")
end
