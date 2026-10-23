#!/usr/bin/env julia
#=
TDS Two-Bus Demo
================

Simple two-bus SMIB system running in phasor TDS mode.

Run with: ./bin/powerful run examples/tds/tds_two_bus_demo.jl
=#

using Powerful
using PowerfulCases
using PowerfulCore.DiffEq: KLUFactorization, FunctionCallingCallback, QNDF
using PowerfulCore: get_series
using PowerfulTDS: prepare_tds_system
using Printf
using Plots

println("=" ^ 60)
println("TDS Two-Bus Demo (Phasor Mode)")
println("=" ^ 60)

# ============================================================================
# 1. Load system and solve power flow
# ============================================================================
case = load("two_bus_branch")
raw_path = file(case, :raw)
dyr_path = file(case, :dyr)

case = parse_network(raw_path)
sys_pf = SystemModel(case)

prob = ACPowerFlowProblem(sys_pf)
pf_sol = solve(prob, Newton())

println("\n[1] Power Flow: ", pf_sol.converged ? "converged" : "FAILED")

# Show power flow results
println("\n[1b] Power Flow Results:")
for i in 1:length(sys_pf.models.Bus)
    bus_num = sys_pf.models.Bus.i[i]
    v = get_var(pf_sol, :Bus, :v, (bus_num,))
    theta = get_var(pf_sol, :Bus, :theta, (bus_num,))
    println("    Bus $bus_num: V=$(round(v, digits=4)) pu, θ=$(round(rad2deg(theta), digits=2))°")
end

# ============================================================================
# 2. Parse dynamic models and build TDS system
# ============================================================================
dyr_data = parse_dynamics(dyr_path, sys_pf)

# Use standard TDS (phasor) system
sys_tds, links = prepare_tds_system(sys_pf, pf_sol, dyr_data)

println("\n[2] TDS System Models:")
for (name, model) in pairs(sys_tds.models)
    println("    $name: $(length(model))")
end

# ============================================================================
# 3. Initialize TDS system
# ============================================================================
init = TDSInitialization(sys_tds, links, pf_sol)

println("\n[3] Initialization:")
for (name, states) in pairs(init.dyn)
    if !isempty(states)
        n = length(first(values(states)))
        println("    $name: $n instances")
    end
end

# Show GENCLS initial states
if haskey(init.dyn, :GENCLS)
    println("\n[3b] GENCLS Initial States:")
    states = init.dyn.GENCLS
    if haskey(states, :delta)
        println("    delta (rad) = $(states.delta)")
        println("    delta (deg) = $(rad2deg.(states.delta))")
    end
    if haskey(states, :omega)
        println("    omega = $(states.omega)")
    end
end

# ============================================================================
# 4. Create TDSProblem
# ============================================================================
tds_prob = TDSProblem(init; tspan=(0.0, 5.0))

println("\n[4] TDSProblem created")
println("    Simulation span: ", tds_prob.tspan)

# ============================================================================
# 5. Solve TDS
# ============================================================================
println("\n[5] Solving TDS (flat run, 5.0s)...")

step_count = Ref(0)
last_t = Ref(0.0)
# Observer called after every step (it never changes the state)
function progress_callback(u, t, integrator)
    step_count[] += 1
    last_t[] = integrator.t
    if step_count[] <= 10 || step_count[] % 100 == 0
        @printf("    Step %d: t=%.6f, dt=%.6e\n", step_count[], integrator.t, integrator.dt)
    end
    return nothing
end

@time sol = solve(tds_prob, QNDF(linsolve=KLUFactorization());
                  callback=FunctionCallingCallback(progress_callback; func_everystep=true),
                  verbose=true, abstol=1e-8, reltol=1e-6)

println("\n[5b] Solver Statistics:")
println("    Return code: $(sol.retcode)")
println("    Steps taken: $(step_count[])")
println("    Last time reached: $(last_t[])")

# Report result
if sol.retcode == :Success
    println("\n[6] Result: SUCCESS")
else
    println("\n[6] Result: $(sol.retcode) - TDS failed at t=$(last_t[])")
end

# Plot trajectory data
if length(sol.t) > 1
    println("\n[7] Plotting trajectory ($(length(sol.t)) time points)...")

    # Use get_series API for all GENCLS variables
    data_delta = get_series(sol, :GENCLS, :delta)
    data_omega = get_series(sol, :GENCLS, :omega)

    status_str = sol.retcode == :Success ? "Complete" : "Failed at t=$(round(last_t[], digits=4))s"

    plots_list = []

    # Plot rotor angles
    p1 = plot(title="Rotor Angles ($status_str)",
              xlabel="Time (s)", ylabel="Delta (rad)",
              legend=:topright)
    for (var, vals) in pairs(data_delta)
        var == :t && continue
        gen_num = split(string(var), "_")[end]
        plot!(p1, data_delta.t, vals, label="Gen $gen_num", linewidth=1.0)
    end
    push!(plots_list, p1)

    # Plot rotor speeds
    p2 = plot(title="Rotor Speeds ($status_str)",
              xlabel="Time (s)", ylabel="Omega (pu)",
              legend=:topright)
    for (var, vals) in pairs(data_omega)
        var == :t && continue
        gen_num = split(string(var), "_")[end]
        plot!(p2, data_omega.t, vals, label="Gen $gen_num", linewidth=1.0)
    end
    push!(plots_list, p2)

    if !isempty(plots_list)
        combined = plot(plots_list..., layout=(length(plots_list), 1), size=(800, 300*length(plots_list)))
        output_file = joinpath(@__DIR__, "tds_two_bus_results.png")
        savefig(combined, output_file)
        println("    Plot saved: $output_file")
    end
else
    println("\n[7] No trajectory data to plot")
end

println("\nTDS Two-Bus Demo Complete!")
