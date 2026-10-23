# # Getting Started with Powerful.jl
#
# Interactive walkthrough of the basic power flow workflow: load a case, solve it,
# inspect the results, change a load and solve again. See the
# [Getting Started guide](../guides/getting_started.md) for detailed explanations.

# ## Setup

using Powerful
using PowerfulCases
using Plots
using DataFrames

# ## 1. Load System and Solve

# Load the IEEE 14-bus test case from PowerfulCases and build the system model
case_bundle = load("ieee14")
case = parse_network(case_bundle.raw)
sys = SystemModel(case)

# Create and solve the power flow problem
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# ## 2. Convergence Check

println("Converged: ", sol.converged)
println("Iterations: ", sol.iterations)
println("Final residual: ", sol.residual_norm)

# ## 3. System Summary

println("System composition:")
for (model, data) in pairs(sys.models)
    println("  $model: $(length(data)) elements")
end

# ## 4. Voltage Profile

# `get_all_vars` returns the values of one variable for all instances of a model,
# in the order of `sys.models.Bus`
bus_numbers = sys.models.Bus.i
voltages = get_all_vars(sol, :Bus, :v)
angles = get_all_vars(sol, :Bus, :theta)

# Plot voltage magnitudes
bar(bus_numbers, voltages,
    xlabel="Bus Number",
    ylabel="Voltage (pu)",
    title="Voltage Profile - IEEE 14 Bus",
    label="",
    ylims=(0.9, 1.1))
hline!([0.95, 1.05], ls=:dash, color=:red, label="0.95 / 1.05 pu")

# ## 5. Voltage Angles

# Convert to degrees for readability
angles_deg = rad2deg.(angles)

bar(bus_numbers, angles_deg,
    xlabel="Bus Number",
    ylabel="Angle (degrees)",
    title="Voltage Angles - IEEE 14 Bus",
    legend=false)

# ## 6. Results Table

# `collect_model_table` gathers the results of one model in physical units
bus_data = collect_model_table(sol, :Bus)
df = DataFrame(bus_data)

# ## 7. Access Individual Variables

# Variables are looked up by external key. Keys are tuples: `(bus_number,)` for buses.
v1 = get_var(sol, :Bus, :v, (1,))
θ1 = get_var(sol, :Bus, :theta, (1,))

println("Bus 1 voltage: $(round(v1, digits=4)) pu")
println("Bus 1 angle: $(round(rad2deg(θ1), digits=2))°")

# ## 8. Modify and Re-solve
#
# Changing system parameters requires creating a new problem from the modified system.
#
# **Note:** Use a load at a bus without a generator to see a voltage change.
# Generator (PV) buses hold their voltage magnitude at the setpoint.

# Increase the load at Bus 4 (a load bus without a generator) by 20%
bus_id = 4
load_idx = findfirst(==(bus_id), sys.models.PQ.i)
original_load = sys.models.PQ.pl[load_idx]
sys.models.PQ.pl[load_idx] *= 1.2
sys.models.PQ.ql[load_idx] *= 1.2

# Create a new problem, warm-started from the previous solution
prob_new = ACPowerFlowProblem(sys; initial_guess=sol.dae.y)
sol_new = solve(prob_new, Newton())

# Check the voltage change at Bus 4
v_base = get_var(sol, :Bus, :v, (bus_id,))
v_new = get_var(sol_new, :Bus, :v, (bus_id,))

println("Original load: $(round(original_load, digits=3)) pu")
println("New load: $(round(sys.models.PQ.pl[load_idx], digits=3)) pu")
println("Converged: ", sol_new.converged)
println("Iterations: ", sol_new.iterations)
println("Bus $bus_id voltage: $(round(v_base, digits=6)) → $(round(v_new, digits=6)) pu")
println("Voltage change: $(round(v_new - v_base, digits=6)) pu ($(round((v_new-v_base)/v_base*100, digits=3))%)")

# ## 9. Compare Results

voltages_new = get_all_vars(sol_new, :Bus, :v)

plot(bus_numbers, voltages,
     label="Base Case", marker=:circle, lw=2)
plot!(bus_numbers, voltages_new,
      label="+20% Load at Bus 4", marker=:square, lw=2)
xlabel!("Bus Number")
ylabel!("Voltage (pu)")
title!("Impact of Load Increase")
hline!([0.95, 1.05], ls=:dash, color=:red, label="")

# ## Next Steps
#
# - **[Power Flow Analysis](../guides/powerflow.md)** - Formulations and algorithms
# - **[Time-Domain Simulation](../guides/tds.md)** - Dynamic analysis
# - **[Accessing Results](../guides/accessing_variables.md)** - Complete API reference
