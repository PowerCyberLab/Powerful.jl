#!/usr/bin/env julia
#=
Y-Bus Two-Bus Demo
==================

Minimal test of Y-bus TDS formulation with a two-bus SMIB system.

This demonstrates the partitioned simulation method where:
1. Network equations are solved on the sparse Y-bus (Newton on the bus voltages)
2. Differential states are integrated using trapezoidal rule
3. No DiffEq.jl - pure partitioned approach

Run with: ./bin/powerful run examples/ybus/ybus_two_bus_demo.jl
=#

using Powerful
using PowerfulCases
using PowerfulYbus
using PowerfulYbus: YBusTDSProblem, solve_network!, solve_trapezoidal
using PowerfulCore: DAE
using PowerfulTDS: init_tds_buffers!, prepare_tds_system
using LinearAlgebra: norm
using SparseArrays: nnz
using Printf

println("=" ^ 60)
println("Y-Bus Two-Bus Demo (Partitioned Simulation)")
println("=" ^ 60)

# ============================================================================
# 1. Load system and solve power flow
# ============================================================================
println("\n[1] Loading two_bus_branch case...")

raw_file = file(load("two_bus_branch"), :raw)
dyr_file = file(load("two_bus_branch"), :dyr)

println("    RAW: $raw_file")
println("    DYR: $dyr_file")

case = parse_network(raw_file)
sys = SystemModel(case)

pf_prob = ACPowerFlowProblem(sys)
pf_sol = solve(pf_prob, Newton())

println("\n[2] Power Flow: ", pf_sol.converged ? "CONVERGED" : "FAILED")

# ============================================================================
# 2. Parse dynamic models and create TDS system
# ============================================================================
dyr_data = parse_dynamics(dyr_file, sys)
sys_tds, links = prepare_tds_system(sys, pf_sol, dyr_data)

# ============================================================================
# 3. Initialize and create TDSProblem
# ============================================================================
println("\n[3] Initializing TDS...")

init = TDSInitialization(sys_tds, links, pf_sol)

tspan = (0.0, 1.0)  # 1 second flat run
tds_prob = TDSProblem(init; tspan=tspan)

# ============================================================================
# 4. Wrap with Y-Bus formulation
# ============================================================================
println("\n[4] Creating YBusTDSProblem...")

ybus_prob = YBusTDSProblem(tds_prob)

# ============================================================================
# 5. Verify initial residuals - KEY DIAGNOSTIC
# ============================================================================
println("\n[5] Checking initial state and residuals...")

am = tds_prob.address_manager
dae = DAE(am; allocate_differential=true)
init_tds_buffers!(dae.x, dae.y, tds_prob)

# Show initial states
println("\n    Initial differential states (x):")
println("      delta_1 = ", @sprintf("%.6f rad (%.2f deg)", dae.x[1], rad2deg(dae.x[1])))
println("      delta_2 = ", @sprintf("%.6f rad (%.2f deg)", dae.x[2], rad2deg(dae.x[2])))
println("      omega_1 = ", @sprintf("%.6f pu", dae.x[3]))
println("      omega_2 = ", @sprintf("%.6f pu", dae.x[4]))

# Show P_m from initialization
if haskey(init.dyn, :GENCLS)
    println("\n    Initial P_m (from init.dyn.GENCLS):")
    println("      P_m_1 = ", @sprintf("%.6f pu", init.dyn.GENCLS.P_m[1]))
    println("      P_m_2 = ", @sprintf("%.6f pu", init.dyn.GENCLS.P_m[2]))
end

# Evaluate at t=0 and trace intermediate values
println("\n    --- DEBUG: Before solve_network! ---")

# Access network cache
cache = ybus_prob.network_cache
E_buf = cache.E_buf

println("    E_internal (before gather): ", E_buf.E_internal)
println("    Z_gen: ", E_buf.Z_gen)

# Manually check what gather_internal_voltages! will produce
using PowerfulYbus: gather_internal_voltages!
gather_internal_voltages!(E_buf, ybus_prob.sys_ws, dae)
println("\n    E_internal (after gather):")
for (i, E) in enumerate(E_buf.E_internal)
    println("      Gen $i: |E| = $(round(abs(E), digits=4)), angle = $(round(rad2deg(angle(E)), digits=2))°")
end

# Now solve the network at these internal voltages
solve_network!(cache, ybus_prob.base_prob.sys)

# Check V_bus from network solution (should match power flow at t=0)
println("\n    V_bus (from solve_network!):")
for (i, V) in enumerate(cache.V_bus)
    println("      Bus $i: |V| = $(round(abs(V), digits=4)), angle = $(round(rad2deg(angle(V)), digits=2))°")
end

# Compare with power flow V_bus
println("\n    V_bus (from power flow - should match at t=0):")
for i in 1:length(init.pf.v)
    v = init.pf.v[i]
    theta = init.pf.theta[i]
    println("      Bus $i: |V| = $(round(v, digits=4)), angle = $(round(rad2deg(theta), digits=2))°")
end

# Show derivatives
println("\n    Initial derivatives (f = dx/dt):")
println("      d(delta_1)/dt = ", @sprintf("%.6f", dae.f[1]))
println("      d(delta_2)/dt = ", @sprintf("%.6f", dae.f[2]))
println("      d(omega_1)/dt = ", @sprintf("%.6f", dae.f[3]), "  <- swing eq: (P_m - P_e - D*(omega-1)) / (2H)")
println("      d(omega_2)/dt = ", @sprintf("%.6f", dae.f[4]), "  <- swing eq: (P_m - P_e - D*(omega-1)) / (2H)")

# Show P_e from network solution
println("\n    P_e from network solution:")
println("      P_e_1 = ", @sprintf("%.6f pu", cache.P_e[1]))
println("      P_e_2 = ", @sprintf("%.6f pu", cache.P_e[2]))

# Show mismatch
if haskey(init.dyn, :GENCLS)
    println("\n    Power mismatch (P_m - P_e):")
    println("      Gen 1: ", @sprintf("%.6f pu", init.dyn.GENCLS.P_m[1] - cache.P_e[1]))
    println("      Gen 2: ", @sprintf("%.6f pu", init.dyn.GENCLS.P_m[2] - cache.P_e[2]))
end

f_norm = norm(dae.f)
g_norm = norm(dae.g)

println("\n    Residual norms:")
# Note: 1e-2 tolerance is reasonable for TDS (small mismatch between PF and TDS formulations)
println("      ||f|| = ", @sprintf("%.2e", f_norm), f_norm < 1e-2 ? " (OK)" : " (PROBLEM)")
println("      ||g|| = ", @sprintf("%.2e", g_norm), g_norm < 1e-2 ? " (OK)" : " (PROBLEM)")

# ============================================================================
# 6. Solve with Trapezoidal (partitioned approach)
# ============================================================================
println("\n[6] Solving with Trapezoidal...")

sol = solve_trapezoidal(ybus_prob; dt=0.005, abstol=1e-6, verbose=false)

println("    Status: ", sol.retcode)
println("    Steps: ", size(sol.u, 2))

# ============================================================================
# Summary (Flat Run)
# ============================================================================
println("\n" * "=" ^ 60)
if f_norm < 1e-2 && g_norm < 1e-2
    println("SUCCESS: Initial residuals are small (||f|| < 1e-2)")
    println("  P_e matches P_m within acceptable TDS tolerance")
    println("  Simulation should be stable")
else
    println("PROBLEM: Initial residuals are large!")
    println("  The P_e computed by solve_network! does not match P_m from initialization.")
    println("  Check: E_internal calculation, Y matrix, P_e formula")
end
println("=" ^ 60)

# ============================================================================
# 7. BusFault Demo (Simple fault event handling)
# ============================================================================
println("\n\n")
println("=" ^ 60)
println("BusFault Demo (Simple Fault Event Handling)")
println("=" ^ 60)

println("\n[7] Creating BusFault event...")

# BusFault: Simple bus fault that auto-creates the Fault model
# No manual create_fault/add_fault! needed!
#
# 2-topology pattern:
#   t < 1.0:  prefault (normal network)
#   1.0 ≤ t < 1.1: fault_on (shunt at bus 1)
#   t ≥ 1.1:  post_fault (shunt removed, network restored)

# BusFault auto-creates Fault model and expands to FaultApply + FaultClear
events = [BusFault(1.0, 1.1, 1, "F1"; xf=1e-4)]

println("    BusFault: bus=1, id=\"F1\"")
println("    time_on=1.0s, time_off=1.1s")
println("    xf=1e-4 (near-bolted fault)")

# Create new YBusTDSProblem with events (topologies pre-computed here)
println("\n[8] Creating YBusTDSProblem with events...")

tspan_fault = (0.0, 2.0)  # 2 second run with fault
tds_prob_fault = TDSProblem(init; tspan=tspan_fault)
ybus_prob_fault = YBusTDSProblem(tds_prob_fault; events=events)

# Show pre-computed topologies
println("\n    Pre-computed topologies:")
topo_cache = ybus_prob_fault.network_cache.topology_cache
for (key, Y) in topo_cache.topologies
    println("      $(key.name): Y size = $(size(Y.Y)), nnz = $(nnz(Y.Y))")
end

# Show event schedule
if !isnothing(ybus_prob_fault.event_schedule)
    sched = ybus_prob_fault.event_schedule
    println("\n    Event schedule:")
    println("      Prefault key: $(sched.prefault_key.name)")
    for (i, (t, key)) in enumerate(zip(sched.times, sched.topology_keys))
        println("      t=$t → $(key.name)")
    end
end

# Solve with events
println("\n[9] Solving with Trapezoidal (fault simulation)...")

sol_fault = solve_trapezoidal(ybus_prob_fault; dt=0.001, abstol=1e-6, verbose=true)

println("\n    Status: ", sol_fault.retcode)
println("    Steps: ", size(sol_fault.u, 2))
println("    Time range: ", sol_fault.t[1], " to ", sol_fault.t[end])

# Show key states at different times
println("\n    State evolution:")
n_points = size(sol_fault.u, 2)
indices_to_show = [1, findfirst(t -> t >= 1.0, sol_fault.t),
                   findfirst(t -> t >= 1.05, sol_fault.t),
                   findfirst(t -> t >= 1.1, sol_fault.t),
                   n_points]
filter!(x -> !isnothing(x), indices_to_show)

for idx in indices_to_show
    t = sol_fault.t[idx]
    delta_1 = sol_fault.u[1, idx]
    delta_2 = sol_fault.u[2, idx]
    omega_1 = sol_fault.u[3, idx]
    omega_2 = sol_fault.u[4, idx]

    phase = if t < 1.0
        "prefault"
    elseif t < 1.1
        "fault_on"
    else
        "post_fault"
    end

    println(@sprintf("      t=%.2f (%s): δ₁=%.3f° δ₂=%.3f° ω₁=%.4f ω₂=%.4f",
            t, phase, rad2deg(delta_1), rad2deg(delta_2), omega_1, omega_2))
end

println("\n" * "=" ^ 60)
println("BusFault Demo Complete")
println("=" ^ 60)

# ============================================================================
# 10. Plot Generator Speed Curves
# ============================================================================
println("\n[10] Plotting generator speed curves...")

using Plots

# Extract time and speed data
t = sol_fault.t
omega_1 = sol_fault.u[3, :]  # Generator 1 speed
omega_2 = sol_fault.u[4, :]  # Generator 2 speed

# Create plot
p = plot(t, omega_1, label="Gen 1 (ω₁)", linewidth=2, color=:blue)
plot!(p, t, omega_2, label="Gen 2 (ω₂)", linewidth=2, color=:red)

# Add event markers
vline!(p, [1.0], label="Fault On (t=1.0s)", linestyle=:dash, color=:gray, linewidth=1)
vline!(p, [1.1], label="Fault Clear (t=1.1s)", linestyle=:dot, color=:gray, linewidth=1)

# Labels and formatting
xlabel!(p, "Time (s)")
ylabel!(p, "Speed (pu)")
title!(p, "Generator Speed Response - BusFault at Bus 1")
plot!(p, legend=:bottomleft, grid=true)

# Save figure
output_path = joinpath(dirname(@__FILE__), "generator_speed_curves.png")
savefig(p, output_path)
println("    Figure saved to: $output_path")
