#!/usr/bin/env julia
"""
TDS fault test - verify tgrad fix works with events.
Uses ieee14 case which is known to work with faults.
"""

using Powerful
using PowerfulCases

println("=" ^ 60)
println("TDS Fault Test: ieee14 (verify tgrad fix)")
println("=" ^ 60)

# Load case
println("\n1. Loading ieee14 case...")

# Parse and solve power flow
sys = SystemModel(parse_network(load("ieee14").raw))

prob = ACPowerFlowProblem(sys)
pf_sol = solve(prob, Newton())
println("   Power flow: ", pf_sol.converged ? "converged" : "FAILED")

# Parse dynamic models
println("\n2. Parsing DYR...")
dyr_data = parse_dynamics(load("ieee14").dyr, sys)

# BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
events = [
    BusFault(0.5, 0.6, 10, "F1"; xf=0.003, rf=0.0),
]

println("\n3. Fault Events:")
println("   t=0.5s: Apply 3-phase fault at bus 10")
println("   t=0.6s: Clear fault")

# 100% constant impedance loads for numerical stability during fault
config = SimConfig(solver = (TDS = TDSConfig(
    convert_pq_to_zip=true, p_zip=(1.0,0.0,0.0), q_zip=(1.0,0.0,0.0)),))

# Create TDS problem
init = setup_tds(pf_sol, dyr_data; config=config)
tds_prob = TDSProblem(init; tspan=(0.0, 2.0), events=events)

println("\n4. TDSProblem: tspan=$(tds_prob.tspan)")

# First solve (includes JIT compilation)
println("\n5. First solve (with compilation)...")
t1 = @elapsed sol1 = solve(tds_prob; initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.01)
println("   Time: $(round(t1, digits=2))s")
println("   Result: $(sol1.retcode)")

if sol1.retcode == :Success
    println("   Steps: $(length(sol1.t))")

    # Second solve (warmed up)
    println("\n6. Second solve (warmed up)...")
    t2 = @elapsed sol2 = solve(tds_prob; initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.01)
    println("   Time: $(round(t2, digits=2))s")
    println("   Result: $(sol2.retcode)")
    println("   Steps: $(length(sol2.t))")

    # Third solve
    println("\n7. Third solve (warmed up)...")
    t3 = @elapsed sol3 = solve(tds_prob; initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.01)
    println("   Time: $(round(t3, digits=2))s")
    println("   Result: $(sol3.retcode)")

    println("\n" * "=" ^ 60)
    println("TIMING SUMMARY")
    println("=" ^ 60)
    println("First solve (with JIT):  $(round(t1, digits=2))s")
    println("Second solve (warm):     $(round(t2, digits=2))s")
    println("Third solve (warm):      $(round(t3, digits=2))s")
else
    println("\n   FAILED at t=$(sol1.t[end])")
end

println("\n" * "=" ^ 60)
