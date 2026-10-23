#!/usr/bin/env julia
"""
DC Power Flow Example

This example demonstrates the DC power flow formulation, which provides a fast
linear approximation of AC power flow by making the following assumptions:
- All voltage magnitudes = 1.0 pu
- Small angle differences: sin(θ_i - θ_j) ≈ θ_i - θ_j
- Branch resistance negligible: r << x
- Reactive power ignored completely

DC power flow is useful for:
- Fast contingency analysis
- Initial guess for AC power flow
- Large-scale transmission planning
- Security-constrained unit commitment

The DC formulation typically converges in 1-2 iterations since the equations
are linear. It approximates AC angles well (5-15% error) but ignores voltage
magnitudes and reactive power completely.
"""

using Powerful
using PowerfulCases
using PowerfulPF: setup_powerflow_newton, solve_powerflow_newton!
using PowerfulCore: get_addresses
using DataFrames
using Printf

println("=" ^ 80)
println("DC Power Flow Example - IEEE 14-bus System")
println("=" ^ 80)

# Load test case
case_path = load("ieee14").raw
case = parse_network(case_path)
println("\nLoaded: IEEE 14-bus system")
println("  Buses: ", length(case.buses))
println("  Generators: ", length(case.generators))
println("  Branches: ", length(case.branches.i))

# Create system model (same for both DC and AC)
sys = SystemModel(case)

println("\n" * "=" ^ 80)
println("Part 1: DC Power Flow Solution")
println("=" ^ 80)

# Create DC power flow problem
prob_dc = ACPowerFlowProblem(sys, DCPowerFlow())

# Setup solver workspace
ws_dc = setup_powerflow_newton(prob_dc)

println("\nDC Problem dimensions:")
println("  Variables: $(length(ws_dc.iter_dae.y)) (angles only)")
println("  Residuals: $(length(ws_dc.iter_dae.g)) (active power only)")

# Solve DC power flow
y0_dc = copy(ws_dc.iter_dae.y)
sol_dc = solve_powerflow_newton!(ws_dc, y0_dc; abstol=1e-10, maxiters=50, verbose=true)

if sol_dc.converged
    println("\n✓ DC converged in $(sol_dc.iterations) iterations (linear problem!)")
    println("  Final residual: $(sol_dc.residual_norm)")
else
    println("\n✗ DC did not converge")
end

println("\n" * "=" ^ 80)
println("Part 2: AC Power Flow Solution (for comparison)")
println("=" ^ 80)

# Create AC power flow problem with PowerInjection formulation
prob_ac = ACPowerFlowProblem(sys, PowerInjection())

# Setup solver workspace
ws_ac = setup_powerflow_newton(prob_ac)

println("\nAC Problem dimensions:")
println("  Variables: $(length(ws_ac.iter_dae.y)) (angles + voltages)")
println("  Residuals: $(length(ws_ac.iter_dae.g)) (P + Q equations)")

# Solve AC power flow
y0_ac = copy(ws_ac.iter_dae.y)
sol_ac = solve_powerflow_newton!(ws_ac, y0_ac; abstol=1e-8, verbose=true)

if sol_ac.converged
    println("\n✓ AC converged in $(sol_ac.iterations) iterations")
    println("  Final residual: $(sol_ac.residual_norm)")
else
    println("\n✗ AC did not converge")
end

println("\n" * "=" ^ 80)
println("Part 3: Comparison of DC vs AC Results")
println("=" ^ 80)

if sol_dc.converged && sol_ac.converged
    # Extract angles from both solutions
    am_dc = prob_dc.address_manager
    am_ac = prob_ac.address_manager

    theta_dc = sol_dc.dae.y[get_addresses(am_dc, :Bus, :theta)]
    theta_ac = sol_ac.dae.y[get_addresses(am_ac, :Bus, :theta)]

    # Convert to degrees for readability
    theta_dc_deg = theta_dc .* (180/π)
    theta_ac_deg = theta_ac .* (180/π)

    # Calculate errors
    abs_error = abs.(theta_dc_deg .- theta_ac_deg)
    rel_error = abs.(abs_error ./ theta_ac_deg) .* 100

    # Create comparison table
    bus_numbers = sys.models.Bus.i

    println("\nBus Angle Comparison (degrees):")
    println("Bus | DC Angle | AC Angle | Abs Error | Rel Error %")
    println("-" ^ 60)

    for i in eachindex(bus_numbers)
        @printf("%3d | %8.3f | %8.3f | %9.4f | %11.2f\n",
                bus_numbers[i], theta_dc_deg[i], theta_ac_deg[i],
                abs_error[i], rel_error[i])
    end

    println("\nError Statistics:")
    println("  Max absolute error: $(maximum(abs_error)) degrees")
    println("  Mean absolute error: $(sum(abs_error)/length(abs_error)) degrees")
    println("  Max relative error: $(maximum(rel_error[2:end])) % (excluding slack)")
    println("  Mean relative error: $(sum(rel_error[2:end])/length(rel_error[2:end])) % (excluding slack)")

    println("\n" * "=" ^ 80)
    println("Performance Comparison")
    println("=" ^ 80)

    println("\nIterations:")
    println("  DC: $(sol_dc.iterations) iterations (linear problem)")
    println("  AC: $(sol_ac.iterations) iterations (nonlinear problem)")
    println("  Speedup: $(sol_ac.iterations / sol_dc.iterations)x fewer iterations")

    println("\nProblem size:")
    println("  DC: $(length(ws_dc.iter_dae.y)) variables")
    println("  AC: $(length(ws_ac.iter_dae.y)) variables")
    println("  Reduction: $(round(100*(1 - length(ws_dc.iter_dae.y)/length(ws_ac.iter_dae.y)), digits=1))% fewer variables")

    println("\n" * "=" ^ 80)
    println("Key Observations")
    println("=" ^ 80)

    println("""
    1. DC converged in $(sol_dc.iterations) iterations vs $(sol_ac.iterations) for AC
       → DC is linear, so Newton converges in 1 iteration (or 2 with roundoff)

    2. DC angles approximate AC angles with $(round(sum(rel_error[2:end])/length(rel_error[2:end]), digits=1))% mean error
       → Good enough for contingency screening and initial estimates

    3. DC ignores voltage magnitudes and reactive power completely
       → Cannot detect voltage stability issues or reactive power limits

    4. DC problem is $(round(100*(1 - length(ws_dc.iter_dae.y)/length(ws_ac.iter_dae.y)), digits=1))% smaller than AC
       → Scales better for very large systems (thousands of buses)

    When to use DC power flow:
    - Fast screening of N-1 contingencies
    - Initial guess for AC power flow
    - Large-scale planning studies where voltage is not critical
    - Unit commitment and economic dispatch

    When NOT to use DC power flow:
    - Voltage stability analysis
    - Reactive power planning
    - Distribution systems with high R/X ratios
    - Systems with significant reactive power issues
    """)
end

println("\n" * "=" ^ 80)
println("Example completed successfully!")
println("=" ^ 80)
