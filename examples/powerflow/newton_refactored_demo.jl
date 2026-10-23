"""
Demonstration of refactored Newton power flow solver.

Shows the new API with explicit setup/iteration separation for:
- Precise benchmarking (excludes setup overhead)
- Multiple solves with different initial conditions
- Workspace reuse for repeated solves
"""

using Powerful
using PowerfulCases
using Printf
using PowerfulPF: setup_powerflow_newton, solve_powerflow_newton!

println("=" ^ 70)
println("Newton Power Flow Solver - Refactored API Demo")
println("=" ^ 70)

# Load test system
println("\nLoading IEEE 14-bus system...")
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

# Create problem (required for setup)
prob = ACPowerFlowProblem(sys, PowerInjection())

# =============================================================================
# New API: Explicit setup + iteration separation
# =============================================================================

println("\n" * "=" ^ 70)
println("1. New API: Explicit Setup + Iterations")
println("=" ^ 70)

# One-time setup (expensive: sparsity detection, coloring, allocation)
println("\n[Setup Phase]")
workspace = setup_powerflow_newton(prob; verbose=true)

# Multiple solves with different initial conditions (fast!)
println("\n[Solve Phase 1: Default initial guess]")
y0 = copy(workspace.iter_dae.y)
sol1 = solve_powerflow_newton!(workspace, y0; abstol=1e-8, verbose=true)

println("\n[Solve Phase 2: Perturbed initial guess]")
y0_perturbed = y0 .* 1.05  # 5% perturbation
sol2 = solve_powerflow_newton!(workspace, y0_perturbed; abstol=1e-8, verbose=true)

# =============================================================================
# High-level API: One-shot solve (for comparison)
# =============================================================================

println("\n" * "=" ^ 70)
println("2. High-level API: One-shot solve")
println("=" ^ 70)

# The high-level API creates problem and solves in one call
sol3 = solve(prob, Newton())
println("High-level solve: $(sol3.converged ? "converged" : "failed") in $(sol3.iterations) iterations")

# =============================================================================
# Benchmarking: Iterations only (no setup overhead)
# =============================================================================

println("\n" * "=" ^ 70)
println("3. Precise Iteration Benchmarking")
println("=" ^ 70)

"""Minimum wall time (s) of `f()` over `n` runs (after one warm-up run)."""
function min_time(f; n = 20)
    f()
    return minimum(@elapsed(f()) for _ in 1:n)
end

println("\nBenchmarking iterations only (excludes setup):")
y_flat = copy(workspace.iter_dae.y)
t_iter = min_time(() -> solve_powerflow_newton!(workspace, copy(y_flat); abstol=1e-8, verbose=false))
@printf("  %.1f μs\n", 1e6 * t_iter)

println("\nBenchmarking full solve (includes setup):")
t_full = min_time(() -> solve(prob, Newton()))
@printf("  %.1f μs\n", 1e6 * t_full)

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)
println("\nNew API Benefits:")
println("  - Explicit setup/iteration separation")
println("  - Workspace reuse for multiple solves")
println("  - Precise benchmarking (excludes setup overhead)")
println("  - Pre-allocated buffers for zero-allocation iterations")
println("\nConvergence:")
println("  Solution 1 (default):    $(sol1.converged ? "OK" : "FAIL") in $(sol1.iterations) iterations")
println("  Solution 2 (perturbed):  $(sol2.converged ? "OK" : "FAIL") in $(sol2.iterations) iterations")
println("  Solution 3 (high-level): $(sol3.converged ? "OK" : "FAIL") in $(sol3.iterations) iterations")
println("\n" * "=" ^ 70)
