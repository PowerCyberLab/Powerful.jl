"""
AMD Addressing Example

Demonstrates the unified API for AMD (Approximate Minimum Degree) address reordering.
AMD reorders variables to minimize fill-in during sparse LU factorization,
improving performance on larger systems.
"""

using Powerful
using PowerfulCases
using PowerfulCore.Foundation: AMDAddressing

# Load a test case
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

# ============================================================================
# Default addressing (natural order)
# ============================================================================
println("Default Addressing (natural order):")
println("=" ^ 60)

prob_natural = ACPowerFlowProblem(sys)
sol_natural = solve(prob_natural, Newton())

println("Converged: ", sol_natural.converged)
println("Iterations: ", sol_natural.iterations)
println("Residual: ", sol_natural.residual_norm)
println()

# ============================================================================
# AMD addressing (optimized order for sparse factorization)
# ============================================================================
println("AMD Addressing (optimized sparse order):")
println("=" ^ 60)

# Just one line with addressing parameter!
prob_amd = ACPowerFlowProblem(sys; addressing=AMDAddressing())
sol_amd = solve(prob_amd, Newton())

println("Converged: ", sol_amd.converged)
println("Iterations: ", sol_amd.iterations)
println("Residual: ", sol_amd.residual_norm)
println()

# ============================================================================
# Verify solutions match (with permutation)
# ============================================================================
println("Verification:")
println("=" ^ 60)

# Get the permutation used
perm = prob_amd.address_manager.applied_permutation
inv_perm = invperm(perm)

# Unpermute AMD solution and compare with natural
y_amd_unpermuted = sol_amd.dae.y[inv_perm]
max_diff = maximum(abs.(sol_natural.dae.y .- y_amd_unpermuted))

println("Maximum difference between solutions: ", max_diff)
println("Solutions match: ", max_diff < 1e-10 ? "✓ Yes" : "✗ No")
println()

# ============================================================================
# Summary
# ============================================================================
println("Summary:")
println("=" ^ 60)
println("Natural addressing: Contiguous variable layout")
println("AMD addressing: Optimized for sparse LU factorization")
println()
println("Both produce identical solutions (just stored differently).")
println("AMD addressing benefits large systems with sparse Jacobians.")
println("=" ^ 60)
