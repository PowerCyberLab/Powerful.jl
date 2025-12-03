# Test what ACTUALLY happens with multi-island systems

using Powerful
using PowerFlowData

# Load multi-island case
case = parse_network(joinpath(pkgdir(Powerful), "cases", "ieee14_island.raw"))
sys = SystemModel(case)

println("Test 1: Create problem with multi-island system directly")
try
    prob = ACPowerFlowProblem(sys)
    println("  ✓ Problem created: $(typeof(prob))")
    println("  sys type: $(typeof(prob.sys))")

    sol = solve(prob, NewtonRaphson())
    println("  ✓ Solution: $(typeof(sol))")
    println("  Converged: $(sol.converged)")
catch e
    println("  ✗ Error: $(typeof(e))")
    println("  $(e)")
end

println("\nTest 2: Split islands manually first")
result = split_islands(sys)
println("  split_islands returned: $(typeof(result))")

if result isa Vector{SystemModel}
    println("  Got $(length(result)) islands")

    # Try creating problem with first island only
    prob1 = ACPowerFlowProblem(result[1])
    sol1 = solve(prob1, NewtonRaphson())
    println("  Island 1 solved: converged=$(sol1.converged), type=$(typeof(sol1))")
end
