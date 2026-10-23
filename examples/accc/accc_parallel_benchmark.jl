"""
Parallel ACCC Benchmark - ACTIVSg2000 (2000 buses)

Compares sequential vs parallel ACCC performance with ~200 branch contingencies.

Usage:
    JULIA_NUM_THREADS=8 julia --project examples/accc/accc_parallel_benchmark.jl
"""

using Powerful
using PowerfulCases
using PowerfulCA
using Printf

# Load the 2000-bus synthetic system
println("=" ^ 70)
println("PARALLEL ACCC BENCHMARK - ACTIVSg2000 (2000 buses)")
println("=" ^ 70)
println()

println("Loading ACTIVSg2000 case...")
t_load = @elapsed begin
    case = load("synthetic/ACTIVSg2000")
    sys = SystemModel(parse_network(case.raw))
end
println("  Load time: $(round(t_load, digits=2))s")

# Get branch info for contingency generation
branches = sys.models.Branch
n_branches = length(branches)
println("  Buses: $(length(sys.models.Bus))")
println("  Branches: $n_branches")
println("  PV generators: $(length(sys.models.PV))")
println("  PQ loads: $(length(sys.models.PQ))")
println()

# Generate contingencies - take first 200 branches
n_contingencies = min(200, n_branches)
println("Generating $n_contingencies branch contingencies...")

contingencies = Vector{Contingency}(undef, n_contingencies)
for i in 1:n_contingencies
    from_bus = Int(branches.i[i])
    to_bus = Int(branches.j[i])
    circuit = String(branches.ckt[i])
    name = "BRANCH_$(from_bus)_$(to_bus)_$(circuit)"
    contingencies[i] = Contingency(name, [BranchOutage(from_bus, to_bus, circuit)])
end
println("  Created $(length(contingencies)) contingencies")
println()

# Report thread count
n_threads = Threads.nthreads()
println("Julia threads: $n_threads")
if n_threads == 1
    println("  WARNING: Running with 1 thread. Set JULIA_NUM_THREADS for parallel speedup.")
end
println()

# Benchmark sequential ACCC
println("=" ^ 70)
println("SEQUENTIAL ACCC")
println("=" ^ 70)
t_seq = @elapsed begin
    result_seq = solve_accc(sys, contingencies; verbose=true)
end
n_converged_seq = count(r -> r.converged, result_seq.contingencies)
println()
println("Sequential total time: $(round(t_seq, digits=3))s")
println("  Per contingency: $(round(t_seq / n_contingencies * 1000, digits=2))ms")
println("  Converged: $n_converged_seq / $n_contingencies")
println()

# Benchmark parallel ACCC
println("=" ^ 70)
println("PARALLEL ACCC ($n_threads threads)")
println("=" ^ 70)
t_par = @elapsed begin
    result_par = solve_accc_parallel(sys, contingencies; verbose=true)
end
n_converged_par = count(r -> r.converged, result_par.contingencies)
println()
println("Parallel total time: $(round(t_par, digits=3))s")
println("  Per contingency: $(round(t_par / n_contingencies * 1000, digits=2))ms")
println("  Converged: $n_converged_par / $n_contingencies")
println()

# Verify results match
println("=" ^ 70)
println("VERIFICATION")
println("=" ^ 70)
max_diff = 0.0
n_matched = 0
n_both_converged = 0
for i in 1:n_contingencies
    r_seq = result_seq.contingencies[i]
    r_par = result_par.contingencies[i]

    if r_seq.converged && r_par.converged
        global n_both_converged += 1
        diff = maximum(abs.(r_seq.y .- r_par.y))
        global max_diff = max(max_diff, diff)
        if diff < 1e-6
            global n_matched += 1
        end
    elseif r_seq.converged == r_par.converged
        global n_matched += 1
    end
end
println("  Both converged: $n_both_converged / $n_contingencies")
println("  Results matched: $n_matched / $n_contingencies")
println("  Max solution difference: $(round(max_diff, sigdigits=3))")
println()

# Summary
println("=" ^ 70)
println("PERFORMANCE SUMMARY")
println("=" ^ 70)
speedup = t_seq / t_par
efficiency = speedup / n_threads * 100
println()
@printf("  Sequential time:  %8.3fs\n", t_seq)
@printf("  Parallel time:    %8.3fs\n", t_par)
@printf("  Speedup:          %8.2fx\n", speedup)
@printf("  Efficiency:       %8.1f%% (speedup / threads)\n", efficiency)
@printf("  Threads:          %8d\n", n_threads)
println()

if speedup > 1.5
    println("  Parallel ACCC provides $(round(speedup, digits=1))x speedup!")
else
    println("  Limited speedup - try increasing JULIA_NUM_THREADS")
end
println()
