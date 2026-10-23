# DC Contingency Calculation (DCCC) Demo - Large Scale
#
# This example demonstrates DC contingency screening on the ACTIVSg70k system
# (70,000 buses, 88,207 branches). DC screening uses PTDF/LODF matrices to
# quickly estimate post-outage flows without solving full power flow.
#
# Key features:
# - LazyPTDF: Computes PTDF columns on-demand (avoids 50+ GB dense matrix)
# - LODF: Line Outage Distribution Factors for flow redistribution
# - Fast N-1 screening: ~8-30ms per contingency on 70k system

using Powerful
using PowerfulCases
using PowerfulCA

# =============================================================================
# Load 70k bus system
# =============================================================================
println("=" ^ 70)
println("DC Contingency Screening Demo - ACTIVSg70k (70,000 buses)")
println("=" ^ 70)

raw_path = load("synthetic/ACTIVSg70k").raw

println("\n[1] Loading system...")
t_load = @elapsed begin
    sys = SystemModel(parse_network(raw_path))
end
println("    Load time: $(round(t_load, digits=2))s")

n_bus = length(sys.models[:Bus])
n_branches = length(sys.models[:Branch])
n_transformers = length(sys.models[:Transformer2W])
n_total_branches = n_branches + n_transformers

println("    Buses: $n_bus")
println("    Branches: $n_branches")
println("    Transformers: $n_transformers")
println("    Total branch elements: $n_total_branches")

# =============================================================================
# Define contingencies
# =============================================================================
println("\n[2] Generating N-1 contingencies...")

# For demo: screen first 1000 branches (full N-1 would be 88k contingencies)
n_contingencies = 1000

# Get branch data to create contingencies
branches = sys.models[:Branch]
contingencies = Contingency[]

for idx in 1:min(n_contingencies, length(branches))
    from_bus = Int(branches.i[idx])
    to_bus = Int(branches.j[idx])
    ckt = String(branches.ckt[idx])
    name = "LINE_$(from_bus)_$(to_bus)_$ckt"
    push!(contingencies, Contingency(name, [BranchOutage(from_bus, to_bus, ckt)]))
end

println("    Generated $(length(contingencies)) contingencies")

# =============================================================================
# Run DC Contingency Screening
# =============================================================================
println("\n[3] Running DC contingency screening...")
println("    Using LazyPTDF (memory-efficient mode)")

t_screen = @elapsed begin
    result = screen_contingencies_dc(
        sys,
        contingencies;
        use_lazy_ptdf=true,  # Essential for large systems!
        loading_threshold=0.8,  # Report branches > 80% loading
        verbose=true
    )
end

println("\n    Screening time: $(round(t_screen, digits=2))s")
println("    Average per contingency: $(round(t_screen/length(contingencies)*1000, digits=1))ms")

# =============================================================================
# Print violation report
# =============================================================================
println("\n[4] Violation Report")
println("-" ^ 70)
print_dc_screening_report(result; max_violations=20)

# =============================================================================
# Programmatic access to results
# =============================================================================
println("\n[5] Programmatic Access")
println("-" ^ 70)

# Summary statistics
n_critical = length(result.critical_contingencies)
n_with_violations = count(r -> !isempty(r.violations), result.results)
max_loading = maximum(r.max_loading for r in result.results)

println("Total contingencies screened: $(length(result.results))")
println("Contingencies with violations: $n_with_violations")
println("Critical contingencies (>100% loading): $n_critical")
println("Maximum loading observed: $(round(max_loading * 100, digits=1))%")

# Show top 5 most severe contingencies
println("\nTop 5 most severe contingencies:")
sorted_results = sort(result.results, by=r -> r.max_loading, rev=true)
for (i, r) in enumerate(sorted_results[1:min(5, length(sorted_results))])
    println("  $i. $(rpad(r.name, 25)) max_loading=$(round(r.max_loading*100, digits=1))%")
    if !isempty(r.violations)
        for v in r.violations[1:min(3, length(r.violations))]
            println("     └─ Branch $(v.branch_idx): $(round(v.loading*100, digits=1))% ($(round(v.flow, digits=1)) MW)")
        end
    end
end

# =============================================================================
# Compare with full LODF (for small subset)
# =============================================================================
println("\n[6] Memory comparison")
println("-" ^ 70)

# Calculate what dense matrices would require
ptdf_size_gb = n_total_branches * n_bus * 8 / 1e9
lodf_size_gb = n_total_branches * n_total_branches * 8 / 1e9

println("Dense PTDF matrix would require: $(round(ptdf_size_gb, digits=1)) GB")
println("Dense LODF matrix would require: $(round(lodf_size_gb, digits=1)) GB")
println("LazyPTDF actual memory: ~50 MB (factorization + cache)")
println("Memory savings: $(round((ptdf_size_gb + lodf_size_gb) / 0.05, digits=0))x")

println("\n" * "=" ^ 70)
println("Demo complete!")
println("=" ^ 70)
