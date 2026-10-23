# Profile DC power flow to find bottleneck

using Powerful
using PowerfulCases
using PowerfulCA
using PowerfulCA: solve_dc_powerflow

raw_path = load("synthetic/ACTIVSg70k").raw
sys = SystemModel(parse_network(raw_path))
println("System loaded: $(length(sys.models[:Bus])) buses, $(length(sys.models[:Branch])) branches\n")

println("=== DC Power Flow Benchmark ===\n")

# Warm up
solve_dc_powerflow(sys)

# Benchmark
t_dc = @elapsed result = solve_dc_powerflow(sys)
println("DC Power Flow: $(round(t_dc*1000, digits=1))ms")

println("\n=== Comparison ===")
println("Before fix: ~34,000ms (34s)")
println("After fix:  $(round(t_dc*1000, digits=1))ms")
println("Speedup:    $(round(34000/t_dc/1000, digits=0))x")
