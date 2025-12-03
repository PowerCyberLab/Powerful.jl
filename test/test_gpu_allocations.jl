"""
Test file to identify and measure allocations in GPU residual computation loop.

This file creates targeted micro-tests to isolate allocation sources in the
evaluate_system_residuals! hot path.
"""

using Test
using Powerful
using Powerful.PowerCore
using Powerful.PowerCore: evaluate_model_residuals!
using Powerful.Models
using PowerFlowData
using CUDA

# Load test case
case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
case = PowerFlowData.parse_network(case_path)

# Parse models
bus = parse_model(Bus, case)
pq = parse_model(PQ, case)

# Create CPU system
am = AddressManager()
sys_cpu = SystemModel(am, (Bus = bus, PQ = pq))
dae_cpu = DAE(Float64, sys_cpu)
sys_ws_cpu = SystemResidualWorkspace(Float64, sys_cpu)

# Transfer to GPU
sys_gpu = to_gpu(sys_cpu)
dae_gpu = DAE(Float64, sys_gpu)
sys_ws_gpu = SystemResidualWorkspace(Float64, sys_gpu)

println("=" ^ 80)
println("GPU Allocation Analysis")
println("=" ^ 80)

# ==============================================================================
# Test 1: Baseline - evaluate_system_residuals! (full call)
# ==============================================================================
println("\n[Test 1] Full evaluate_system_residuals! call (with CUDA.synchronize)")

# Warm-up
evaluate_system_residuals!(sys_ws_gpu, dae_gpu)
CUDA.synchronize()

# Measure
allocs_full = @allocated begin
    evaluate_system_residuals!(sys_ws_gpu, dae_gpu)
    CUDA.synchronize()
end

println("  Allocations: ", allocs_full, " bytes")

# ==============================================================================
# Test 2: Just the top-level function (no synchronize)
# ==============================================================================
println("\n[Test 2] evaluate_system_residuals! WITHOUT CUDA.synchronize()")

allocs_no_sync = @allocated begin
    evaluate_system_residuals!(sys_ws_gpu, dae_gpu)
end

println("  Allocations: ", allocs_no_sync, " bytes")
println("  Sync cost: ", allocs_full - allocs_no_sync, " bytes")

# ==============================================================================
# Test 3: Manual iteration over model workspaces (mimics internal loop)
# ==============================================================================
println("\n[Test 3] Manual loop over model workspaces")

allocs_manual = @allocated begin
    fill!(dae_gpu.g, 0.0)

    for (model_name, workspace) in pairs(sys_ws_gpu.model_workspaces)
        model_vec = getfield(sys_ws_gpu.model_vectors, model_name)
        metadata = getfield(sys_ws_gpu.model_metadata, model_name)
        evaluate_model_residuals!(workspace, model_vec, dae_gpu, metadata)
    end
end

println("  Allocations: ", allocs_manual, " bytes")

# ==============================================================================
# Test 4: Isolated operations - fill!
# ==============================================================================
println("\n[Test 4] Just fill!(dae_gpu.g, 0.0)")

allocs_fill = @allocated begin
    fill!(dae_gpu.g, 0.0)
end

println("  Allocations: ", allocs_fill, " bytes")

# ==============================================================================
# Test 5: Isolated operations - pairs() iterator
# ==============================================================================
println("\n[Test 5] Just iterating pairs(sys_ws_gpu.model_workspaces)")

allocs_pairs = @allocated begin
    for (model_name, workspace) in pairs(sys_ws_gpu.model_workspaces)
        # Do nothing, just iterate
    end
end

println("  Allocations: ", allocs_pairs, " bytes")

# ==============================================================================
# Test 6: Isolated operations - getfield on NamedTuple
# ==============================================================================
println("\n[Test 6] getfield on model_vectors NamedTuple")

allocs_getfield = @allocated begin
    for (model_name, _) in pairs(sys_ws_gpu.model_workspaces)
        model_vec = getfield(sys_ws_gpu.model_vectors, model_name)
        metadata = getfield(sys_ws_gpu.model_metadata, model_name)
    end
end

println("  Allocations: ", allocs_getfield, " bytes")

# ==============================================================================
# Test 7: Single model evaluation (Bus)
# ==============================================================================
println("\n[Test 7] Evaluate single model (Bus) residuals")

bus_ws = sys_ws_gpu.model_workspaces.Bus
bus_vec = sys_ws_gpu.model_vectors.Bus
bus_meta = sys_ws_gpu.model_metadata.Bus

allocs_bus = @allocated begin
    evaluate_model_residuals!(bus_ws, bus_vec, dae_gpu, bus_meta)
end

println("  Allocations: ", allocs_bus, " bytes")

# ==============================================================================
# Test 8: Single model evaluation (PQ)
# ==============================================================================
println("\n[Test 8] Evaluate single model (PQ) residuals")

pq_ws = sys_ws_gpu.model_workspaces.PQ
pq_vec = sys_ws_gpu.model_vectors.PQ
pq_meta = sys_ws_gpu.model_metadata.PQ

allocs_pq = @allocated begin
    evaluate_model_residuals!(pq_ws, pq_vec, dae_gpu, pq_meta)
end

println("  Allocations: ", allocs_pq, " bytes")

# ==============================================================================
# Test 9: Breakdown of evaluate_model_residuals! steps (PQ)
# ==============================================================================
println("\n[Test 9] Step-by-step PQ model evaluation")

# Warmup for accurate measurements
gather_variables!(pq_ws, dae_gpu)
zero_residuals!(pq_ws)
calc_residuals!(pq_ws, pq_vec)
scatter_residuals!(pq_ws, dae_gpu, pq_meta)
CUDA.synchronize()

println("  [9a] gather_variables!")
allocs_gather = @allocated begin
    gather_variables!(pq_ws, dae_gpu)
end
println("      Allocations: ", allocs_gather, " bytes")

println("  [9b] zero_residuals!")
allocs_zero = @allocated begin
    zero_residuals!(pq_ws)
end
println("      Allocations: ", allocs_zero, " bytes")

println("  [9c] calc_residuals!")
allocs_calc = @allocated begin
    calc_residuals!(pq_ws, pq_vec)
end
println("      Allocations: ", allocs_calc, " bytes")

println("  [9d] scatter_residuals!")
allocs_scatter = @allocated begin
    scatter_residuals!(pq_ws, dae_gpu, pq_meta)
end
println("      Allocations: ", allocs_scatter, " bytes")

# ==============================================================================
# Test 10: Analyze gather_variables! in detail (with new kernel implementation)
# ==============================================================================
println("\n[Test 10] Detailed gather_variables! analysis (optimized)")

println("  [10a] pairs(pq_ws.var_addrs)")
allocs_var_pairs = @allocated begin
    for (var_name, addrs) in pairs(pq_ws.var_addrs)
        # Do nothing
    end
end
println("      Allocations: ", allocs_var_pairs, " bytes")

println("  [10b] getfield(pq_ws.var_vals, var_name)")
allocs_var_getfield = @allocated begin
    for (var_name, addrs) in pairs(pq_ws.var_addrs)
        vals = getfield(pq_ws.var_vals, var_name)
    end
end
println("      Allocations: ", allocs_var_getfield, " bytes")

println("  [10c] Custom gather kernel (current implementation)")
allocs_kernel_gather = @allocated begin
    gather_variables!(pq_ws, dae_gpu)
end
println("      Allocations: ", allocs_kernel_gather, " bytes")

# ==============================================================================
# Test 11: CUDA kernel launch overhead
# ==============================================================================
println("\n[Test 11] CUDA kernel launch overhead (calc_residuals!)")

# First check if this is actually launching a CUDA kernel
println("  PQ model parameters:")
println("    Length: ", length(pq_vec))
println("    Type of pq.pl: ", typeof(pq_vec.pl))
println("    Type of var_vals.v: ", typeof(pq_ws.var_vals.v))
println("    Type of active_mask: ", typeof(pq_ws.active_mask))

# Multiple calls to measure consistency
allocs_kernel_1 = @allocated calc_residuals!(pq_ws, pq_vec)
allocs_kernel_2 = @allocated calc_residuals!(pq_ws, pq_vec)
allocs_kernel_3 = @allocated calc_residuals!(pq_ws, pq_vec)

println("  Kernel call 1: ", allocs_kernel_1, " bytes")
println("  Kernel call 2: ", allocs_kernel_2, " bytes")
println("  Kernel call 3: ", allocs_kernel_3, " bytes")

# ==============================================================================
# Test 12: scatter_residuals! CUDA kernel analysis
# ==============================================================================
println("\n[Test 12] scatter_residuals! CUDA kernel analysis")

# Check residual access patterns
println("  PQ residuals:")
for res in pq_meta.residuals
    println("    ", res.name, " - ", res.access)
end

# Multiple calls
allocs_scatter_1 = @allocated scatter_residuals!(pq_ws, dae_gpu, pq_meta)
allocs_scatter_2 = @allocated scatter_residuals!(pq_ws, dae_gpu, pq_meta)
allocs_scatter_3 = @allocated scatter_residuals!(pq_ws, dae_gpu, pq_meta)

println("  Scatter call 1: ", allocs_scatter_1, " bytes")
println("  Scatter call 2: ", allocs_scatter_2, " bytes")
println("  Scatter call 3: ", allocs_scatter_3, " bytes")

# ==============================================================================
# Summary
# ==============================================================================
println("\n" * "=" ^ 80)
println("ALLOCATION SUMMARY")
println("=" ^ 80)
println("Full call (with sync):        ", allocs_full, " bytes")
println("Without sync:                  ", allocs_no_sync, " bytes")
println("Manual loop:                   ", allocs_manual, " bytes")
println("  ├─ fill!:                    ", allocs_fill, " bytes")
println("  ├─ pairs() iteration:        ", allocs_pairs, " bytes")
println("  ├─ getfield calls:           ", allocs_getfield, " bytes")
println("  ├─ Bus model:                ", allocs_bus, " bytes")
println("  └─ PQ model:                 ", allocs_pq, " bytes")
println()
println("PQ model breakdown:")
println("  ├─ gather_variables!:        ", allocs_gather, " bytes")
println("  ├─ zero_residuals!:          ", allocs_zero, " bytes")
println("  ├─ calc_residuals!:          ", allocs_calc, " bytes")
println("  └─ scatter_residuals!:       ", allocs_scatter, " bytes")
println()
println("gather_variables! breakdown (optimized):")
println("  ├─ pairs() iteration:        ", allocs_var_pairs, " bytes")
println("  ├─ getfield calls:           ", allocs_var_getfield, " bytes")
println("  └─ kernel gather (total):    ", allocs_kernel_gather, " bytes")
println("=" ^ 80)

# ==============================================================================
# Results Analysis
# ==============================================================================
println("\nOPTIMIZATION RESULTS:")
println()
println("BEFORE (old broadcast method):")
println("  gather_variables!:              ~80 KB")
println("  evaluate_system_residuals!:     ~34 KB")
println()
println("AFTER (custom gather kernel):")
println("  gather_variables!:              ", allocs_gather, " bytes (~", round(allocs_gather/1024, digits=1), " KB)")
println("  evaluate_system_residuals!:     ", allocs_full, " bytes (~", round(allocs_full/1024, digits=1), " KB)")
println()

reduction_gather = (1 - allocs_gather/82656) * 100
reduction_full = (1 - allocs_full/34688) * 100

println("IMPROVEMENT:")
println("  gather_variables! reduction:    ", round(reduction_gather, digits=1), "%")
println("  Full call reduction:            ", round(reduction_full, digits=1), "%")
println()

if allocs_gather < 5000 && allocs_full < 15000
    println("✓ SUCCESS: Custom gather kernel eliminated GPU indexing allocations!")
    println("  Remaining allocations are minimal CUDA kernel launch overhead.")
    println("  Expected ~1 KB per kernel launch, which is unavoidable.")
else
    println("⚠ WARNING: Allocations higher than expected")
    if allocs_gather > 5000
        println("  gather_variables! still allocating: ", allocs_gather, " bytes")
    end
    if allocs_full > 15000
        println("  Full call still allocating: ", allocs_full, " bytes")
    end
end
