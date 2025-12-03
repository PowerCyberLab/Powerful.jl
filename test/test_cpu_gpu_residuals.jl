"""
Test that CPU and GPU residual calculations produce identical results.
"""

using TestItems

@testitem "CPU vs GPU Residual Equivalence" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData
    using CUDA

    # Skip if no GPU available
    if !CUDA.functional()
        @info "CUDA not functional - skipping GPU tests"
        return
    end

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

    # Evaluate CPU residuals
    evaluate_system_residuals!(sys_ws_cpu, dae_cpu)
    g_cpu = copy(dae_cpu.g)

    # Create GPU system
    sys_gpu = to_gpu(sys_cpu)
    dae_gpu = DAE(Float64, sys_gpu)
    sys_ws_gpu = SystemResidualWorkspace(Float64, sys_gpu)

    # Evaluate GPU residuals
    evaluate_system_residuals!(sys_ws_gpu, dae_gpu)
    CUDA.synchronize()
    g_gpu = Array(dae_gpu.g)  # Transfer back to CPU

    # Test: Residual values match to machine precision
    max_diff = maximum(abs.(g_cpu .- g_gpu))
    @test max_diff < 1e-10

    rms_diff = sqrt(sum((g_cpu .- g_gpu).^2) / length(g_cpu))
    @test rms_diff < 1e-10

    # Test: Relative error for non-zero values
    nonzero_mask = abs.(g_cpu) .> 1e-12
    if any(nonzero_mask)
        rel_errors = abs.((g_cpu[nonzero_mask] .- g_gpu[nonzero_mask]) ./ g_cpu[nonzero_mask])
        max_rel_err = maximum(rel_errors)
        @test max_rel_err < 1e-10
    end

    # Test: Same non-zero pattern
    cpu_nonzero = count(x -> abs(x) > 1e-12, g_cpu)
    gpu_nonzero = count(x -> abs(x) > 1e-12, g_gpu)
    @test cpu_nonzero == gpu_nonzero

    cpu_nz_indices = findall(x -> abs(x) > 1e-12, g_cpu)
    gpu_nz_indices = findall(x -> abs(x) > 1e-12, g_gpu)
    @test cpu_nz_indices == gpu_nz_indices

    # Test: Minimal allocations (warm up first to avoid compilation)
    for _ in 1:3
        evaluate_system_residuals!(sys_ws_cpu, dae_cpu)
        evaluate_system_residuals!(sys_ws_gpu, dae_gpu)
        CUDA.synchronize()
    end

    alloc_cpu = @allocated evaluate_system_residuals!(sys_ws_cpu, dae_cpu)
    @test alloc_cpu < 5000  # Should be ~1.2KB, but allow some margin

    alloc_gpu = @allocated begin
        evaluate_system_residuals!(sys_ws_gpu, dae_gpu)
        CUDA.synchronize()
    end
    @test alloc_gpu < 20000  # Should be ~6KB, but kernel launch has overhead
end
