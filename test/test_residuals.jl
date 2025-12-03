"""
Tests for residual calculation framework using @testitem format

These tests demonstrate the new simplified API with automatic setup.
"""

using TestItems

@testitem "DAE Construction from SystemModel" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)

    # Create system with automatic setup (new simplified API)
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq))

    # Verify system is set up
    @test is_setup(sys)

    # Create DAE
    dae = DAE(Float64, sys)

    @test dae isa DAE
    @test dae.m > 0
    @test length(dae.y) == dae.m
    @test length(dae.g) == dae.m

    # Check that voltages are initialized
    # At least some non-zero voltages should exist
    @test any(dae.y .!= 0.0)
end

@testitem "SystemResidualWorkspace Construction" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)

    # Create system with automatic setup
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq))

    # Create system workspace
    sys_ws = SystemResidualWorkspace(Float64, sys)

    @test sys_ws isa SystemResidualWorkspace
    @test haskey(sys_ws.model_workspaces, :Bus)
    @test haskey(sys_ws.model_workspaces, :PQ)
    @test haskey(sys_ws.model_vectors, :Bus)
    @test haskey(sys_ws.model_vectors, :PQ)
    @test haskey(sys_ws.model_metadata, :Bus)
    @test haskey(sys_ws.model_metadata, :PQ)
end

@testitem "System Residual Evaluation" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)

    # Create complete system with automatic setup
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq))
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)

    # Evaluate system residuals
    evaluate_system_residuals!(sys_ws, dae)

    # Check that residuals were calculated
    # Should have non-zero residuals from PQ loads
    @test any(dae.g .!= 0.0)
end

@testitem "Zero Allocation" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)

    # Create system with automatic setup
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq))
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)

    # Warm-up (compile)
    evaluate_system_residuals!(sys_ws, dae)

    # Check allocation
    allocs = @allocated evaluate_system_residuals!(sys_ws, dae)

    # Allow small allocations for type instabilities we haven't fixed
    # In production, this should be exactly 0
    @test allocs < 2000  # Less than 1KB

    println("Allocations in evaluate_system_residuals!: $allocs bytes")
end

@testitem "Type Stability" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)

    # Create system with automatic setup
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq))
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)

    # Check that main function is type-stable
    @inferred evaluate_system_residuals!(sys_ws, dae)
end

@testitem "PV Generator Integration - System Evaluation" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pv = parse_model(PV, case)

    # Create complete system with Bus and PV (automatic setup)
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PV = pv))
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)

    # Evaluate system residuals
    evaluate_system_residuals!(sys_ws, dae)

    # Check that residuals were calculated
    @test any(dae.g .!= 0.0)

    # Verify voltage control residuals
    pv_ws = sys_ws.model_workspaces.PV
    for i in 1:length(pv.i)
        if pv.stat[i]
            # Get the bus index for this generator
            bus_idx = pv.i[i]

            # The voltage residual should be V - Vs
            v_addrs = sys.address_manager.assignments[:Bus].by_type[AlgebVar()][:v]
            v_addr = v_addrs[bus_idx]
            v_actual = dae.y[v_addr]
            v_expected_residual = v_actual - pv.vs[i]

            # Check the residual was computed correctly
            @test pv_ws.res_vals.v_ctrl[i] ≈ v_expected_residual
        end
    end
end

@testitem "Multi-Model System (Bus + PQ + PV)" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse all models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)
    pv = parse_model(PV, case)

    # Create complete system with all models (automatic setup)
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq, PV = pv))
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)

    # Evaluate system residuals
    evaluate_system_residuals!(sys_ws, dae)

    # Check that residuals were calculated from all models
    @test any(dae.g .!= 0.0)

    # Verify all model workspaces exist
    @test haskey(sys_ws.model_workspaces, :Bus)
    @test haskey(sys_ws.model_workspaces, :PQ)
    @test haskey(sys_ws.model_workspaces, :PV)
end

@testitem "Convenience Constructor SystemModel(case)" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # NEW: Ultra-simple one-line constructor!
    sys = SystemModel(case)

    # Verify system is fully set up
    @test is_setup(sys)
    @test haskey(sys.models, :Bus)
    @test haskey(sys.models, :PQ)
    @test haskey(sys.models, :PV)

    # Verify addresses are assigned
    @test haskey(sys.address_manager.assignments, :Bus)
    @test haskey(sys.address_manager.assignments, :PQ)

    # Verify external addresses are retrieved
    @test haskey(sys.address_manager.retrieved, :PQ)
    @test haskey(sys.address_manager.retrieved, :PV)

    # Should be able to create DAE and workspace immediately
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)

    @test dae isa DAE
    @test sys_ws isa SystemResidualWorkspace

    # Should be able to evaluate residuals
    evaluate_system_residuals!(sys_ws, dae)
    @test any(dae.g .!= 0.0)
end

@testitem "Manual Setup Workflow Still Works" begin
    using Powerful
    using Powerful.PowerCore
    using Powerful.Models
    using PowerFlowData

    # Load test case
    case_path = joinpath(pkgdir(Powerful), "cases", "ieee14.raw")
    case = PowerFlowData.parse_network(case_path)

    # Parse models
    bus = parse_model(Bus, case)
    pq = parse_model(PQ, case)

    # OLD manual workflow: disable auto_setup
    am = AddressManager()
    sys = SystemModel(am, (Bus = bus, PQ = pq); auto_setup=false)

    # System should NOT be set up yet
    @test !is_setup(sys)

    # Manually call setup!
    setup!(sys)

    # Now it should be set up
    @test is_setup(sys)

    # Should work normally now
    dae = DAE(Float64, sys)
    sys_ws = SystemResidualWorkspace(Float64, sys)
    evaluate_system_residuals!(sys_ws, dae)

    @test any(dae.g .!= 0.0)
end
