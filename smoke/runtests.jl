# Smoke tests for the public release: completion and coarse sanity bounds only.
# No precise reference values, tolerances on trajectories, event timing or edge
# cases, and nothing that specifies the automaton runtime's behavior.
using Test
using Powerful
using PowerfulCases
using PowerfulCA
using PowerfulYbus

@testset "Smoke" begin

    @testset "Power flow" begin
        sys = SystemModel(parse_network(load("ieee14").raw))
        sol = solve(ACPowerFlowProblem(sys), Newton())

        @test sol.converged
        @test sol.iterations <= 10

        v = get_all_vars(sol, :Bus, :v)
        @test all(0.9 .<= v .<= 1.1)
    end

    @testset "Phasor TDS" begin
        case = load("ieee14")
        sys = SystemModel(parse_network(case.raw))
        pf_sol = solve(ACPowerFlowProblem(sys), Newton())
        config = SimConfig(solver = (TDS = TDSConfig(
            convert_pq_to_zip = true, p_zip = (1.0, 0.0, 0.0), q_zip = (1.0, 0.0, 0.0)),))
        dyr_data = parse_dynamics(case.dyr, sys)
        init = setup_tds(pf_sol, dyr_data; config = config)
        prob = TDSProblem(init; tspan = (0.0, 2.0), events = [BusFault(0.5, 0.6, 7)])
        sol = solve(prob)

        @test sol.retcode == :Success

        v = get_series(sol, :Bus, :v)
        @test all(v[k][end] > 0.9 for k in keys(v) if k != :t)

        omega = get_series(sol, :GENCLS, :omega)
        @test all(abs(omega[k][end] - 1.0) <= 0.01 for k in keys(omega) if k != :t)
    end

    @testset "Contingency analysis" begin
        sys = SystemModel(parse_network(load("ieee14").raw))
        contingencies = [
            Contingency("LINE_1_5", [BranchOutage(1, 5, "1")]),
            Contingency("LINE_2_3", [BranchOutage(2, 3, "1")]),
            Contingency("GEN_3", [GeneratorOutage(3, "1")]),
        ]

        result = solve_accc(sys, contingencies)

        @test length(result.contingencies) == length(contingencies)
    end

    @testset "Y-bus solver" begin
        case = load("ieee14")
        sys = SystemModel(parse_network(case.raw))
        pf_sol = solve(ACPowerFlowProblem(sys), Newton())
        config = SimConfig(solver = (TDS = TDSConfig(
            convert_pq_to_zip = true, p_zip = (1.0, 0.0, 0.0), q_zip = (1.0, 0.0, 0.0)),))
        init = setup_tds(pf_sol, parse_dynamics(case.dyr, sys); config = config)
        prob = TDSProblem(init; tspan = (0.0, 2.0), events = [BusFault(0.5, 0.6, 7)])
        sol = solve(YBusTDSProblem(prob))

        @test sol.retcode == :Success
    end

end
