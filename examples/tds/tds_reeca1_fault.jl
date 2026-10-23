#!/usr/bin/env julia
#= REECA1 Fault Test - Mixed synchronous and renewable generators with electrical controls =#

using Powerful
using PowerfulCases
using PowerfulCore.DiffEq: FBDF, KLUFactorization
using PowerfulCore: get_series
using PowerfulTDS: validate_tds_initialization
using Plots

"""
Test REECA1 with specified control mode flags.

REECA1 Control Modes:
- VFLAG=0, QFLAG=1: V control mode (Vsel = Vref1 + PFsel -> PIV -> Iqcmd)
- VFLAG=1, QFLAG=1: Q control mode (Qerr -> PIQ -> Vsel -> PIV -> Iqcmd)
- VFLAG=x, QFLAG=0: Constant Q/PF mode (s4_y = PFsel/vp directly)
- PFFLAG=0: Use Qref for PFsel
- PFFLAG=1: Use Qcpf (Q from power factor) for PFsel
"""
function test_reeca1_mode(; VFLAG::Float64, QFLAG::Float64, PFFLAG::Float64=0.0,
                          tspan=(0.0, 3.0), plot_results::Bool=false)
    mode_name = if QFLAG == 0
        "Constant Q/PF (QFLAG=0)"
    elseif VFLAG == 0
        "V control (VFLAG=0, QFLAG=1)"
    else
        "Q control (VFLAG=1, QFLAG=1)"
    end
    println("\n" * "=" ^ 60)
    println("Testing REECA1 Mode: $mode_name, PFFLAG=$PFFLAG")
    println("=" ^ 60)

    # Load system
    case = load("ieee14")
    sys = SystemModel(parse_network(case.raw))

    # Solve power flow
    prob = ACPowerFlowProblem(sys)
    pf_sol = solve(prob, Newton())
    if !pf_sol.converged
        println("Power flow FAILED")
        return (success=false, mode=mode_name, reason="Power flow failed")
    end

    # Parse dynamic models
    dyr_file = file(case, :dyr, variant="reeca1")
    dyr_data = parse_dynamics(dyr_file, sys)

    if !hasproperty(dyr_data, :REECA1) || length(dyr_data.REECA1) == 0
        println("No REECA1 models found")
        return (success=false, mode=mode_name, reason="No REECA1 models")
    end

    # Modify REECA1 flags for this test
    for i in 1:length(dyr_data.REECA1)
        dyr_data.REECA1.VFLAG[i] = VFLAG
        dyr_data.REECA1.QFLAG[i] = QFLAG
        dyr_data.REECA1.PFFLAG[i] = PFFLAG
    end

    println("REECA1 flags: VFLAG=$VFLAG, QFLAG=$QFLAG, PFFLAG=$PFFLAG")

    # BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
    events = [BusFault(0.5, 0.6, 4, "F1"; xf=0.0001, rf=0.0)]
    config = SimConfig(solver = (TDS = TDSConfig(
        convert_pq_to_zip=true, p_zip=(1.0,0.0,0.0), q_zip=(1.0,0.0,0.0)),))

    # Create TDS problem
    init = setup_tds(pf_sol, dyr_data; config=config)
    tds_prob = TDSProblem(init; tspan=tspan, events=events)

    # Validate initialization
    validation = validate_tds_initialization(tds_prob; atol=1e-6)
    println("Initialization: $(validation.is_valid ? "PASSED" : "FAILED")")
    println("  Max diff residual: $(validation.max_diff_residual)")
    println("  Max alg residual: $(validation.max_alg_residual)")

    if !validation.is_valid
        for (addr, val, model, eq) in validation.non_zero_diff
            println("    [$addr] $model.$eq = $val (diff)")
        end
        for (addr, val, model, eq) in validation.non_zero_alg
            println("    [$addr] $model.$eq = $val (alg)")
        end
        return (success=false, mode=mode_name, reason="Initialization failed",
                max_diff=validation.max_diff_residual, max_alg=validation.max_alg_residual)
    end

    # Solve TDS
    println("Solving TDS...")
    sol = solve(tds_prob, FBDF(linsolve=KLUFactorization());
                initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.1)

    if sol.retcode != :Success
        println("TDS FAILED: $(sol.retcode) at t=$(sol.t[end])")
        return (success=false, mode=mode_name, reason="TDS failed: $(sol.retcode)",
                final_t=sol.t[end])
    end

    println("TDS SUCCESS: $(length(sol) - 1) accepted steps")

    # Collect data for plotting/analysis
    if plot_results
        # Get REECA1 variables using get_series
        vars_to_check = [:s0_y, :S1_y, :PIQ_xi, :PIV_xi, :s4_y, :pfilt_y, :s5_y, :Ipcmd_out, :Iqcmd_out]
        println("\nREECA1 variable ranges:")
        for v in vars_to_check
            data_v = get_series(sol, :REECA1, v)
            for (var, vals) in pairs(data_v)
                var == :t && continue
                println("  $var: min=$(round(minimum(vals), digits=4)), max=$(round(maximum(vals), digits=4))")
            end
        end
    end

    return (success=true, mode=mode_name, steps=length(sol) - 1, tds_prob=tds_prob, sol=sol)
end

# ============================================================================
# Main test: Run all REECA1 control mode combinations
# ============================================================================
println("REECA1 Control Mode Tests")
println("=" ^ 60)

results = []

# Test all mode combinations
# Mode 1: Q control (VFLAG=1, QFLAG=1) - original DYR setting
push!(results, test_reeca1_mode(VFLAG=1.0, QFLAG=1.0))

# Mode 2: V control (VFLAG=0, QFLAG=1) - had sign bug in Qref0 init
push!(results, test_reeca1_mode(VFLAG=0.0, QFLAG=1.0))

# Mode 3: Constant Q mode (QFLAG=0) - bypasses PI controllers
push!(results, test_reeca1_mode(VFLAG=0.0, QFLAG=0.0))

# Mode 4: PF control mode (PFFLAG=1)
push!(results, test_reeca1_mode(VFLAG=1.0, QFLAG=1.0, PFFLAG=1.0))

# Summary
println("\n" * "=" ^ 60)
println("SUMMARY")
println("=" ^ 60)
n_passed = count(r -> r.success, results)
n_total = length(results)
println("Passed: $n_passed / $n_total")
for r in results
    status = r.success ? "PASS" : "FAIL"
    info = r.success ? "$(r.steps) steps" : r.reason
    println("  [$status] $(r.mode): $info")
end

# For detailed plotting, run the last successful test with plots
last_success = filter(r -> r.success, results)
if !isempty(last_success) && length(ARGS) > 0 && ARGS[1] == "--plot"
    r = last_success[end]
    sol = r.sol

    data_gen = get_series(sol, :GENROU, :omega)
    data_Ip = get_series(sol, :REGCA1, :Ip)
    data_Iq = get_series(sol, :REGCA1, :Iq)
    data_Ipcmd = get_series(sol, :REECA1, :Ipcmd_out)
    data_Iqcmd = get_series(sol, :REECA1, :Iqcmd_out)

    p1 = plot(title="Generator Speeds", xlabel="Time (s)", ylabel="ω (pu)", legend=:topright)
    for (var, vals) in pairs(data_gen)
        var == :t && continue
        plot!(p1, data_gen.t, vals, label="Gen $(split(string(var), "_")[3])", linewidth=1.5)
    end
    vline!(p1, [0.5, 0.6], linestyle=:dash, color=:red, label="Fault", linewidth=1)

    p2 = plot(title="REGCA1 Currents", xlabel="Time (s)", ylabel="I (pu)", legend=:topright)
    for (var, vals) in pairs(data_Ip)
        var == :t && continue
        plot!(p2, data_Ip.t, vals, label="Ip Bus $(split(string(var), "_")[3])", linewidth=1.5)
    end
    for (var, vals) in pairs(data_Iq)
        var == :t && continue
        plot!(p2, data_Iq.t, vals, label="Iq Bus $(split(string(var), "_")[3])", linewidth=1.5)
    end
    vline!(p2, [0.5, 0.6], linestyle=:dash, color=:red, label="Fault", linewidth=1)

    p3 = plot(title="REECA1 Current Commands ($(r.mode))", xlabel="Time (s)", ylabel="I (pu)", legend=:topright)
    for (var, vals) in pairs(data_Ipcmd)
        var == :t && continue
        plot!(p3, data_Ipcmd.t, vals, label="Ipcmd Bus $(split(string(var), "_")[4])", linewidth=1.5)
    end
    for (var, vals) in pairs(data_Iqcmd)
        var == :t && continue
        plot!(p3, data_Iqcmd.t, vals, label="Iqcmd Bus $(split(string(var), "_")[4])", linewidth=1.5)
    end
    vline!(p3, [0.5, 0.6], linestyle=:dash, color=:red, label="Fault", linewidth=1)
    display(plot(p1, p2, p3, layout=(3,1), size=(900, 800)))
end

# Exit with error if any test failed
if n_passed < n_total
    exit(1)
end
