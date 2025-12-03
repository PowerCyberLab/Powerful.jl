"""
TDS TGOV1 test - test GENCLS + TGOV1 initialization and flat run
"""

using Powerful
using PowerFlowData
using DifferentialEquations

# Load case files
raw_file = joinpath(pkgdir(Powerful), "cases", "two_bus_transformer.raw")
dyr_file = joinpath(pkgdir(Powerful), "cases", "two_bus_tgov1.dyr")

println("="^60)
println("TGOV1 Integration Test - Two Bus Transformer Case")
println("="^60)

# Step 1: Solve power flow
println("\n[1] Solving power flow...")
case = PowerFlowData.parse_network(raw_file)
sys_pf = SystemModel(case)
pf_prob = ACPowerFlowProblem(sys_pf)
pf_sol = solve(pf_prob, NewtonRaphson())
println("    Power flow converged: $(pf_sol.converged)")

# Show generator power
if !isempty(sys_pf.models.Slack)
    slack = sys_pf.models.Slack
    println("    Slack bus $(slack.i[1]): P_gen will be computed from power balance")
end

# Step 2: Parse DYR data
println("\n[2] Parsing DYR file...")
dyr_data = parse_dyr(dyr_file, sys_pf)
println("    Models found: ", keys(dyr_data))

if haskey(dyr_data, :GENCLS)
    gencls = dyr_data.GENCLS
    println("    GENCLS instances: $(length(gencls))")
    for i in 1:length(gencls)
        println("      [$i] Bus=$(gencls.BUS[i]), ID=$(gencls.ID[i]), H=$(gencls.H[i]), D=$(gencls.D[i])")
    end
end

if haskey(dyr_data, :TGOV1)
    tgov1 = dyr_data.TGOV1
    println("    TGOV1 instances: $(length(tgov1))")
    for i in 1:length(tgov1)
        println("      [$i] Bus=$(tgov1.BUS[i]), ID=$(tgov1.ID[i]), R=$(tgov1.R[i]), T1=$(tgov1.T1[i]), T2=$(tgov1.T2[i]), T3=$(tgov1.T3[i])")
    end
end

# Step 3: Create TDS problem with both models
println("\n[3] Creating TDS problem...")
dynamic_models = (GENCLS=dyr_data.GENCLS,)  # Start with just GENCLS for now

# Add TGOV1 if present
if haskey(dyr_data, :TGOV1)
    dynamic_models = (GENCLS=dyr_data.GENCLS, TGOV1=dyr_data.TGOV1)
    println("    Including TGOV1 in dynamic models")
end

prob = TDSProblem(sys_pf, pf_sol, dynamic_models; tspan=(0.0, 1.0))
println("    TDS problem created successfully")

# Step 4: Check initialization
println("\n[4] Checking initialization...")
# The initial DAE vector is in prob
y0 = prob.u0
println("    Initial state vector length: $(length(y0))")

# Step 5: Validate initial residuals
println("\n[5] Validating initial residuals...")
try
    validation = validate_tds_initialization(prob)
    println("    Validation completed")
    print_validation_report(validation)
catch e
    println("    Validation error: $e")
end

# Step 6: Run flat simulation
println("\n[6] Running flat simulation (1.0 sec)...")
sol = solve(prob)
println("    Return code: $(sol.retcode)")
println("    Time steps: $(length(sol.t))")
println("    Final time: $(sol.t[end])")

# Step 7: Check state drift
println("\n[7] Checking state drift (should be ~0 for flat run)...")
y_init = sol.u[1]
y_final = sol.u[end]
max_drift = maximum(abs.(y_final - y_init))
println("    Max state drift: $max_drift")

if max_drift < 1e-6
    println("    ✓ Flat run successful - states stable")
else
    println("    ⚠ States drifted - check initialization")

    # Show which states drifted
    drift = abs.(y_final - y_init)
    drift_idx = findall(drift .> 1e-8)
    if !isempty(drift_idx)
        println("    Drifted states (index, drift):")
        for idx in drift_idx[1:min(10, length(drift_idx))]
            println("      [$idx]: $(drift[idx])")
        end
    end
end

println("\n" * "="^60)
println("Test completed")
println("="^60)
