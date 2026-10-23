#!/usr/bin/env julia
"""
TDS fault run for ACTIVSg2000 (singlegen variant) with voltage plotting.
"""

using Powerful
using PowerfulCases
using PowerfulCore: create_adcache, ForwardDiffBackend
using PowerfulTDS: TDSAdapter, ODEInterface, evaluate!
using PowerfulCore: AlgebRes, get_series
using Plots

println("=" ^ 60)
println("TDS Fault Run: ACTIVSg2000 (singlegen variant)")
println("=" ^ 60)

# Load case using collection/case format to avoid ambiguity
case = load("synthetic/ACTIVSg2000")
raw_path = file(case, :raw, variant="singlegen")
dyr_path = file(case, :dyr, variant="singlegen")
println("\n1. Loading case...")
println("   RAW: $(basename(raw_path))")
println("   DYR: $(basename(dyr_path))")

# Create system
println("\n2. Setting up system...")
sys_pf = SystemModel(parse_network(raw_path))

# Solve power flow
println("\n3. Solving power flow...")
pf_sol = solve(ACPowerFlowProblem(sys_pf), Newton())
println("   Converged: $(pf_sol.converged)")
println("   Iterations: $(pf_sol.iterations)")

# Parse dynamics
println("\n4. Parsing dynamic models...")
dyr_data = parse_dynamics(dyr_path, sys_pf)
n_gencls = haskey(dyr_data, :GENCLS) ? length(dyr_data.GENCLS) : 0
n_sexs = haskey(dyr_data, :SEXS) ? length(dyr_data.SEXS) : 0
n_tgov1 = haskey(dyr_data, :TGOV1) ? length(dyr_data.TGOV1) : 0
println("   GENCLS: $n_gencls, SEXS: $n_sexs, TGOV1: $n_tgov1")

# BusFault auto-creates the Fault model and expands to FaultApply/FaultClear
events = [BusFault(0.1, 0.2, 5001, "F1"; xf=0.0001, rf=0.0)]
println("\n5. Fault events:")
println("   t=0.1s: Apply fault at bus 5001")
println("   t=0.2s: Clear fault")

# 100% constant impedance loads for numerical stability during fault
config = SimConfig(
    solver = (
        TDS = TDSConfig(
            convert_pq_to_zip = true,
            p_zip = (1.0, 0.0, 0.0),
            q_zip = (1.0, 0.0, 0.0),
        ),
    ),
)

# Create TDS problem
println("\n6. Creating TDS problem...")
init = setup_tds(pf_sol, dyr_data; config=config)
prob = TDSProblem(init; tspan=(0.0, 0.5), events=events)

# ============================================================================
# VALIDATION: Check steady-state initialization
# ============================================================================
println("\n7. Validating steady-state initialization...")

adcache = create_adcache(prob.address_manager, prob.sys, prob.formulation, ForwardDiffBackend())
adapter = TDSAdapter(adcache, prob)
u0 = vcat(adapter.x_buf, adapter.y_buf)
du = similar(u0)
evaluate!(du, adapter, u0, 0.0, ODEInterface())

n_diff = adapter.n_diff
dx = du[1:n_diff]
g = du[n_diff+1:end]

max_dx = maximum(abs, dx)
max_g = maximum(abs, g)

# Get bus residuals
am = prob.address_manager
p_addrs = am.addresses[(AlgebRes(), :Bus, :p)]
q_addrs = am.addresses[(AlgebRes(), :Bus, :q)]

max_p = maximum(abs, g[p_addrs])
max_q = maximum(abs, g[q_addrs])

println("   Max |dx/dt|: $max_dx")
println("   Max |g|: $max_g")
println("   Max |P residual|: $max_p")
println("   Max |Q residual|: $max_q")

# Check tolerance
tol = 1e-6
if max_g > tol
    n_p_issues = sum(abs(g[addr]) > 1e-3 for addr in p_addrs)
    n_q_issues = sum(abs(g[addr]) > 1e-3 for addr in q_addrs)
    println("\n   ✗ NOT AT STEADY STATE!")
    println("   P issues (>1e-3): $n_p_issues")
    println("   Q issues (>1e-3): $n_q_issues")
    error("TDS initialization failed: max algebraic residual = $max_g > $tol")
end

println("   ✓ Steady-state validation PASSED")

# ============================================================================
# Run TDS simulation
# ============================================================================
println("\n8. Running TDS simulation (0 to 0.5 sec)...")

sol = solve(
    prob;
    initializealg = PowerFlowInit(; initial_guess = :pre_fault),
    dtmax = 0.01,
)

println("\n" * "=" ^ 60)
println("RESULTS")
println("=" ^ 60)
println("Return code: $(sol.retcode)")
println("Time steps: $(length(sol.t))")
println("Final time: $(sol.t[end]) s")

if sol.retcode == :Success
    println("\n✓ TDS FAULT RUN SUCCESS!")

    # ============================================================================
    # Plot bus voltages around faulted bus
    # ============================================================================
    println("\n9. Plotting bus voltages around faulted bus...")

    # Buses around the faulted bus (5001)
    buses_to_plot = [5001, 5002, 5003, 5004]

    # Create voltage plot
    p = plot(
        title = "Bus Voltages Around Faulted Bus 5001",
        xlabel = "Time (s)",
        ylabel = "Voltage (pu)",
        legend = :bottomright,
        size = (900, 500),
    )

    for bus in buses_to_plot
        t, v = get_series(sol, :Bus, :v, bus)
        plot!(p, t, v, label = "Bus $bus", linewidth = 1.5)
    end

    # Add fault timing markers
    vline!(p, [0.1, 0.2], linestyle = :dash, color = :red, label = "Fault", linewidth = 1)

    display(p)
else
    println("\n✗ TDS FAULT RUN FAILED: $(sol.retcode)")
end
