#!/usr/bin/env julia
#= REECA1 Current Command Plot - Compare TDS vs Ybus formulations

This script compares REECA1 (Renewable Electrical Controller) outputs between
the PowerfulTDS (DifferentialEquations.jl) and PowerfulYbus (trapezoidal) solvers.

KEY INSIGHT: Address semantics
- StateVar addresses are direct DAE positions (1 to n_states)
- AlgebVar addresses are relative to algebraic space, need offset: dae_pos = n_states + alg_addr
- Solution vector layout: [states; algebraics]
=#

using Powerful
using PowerfulCases
using PowerfulYbus
using PowerfulYbus: solve_trapezoidal
using PowerfulCore.DiffEq: FBDF, KLUFactorization
using PowerfulCore: get_addresses, StateVar, AlgebVar
using Plots

# ============================================================================
# Helper: Convert address to DAE position
# ============================================================================
"""
    to_dae_position(am, model::Symbol, var::Symbol) -> Vector{Int}

Convert addresses to DAE vector positions.
State variables: direct position
Algebraic variables: n_states + algebraic_addr
"""
function to_dae_position(am, model::Symbol, var::Symbol)
    addrs = get_addresses(am, model, var)

    # Check if this is a state or algebraic variable
    # by looking up in the address dictionaries
    key_state = (StateVar(), model, var)
    key_alg = (AlgebVar(), model, var)

    if haskey(am.addresses, key_state)
        # State variable - address is direct DAE position
        return Int.(addrs)
    elseif haskey(am.addresses, key_alg)
        # Algebraic variable - need offset
        n_states = length(am.reverse_map.state_vars)
        return [n_states + Int(a) for a in addrs]
    else
        error("Variable $model.$var not found in address manager")
    end
end

# ============================================================================
# 1. Setup
# ============================================================================
case = load("ieee14")
tspan = (0.0, 3.0)
fault_bus = 4
fault_xf = 0.01

config = SimConfig(solver = (TDS = TDSConfig(
    convert_pq_to_zip=true, p_zip=(1.0, 0.0, 0.0), q_zip=(1.0, 0.0, 0.0)),))

# ============================================================================
# 2. Solve with TDS (DifferentialEquations.jl)
# ============================================================================
println("[1] Solving with PowerfulTDS...")

sys_tds = SystemModel(parse_network(case.raw))
pf_sol_tds = solve(ACPowerFlowProblem(sys_tds), Newton())
@assert pf_sol_tds.converged "Power flow failed"

dyr_file = file(case, :dyr, variant="reeca1")
dyr_data_tds = parse_dynamics(dyr_file, sys_tds)

# Set REECA1 to Q control mode
for i in 1:length(dyr_data_tds.REECA1)
    dyr_data_tds.REECA1.VFLAG[i] = 1.0
    dyr_data_tds.REECA1.QFLAG[i] = 1.0
    dyr_data_tds.REECA1.PFFLAG[i] = 0.0
end

events_tds = [BusFault(0.5, 0.6, fault_bus, "F1"; xf=fault_xf, rf=0.0)]
init_tds = setup_tds(pf_sol_tds, dyr_data_tds; config=config)
prob_tds = TDSProblem(init_tds; tspan=tspan, events=events_tds)

sol_tds = solve(prob_tds, FBDF(linsolve=KLUFactorization());
                initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.1)
println("    TDS: $(sol_tds.retcode), $(sol_tds.stats.naccept) steps")

# ============================================================================
# 3. Solve with Ybus (Trapezoidal)
# ============================================================================
println("[2] Solving with PowerfulYbus...")

sys_ybus = SystemModel(parse_network(case.raw))
pf_sol_ybus = solve(ACPowerFlowProblem(sys_ybus), Newton())
dyr_data_ybus = parse_dynamics(dyr_file, sys_ybus)

# Set same REECA1 flags
for i in 1:length(dyr_data_ybus.REECA1)
    dyr_data_ybus.REECA1.VFLAG[i] = 1.0
    dyr_data_ybus.REECA1.QFLAG[i] = 1.0
    dyr_data_ybus.REECA1.PFFLAG[i] = 0.0
end

events_ybus = [BusFault(0.5, 0.6, fault_bus, "F1"; xf=fault_xf, rf=0.0)]
init_ybus = setup_tds(pf_sol_ybus, dyr_data_ybus; config=config)
prob_ybus_base = TDSProblem(init_ybus; tspan=tspan)
prob_ybus = YBusTDSProblem(prob_ybus_base; events=events_ybus)

sol_ybus = solve_trapezoidal(prob_ybus; dt=0.001, abstol=1e-6, verbose=false)
println("    Ybus: $(sol_ybus.retcode)")

# ============================================================================
# 4. Extract data using correct DAE positions
# ============================================================================
println("[3] Extracting data...")

t_tds = sol_tds.t
t_ybus = sol_ybus.t

am = prob_tds.address_manager

# Get DAE positions (handles state vs algebraic automatically)
reeca_Ipcmd_pos = to_dae_position(am, :REECA1, :Ipcmd_out)[1]
reeca_Iqcmd_pos = to_dae_position(am, :REECA1, :Iqcmd_out)[1]
regca_Ip_pos = to_dae_position(am, :REGCA1, :Ip)[1]
regca_Iq_pos = to_dae_position(am, :REGCA1, :Iq_lag)[1]

println("    DAE positions:")
println("      REECA1 Ipcmd_out: $reeca_Ipcmd_pos")
println("      REECA1 Iqcmd_out: $reeca_Iqcmd_pos")
println("      REGCA1 Ip: $regca_Ip_pos")
println("      REGCA1 Iq_lag: $regca_Iq_pos")

# Extract time series
tds_matrix = hcat(sol_tds.u...)  # (n_vars, n_times)

Ipcmd_tds = tds_matrix[reeca_Ipcmd_pos, :]
Iqcmd_tds = tds_matrix[reeca_Iqcmd_pos, :]
Ip_tds = tds_matrix[regca_Ip_pos, :]
Iq_tds = tds_matrix[regca_Iq_pos, :]

Ipcmd_ybus = sol_ybus.u[reeca_Ipcmd_pos, :]
Iqcmd_ybus = sol_ybus.u[reeca_Iqcmd_pos, :]
Ip_ybus = sol_ybus.u[regca_Ip_pos, :]
Iq_ybus = sol_ybus.u[regca_Iq_pos, :]

# Verify values match
println("\n    Initial values (t=0):")
println("      REECA1 Ipcmd_out: TDS=$(round(Ipcmd_tds[1], digits=4)), Ybus=$(round(Ipcmd_ybus[1], digits=4))")
println("      REGCA1 Ip:        TDS=$(round(Ip_tds[1], digits=4)), Ybus=$(round(Ip_ybus[1], digits=4))")

# ============================================================================
# 5. Plot comparison
# ============================================================================
println("[4] Plotting...")

ybus_step = 5  # Downsample Ybus markers for cleaner plot

function plot_comparison!(sp, t_tds, y_tds, t_ybus, y_ybus, title_str; color=:blue)
    scatter!(sp, t_ybus[1:ybus_step:end], y_ybus[1:ybus_step:end],
             marker=:x, ms=2, alpha=0.5, label="Ybus", color=color)
    plot!(sp, t_tds, y_tds, lw=2, alpha=0.8, label="TDS", color=color)
    vline!(sp, [0.5, 0.6], label="", ls=:dot, color=:gray, lw=1)
    title!(sp, title_str)
end

p = plot(layout=(2, 2), size=(1000, 700), dpi=150,
         left_margin=5Plots.mm, bottom_margin=5Plots.mm)

plot_comparison!(p[1], t_tds, Ipcmd_tds, t_ybus, Ipcmd_ybus,
                 "REECA1 Ipcmd_out (Bus 2)"; color=:blue)
xlabel!(p[1], "Time (s)")
ylabel!(p[1], "Current (pu)")

plot_comparison!(p[2], t_tds, Iqcmd_tds, t_ybus, Iqcmd_ybus,
                 "REECA1 Iqcmd_out (Bus 2)"; color=:red)
xlabel!(p[2], "Time (s)")
ylabel!(p[2], "Current (pu)")

plot_comparison!(p[3], t_tds, Ip_tds, t_ybus, Ip_ybus,
                 "REGCA1 Ip (Bus 2)"; color=:blue)
xlabel!(p[3], "Time (s)")
ylabel!(p[3], "Current (pu)")

plot_comparison!(p[4], t_tds, Iq_tds, t_ybus, Iq_ybus,
                 "REGCA1 Iq_lag (Bus 2)"; color=:red)
xlabel!(p[4], "Time (s)")
ylabel!(p[4], "Current (pu)")

display(p)
savefig(p, joinpath(@__DIR__, "reeca1_tds_vs_ybus.png"))
println("\nPlot saved to reeca1_tds_vs_ybus.png")
