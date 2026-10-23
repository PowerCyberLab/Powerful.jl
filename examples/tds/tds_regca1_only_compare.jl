#!/usr/bin/env julia
#= REGCA1 Only Comparison - Compare TDS vs Ybus without REECA1 controller =#

using Powerful
using PowerfulCases
using PowerfulYbus
using PowerfulYbus: solve_trapezoidal
using PowerfulCore.DiffEq: FBDF, KLUFactorization, ReturnCode
using PowerfulCore: get_addresses, get_series
using Plots

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
# 2. Solve with TDS
# ============================================================================
println("[1] Solving with PowerfulTDS...")

sys_tds = SystemModel(parse_network(case.raw))
pf_sol_tds = solve(ACPowerFlowProblem(sys_tds), Newton())
@assert pf_sol_tds.converged "Power flow failed"

# Use reeca1 variant (REGCA1 with REECA1 controller)
dyr_file = file(case, :dyr, variant="reeca1")
dyr_data_tds = parse_dynamics(dyr_file, sys_tds)

println("    REGCA1 buses: $(dyr_data_tds.REGCA1.BUS)")

events_tds = [BusFault(0.5, 0.6, fault_bus, "F1"; xf=fault_xf, rf=0.0)]
init_tds = setup_tds(pf_sol_tds, dyr_data_tds; config=config)
prob_tds = TDSProblem(init_tds; tspan=tspan, events=events_tds)

sol_tds = solve(prob_tds, FBDF(linsolve=KLUFactorization());
                initializealg=PowerFlowInit(; initial_guess=:pre_fault), dtmax=0.1)
println("    TDS: $(sol_tds.retcode), $(sol_tds.stats.naccept) steps")

# ============================================================================
# 3. Solve with Ybus
# ============================================================================
println("[2] Solving with PowerfulYbus...")

sys_ybus = SystemModel(parse_network(case.raw))
pf_sol_ybus = solve(ACPowerFlowProblem(sys_ybus), Newton())
dyr_data_ybus = parse_dynamics(dyr_file, sys_ybus)

events_ybus = [BusFault(0.5, 0.6, fault_bus, "F1"; xf=fault_xf, rf=0.0)]
init_ybus = setup_tds(pf_sol_ybus, dyr_data_ybus; config=config)
prob_ybus_base = TDSProblem(init_ybus; tspan=tspan)
prob_ybus = YBusTDSProblem(prob_ybus_base; events=events_ybus)

sol_ybus = solve_trapezoidal(prob_ybus; dt=0.001, abstol=1e-6, verbose=false)
println("    Ybus: $(sol_ybus.retcode)")

# ============================================================================
# 4. Extract data
# ============================================================================
println("[3] Extracting data...")

t_tds = sol_tds.t
t_ybus = sol_ybus.t

# TDS data - REGCA1 and GENROU
data_tds_regca = merge((get_series(sol_tds, :REGCA1, var) for var in (:Ip, :Iq_lag, :Ipcmd, :Iqcmd))...)
data_tds_genrou = get_series(sol_tds, :GENROU, :omega)

# Ybus addresses
am_ybus = prob_ybus_base.address_manager
regca_Ip_addrs = get_addresses(am_ybus, :REGCA1, :Ip)
regca_Iq_addrs = get_addresses(am_ybus, :REGCA1, :Iq_lag)
regca_Ipcmd_addrs = get_addresses(am_ybus, :REGCA1, :Ipcmd)
regca_Iqcmd_addrs = get_addresses(am_ybus, :REGCA1, :Iqcmd)
genrou_omega_addrs = get_addresses(am_ybus, :GENROU, :omega)

# Print variable names
println("\nTDS REGCA1 variables:")
for var in propertynames(data_tds_regca)
    var == :t && continue
    println("  $var")
end

println("\nYbus addresses:")
println("  REGCA1 Ip: $regca_Ip_addrs")
println("  REGCA1 Iq_lag: $regca_Iq_addrs")
println("  REGCA1 Ipcmd: $regca_Ipcmd_addrs")
println("  REGCA1 Iqcmd: $regca_Iqcmd_addrs")

# Build bus → model index mapping
regca_bus_to_idx = Dict{Int, Int}()
for (idx, bus) in enumerate(dyr_data_ybus.REGCA1.BUS)
    regca_bus_to_idx[bus] = idx
end
genrou_bus_to_idx = Dict{Int, Int}()
for (idx, bus) in enumerate(dyr_data_ybus.GENROU.BUS)
    genrou_bus_to_idx[bus] = idx
end
println("\nREGCA1 bus → idx: $regca_bus_to_idx")
println("GENROU bus → idx: $genrou_bus_to_idx")

# ============================================================================
# 5. Plot comparison - Bus 6 and Bus 8
# ============================================================================
println("[4] Plotting...")

ybus_step = 5

function plot_cmp!(sp, t_tds, y_tds, t_ybus, y_ybus, title_str; kwargs...)
    scatter!(sp, t_ybus[1:ybus_step:end], y_ybus[1:ybus_step:end],
             marker=:x, ms=2, alpha=0.5, label="Ybus"; kwargs...)
    plot!(sp, t_tds, y_tds, lw=2, alpha=0.8, label="TDS"; kwargs...)
    vline!(sp, [0.5, 0.6], label="", ls=:dot, color=:gray, lw=1)
    title!(sp, title_str)
end

p = plot(layout=(3, 2), size=(1000, 900), dpi=150,
         left_margin=5Plots.mm, bottom_margin=5Plots.mm)

# Get first two buses from the model
bus1 = Int(dyr_data_ybus.REGCA1.BUS[1])
bus2 = Int(dyr_data_ybus.REGCA1.BUS[2])

# Row 1: REGCA1 Ip for both buses
idx_1 = regca_bus_to_idx[bus1]
Ip_1_tds = getproperty(data_tds_regca, Symbol("REGCA1_Ip_$(bus1)_1"))
Ip_1_ybus = sol_ybus.u[Int(regca_Ip_addrs[idx_1]), :]
plot_cmp!(p[1], t_tds, Ip_1_tds, t_ybus, Ip_1_ybus, "REGCA1 Ip (Bus $bus1)"; color=:blue)
xlabel!(p[1], "Time (s)")
ylabel!(p[1], "Current (pu)")

idx_2 = regca_bus_to_idx[bus2]
Ip_2_tds = getproperty(data_tds_regca, Symbol("REGCA1_Ip_$(bus2)_1"))
Ip_2_ybus = sol_ybus.u[Int(regca_Ip_addrs[idx_2]), :]
plot_cmp!(p[2], t_tds, Ip_2_tds, t_ybus, Ip_2_ybus, "REGCA1 Ip (Bus $bus2)"; color=:blue)
xlabel!(p[2], "Time (s)")
ylabel!(p[2], "Current (pu)")

# Row 2: REGCA1 Iq_lag for both buses
Iq_1_tds = getproperty(data_tds_regca, Symbol("REGCA1_Iq_lag_$(bus1)_1"))
Iq_1_ybus = sol_ybus.u[Int(regca_Iq_addrs[idx_1]), :]
plot_cmp!(p[3], t_tds, Iq_1_tds, t_ybus, Iq_1_ybus, "REGCA1 Iq_lag (Bus $bus1)"; color=:red)
xlabel!(p[3], "Time (s)")
ylabel!(p[3], "Current (pu)")

Iq_2_tds = getproperty(data_tds_regca, Symbol("REGCA1_Iq_lag_$(bus2)_1"))
Iq_2_ybus = sol_ybus.u[Int(regca_Iq_addrs[idx_2]), :]
plot_cmp!(p[4], t_tds, Iq_2_tds, t_ybus, Iq_2_ybus, "REGCA1 Iq_lag (Bus $bus2)"; color=:red)
xlabel!(p[4], "Time (s)")
ylabel!(p[4], "Current (pu)")

# Row 3: GENROU omega (pick first two generators)
genrou_buses = collect(keys(genrou_bus_to_idx))
sort!(genrou_buses)
gen_bus1 = genrou_buses[1]
gen_bus2 = genrou_buses[2]

gen_idx_1 = genrou_bus_to_idx[gen_bus1]
omega_1_tds = getproperty(data_tds_genrou, Symbol("GENROU_omega_$(gen_bus1)_1"))
omega_1_ybus = sol_ybus.u[Int(genrou_omega_addrs[gen_idx_1]), :]
plot_cmp!(p[5], t_tds, omega_1_tds, t_ybus, omega_1_ybus, "GENROU omega (Bus $gen_bus1)"; color=:green)
xlabel!(p[5], "Time (s)")
ylabel!(p[5], "Speed (pu)")

gen_idx_2 = genrou_bus_to_idx[gen_bus2]
omega_2_tds = getproperty(data_tds_genrou, Symbol("GENROU_omega_$(gen_bus2)_1"))
omega_2_ybus = sol_ybus.u[Int(genrou_omega_addrs[gen_idx_2]), :]
plot_cmp!(p[6], t_tds, omega_2_tds, t_ybus, omega_2_ybus, "GENROU omega (Bus $gen_bus2)"; color=:green)
xlabel!(p[6], "Time (s)")
ylabel!(p[6], "Speed (pu)")

display(p)
savefig(p, joinpath(@__DIR__, "regca1_only_tds_vs_ybus.png"))
println("\nPlot saved to regca1_only_tds_vs_ybus.png")
