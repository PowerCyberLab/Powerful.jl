"""
REECA1 Bode Plot Example: IBR Active Power Control Response

This example demonstrates frequency response analysis for a grid-following
inverter-based resource (IBR) using the REECA1 electrical control model.

Transfer Functions Analyzed:
- Pref0 → Ip: Active power reference to active current output
- Pref0 → pfilt_y: Active power reference to filtered power (internal dynamics)

These transfer functions reveal the IBR's control bandwidth and phase margins,
which are critical for grid stability assessment and controller tuning.
"""

using Powerful
using PowerfulCases
using PowerfulTDS: linearize
using ControlSystemsBase: bode
using Plots

# =============================================================================
# Load System with IBR (IEEE 14-bus with REECA1/REGCA1)
# =============================================================================

println("Loading IEEE 14-bus system with REECA1/REGCA1 IBR models...")

case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# Solve power flow
prob = ACPowerFlowProblem(sys)
pf_sol = solve(prob, Newton())
if !pf_sol.converged
    error("Power flow failed to converge")
end
println("Power flow converged.")

# Parse dynamic models (REECA1 variant)
dyr_file = file(case, :dyr, variant="reeca1")
dyr_data = parse_dynamics(dyr_file, sys)

println("Dynamic models loaded:")
for name in propertynames(dyr_data)
    data = getproperty(dyr_data, name)
    if length(data) > 0
        println("  $name: $(length(data)) instance(s)")
    end
end

# Setup TDS - no events needed for linearization (steady-state analysis)
init = setup_tds(pf_sol, dyr_data)
tds_prob = TDSProblem(init; tspan=(0.0, 1.0))

# =============================================================================
# Linearization 1: Pref0 → Ip (Power Reference to Active Current)
# =============================================================================

println("\n" * "="^60)
println("Linearizing: REECA1.Pref0 → REGCA1.Ip")
println("="^60)

# For simplicity, analyze first IBR instance
# 1:1 means "first instance only" (single input/output)
inputs_p = [(:REECA1, :Pref0, 1:1)]    # Active power reference setpoint
outputs_p = [(:REGCA1, :Ip, 1:1)]       # Active current output

# Linearize at operating point
sys_p = linearize(tds_prob, inputs_p, outputs_p; verbose=true)

println("\nState-space model: A=$(size(sys_p.A)), B=$(size(sys_p.B)), C=$(size(sys_p.C)), D=$(size(sys_p.D))")

# =============================================================================
# Linearization 2: Pref0 → pfilt_y (Internal Power Filter Dynamics)
# =============================================================================

println("\n" * "="^60)
println("Linearizing: REECA1.Pref0 → REECA1.pfilt_y")
println("="^60)

inputs_pfilt = [(:REECA1, :Pref0, 1:1)]
outputs_pfilt = [(:REECA1, :pfilt_y, 1:1)]

sys_pfilt = linearize(tds_prob, inputs_pfilt, outputs_pfilt; verbose=true)

println("\nState-space model: A=$(size(sys_pfilt.A)), B=$(size(sys_pfilt.B)), C=$(size(sys_pfilt.C)), D=$(size(sys_pfilt.D))")

# =============================================================================
# Compute Frequency Response using ControlSystemsBase
# =============================================================================

frequencies = 10.0 .^ range(-2, 2, length=100)

# Pref0 → Ip
mag_p, phase_p, w = bode(sys_p, frequencies)
mag_p_db = 20 .* log10.(mag_p[1, 1, :])
phase_p_deg = phase_p[1, 1, :]

# Pref0 → pfilt_y
mag_pfilt, phase_pfilt, _ = bode(sys_pfilt, frequencies)
mag_pfilt_db = 20 .* log10.(mag_pfilt[1, 1, :])
phase_pfilt_deg = phase_pfilt[1, 1, :]

# =============================================================================
# Plot Results
# =============================================================================

# Plot 1: Pref0 → Ip
p1_mag = plot(w, mag_p_db,
    ylabel="Magnitude (dB)",
    xscale=:log10,
    legend=false,
    title="REECA1: Pref0 → Ip",
    linewidth=2,
    grid=true,
    minorgrid=true,
    xformatter=x -> ""
)

p1_phase = plot(w, phase_p_deg,
    xlabel="Frequency (rad/s)",
    ylabel="Phase (deg)",
    xscale=:log10,
    legend=false,
    linewidth=2,
    grid=true,
    minorgrid=true,
    ylims=(-360, 0)
)

fig1 = plot(p1_mag, p1_phase, layout=(2, 1))

# Plot 2: Pref0 → pfilt_y
p2_mag = plot(w, mag_pfilt_db,
    ylabel="Magnitude (dB)",
    xscale=:log10,
    legend=false,
    title="REECA1: Pref0 → pfilt_y",
    linewidth=2,
    grid=true,
    minorgrid=true,
    color=:orange,
    xformatter=x -> ""
)

p2_phase = plot(w, phase_pfilt_deg,
    xlabel="Frequency (rad/s)",
    ylabel="Phase (deg)",
    xscale=:log10,
    legend=false,
    linewidth=2,
    grid=true,
    minorgrid=true,
    color=:orange,
    ylims=(-180, 0)
)

fig2 = plot(p2_mag, p2_phase, layout=(2, 1))

# Combine side by side
fig = plot(fig1, fig2, layout=(1, 2), size=(600, 350), dpi=200)

# Save to file
savefig(fig, joinpath(@__DIR__, "bode_reeca1_plot.png"))
println("\nSaved REECA1 Bode plot to examples/bode/bode_reeca1_plot.png")

# =============================================================================
# Frequency Response Metrics
# =============================================================================

println("\n" * "="^60)
println("Frequency Response Metrics")
println("="^60)

"""
    find_crossover_frequency(mag, freq)

Find the 0 dB crossover frequency (where magnitude = 1).
For monotonically decreasing systems, this is more meaningful than -3dB bandwidth.
"""
function find_crossover_frequency(mag, freq)
    mag_db = 20 .* log10.(mag)
    # Find where magnitude crosses 0 dB (interpolate between points)
    for i in 2:length(mag_db)
        if mag_db[i-1] >= 0 && mag_db[i] < 0
            # Linear interpolation to find exact crossing
            t = -mag_db[i-1] / (mag_db[i] - mag_db[i-1])
            return freq[i-1] + t * (freq[i] - freq[i-1]), i
        end
    end
    return nothing, nothing
end

# Pref0 → Ip metrics
dc_gain_db_p = mag_p_db[1]
crossover_p, crossover_idx_p = find_crossover_frequency(mag_p[1, 1, :], w)

println("Pref0 → Ip (Power Reference to Active Current):")
println("  DC gain: $(round(dc_gain_db_p, digits=1)) dB")
if crossover_p !== nothing
    println("  Crossover frequency (0 dB): $(round(crossover_p, digits=2)) rad/s ($(round(crossover_p/(2π), digits=2)) Hz)")
    phase_at_crossover = phase_p_deg[crossover_idx_p]
    phase_margin = 180 + phase_at_crossover
    println("  Phase at crossover: $(round(phase_at_crossover, digits=1))°")
    println("  Phase margin: $(round(phase_margin, digits=1))°")
else
    println("  Crossover frequency: Not found in range (gain stays above 0 dB)")
end

# Pref0 → pfilt_y metrics
dc_gain_db_pfilt = mag_pfilt_db[1]
crossover_pfilt, crossover_idx_pfilt = find_crossover_frequency(mag_pfilt[1, 1, :], w)

println("\nPref0 → pfilt_y (Power Reference Filter):")
println("  DC gain: $(round(dc_gain_db_pfilt, digits=1)) dB")
if crossover_pfilt !== nothing
    println("  Crossover frequency (0 dB): $(round(crossover_pfilt, digits=2)) rad/s ($(round(crossover_pfilt/(2π), digits=2)) Hz)")
    phase_at_crossover = phase_pfilt_deg[crossover_idx_pfilt]
    phase_margin = 180 + phase_at_crossover
    println("  Phase at crossover: $(round(phase_at_crossover, digits=1))°")
    println("  Phase margin: $(round(phase_margin, digits=1))°")
else
    println("  Crossover frequency: Not found in range")
end

# Interpretation
println("\n--- Interpretation ---")
println("For monotonically decreasing Bode plots (integrating systems):")
println("• Crossover frequency = frequency where gain = 1 (0 dB)")
println("• Phase margin = 180° + phase at crossover (stability metric)")
println("• Below crossover: system tracks reference (gain > 1)")
println("• Above crossover: system attenuates disturbances (gain < 1)")

println("\n" * "="^60)
println("Analysis complete.")
println("="^60)
