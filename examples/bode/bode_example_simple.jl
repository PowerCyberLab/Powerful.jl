"""
Bode Plot Example: Governor-to-Rotor Speed Transfer Function

This example demonstrates how to compute and plot the frequency response
from governor power reference (TGOV1.pref0) to rotor speed (GENCLS.omega).

The transfer function G(jω) = C(jωI - A)⁻¹B + D is computed from the
linearized DAE system at steady state.
"""

using Powerful
using PowerfulCases
using PowerfulTDS: linearize
using ControlSystemsBase: bode
using Plots

# =============================================================================
# Load System
# =============================================================================

bundle = load("two_bus_branch")
raw = bundle.raw
dyr = file(bundle, :dyr, variant = "tgov1")  # GENCLS with a TGOV1 governor

# Parse network and solve power flow
case = parse_network(raw)
sys = SystemModel(case)
pf_sol = solve(ACPowerFlowProblem(sys), Newton())

# Parse dynamic models and setup TDS
dyr_data = parse_dynamics(dyr, sys)
init = setup_tds(pf_sol, (GENCLS=dyr_data.GENCLS, TGOV1=dyr_data.TGOV1))
prob = TDSProblem(init; tspan=(0.0, 1.0))

# =============================================================================
# Compute Linearized State-Space Model
# =============================================================================

# Define inputs and outputs
# - inputs: (model_name, field_name, instances) where instances is 1-indexed
#   1:1 means "first instance only" (single input)
# - outputs: same format
inputs = [(:TGOV1, :pref0, 1:1)]   # Governor power reference setpoint
outputs = [(:GENCLS, :omega, 1:1)] # Rotor speed

# Linearize at operating point - returns ControlSystemsBase.StateSpace
state_space = linearize(prob, inputs, outputs; verbose=true)

# Print state-space summary
println("\nState-Space Model:")
println("  A: $(size(state_space.A))")
println("  B: $(size(state_space.B))")
println("  C: $(size(state_space.C))")
println("  D: $(size(state_space.D))")

# =============================================================================
# Compute Frequency Response using ControlSystemsBase
# =============================================================================

# Define frequency range
frequencies = 10.0 .^ range(-2, 2, length=100)

# Compute Bode data using ControlSystemsBase
mag, phase, w = bode(state_space, frequencies)

# Extract data for plotting (ControlSystemsBase returns (ny, nu, nω) arrays)
mag_db = 20 .* log10.(mag[1, 1, :])
phase_deg = phase[1, 1, :]

# Print DC gain
dc_gain_db = mag_db[1]
println("\nDC Gain: $(round(dc_gain_db, digits=2)) dB")

# =============================================================================
# Plot Results
# =============================================================================

# Magnitude plot
p1 = plot(w, mag_db,
    ylabel="Magnitude (dB)",
    xscale=:log10,
    legend=false,
    title="Bode Plot: TGOV1.pref0 → GENCLS.omega",
    linewidth=2,
    grid=true,
    minorgrid=true,
    xformatter=x -> ""
)

# Phase plot
p2 = plot(w, phase_deg,
    xlabel="Frequency (rad/s)",
    ylabel="Phase (deg)",
    xscale=:log10,
    legend=false,
    linewidth=2,
    grid=true,
    minorgrid=true
)

# Combine into single figure
fig = plot(p1, p2, layout=(2, 1), size=(400, 300), dpi=200)

# Save to file
savefig(fig, joinpath(@__DIR__, "bode_plot.png"))
println("\nSaved Bode plot to examples/bode/bode_plot.png")
