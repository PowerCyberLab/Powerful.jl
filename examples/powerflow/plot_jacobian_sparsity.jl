#!/usr/bin/env julia
"""
Jacobian sparsity pattern visualization for ACTIVSg 70k bus system.

Produces a spy plot of the power flow Jacobian sparsity structure
suitable for presentations.
"""

using Powerful
using PowerfulCore: jacobian_sparsity_pattern
using PowerfulCases
using SparseArrays
using Plots

println("Loading ACTIVSg 70k bus system...")
case_bundle = load("synthetic/ACTIVSg70k")
case = parse_network(case_bundle.raw)
sys = SystemModel(case)
prob = ACPowerFlowProblem(sys, PowerInjection())

n_buses = length(sys.models.Bus)
println("  Buses: $n_buses")
println("  PQ loads: $(length(sys.models.PQ))")
println("  PV generators: $(length(sys.models.PV))")
println("  Branches: $(length(sys.models.Branch))")

println("Computing Jacobian sparsity pattern...")
am = prob.address_manager
@time sparsity = jacobian_sparsity_pattern(am, sys, PowerInjection())

n = size(sparsity, 1)
nz = nnz(sparsity)
density = nz / n^2 * 100

println("  Matrix size: $n × $n")
println("  Nonzeros: $nz")
println("  Density: $(round(density, digits=4))%")
println("  Sparsity: $(round(100 - density, digits=4))%")

# ============================================================
# Figure: spy plot
# ============================================================
println("Rendering spy plot...")

rows, cols, _ = findnz(sparsity)

p = scatter(cols, rows;
    markersize        = 0.3,
    markerstrokewidth = 0,
    color             = :black,
    alpha             = 0.6,
    legend            = false,
    yflip             = true,
    aspect_ratio      = 1,
    xlims             = (1, n),
    ylims             = (1, n),
    xlabel            = "Variable index",
    ylabel            = "Equation index",
    title             = "ACTIVSg 70k — Power Flow Jacobian Sparsity Pattern\n" *
                        "$n × $n   |   $nz nonzeros   |   $(round(density, digits=4))% density",
    titlefontsize     = 10,
    guidefontsize     = 10,
    tickfontsize      = 9,
    size              = (620, 600),
    dpi               = 300,
    top_margin        = 6Plots.mm,
    bottom_margin     = 4Plots.mm,
    left_margin       = 6Plots.mm,
    right_margin      = 4Plots.mm,
)

out_path = joinpath(@__DIR__, "jacobian_sparsity_70k.png")
savefig(p, out_path)
println("Saved to $out_path")

display(p)
