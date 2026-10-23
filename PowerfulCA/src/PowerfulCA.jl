"""
PowerfulCA - Contingency Analysis for Power Systems

This package provides AC and DC contingency analysis capabilities including:
- Sequential and parallel AC contingency calculation (ACCC)
- N-1 and N-k outage analysis
- Violation reporting (voltage limits, branch overloads)

# Usage
```julia
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))

contingencies = [
    Contingency("LINE_1_2", [BranchOutage(1, 2, "1")]),
    Contingency("GEN_3", [GeneratorOutage(3, "1")]),
]

# Sequential solve
result = solve_accc(sys, contingencies)

# Parallel solve (uses all available threads)
result = solve_accc_parallel(sys, contingencies; verbose=true)

# Print violation report
print_accc_report(result)
```
"""
module PowerfulCA

using DocStringExtensions
using LinearAlgebra
using SparseArrays

using Printf: Printf, @sprintf

# Import CommonSolve to use solve
import CommonSolve: solve

# PowerfulCore - Foundation types and DAE infrastructure
using PowerfulCore
using PowerfulCore: SystemModel
using PowerfulCore: build_status_dict, copy_status_dict, compute_observations, get_addresses
using PowerfulCore: analyze_topology, build_reverse_terminal_map, build_incidence_matrix

# PowerfulPF - Power flow solvers and models
using PowerfulPF
using PowerfulPF: ACPowerFlowProblem, Newton, setup_powerflow_newton, solve_powerflow_newton!
using PowerfulPF: Branch, calc_branch_power_flows

# ============================================================================
# Contingency Types (defined in PowerfulCore, shared with GPU screening)
# ============================================================================
export OutageAction, BranchOutage, GeneratorOutage, Transformer2WOutage, Contingency

# ============================================================================
# DC Contingency Analysis (PTDF/LODF)
# ============================================================================
include("ptdf.jl")
export PTDF, compute_ptdf, LazyPTDF
export get_ptdf_row, get_ptdf_column, get_ptdf_element

include("lodf.jl")
export LODF, compute_lodf, compute_lodf_column
export estimate_post_outage_flows

# ============================================================================
# AC Contingency Calculation
# ============================================================================
include("accc.jl")
export solve_accc, solve_accc_parallel

# ============================================================================
# DC Contingency Screening
# ============================================================================
include("dccc.jl")
export screen_contingencies_dc, DCContingencyScreenResult, DCScreeningResult
export BranchViolation, print_dc_screening_report

# ============================================================================
# Reporting
# ============================================================================
include("report.jl")
export print_accc_report

end # module
