"""
    Powerful

Umbrella package of Powerful.jl, a power system simulation framework. It re-exports
the curated public API for power flow, time-domain simulation, eigenvalue analysis
and reporting; model types live in the [`Powerful.Models`](@ref) submodule.

Test cases come from PowerfulCases, which is loaded explicitly. Functionality beyond
the public API is imported from the sub-packages:
- `PowerfulCore`: foundation types, DAE infrastructure, address management
- `PowerfulPF`: power flow models and solvers (advanced API)
- `PowerfulTDS`: time-domain simulation models and solvers
- `PowerfulYbus`: Y-bus time-domain simulation, standalone
- `PowerfulCases`: test case data

# Quick Start
```julia
using Powerful
using PowerfulCases

# Load a system and solve the power flow
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# Load dynamic models and run a time-domain simulation
dyr_data = parse_dynamics(case.dyr, sys)
init = setup_tds(sol, dyr_data)
tds_prob = TDSProblem(init; tspan=(0.0, 10.0))
tds_sol = solve(tds_prob)
```

# Advanced Usage
Internal functions are imported from the sub-package that defines them:
```julia
using PowerfulCore: AddressManager, setup_with_formulation!, create_adcache
using PowerfulPF: setup_powerflow_newton, NewtonWorkspace
using PowerfulTDS: InitContext
```
"""
module Powerful

# Import all sub-packages (not @reexport - we curate exports explicitly)
using PowerfulCore
using PowerfulPF
using PowerfulTDS
using PowerfulCases

# ============================================================================
# PUBLIC API - Core Types
# ============================================================================
export SystemModel, ModelLinks
export DAE  # For advanced users
export TopologyInfo, analyze_topology, analyze_topology!, split_islands
export parse_network, Network  # Re-exported from PowerFlowData via PowerfulCore

# Formulations
export PowerInjection, CurrentInjection, DCPowerFlow

# Configuration
export SimConfig, SystemConfig, load_config, load_config_from_string
export TDSConfig
export get_solver_config

# Bus type constants
export BUS_PQ, BUS_PV, BUS_SLACK, BUS_ISOLATED

# ============================================================================
# PUBLIC API - Power Flow
# ============================================================================
# Problems
export ACPowerFlowProblem
export PowerFlowProblem  # Abstract type

# Algorithms
export Newton, RobustNR
export FastDecoupledXB, FastDecoupledBX

# Solutions
export solve  # From CommonSolve
export PowerFlowSolution, MultiIslandSolution

# Variable access
export get_var, get_var_by_index
export get_all_vars, get_vars_dict

# Reporting
export collect_model_table, print_model_table
export get_available_outputs, get_default_outputs
export ReportConfig, UnitPreferences, FormattingPreferences, ColumnPreferences
export Default, AllPerUnit, AllPhysical, Minimal, Compact
export PowerFlowReport, build_powerflow_report, format_text

# ============================================================================
# PUBLIC API - Events (used in examples)
# ============================================================================
# Phasor domain events
export BusFault
export BranchTrip, BranchClose
export FieldModification, SetpointChange

# ============================================================================
# PUBLIC API - TDS (Time Domain Simulation)
# ============================================================================
export TDSProblem, TDSInitialization, setup_tds
export TDSSolution, get_series
export parse_dynamics
export validate_tds_initialization, print_validation_report
export PowerFlowInit  # DAE initialization algorithm for DifferentialEquations.jl
export compute_eigenvalues, EigenvalueResult
export print_eigenvalue_report, write_eigenvalue_report
export eigenvalue_sensitivity, EigenvalueSensitivity, damping_ratio_sensitivity
export with_parameters
export linearize
export playback_voltage  # The recorded voltage of a VoltagePlayback at a time (its start)

# ============================================================================
# PUBLIC API - Test Cases
# ============================================================================
export CaseBundle, list_files

# Options of the public constructors and solvers
using PowerfulCore.DiffEq: QNDF, FBDF, KLUFactorization, NoInit, BrownFullBasicInit
export NaturalAddressing, AMDAddressing, CustomPermutation  # `addressing` of problems
export KLUSolver  # `linear_solver` of Newton
export ForwardDiffBackend, EnzymeBackend  # `backend` of problems
export QNDF, FBDF, KLUFactorization  # integrators of solve(::TDSProblem, alg)
export NoInit, BrownFullBasicInit  # `initializealg` of solve(::TDSProblem), with PowerFlowInit

# Injection recorder (the `recorder` of solve(::TDSProblem))
export InjectionMonitor, InjectionRecorder, contributions_table, get_contributions,
       get_timeseries, list_contributors, list_monitored_addresses

# ============================================================================
# MODELS SUBMODULE - All model types
# ============================================================================
# Models are organized in a submodule to keep the main namespace clean.
# Users can access models via:
#   using Powerful.Models: Bus, GENCLS, TGOV1
# or import all:
#   using Powerful.Models
"""
    Powerful.Models

Submodule containing all model types.

Model types are auto-registered from sub-packages (PowerfulPF, PowerfulTDS).

# Usage
```julia
using Powerful
using Powerful.Models  # Import all model types

# Or import specific models
using Powerful.Models: Bus, GENCLS
```

# Listing Available Models
```julia
names(Powerful.Models)  # List all exported symbols
```
"""
module Models
    using PowerfulPF
    using PowerfulTDS

    # Dynamically export all models from sub-packages
    for m in PowerfulPF.EXPORTED_MODELS
        @eval using PowerfulPF: $m
        @eval export $m
    end

    for m in PowerfulTDS.EXPORTED_MODELS
        @eval using PowerfulTDS: $m
        @eval export $m
    end
end

end # module Powerful
