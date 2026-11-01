"""
PowerfulTDS - Time Domain Simulation Models and Solvers

This package provides dynamic models (synchronous generators, exciters, governors,
renewable converters and controllers) and the phasor time-domain simulation,
initialization, event, eigenvalue and frequency response tools built on them.

# Usage
```julia
using Powerful
using PowerfulCases

# Solve the power flow first
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))
sol = solve(ACPowerFlowProblem(sys), Newton())

# Parse dynamic data, initialize and simulate
dyr_data = parse_dynamics(case.dyr, sys)
init = setup_tds(sol, dyr_data)
tds_sol = solve(TDSProblem(init; tspan=(0.0, 10.0)))
```
"""
module PowerfulTDS

using DocStringExtensions
using StructArrays
using SparseArrays
using LinearAlgebra

# External dependencies
using InlineStrings
using OrderedCollections
using PowerDynData

# Import CommonSolve to extend solve
import CommonSolve: solve

# PowerfulPF and PowerfulCore
using PowerfulPF
using PowerfulPF: calc_zip_load  # ZIP load physics (ZIP conversion at the operating voltage)
using PowerfulCore

# Import traits to extend from PowerfulCore
import PowerfulCore: model_metadata, calc_residuals!, calc_automaton_values!, get_terminals, get_key_name, compute_observations
import PowerfulCore: evaluate_system_residuals!, build_status_dict  # For zero-allocation hot path in TDSAdapter
using PowerfulCore: set_hold!, check_formulation_support
import PowerfulCore: implements_group, group_var_mapping, group_res_mapping
import PowerfulCore: to_callback, register_solver_config!, create_adcache
import PowerfulCore: get_replacement_map, get_data_inheritance
import PowerfulCore: get_controlled_group, parse_dynamic_model, get_dyr_model_name
import PowerfulCore: requires_initialization, find_control_terminal, check_parameter_domains
using PowerfulCore: is_zip_load, declares_zip_load, declares_requires_initialization, declares_check_parameter_domains
import PowerfulCore: base_quantity_role, supports_formulation, supports_format

# Initialize model registration (must be before model includes)
@init_model_registry

# ============================================================================
# Initialization Types (included early - models need InitContext)
# ============================================================================
include("solvers/initialization.jl")

# ============================================================================
# Generator Models
# ============================================================================
include("models/generators/gencls.jl")
include("models/generators/genrou.jl")
include("models/generators/gentpj1.jl")
include("models/generators/gensal.jl")


# ============================================================================
# Governor Models
# ============================================================================
include("models/governors/tgov1.jl")


# ============================================================================
# Exciter Models
# ============================================================================
include("models/exciters/sexs.jl")
include("models/exciters/ieeex1.jl")


# ============================================================================
# Renewable Generator Models
# ============================================================================
include("models/renewable/regca1.jl")
include("models/renewable/regfm_a1.jl")


# ============================================================================
# Renewable Electrical Control Models
# ============================================================================
include("models/renewable/reeca1.jl")


# ============================================================================
# Network Models (Fault for TDS)
# ============================================================================
include("models/fault.jl")

# ============================================================================
# Sources (recorded voltage played back at a bus)
# ============================================================================
include("models/sources/voltage_playback.jl")

# GPU support of the dynamic models (none has a GPU kernel)
include("models/gpu_support.jl")

# Phase representation of the dynamic models (positive sequence)
include("models/phase_representation.jl")

# Whether the dynamic models are ZIP loads (none is)
include("models/zip_load.jl")

# Parameter domains of the hand-maintained dynamic models (none has constraints of its own;
# the DSL-generated models declare theirs with the model)
include("models/check_parameter_domains.jl")

# Docstrings of the dynamic model structs (kept out of the generated model files)
include("models/docstrings.jl")


# ============================================================================
# TDS Solver
# ============================================================================
include("solvers/linking.jl")
include("solvers/tds.jl")         # TDSProblem, TDSInitialization definitions
include("solvers/value_degree.jl")  # The polynomial degree of the values of automata (traced)
include("solvers/rounding_bound.jl")  # The rounding of the values of automata (a running error bound)
include("solvers/state_projection.jl")  # The correction of states onto linear constraints (index reduction, governed states)
include("solvers/automaton_runtime.jl")  # Automata: start and transitions (held by TDSAdapter)
include("solvers/governed_states.jl")  # States a mode governs (Mode.states): admission, entry, start
include("solvers/adapter.jl")      # TDSAdapter (uses TDSProblem)
include("solvers/reinit.jl")       # Algebraic reinitialization (uses TDSAdapter)
include("solvers/functors.jl")     # Precompilable functors for ODE callbacks
include("solvers/erased_callbacks.jl")  # Callbacks with small types (ErasedCall, erase_callback, erase_init)
include("solvers/index_reduction.jl")  # Reads/uses of a state-only constraint reduction, generic over its algorithm
include("solvers/initial_limits.jl")  # Initial values against their declared limits
include("solvers/events.jl")       # Event callbacks
include("solvers/automaton_rates.jl")  # Rates of the values of automata along the trajectory
include("solvers/step_extremum.jl")  # Extrema of the automaton callback's components inside one step
include("solvers/automata.jl")     # Automaton callback and settling after events
include("solvers/validation.jl")       # Validation (uses TDSProblem, TDSAdapter)
include("solvers/linearization_utils.jl")  # Shared utilities for eigenvalue and bode
include("solvers/eigenvalue.jl")       # Eigenvalue analysis for small-signal stability
include("solvers/bode.jl")             # Linearization for frequency response analysis
include("solvers/parameter_override.jl")  # with_parameters: a copy of dynamic models with some parameters overridden
include("solvers/parameter_seeding.jl")  # Model parameters seeded with Dual numbers (sensitivities)
include("solvers/eigen_sensitivity.jl")  # Eigenvalue derivatives with respect to model parameters
include("solvers/solve.jl")            # solve() method

# ============================================================================
# Formats - DYR file parsing
# ============================================================================
include("formats/psse_dyr.jl")
# exports: parse_dynamics, create_model_links_for_type

# ============================================================================
# Exported Model Names (for umbrella package)
# ============================================================================
# This must come after all @register_model calls (model includes)
const EXPORTED_MODELS = copy(_EXPORTED_MODEL_NAMES)

# ============================================================================
# Runtime Initialization
# ============================================================================
function __init__()
    CHECK_MEMORY_SYNC[] = haskey(ENV, "POWERFUL_CHECK_MEMORY_SYNC")
    # Register solver config (must be in __init__ to survive precompilation)
    register_solver_config!(:TDS, TDSConfig)

    # Register models
    register_pending_models!(_PENDING_MODELS)
end

# ============================================================================
# Interface for packages that build on PowerfulTDS (public, not exported)
# ============================================================================

# Hooks for bus models of other packages (bus faults)
public accepts_bus_fault, bus_fault

# Initialization and DAE problem construction for systems built around the TDS models
public initialize_dynamic_models, create_diffeq_problem

# ============================================================================
# Precompilation Workload
# ============================================================================
include("precompile.jl")

end # module
