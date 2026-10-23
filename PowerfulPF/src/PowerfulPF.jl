"""
PowerfulPF - Power Flow Models and Solvers

This package provides network models (Bus, PQ, PV, Branch, etc.) and power flow
solvers (Newton-Raphson, Fast Decoupled) for steady-state power system analysis.

# Usage
```julia
using PowerfulCore, PowerfulPF
using PowerfulCases

# Load and solve (parse_network and SystemModel come from PowerfulCore)
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())  # or RobustNR()
```
"""
module PowerfulPF

using DocStringExtensions
using SparseArrays
using LinearAlgebra

using StructArrays: StructArrays, StructArray

# External dependencies
using InlineStrings: InlineStrings, String3, String15
using ForwardDiff: ForwardDiff
using DifferentiationInterface: DifferentiationInterface
using SparseMatrixColorings: SparseMatrixColorings
using KLU: KLU, klu, klu!
using ADTypes: ADTypes
using OrderedCollections: OrderedCollections, OrderedDict
using PrecompileTools: PrecompileTools, @compile_workload, @setup_workload
using PrettyTables: PrettyTables
using Printf: Printf, @sprintf
using Statistics: Statistics

# Import CommonSolve to extend solve
import CommonSolve: solve

# PowerfulCore

using PowerfulCore
using PowerfulCore: check_formulation_support
using PowerfulCore: is_gpu_array_type, gpu_linear_solver  # GPU hooks (a GPU package extends them)

# Import traits to extend
import PowerfulCore: model_metadata, calc_residuals!, compute_observations, get_terminals
import PowerfulCore: supports_format, format_records, parse_model, get_key_name, to_callback, set_algebraic_guess!
import PowerfulCore: get_replacement_map, get_data_inheritance
import PowerfulCore: parse_dynamic_model, get_dyr_model_name
import PowerfulCore: base_quantity_role, supports_formulation
import PowerfulCore: check_parameter_domains

# Initialize model registration (must be before model includes)
@init_model_registry

# ============================================================================
# Solution Types (included early for model compute_observations methods)
# ============================================================================
include("solvers/solution_types.jl")
# Exports: PowerFlowSolution, MultiIslandSolution, get_var/get_all_vars/get_vars_dict for solutions

# Model traits (get_replacement_map, requires_initialization, etc.) are now in PowerfulCore

# ============================================================================
# Group Definitions
# ============================================================================
include("models/group_definitions.jl")
register_standard_groups!()
export register_standard_groups!

# ============================================================================
# Topological Models (Areas, Zones, Owners)
# ============================================================================
include("models/topological/topological.jl")

# ============================================================================
# Power Flow Models
# ============================================================================
include("models/dc_physics.jl")  # DC power flow physics shared by kernels, B-matrix and DC screening
include("models/bus.jl")
include("models/pq.jl")
include("models/pv.jl")
include("models/slack.jl")
include("models/bus_types.jl")  # Bus classification from in-service devices
include("models/branch.jl")
include("models/transformer.jl")
include("models/transformer_conversion.jl")
include("models/shunt.jl")
include("models/switched_shunt.jl")  # held at binit, with the physics of shunt.jl (AbstractShunt)

# Phase representation of the models (positive sequence)
include("models/phase_representation.jl")

# Whether the models are ZIP loads (PQ is)
include("models/zip_load.jl")

# Whether the models require initialization (none does: they start from the power flow)
include("models/requires_initialization.jl")

# Parameter domains of the power flow models (none has constraints of its own)
include("models/check_parameter_domains.jl")

# Post-processing transformations
include("models/post_processing.jl")

# ============================================================================
# Solvers
# ============================================================================
# NOTE: sparsity_pattern.jl and autodiff_cache.jl live in PowerfulCore; `using PowerfulCore`
# above makes their exports available here:
# - jacobian_sparsity_pattern, matrix_colors, PrecomputedSparsityDetector
# - ForwardDiffCache, create_adcache, get_tds_dual_buffers

include("solvers/power_flow_api.jl")
# Exports from power_flow_api.jl

include("solvers/powerflow.jl")
# PowerFlowSolution is exported from solution_types.jl
export solve_powerflow_newton!, setup_powerflow_newton, NewtonWorkspace, ResidualFunctor
export setup_powerflow, solve_powerflow!

include("solvers/fastdecoupled.jl")
export solve_powerflow_fastdecoupled, solve_powerflow_fastdecoupled!, setup_powerflow_fastdecoupled

include("solvers/bmatrix.jl")
export gather_branch_susceptances, build_B_from_incidence, extract_from_to_idx, dc_network, dc_bus_injections

# ============================================================================
# Variable Accessors
# ============================================================================
include("access/variable_accessors.jl")
export get_var, get_var_by_index
export get_all_vars, get_vars_dict
export bus_voltages, terminal_voltages, generator_power

# ============================================================================
# I/O - Reports and output formatting
# ============================================================================
include("io/reports.jl")
export get_available_outputs, get_default_outputs

include("io/text_format.jl")

# ============================================================================
# Formats - RAW file parsing
# ============================================================================
include("formats/psse_raw.jl")
export load_models

# ============================================================================
# Exported Model Names (for umbrella package)
# ============================================================================
# This must come after all @register_model calls (model includes)
const EXPORTED_MODELS = copy(_EXPORTED_MODEL_NAMES)

# ============================================================================
# Model Registration (runtime)
# ============================================================================
# Register models collected during precompilation to the shared registry.
# Individual model files use @register_model which populates _PENDING_MODELS.
function __init__()
    register_standard_groups!()  # Register group interfaces (must be first)
    register_pending_models!(_PENDING_MODELS)

    # Provide helpful error when NonlinearSolve extension is not loaded
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f === _solve && length(argtypes) >= 2 &&
                argtypes[2] <: AlgorithmConfig{<:NewtonMethod, NonlinearSolveBackend}
            print(io, "\n\nNonlinearSolve.jl backend requested but not loaded.\n")
            print(io, "To use RobustNR(), first load NonlinearSolve:\n\n")
            print(io, "    using NonlinearSolve\n\n")
            print(io, "Or use the native Newton solver (no extra dependency):\n\n")
            print(io, "    solve(prob, Newton())\n")
        end
    end
end

# ============================================================================
# Interface for packages that build on PowerfulPF (public, not exported)
# ============================================================================

# Residual kernels of the network models, for other array types and formulations
public bus_residual_kernel!, pq_residual_kernel!, pv_residual_kernel!, branch_residual_kernel!,
       shunt_residual_kernel!, slack_residual_kernel!, slack_current_kernel!, slack_dc_kernel!

# Load, shunt and transformer physics shared with other network representations
public calc_zip_load, calc_shunt_power, calc_transformer_power_flows

# The fixed shunt models (AbstractShunt) and their admittance, for the network builders
public fixed_shunt_models, shunt_admittance_params, shunt_admittance

# Newton linear solver selector implemented by a GPU package
public CUDSSSolver

# ============================================================================
# Precompilation Workload
# ============================================================================
include("precompile.jl")

end # module
