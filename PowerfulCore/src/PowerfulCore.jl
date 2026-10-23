"""
PowerfulCore - Foundation module for power systems simulation

This package provides the core type system, DAE infrastructure, and system utilities
that form the foundation of the Powerful ecosystem.

# Contents
- Foundation: Core types (AbstractModel, ModelVar, ModelResidual, etc.)
- DAEInfra: Address management, gather/scatter, residual workspaces
- System: SystemBuffers, per-unit conversion, topology analysis
- Registry: Model registration system

# Usage
```julia
using PowerfulCore

# Access core types
abstract type MyModel <: AbstractModel end

# Use address management
am = AddressManager()
assign_address!(am, :MyModel, :x, AlgebVar, 5)

# Use DAE vectors
dae = DAE(Float64, Vector, 10, 5)
```
"""
module PowerfulCore

using DocStringExtensions
using StructArrays
using SparseArrays
using AMD
using TOML
using Logging

# ============================================================================
# Verbose Printing Macro - Zero-allocation conditional output
# ============================================================================
"""
    @vprint "message with \$interpolation"

Print message only if `verbose` variable in scope is true.
Zero allocation when `verbose=false` because string interpolation is deferred.

Requires a `verbose::Bool` variable to exist in the calling scope.

# Example
```julia
function solve(prob; verbose=false)
    @vprint "Setting up problem..."
    # work...
    @vprint "Solving with \$(typeof(alg))..."
end
```
"""
macro vprint(expr)
    quote
        if $(esc(:verbose))
            println($(esc(expr)))
        end
    end
end

export @vprint

# Re-export PowerFlowData's main API so downstream packages don't need direct dependency
using PowerFlowData: parse_network, Network
export parse_network, Network

# ============================================================================
# DiffEq namespace - Unified OrdinaryDiffEq exports for all packages
# ============================================================================
# This replaces the umbrella DifferentialEquations.jl dependency with specific
# modular OrdinaryDiffEq packages for smaller sysimage and faster compilation.

"""
    DiffEq

Namespace containing all OrdinaryDiffEq symbols needed for time-domain simulation.
This provides a unified interface that replaces the DifferentialEquations.jl umbrella package.

Includes:
- ODE infrastructure: ODEFunction, ODEProblem
- BDF solvers: QNDF, FBDF (for stiff DAE systems)
- Callbacks: CallbackSet, DiscreteCallback, SavingCallback, PresetTimeCallback, FunctionCallingCallback
- Utilities: NoInit, KLUFactorization, u_modified!, solve
- Symbolic indexing: SymbolCache (variable names for `ODEFunction(...; sys=...)`)
- Return codes: ReturnCode
"""
module DiffEq
    using OrdinaryDiffEqCore: ODEFunction, ODEProblem
    using OrdinaryDiffEqBDF: QNDF, FBDF
    using OrdinaryDiffEqSDIRK: Trapezoid
    using DiffEqCallbacks: CallbackSet, DiscreteCallback, SavingCallback, SavedValues, FunctionCallingCallback, PresetTimeCallback
    using DiffEqBase: NoInit, BrownFullBasicInit, u_modified!
    using SciMLBase: ReturnCode, set_proposed_dt!, VectorContinuousCallback, RightRootFind,
                     derivative_discontinuity!, terminate!
    using SymbolicIndexingInterface: SymbolCache
    using LinearSolve: KLUFactorization
    using CommonSolve: solve, init, solve!

    export ODEFunction, ODEProblem, QNDF, FBDF
    export CallbackSet, DiscreteCallback, SavingCallback, SavedValues, FunctionCallingCallback, PresetTimeCallback
    export NoInit, BrownFullBasicInit, KLUFactorization, solve, init, solve!, u_modified!
    export ReturnCode, set_proposed_dt!, VectorContinuousCallback, RightRootFind,
           derivative_discontinuity!, terminate!
    export SymbolCache
end

export DiffEq

# ============================================================================
# Foundation - Core type system (rarely changes)
# ============================================================================
include("foundation/Foundation.jl")
using .Foundation

# Re-export all Foundation exports
export VIRTUAL_BUS_BASE
export AbstractModel
export AddressableType, AddressableVar, AddressableRes
export AlgebVar, StateVar, ObservableVar, AlgebRes, StateRes
export VarTrait, Internal, External
export VarProperty, Description, Units, MassCoefficient, PropertyDict, eval_mass_coeff
export ModelVar, ModelResidual
export SharedRes, PrivateRes
export AddressManager, ModelAddressAssignment, ModelRetrievedAddresses, ReverseAddressMap
export AddressingStrategy, NaturalAddressing, AMDAddressing, CustomPermutation
export GroupLinkRef, GroupRefKey, IndexOrder, is_heterogeneous, build_terminal_refs
export PUType, NoPU
export PowerNominal, ActivePower, ReactivePower, VoltageNominal, CurrentNominal
export ImpedanceNominal, AdmittanceNominal, AngleDegrees
export PowerPU, VoltagePU, CurrentPU, ImpedancePU, AdmittancePU
export RatingType, NominalPower, NominalVoltage, NominalCurrent
export RatingSource, OwnField, LinkedDevice, RatingRatio
export ReferenceType, BaseVoltage, BasePower, BaseCurrent, Dimensionless
export RetrievalPhase, PrePerUnit, PostPerUnit
export ModelParam, ModelObservable, ModelRating, Terminal, ExtParam, InitParam, InitInput
export VariableLimit, InjectionForms
export Disturbance, to_callback
export ModelMetadata
export SystemModel, ModelLinks
export DataFormat, PSSE, PSSEDYR

# Constants
export BUS_PQ, BUS_PV, BUS_SLACK, BUS_ISOLATED

# Format traits
export FormatSupport, NoFormatSupport
export supports_format, parse_model, is_format_supported, validate_raw
export format_records, declares_format_records, record_counts
export load_models

# Formulation traits
export FormulationTrait, PowerInjection, CurrentInjection, DCPowerFlow
export DEFAULT_FORMULATION, supports_formulation

# Model variable/residual utilities
export is_internal, is_external
export get_property, get_description, get_units, get_mass_coeff
export external_vars

# Model metadata interface
export get_base_model_type
export model_metadata, calc_residuals!, compute_observations, get_status_field
export has_state_vars
export get_key_name

# Metadata query predicates
export is_state_var, is_algeb_var, is_state_res, is_algeb_res
export is_private_res, is_shared_res
export variable_limits
export get_metadata_for_formulation

# Base quantity traits
export BaseQuantityRole, OwnsBase, InheritsBase, NoBase
export base_quantity_role, owns_base, inherits_base, needs_base
export get_base_quantity_vector

# Group system
export GroupVar, GroupRes, GroupInterface
export register_group!, get_group_interface, has_group_interface, all_groups
export implements_group, get_controlled_group, group_var_mapping, group_res_mapping
export validate_group_implementation
export resolve_group_var_name, resolve_group_res_name
export is_group_reference

# Terminal utilities (from model_metadata.jl)
export get_terminals, find_terminal

# Model traits (TDS/dynamic model support)
export PowerSystemComponent, PowerFlowModel, DynamicModel
export GeneratorModel, LoadModel, NetworkElement
export requires_initialization, get_dyr_model_name
export check_parameter_domains
export get_replacement_map, get_data_inheritance
export find_control_terminal, parse_dynamic_model

# Solution interface (abstract interface for steady-state solutions)
export AbstractSteadyStateSolution
export get_system, get_address_manager, get_state_vector, is_converged, get_residual_norm

# ============================================================================
# Model Registry (must come before System for all_models availability)
# ============================================================================
include("registry.jl")
export @register_model, @init_model_registry, MODEL_REGISTRY
export get_registered_models, get_model_type
export all_models
export register_pending_models!

# ============================================================================
# DAE Infrastructure - Address management, vectors, workspaces
# ============================================================================
# NOTE: config_types.jl (SimConfig) is included inside DAEInfra; SystemConfig is
# defined in Foundation (owned by SystemModel)
# because ResidualBuffers needs it, and DAEInfra is a submodule.
include("dae/DAEInfra.jl")
using .DAEInfra
using .DAEInfra: is_gpu_array_type, device_sparse  # GPU hooks (not exported)

# Re-export DAEInfra exports
export DAE
export assign_address!, get_addresses, get_var_address_type
export remap_addresses!
export build_reverse_map!, lookup_address, variable_name
export setup_with_formulation!
export retrieve_addresses, injection_routes, injection_form, UNROUTED
export get_key
export ext2int
export get_index_by_key, find_index, normalize_key, get_field_by_key, get_foreign_indices, get_foreign_link_ref
export normalize_ref_key, instance_key, build_index_map, get_group_ref, has_group_ref, build_group_lookup
export gather_from_group
export max_address_in_space
export gather!, gather_cpu!
export scatter!, scatter_cpu!, scatter_unmasked!, PartialAddresses
export ResidualBuffers
export AutomatonState, calc_automaton_values!, evaluate_automaton_values!, enforce_governed!
export AffineRelation, refresh_affine_governed!, enforce_governed_affine!
export gather_variables!, variable_gather, gather_residuals!, zero_residuals!, scatter_residuals!, scatter_residuals_accumulate!
export SystemConfig, SimConfig  # SystemConfig from Foundation, SimConfig from config_types.jl
export jacobian_sparsity_pattern, instance_residual_rows, matrix_colors, PrecomputedSparsityDetector
export TDSSolution, get_series, state_address

# ============================================================================
# System - High-level system operations
# ============================================================================
include("system/System.jl")
using .System
# Extended by the AD caches, which copy the memory of the automata into their workspaces
import .System: sync_automaton_memory!

# Re-export System exports
export set_algebraic_guess!
export filter_instances, add_instances!, derive_system
export SystemBuffers
export evaluate_system_residuals!, evaluate_generator_contributions!
export retrieve_external_params, merge_external, external_dependents, refresh_external!
export extract_device_ratings, extract_external_bases
export calc_pu_coefficients, convert_to_per_unit
export convert_system_to_per_unit, ParsedSystemData
export TopologyInfo, analyze_topology, analyze_topology!, is_bridge_device
export find_devices_on_vertex, build_reverse_terminal_map, disable_devices!, get_partition_after_outage
export terminal_vertex_indices, in_service
export split_islands
export OutageAction, BranchOutage, GeneratorOutage, Transformer2WOutage, Contingency
export resolve_outage, apply_contingency!, deenergize_islands!
export compute_amd_ordering

# Configuration system exports
export load_config, load_config_from_string
export get_solver_config
export register_solver_config!

# ============================================================================
# Injection Recorder - Convenience constructor (needs SystemBuffers from System)
# ============================================================================
# Core types and functions are in DAEInfra. This constructor is defined here
# because it needs SystemBuffers which is loaded after DAEInfra.
export InjectionMonitor, InjectionRecorder
export record_model_contributions!, commit_contributions!
export get_contributions, get_timeseries
export contributions_table
export list_monitored_addresses, list_contributors

"""
    InjectionRecorder(monitor::InjectionMonitor{K}, prob::TDSProblem) where K

Convenience constructor that creates recorder from TDSProblem.
Creates a temporary SystemBuffers internally for address mapping setup.

# Arguments
- `monitor`: InjectionMonitor specifying which addresses to track
- `prob`: TDSProblem (provides sys, address_manager, and formulation)

# Example
```julia
bus_keys = [(i,) for i in tds_prob.sys.models.Bus.i]
monitor = InjectionMonitor(:Bus => [:p => bus_keys, :q => bus_keys])
recorder = InjectionRecorder(monitor, tds_prob)
sol = solve(tds_prob; recorder=recorder)
```
"""
function DAEInfra.InjectionRecorder(
    monitor::DAEInfra.InjectionMonitor{K},
    prob  # TDSProblem - duck typed to avoid circular dependency
) where K
    sys_ws = SystemBuffers(Float64, prob.address_manager, prob.sys, prob.formulation)
    return DAEInfra.InjectionRecorder(monitor, sys_ws, prob.sys, prob.address_manager)
end

# ============================================================================
# ADCache - Autodiff workspace caching (depends on both DAEInfra and System)
# ============================================================================
# Included here (not in DAEInfra) because it needs SystemBuffers from System
# ADCache accessors extend the Foundation solution-interface function
import .Foundation: get_address_manager, get_system
include("dae/autodiff_cache.jl")

# ============================================================================
# Newton-Raphson engine (power flow, contingency analysis, TDS reinitialization)
# ============================================================================
include("solvers/newton.jl")
export NewtonPolicy, NewtonBuffers, newton_solve!, BUS_VOLTAGE_STEP_LIMIT
export AbstractNewtonLinearSolver, KLULinearSolver
export factorize!, linear_solve!, reset_factorization!, is_factorized

# ============================================================================
# Interface for packages that build on PowerfulCore (public, not exported)
# ============================================================================

# GPU hooks, extended by a GPU package
using .System: synchronize_if_gpu, to_gpu, check_gpu_support
using .Foundation: GPUSupport, SupportedOnGPU, UnsupportedOnGPU, gpu_support, declares_gpu_support
public is_gpu_array_type, device_sparse, synchronize_if_gpu, to_gpu, gpu_linear_solver, check_gpu_support,
    GPUSupport, SupportedOnGPU, UnsupportedOnGPU, gpu_support, declares_gpu_support

# Phase representation of the models (positive sequence or three-phase)
using .Foundation: PhaseRepresentation, PositiveSequence, ThreePhase, phase_representation,
    declares_phase_representation
public PhaseRepresentation, PositiveSequence, ThreePhase, phase_representation, declares_phase_representation

# Whether a model is a ZIP load the time-domain simulation converts
using .Foundation: is_zip_load, declares_zip_load
public is_zip_load, declares_zip_load

# Whether a model is declared to require initialization (requires_initialization)
using .Foundation: declares_requires_initialization
public declares_requires_initialization

# Whether a model declares its parameter domains (check_parameter_domains)
using .Foundation: declares_check_parameter_domains
public declares_check_parameter_domains

# System status and topology helpers for analyses built on a SystemModel
public build_status_dict, copy_status_dict, build_incidence_matrix, check_formulation_support

# Dual workspace invalidation for packages that change setpoints (e.g. a parameter system)
public invalidate_dual_workspaces!

# Automata: metadata of a model's discrete behaviour (Foundation/automata.jl) and its
# state in the workspaces (dae/automaton_state.jl), for models and the simulation
using .Foundation: ValueRef, Condition, Atom, AllOf, AnyOf, Negation, atoms, Mode, Transition,
    InitialMode, Held, Automaton, MemoryQuantity, discrete_quantities, automaton_values,
    check_automata, has_automata, is_condition, governed_states, governs_states
public ValueRef, Condition, Atom, AllOf, AnyOf, Negation, atoms, Mode, Transition, InitialMode,
    Held, Automaton, MemoryQuantity, discrete_quantities, automaton_values, check_automata, governed_states, governs_states,
    has_automata, is_condition, automaton_states, sync_automaton_memory!, automaton_memory_synced, with_float_columns

# Parameters: model fields and setpoints written by element type (parameter systems,
# linearization inputs)
public ParameterSlice, parameter_array, register_parameter_fields!, register_parameters!
public gather_parameters, scatter_parameters!

end # module
