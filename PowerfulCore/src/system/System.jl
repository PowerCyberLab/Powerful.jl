"""
System module - High-level system operations for Powerful.jl

This module provides:
- SystemBuffers: Workspace for residual evaluation
- Per-unit conversion utilities
- Topology analysis (island detection, AMD ordering)
- Configuration system (SimulationConfig, settings loading)
"""
module System

using DocStringExtensions
using StructArrays
using SparseArrays
using AMD
using TOML

using ..Foundation
using ..Foundation: AbstractModel, ModelMetadata
using ..Foundation: SystemModel, FormulationTrait, IndexOrder
using ..Foundation: model_metadata, get_base_model_type, get_terminals
using ..Foundation: PrePerUnit, PostPerUnit
using ..Foundation: has_automata
using ..Foundation: gpu_support, declares_gpu_support, UnsupportedOnGPU

using ..DAEInfra

import ..all_models
using ..DAEInfra: AddressManager, ResidualBuffers, DAE
using ..DAEInfra: gather_variables!, gather_residuals!, scatter_residuals!, scatter_residuals_accumulate!
using ..DAEInfra: record_model_contributions!  # Injection recorder
using ..DAEInfra: instance_key  # For reverse Terminal lookups

# ============================================================================
# Model Utilities
# ============================================================================
include("model_utils.jl")
export set_algebraic_guess!

# ============================================================================
# Configuration System
# ============================================================================
# NOTE: SimConfig is in config_types.jl (SystemConfig is in Foundation), included
# early in PowerfulCore.jl (before DAEInfra) for use in ResidualBuffers.
# This file contains the config access functions, TOML loading, and registry.
include("config.jl")
# SystemConfig, SimConfig already exported from PowerfulCore.jl

# ============================================================================
# System Buffers and Residual Evaluation
# ============================================================================
include("system_buffers.jl")

# ============================================================================
# Per-Unit Conversion
# ============================================================================
include("per_unit.jl")

# ============================================================================
# System Assembly (removing / adding instances)
# ============================================================================
include("assembly.jl")
export filter_instances, add_instances!, derive_system

# ============================================================================
# Topology Analysis
# ============================================================================
include("topology.jl")
export TopologyInfo, analyze_topology, analyze_topology!, is_bridge_device
export find_devices_on_vertex, build_reverse_terminal_map, disable_devices!, get_partition_after_outage
export compute_amd_ordering
export build_adjacency_matrix, build_incidence_matrix, extract_from_to_idx
export terminal_vertex_indices, in_service

include("outages.jl")
export OutageAction, BranchOutage, GeneratorOutage, Transformer2WOutage, Contingency
export resolve_outage, apply_contingency!, deenergize_islands!

include("multi_island.jl")
export split_islands

# ============================================================================
# SystemModel Convenience Constructors
# ============================================================================
# Must come after per_unit.jl (ParsedSystemData) and system_buffers.jl (retrieve_external_params)
include("system_builder.jl")
# No exports needed - constructors extend SystemModel from Foundation

end # module
