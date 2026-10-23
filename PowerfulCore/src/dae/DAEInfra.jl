"""
DAEInfra module - Differential-Algebraic Equation infrastructure for Powerful.jl

This module provides the core DAE data structures and operations:
- Address management for mapping model variables to DAE vectors
- External-to-internal index conversion
- DAE vector storage (variables and residuals)
- Gather/scatter kernels for efficient data movement
- Residual workspace for model evaluation
"""
module DAEInfra

using StructArrays
using SparseArrays: SparseMatrixCSC
using DocStringExtensions

using ..Foundation
using ..Foundation: AbstractModel, ModelMetadata
using ..Foundation: AddressManager, ModelAddressAssignment, ModelRetrievedAddresses
using ..Foundation: AlgebVar, AlgebRes, StateVar, StateRes, ObservableVar
using ..Foundation: GroupLinkRef, GroupRefKey, IndexOrder, FormulationTrait
using ..Foundation: is_internal, model_metadata, get_status_field, get_key_name
using ..Foundation: is_state_res, is_shared_res
using ..Foundation: check_automata, discrete_quantities, automaton_values, has_automata, governs_states
# Solution interface: TDSSolution methods extend the Foundation functions
import ..Foundation: get_system, get_address_manager, get_state_vector, is_converged

# ============================================================================
# Address Manager - Variable/residual address assignment
# ============================================================================
include("address_manager.jl")
export assign_address!, get_addresses, get_var_address_type
export remap_addresses!
export build_reverse_map!, lookup_address, variable_name
export setup_with_formulation!
export retrieve_addresses, injection_routes, injection_form, UNROUTED
export get_key, get_key_name

# ============================================================================
# External-to-Internal Indexing
# ============================================================================
include("ext2int.jl")
export GroupLinkRef, GroupRefKey
export build_index_map, get_group_ref, has_group_ref, build_group_lookup

# Unified helpers for group references

"""
    get_group_ref(refs::Dict, model::Symbol, key) -> GroupLinkRef

Retrieve a GroupLinkRef for the given model/key from a group reference map
(`IndexOrder.foreign_refs`). Normalizes the key and throws a clear error if
missing.

# Arguments
- `refs::Dict` - Group reference dictionary (e.g., `sys.order.foreign_refs`)
- `model::Symbol` - Model name (e.g., `:TGOV1`, `:PQ`)
- `key` - Field key (Symbol for single-field, Vector/Tuple for composite keys)

# Returns
`GroupLinkRef` containing target types and indices

# Example
```julia
# Single-field terminal (bus connection)
ref = get_group_ref(sys.order.foreign_refs, :PQ, :i)

# Composite-key terminal (controller link)
ref = get_group_ref(sys.order.foreign_refs, :TGOV1, [:BUS, :ID])
ref = get_group_ref(sys.order.foreign_refs, :TGOV1, (:BUS, :ID))  # Also works
```

# See Also
- `has_group_ref` - Check if reference exists before retrieval
"""
function get_group_ref(refs::Dict, model::Symbol, key)
    field_key = normalize_ref_key(key)
    if !haskey(refs, model)
        available = sort(collect(keys(refs)))
        error("Group references for model :$model not found. Available models: $available")
    end
    model_refs = refs[model]
    if !haskey(model_refs, field_key)
        available_keys = sort([string(k) for k in keys(model_refs)])
        error("Group reference :$model.$field_key not found. Available keys: $available_keys")
    end
    return model_refs[field_key]
end

"""
    get_group_ref(order::IndexOrder, model::Symbol, key) -> GroupLinkRef

Convenience wrapper that extracts `foreign_refs` from IndexOrder.
"""
get_group_ref(order::IndexOrder, model::Symbol, key) = get_group_ref(order.foreign_refs, model, key)

"""
    has_group_ref(refs::Dict, model::Symbol, key) -> Bool

Check if a GroupLinkRef exists for model/key. Normalizes the key automatically.

# Arguments
- `refs::Dict` - Group reference dictionary
- `model::Symbol` - Model name
- `key` - Field key (Symbol, Vector{Symbol}, or Tuple{Vararg{Symbol}})

# Returns
`true` if the reference exists, `false` otherwise

# Example
```julia
if has_group_ref(sys.order.foreign_refs, :TGOV1, (:BUS, :ID))
    ref = get_group_ref(sys.order.foreign_refs, :TGOV1, (:BUS, :ID))
end
```

# See Also
- `get_group_ref` - Retrieve the reference (throws if missing)
"""
function has_group_ref(refs::Dict, model::Symbol, key)
    field_key = normalize_ref_key(key)
    return haskey(refs, model) && haskey(refs[model], field_key)
end
has_group_ref(order::IndexOrder, model::Symbol, key) = has_group_ref(order.foreign_refs, model, key)

# ============================================================================
# DAE Vectors - Storage for variables and residuals
# ============================================================================
include("vectors.jl")
export DAE  # The struct for DAE vectors
export max_address_in_space

# ============================================================================
# Gather/Scatter Kernels - CPU/GPU data movement
# ============================================================================
include("gather_scatter.jl")
export gather!, gather_cpu!
export scatter!, scatter_cpu!, scatter_unmasked!, PartialAddresses

# ============================================================================
# Configuration Types (needed by ResidualBuffers)
# ============================================================================
# SimConfig struct - must be defined before residual_workspace.jl
include("../system/config_types.jl")

# ============================================================================
# Automaton State - Modes, discrete quantities and memory of a model's automata
# (held by the residual workspace)
# ============================================================================
include("automaton_state.jl")
export AutomatonState, calc_automaton_values!, evaluate_automaton_values!, enforce_governed!
export AffineRelation, refresh_affine_governed!, enforce_governed_affine!

# ============================================================================
# Residual Workspace - Pre-allocated buffers for model evaluation
# ============================================================================
include("residual_workspace.jl")
export ResidualBuffers
export gather_variables!, variable_gather, gather_residuals!, zero_residuals!, scatter_residuals!, scatter_residuals_accumulate!
# Low-level tuple builders for custom scatter patterns (used by Y-bus formulation)
export build_res_scatter_tuple, scatter_res_add!, scatter_res_set!

# ============================================================================
# Injection Recorder - Track per-model contributions to SharedRes addresses
# ============================================================================
include("injection_recorder.jl")
export InjectionMonitor, InjectionRecorder
export record_model_contributions!, commit_contributions!
export get_contributions, get_timeseries
export contributions_table
export list_monitored_addresses, list_contributors

# ============================================================================
# Sparsity Pattern - Jacobian sparsity detection from model structure
# ============================================================================
include("sparsity_pattern.jl")

# ============================================================================
# TDS Solution - Unified solution type for TDS and Ybus results
# ============================================================================
include("tds_solution.jl")
export TDSSolution, get_series, state_address

end # module
