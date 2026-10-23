
#==============================================================================
# Island Splitting
==============================================================================#

"""
    split_islands(sys::SystemModel; auto_analyze::Bool=true) -> Union{SystemModel, Vector{SystemModel}}

Extract disconnected islands as separate SystemModel instances.

Returns the original system if connected, or a vector of independent island systems.
Each island has its own AddressManager and IndexOrder, and can be solved independently.

**Note**: `ACPowerFlowProblem` calls this automatically during construction. For most
use cases, simply call `solve(ACPowerFlowProblem(sys), Newton())` which handles
islands transparently.

# Arguments
- `sys::SystemModel`: Original system
- `auto_analyze::Bool=true`: Run topology analysis if not cached

# Returns
- `SystemModel`: If single island (returns original, zero overhead)
- `Vector{SystemModel}`: If multiple disconnected islands

# Example
```julia
sys = SystemModel(case)
result = split_islands(sys)

if result isa Vector
    println("System has \$(length(result)) islands")
    for island in result
        sol = solve(ACPowerFlowProblem(island), Newton())
    end
else
    sol = solve(ACPowerFlowProblem(result), Newton())
end
```

# See Also
- [`analyze_topology`](@ref): Analyze connectivity without splitting
- [`ACPowerFlowProblem`](@ref PowerfulPF.ACPowerFlowProblem): Auto-handles islands during construction
"""
function split_islands(sys::SystemModel; auto_analyze::Bool=true)
    # Get or compute topology
    topo = if isnothing(sys.topology) && auto_analyze
        analyze_topology(sys)
    elseif !isnothing(sys.topology)
        sys.topology
    else
        error("No topology information. Call analyze_topology!(sys) first or use auto_analyze=true")
    end

    # Single island? Return original system (no need to wrap)
    if topo.n_islands == 1
        return sys
    end

    @debug "Splitting system into $(topo.n_islands) islands"

    # Create separate SystemModel for each island
    islands = Vector{SystemModel}(undef, topo.n_islands)

    for i in 1:topo.n_islands
        islands[i] = _extract_island(sys, topo, i)
    end

    # Return vector directly (no wrapper needed)
    return islands
end


"""
    _extract_island(sys::SystemModel, topo::TopologyInfo, island_idx::Int) -> SystemModel

Extract a single island as a new SystemModel (internal function).

Creates a complete independent SystemModel containing only the buses and devices
belonging to the specified island.
"""
function _extract_island(sys::SystemModel, topo::TopologyInfo, island_idx::Int)
    # Filter models to this island using Terminal API
    keep = Dict(name => (isempty(m) ? falses(0) : _get_island_mask(name, m, island_idx, topo.island_map, sys.order))
                for (name, m) in pairs(sys.models))
    island_models = filter_instances(sys.models, keep; drop_empty=true)

    # Filter per-instance original data and pu_coeffs with the same masks
    island_parsed_data = filter_instances(sys.parsed_data, keep; drop_empty=true)
    island_pu_coeffs = filter_instances(sys.pu_coeffs, keep; drop_empty=true)

    # Rebuild index mappings for island (creates new local IndexOrder with foreign_refs)
    models_indexed, index_order = ext2int(island_models)

    # Models are already in per-unit from parent filtering, so both phases retrieve at once (the
    # PrePerUnit ones recovered to the physical value through the coefficients)
    external = retrieve_external_params(models_indexed, models_indexed, index_order; pu_coeffs = island_pu_coeffs)

    # Create new SystemModel with independent address space
    # Note: AddressManager is created during problem construction, not stored in SystemModel
    island_sys = SystemModel(
        models_indexed,
        external,
        index_order;
        settings=sys.settings,
        array_type=sys.array_type,
        parsed_data=island_parsed_data,  # Filtered original data matching island instances
        pu_coeffs=island_pu_coeffs,  # Filtered coefficients matching island instances
        topology=nothing      # Island is connected by definition (single component)
    )

    return island_sys
end


"""
    _get_island_mask(model_name, model_vec, island_idx, island_map, order) -> BitVector

Determine which instances of a model belong to specified island.

Uses Terminal API to generically handle all model types:
- Bus: Direct lookup in island_map
- Multi-terminal devices (Branch, Transformer): All terminals must be in island
- Single-attachment devices (PQ, PV, Shunt): Attachment bus must be in island
- Topological models (Area, Zone): Skipped (no spatial location)
"""
function _get_island_mask(
    model_name::Symbol,
    model_vec::StructArray,
    island_idx::Int,
    island_map::Vector{Int},
    order::IndexOrder
)
    n = length(model_vec)
    mask = falses(n)
    model_type = get_base_model_type(model_vec)
    terminals = get_terminals(model_type)

    # Special case: Bus - direct mapping via island_map
    # Check by type name instead of subtype to avoid import issues
    type_name = nameof(model_type)
    if type_name == :Bus || type_name == :BusNum
        for i in 1:n
            mask[i] = (island_map[i] == island_idx)
        end
        return mask
    end

    # Models with terminals - ALL terminals must be in this island
    # This works for 1-terminal (PQ, PV, Shunt), 2-terminal (Branch),
    # or N-terminal devices (future 3-winding transformers, etc.)
    if !isempty(terminals)
        # Only :Bus vertices define islands for now (future: :Node for DC); devices on
        # a bus through a group terminal (e.g. generators via :bus_group) count too.
        # An instance attached to no bus keeps the previous behavior (every island).
        bus_vertices = [terminal_vertex_indices(order, model_name, t, :Bus) for t in _vertex_terminals(model_type, :Bus)]
        for i in 1:n
            mask[i] = all(bv -> bv[i] == 0 || island_map[bv[i]] == island_idx, bus_vertices)
        end
        return mask
    end

    # Topological models (Area, Zone, Owner) - no terminals, no clear island membership
    # Skip with debug message (could assign to all islands in future)
    @debug "Model $(nameof(model_type)) has no terminals - skipping in island extraction"

    return mask  # All false
end

