"""
    topology.jl

Network topology analysis for identifying disconnected islands and connectivity.

Uses the Terminal API (get_terminals) to generically analyze network connectivity
without hardcoded knowledge of specific model types.

Dependencies: SparseArrays, get_terminals (provided by System module)
"""

"""
    TopologyInfo

Network topology analysis results.

# Fields
- `n_islands::Int`: Number of disconnected islands
- `island_map::Vector{Int}`: Maps bus internal index → island number (1-based)
- `slack_buses::Vector{Int}`: Slack bus internal index per island
- `adjacency::SparseMatrixCSC{Bool, Int}`: Bus connectivity matrix
- `components::Vector{Vector{Int}}`: Bus internal indices per island
- `incidence::Union{Nothing, SparseMatrixCSC{Int8, Int}}`: Branch-bus incidence matrix (L×N)
  - Row l has +1 at from-bus, -1 at to-bus for branch l
  - Only includes active (in-service) 2-terminal devices
  - Used for PTDF/LODF computation and deriving B matrices
- `is_bridge::BitVector`: Which rows of incidence matrix are bridges (cut-edges)
  - `is_bridge[l]` = true if removing edge l disconnects the graph
  - Used for N-1 contingency pre-screening
- `branch_indices::Vector{Tuple{Symbol, Int}}`: Maps incidence row → (model_name, instance_idx)
  - Includes ALL 2-terminal devices (Branch, Transformer, TCSC, etc.) via Terminal trait
  - `branch_indices[l]` = `(:Branch, 5)` means row l is Branch instance 5
  - `branch_indices[l]` = `(:Transformer, 3)` means row l is Transformer instance 3

# Example
```julia
topo = analyze_topology(sys)

# Find all bridge outages (any 2-terminal device type)
for l in findall(topo.is_bridge)
    model_name, instance_idx = topo.branch_indices[l]
    println("Bridge: \$model_name instance \$instance_idx")
    # model_name ∈ {:Branch, :Transformer, :TCSC, ...}
end

# Incidence matrix for PTDF computation
A = topo.incidence  # L×N sparse matrix
```
"""
struct TopologyInfo
    n_islands::Int
    island_map::Vector{Int}
    slack_buses::Vector{Int}
    adjacency::SparseMatrixCSC{Bool, Int}
    components::Vector{Vector{Int}}
    incidence::Union{Nothing, SparseMatrixCSC{Int8, Int}}
    is_bridge::BitVector
    branch_indices::Vector{Tuple{Symbol, Int}}
end

"""
    analyze_topology(sys::SystemModel; vertex_model::Symbol=:Bus) -> TopologyInfo

Analyze power system topology to identify islands and connectivity.

Uses the Terminal API (get_terminals) to generically discover connectivity
from all models in the system. The vertex model must carry a bus type code `ide`
(the slack vertex of an island has `ide == BUS_SLACK`).

# Algorithm
1. Build adjacency matrix by querying all models for their terminals
2. Filter by device status (only active devices)
3. Run connected components algorithm (BFS)
4. Identify slack/reference vertex per island (or warn if missing)

# Arguments
- `sys::SystemModel`: System to analyze
- `vertex_model::Symbol=:Bus`: Vertex model to analyze

# Returns
`TopologyInfo` containing islands, connectivity, and slack buses

# Example
```julia
sys = SystemModel(case)

# Analyze AC topology
topo = analyze_topology(sys)  # Defaults to :Bus
println("Found \$(topo.n_islands) AC islands")
```

# Caching
```julia
# Cache in SystemModel for reuse
sys.topology = analyze_topology(sys)

# Or use analyze_topology! to cache automatically
analyze_topology!(sys)  # Sets sys.topology
```
"""
function analyze_topology(sys::SystemModel; vertex_model::Symbol=:Bus, build_incidence::Bool=true)
    # Validate vertex model exists
    if !haskey(sys.models, vertex_model)
        error("Cannot analyze topology: no $vertex_model model found in system")
    end
    n_vertices = length(sys.models[vertex_model])

    # Build adjacency matrix using generic Terminal API
    adj = build_adjacency_matrix(sys, vertex_model)

    # Find connected components
    island_map, n_islands = _find_connected_components(adj, n_vertices)
    components = _group_by_island(island_map, n_islands)

    # Identify slack/reference vertices per island
    slack_vertices = _find_slack_per_island(sys.models[vertex_model], island_map, n_islands)

    # Build incidence matrix and find bridges
    if build_incidence
        incidence, branch_indices = build_incidence_matrix(sys, vertex_model)
        is_bridge = _find_bridges(incidence, n_vertices)
        return TopologyInfo(n_islands, island_map, slack_vertices, adj, components,
                           incidence, is_bridge, branch_indices)
    else
        # Empty bridge info when no incidence matrix
        return TopologyInfo(n_islands, island_map, slack_vertices, adj, components,
                           nothing, BitVector(), Tuple{Symbol, Int}[])
    end
end

"""
    analyze_topology!(sys::SystemModel) -> TopologyInfo

Analyze topology and cache result in `sys.topology`.

If topology is already cached, returns the cached result without re-analyzing.
To force re-analysis, set `sys.topology = nothing` first.

Returns the TopologyInfo and stores it in the SystemModel.

# Example
```julia
sys = SystemModel(case)
topo = analyze_topology!(sys)
# Now sys.topology contains the cached result

# Calling again returns same cached object
topo2 = analyze_topology!(sys)
@assert topo2 === topo  # Same object
```
"""
function analyze_topology!(sys::SystemModel)
    # Return cached result if available
    if !isnothing(sys.topology)
        return sys.topology
    end

    # Analyze and cache
    sys.topology = analyze_topology(sys)
    return sys.topology
end

#==============================================================================
# Generic Adjacency Matrix Construction
==============================================================================#

"""
    build_adjacency_matrix(sys::SystemModel, vertex_model::Symbol) -> SparseMatrixCSC{Bool}

Build adjacency matrix by querying all models for their terminals.

Generic implementation that:
- Iterates over ALL models in the system
- Queries each model type for its terminals via `get_terminals()`
- Filters terminals for the requested vertex type
- Adds edges based on terminal declarations
- Handles device status (active/inactive)

No hardcoded knowledge of Branch, Transformer, or any specific model type.

# Arguments
- `sys::SystemModel`: System containing models
- `vertex_model::Symbol`: Vertex type to build graph for (`:Bus`, `:Node`, etc.)

# Returns
Symmetric sparse adjacency matrix (undirected graph)

# Example
```julia
# Build AC bus connectivity
adj = build_adjacency_matrix(sys, :Bus)

# Future: Build DC node connectivity
adj_dc = build_adjacency_matrix(sys, :Node)
```
"""
function build_adjacency_matrix(sys::SystemModel, vertex_model::Symbol)
    n_vertices = length(sys.models[vertex_model])
    adj = spzeros(Bool, n_vertices, n_vertices)

    # Iterate over ALL models in the system
    for (model_name, model_vec) in pairs(sys.models)
        # Skip empty models
        if length(model_vec) == 0
            continue
        end

        # Get terminal declaration for this model type
        model_type = get_base_model_type(model_vec)
        terminals = get_terminals(model_type)

        # Skip models with no terminals (e.g., loads, generators)
        if isempty(terminals)
            continue
        end

        # Filter terminals for the requested vertex type
        relevant_terminals = filter(t -> t.vertex_model == vertex_model, terminals)

        # Skip if this model doesn't connect vertices of the requested type
        if isempty(relevant_terminals)
            continue
        end

        # Add edges for this model
        _add_model_connections!(adj, model_name, model_vec, relevant_terminals, sys.order)
    end

    # Make symmetric (undirected graph)
    return adj .| adj'
end

"""
    terminal_vertex_indices(order::IndexOrder, model_name::Symbol, terminal) -> Vector{Int}

Internal index of the vertex that each instance of `model_name` connects to through
`terminal` (a direct terminal such as `Terminal(:Bus, :i)`), from `foreign_refs`.
"""
terminal_vertex_indices(order::IndexOrder, model_name::Symbol, terminal) =
    order.foreign_refs[model_name][normalize_ref_key(terminal.ref_key)].target_idx

"""
    terminal_vertex_indices(order, model_name, terminal, vertex_model) -> Vector{Int}

Per instance of `model_name`, the internal index of the `vertex_model` instance it
attaches to through `terminal`, or 0 if it attaches elsewhere. A direct terminal on
`vertex_model` (e.g. `Terminal(:Bus, :i)`) attaches every instance; a group terminal
(e.g. `Terminal(:bus_group, :BUS)`) attaches the instances whose group target is a
`vertex_model` instance, so a generator on a phasor bus counts as a device of that Bus.
"""
function terminal_vertex_indices(order::IndexOrder, model_name::Symbol, terminal, vertex_model::Symbol)
    ref = order.foreign_refs[model_name][normalize_ref_key(terminal.ref_key)]
    terminal.vertex_model === vertex_model && return ref.target_idx
    return [ref.target_type[i] === vertex_model ? ref.target_idx[i] : 0 for i in eachindex(ref.target_idx)]
end

# Terminals through which instances of `model_type` can attach to `vertex_model`:
# direct terminals on it, and group terminals (resolved per instance)
_vertex_terminals(model_type, vertex_model::Symbol) =
    [t for t in get_terminals(model_type) if t.vertex_model === vertex_model || is_group_reference(t.vertex_model)]

"""
    in_service(model_vec::StructArray) -> BitVector

Whether each instance is in service: its status field (per `get_status_field`) is
nonzero. Models without a status field are always in service. The status is
compared where it lives (CPU or GPU) and the result is a host `BitVector`.
"""
function in_service(model_vec::StructArray)
    status_field = get_status_field(model_vec)
    isnothing(status_field) && return trues(length(model_vec))
    return BitVector(Array(.!iszero.(getproperty(model_vec, status_field))))
end

"""
    _add_model_connections!(adj, model_name, model_vec, terminals, order)

Add edges from a model's instances based on its terminal declarations.

Generic implementation that works for any model type:
- Handles 2-terminal (Branch), 3-terminal (Transformer3W), N-terminal devices
- Respects device status (st, stat, status fields)
- Vertex indices come from the pre-computed `foreign_refs`

# Arguments
- `adj::SparseMatrixCSC{Bool}`: Adjacency matrix to modify
- `model_name::Symbol`: Model name in the system
- `model_vec`: StructArray of model instances
- `terminals::Vector{Terminal}`: Filtered terminals for this vertex type
- `order::IndexOrder`: Index order holding `foreign_refs`
"""
function _add_model_connections!(adj, model_name, model_vec, terminals, order)
    active = in_service(model_vec)
    terminal_vertices = [terminal_vertex_indices(order, model_name, t) for t in terminals]

    for i in 1:length(model_vec)
        active[i] || continue

        # Vertex indices for all terminals of this device instance
        vertex_indices = [tv[i] for tv in terminal_vertices]

        # Add edges between all pairs of vertices
        # For 2-terminal: adds one edge (i,j)
        # For 3-terminal: adds three edges (i,j), (i,k), (j,k)
        # For N-terminal: adds N*(N-1)/2 edges
        for (idx1, idx2) in Iterators.product(vertex_indices, vertex_indices)
            if idx1 != idx2
                adj[idx1, idx2] = true
            end
        end
    end
end


"""
    compute_amd_ordering(pattern::SparseMatrixCSC{Bool}) -> Vector{Int}

Compute AMD (Approximate Minimum Degree) ordering for sparse matrix.

Takes a Jacobian sparsity pattern and returns a permutation vector that
minimizes fill-in during sparse factorization.

# Arguments
- `pattern`: Boolean sparsity pattern (typically from build_jacobian_sparsity_pattern)

# Returns
- `perm::Vector{Int}`: Permutation where new position `k` holds original variable
  `perm[k]` (the reordered matrix is `pattern[perm, perm]`)

# Example
```julia
perm = compute_amd_ordering(pattern)
prob = ACPowerFlowProblem(sys; addressing=CustomPermutation(perm))
```

# Notes
- AMD requires a symmetric pattern - this function symmetrizes if needed
"""
function compute_amd_ordering(pattern::SparseMatrixCSC{Bool})
    # Make symmetric for AMD (it expects symmetric pattern)
    sym_pattern = pattern .| pattern'

    # Compute AMD permutation (AMD is imported at module level)
    perm = AMD.amd(sym_pattern)

    return perm
end

#==============================================================================
# Incidence Matrix Construction
==============================================================================#

"""
    build_incidence_matrix(sys::SystemModel, vertex_model::Symbol=:Bus) -> SparseMatrixCSC{Int8, Int}

Build the branch-bus incidence matrix A (L×N) for all active 2-terminal devices.

# Structure
- L rows = number of active 2-terminal devices (branches, transformers, etc.)
- N columns = number of vertices (buses)
- Row l has +1 at from-bus column, -1 at to-bus column

# Usage
The incidence matrix enables efficient computation of:
- Bus susceptance matrix: `B = A' D A` where `D = diag(1/x_l)`
- PTDF: `PTDF = D · A · B⁻¹` (via sparse solves, never explicit inverse)
- LODF: Derived from PTDF
- Adjacency matrix: `adj = (A' * A) .!= 0` (derived, not primary)

# Returns
- `SparseMatrixCSC{Int8, Int}`: Incidence matrix with values {-1, 0, +1}
- Also returns `branch_indices::Vector{Tuple{Symbol, Int}}` mapping row l to (model_name, instance_idx)

# Example
```julia
A, branch_idx = build_incidence_matrix(sys)
# A[5, :] has +1 at from-bus, -1 at to-bus for branch 5
# branch_idx[5] = (:Branch, 3) means row 5 corresponds to Branch instance 3
```
"""
function build_incidence_matrix(sys::SystemModel, vertex_model::Symbol=:Bus)
    n_vertices = length(sys.models[vertex_model])

    # First pass: count active 2-terminal devices and collect their endpoints
    endpoints = Vector{Tuple{Symbol, Int, Int, Int}}()  # (model_name, instance_idx, from_idx, to_idx)

    for (model_name, model_vec) in pairs(sys.models)
        length(model_vec) == 0 && continue

        # Get terminal declaration for this model type
        model_type = get_base_model_type(model_vec)
        terminals = get_terminals(model_type)

        # Only process 2-terminal devices
        relevant_terminals = filter(t -> t.vertex_model == vertex_model, terminals)
        length(relevant_terminals) != 2 && continue

        active = in_service(model_vec)
        from_vertices = terminal_vertex_indices(sys.order, model_name, relevant_terminals[1])
        to_vertices = terminal_vertex_indices(sys.order, model_name, relevant_terminals[2])

        for i in 1:length(model_vec)
            active[i] || continue
            push!(endpoints, (model_name, i, from_vertices[i], to_vertices[i]))
        end
    end

    # Build sparse matrix using COO format
    n_branches = length(endpoints)
    if n_branches == 0
        # Return empty matrix if no branches
        branch_indices = Tuple{Symbol, Int}[]
        return spzeros(Int8, 0, n_vertices), branch_indices
    end

    # COO triplets: 2 entries per branch (+1 at from, -1 at to)
    I_row = Vector{Int}(undef, 2 * n_branches)
    J_col = Vector{Int}(undef, 2 * n_branches)
    V_val = Vector{Int8}(undef, 2 * n_branches)

    branch_indices = Vector{Tuple{Symbol, Int}}(undef, n_branches)

    for (l, (model_name, instance_idx, from_idx, to_idx)) in enumerate(endpoints)
        # +1 at from-bus
        I_row[2l - 1] = l
        J_col[2l - 1] = from_idx
        V_val[2l - 1] = Int8(1)

        # -1 at to-bus
        I_row[2l] = l
        J_col[2l] = to_idx
        V_val[2l] = Int8(-1)

        branch_indices[l] = (model_name, instance_idx)
    end

    A = sparse(I_row, J_col, V_val, n_branches, n_vertices)
    return A, branch_indices
end

#==============================================================================
# Connected Components Algorithm
==============================================================================#

"""
    _find_connected_components(adj::SparseMatrixCSC, n::Int; excluded_edge=(0, 0)) -> (island_map, n_islands)

Find connected components using breadth-first search (BFS).

# Arguments
- `adj::SparseMatrixCSC`: Adjacency matrix (assumed symmetric/undirected)
- `n::Int`: Number of vertices
- `excluded_edge`: `(from, to)` vertex pair treated as disconnected (in both
  directions); the default `(0, 0)` excludes nothing

# Returns
- `island_map::Vector{Int}`: Maps vertex internal index → island number (1-based)
- `n_islands::Int`: Number of islands found

Island numbers follow the lowest vertex index in each island. Neighbors are read
directly from the CSC arrays and the queue is preallocated, so the search makes
no per-vertex allocations.
"""
function _find_connected_components(adj::SparseMatrixCSC, n::Int; excluded_edge::Tuple{Int, Int}=(0, 0))
    excluded_from, excluded_to = excluded_edge
    visited = falses(n)
    island_map = zeros(Int, n)
    island_num = 0

    colptr = adj.colptr
    rowval = adj.rowval
    queue = Vector{Int}(undef, n)

    for start in 1:n
        visited[start] && continue

        island_num += 1
        visited[start] = true
        island_map[start] = island_num

        queue[1] = start
        head = 1
        tail = 1

        while head <= tail
            node = queue[head]
            head += 1

            for k in colptr[node]:(colptr[node+1]-1)
                neighbor = rowval[k]

                # Skip the excluded edge (in either direction)
                if (node == excluded_from && neighbor == excluded_to) ||
                   (node == excluded_to && neighbor == excluded_from)
                    continue
                end

                if !visited[neighbor]
                    visited[neighbor] = true
                    island_map[neighbor] = island_num
                    tail += 1
                    queue[tail] = neighbor
                end
            end
        end
    end

    return island_map, island_num
end

"""
    _group_by_island(island_map::Vector{Int}, n_islands::Int) -> Vector{Vector{Int}}

Group vertex internal indices by island number.

# Returns
Vector of vectors, where `components[i]` contains all vertex internal indices
belonging to island `i`.

# Example
```julia
island_map = [1, 1, 2, 1, 2]  # 5 vertices, 2 islands
components = _group_by_island(island_map, 2)
# → [[1, 2, 4], [3, 5]]  # Island 1 has vertices {1,2,4}, island 2 has {3,5}
```
"""
function _group_by_island(island_map::Vector{Int}, n_islands::Int)
    components = [Int[] for _ in 1:n_islands]
    for (vertex_idx, island_idx) in enumerate(island_map)
        push!(components[island_idx], vertex_idx)
    end
    return components
end

#==============================================================================
# Bridge (Cut-Edge) Detection
==============================================================================#

"""
    is_bridge_device(topo::TopologyInfo, sys::SystemModel, model_name::Symbol, key) -> Bool

Check if a 2-terminal device outage creates an island.

Works for any model type with 2 terminals (Branch, Transformer2W, TCSC, etc.).
Throws `KeyError` if model or device key not found.

# Arguments
- `topo`: TopologyInfo from `analyze_topology(sys)`
- `sys`: SystemModel (for key lookup)
- `model_name`: Model type symbol (`:Branch`, `:Transformer2W`, etc.)
- `key`: External key tuple (e.g., `(from_bus, to_bus, ckt)` for Branch)

# Returns
`true` if device outage creates an island, `false` otherwise.

# Throws
- `KeyError`: If model type or device key not found

# Example
```julia
topo = analyze_topology(sys)

# Check Branch outage
is_bridge_device(topo, sys, :Branch, (1, 2, "1"))

# Check Transformer outage
is_bridge_device(topo, sys, :Transformer2W, (100, 101, "1"))
```
"""
function is_bridge_device(topo::TopologyInfo, sys::SystemModel, model_name::Symbol, key)
    l = _find_edge_row(topo, model_name, get_index_by_key(sys, model_name, key))

    # Device not in incidence matrix (inactive or not 2-terminal)
    isnothing(l) && return false
    return topo.is_bridge[l]
end

"""
    _find_edge_row(topo, model_name, instance_idx) -> Union{Int, Nothing}

Incidence-matrix row of a 2-terminal device instance, or `nothing` when the
device is not in the incidence matrix (inactive or not 2-terminal). O(L) search.
"""
_find_edge_row(topo::TopologyInfo, model_name::Symbol, instance_idx::Int) =
    findfirst(==((model_name, instance_idx)), topo.branch_indices)

# Partition result for an outage that does not change connectivity
_unsplit_partition(topo::TopologyInfo) = (
    is_split = false,
    new_islands = topo.components,
    n_new_islands = topo.n_islands,
    affected_islands = Int[],
    smaller_affected = 0
)

"""
    get_partition_after_outage(topo, edge_idx) -> NamedTuple

Compute connected components after removing an edge.

For contingency screening, identifies which buses become isolated or form
new islands when an edge is disconnected. Works correctly even if the
original system has multiple islands.

# Arguments
- `topo::TopologyInfo`: Topology from `analyze_topology(sys; build_incidence=true)`
- `edge_idx::Int`: Row index in incidence matrix (1-based)

# Returns
NamedTuple with:
- `is_split::Bool`: True if removal increases the number of islands
- `new_islands::Vector{Vector{Int}}`: All islands after edge removal
- `n_new_islands::Int`: Number of islands after removal
- `affected_islands::Vector{Int}`: Indices of islands that were split (in new_islands)
- `smaller_affected::Int`: Index of smallest affected island, or 0 if not split

# Example
```julia
topo = analyze_topology(sys; build_incidence=true)

result = get_partition_after_outage(topo, 5)
if result.is_split
    # Get the smallest newly-created island
    isolated_buses = result.new_islands[result.smaller_affected]
    println("Removing edge 5 isolates buses: \$isolated_buses")
end
```

# Algorithm
Recomputes connected components on the adjacency matrix with the specified
edge virtually removed. O(V+E) time complexity.
"""
function get_partition_after_outage(topo::TopologyInfo, edge_idx::Int)
    # Validate input
    if isnothing(topo.incidence)
        error("TopologyInfo has no incidence matrix. Use analyze_topology(sys; build_incidence=true)")
    end

    n_edges = size(topo.incidence, 1)
    if edge_idx < 1 || edge_idx > n_edges
        error("Edge index $edge_idx out of range [1, $n_edges]")
    end

    # Non-bridge edges don't change connectivity
    topo.is_bridge[edge_idx] || return _unsplit_partition(topo)

    # Extract endpoints of this edge from incidence matrix
    n_vertices = size(topo.incidence, 2)
    from_idx, to_idx = extract_from_to_idx(topo.incidence)
    from_bus, to_bus = from_idx[edge_idx], to_idx[edge_idx]

    # Recompute connected components with edge virtually removed
    island_map, n_islands = _find_connected_components(
        topo.adjacency, n_vertices; excluded_edge=(from_bus, to_bus)
    )
    new_islands = _group_by_island(island_map, n_islands)

    # Find which islands were affected (the two endpoints will be in different islands now)
    island_from = island_map[from_bus]
    island_to = island_map[to_bus]
    affected = [island_from, island_to]

    # Find smaller of the affected islands
    smaller = length(new_islands[island_from]) <= length(new_islands[island_to]) ? island_from : island_to

    return (
        is_split = true,
        new_islands = new_islands,
        n_new_islands = n_islands,
        affected_islands = affected,
        smaller_affected = smaller
    )
end

"""
    get_partition_after_outage(topo, sys, model_name, key) -> NamedTuple

Compute partition after removing a specific device by its external key.

Convenience overload that looks up the edge index from model/key.

# Example
```julia
topo = analyze_topology(sys; build_incidence=true)
result = get_partition_after_outage(topo, sys, :Branch, (1, 5, "1"))
```
"""
function get_partition_after_outage(
    topo::TopologyInfo,
    sys::SystemModel,
    model_name::Symbol,
    key
)
    l = _find_edge_row(topo, model_name, get_index_by_key(sys, model_name, key))

    # Device not in incidence matrix (inactive or not 2-terminal)
    isnothing(l) && return _unsplit_partition(topo)
    return get_partition_after_outage(topo, l)
end

"""
    extract_from_to_idx(A::SparseMatrixCSC{Int8}) -> (from_indices, to_indices)

Endpoints of every edge of an L×N incidence matrix in a single O(nnz) pass:
`from_indices[l]` is the vertex with +1 in row `l`, `to_indices[l]` the vertex with -1.
"""
function extract_from_to_idx(A::SparseMatrixCSC{Int8})
    n_edges = size(A, 1)
    from_indices = zeros(Int, n_edges)
    to_indices = zeros(Int, n_edges)

    # CSC is column (vertex) major, so walk each vertex's column of edges
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    for vertex in 1:size(A, 2)
        for k in colptr[vertex]:(colptr[vertex+1]-1)
            if nzval[k] == 1
                from_indices[rowval[k]] = vertex
            else  # -1
                to_indices[rowval[k]] = vertex
            end
        end
    end

    return from_indices, to_indices
end

#==============================================================================
# Reverse Terminal Lookup (for contingency cascading)
==============================================================================#

"""
    find_devices_on_vertex(sys, vertex_model, vertex_idx) -> Dict{Symbol, Vector{Int}}

Find all model instances with any terminal connecting to a specific vertex.

Uses the Terminal API generically for any vertex model that devices reference
directly (e.g. `:Bus`). Group names such as `:generators` are not
vertex models and raise an error. Useful for bridge contingencies (all devices
on an isolated bus).

# Arguments
- `sys::SystemModel`: System with models and index order
- `vertex_model::Symbol`: Vertex model in the system (e.g. `:Bus`)
- `vertex_idx::Int`: Internal index of the vertex

# Returns
`Dict{Symbol, Vector{Int}}` mapping model names to instance indices.

# Example
```julia
# Find all devices on bus 5 (internal index)
devices = find_devices_on_vertex(sys, :Bus, 5)
# => Dict(:PQ => [1, 3], :PV => [2], :Branch => [7, 12], ...)
```
"""
function find_devices_on_vertex(
    sys::SystemModel,
    vertex_model::Symbol,
    vertex_idx::Int
)
    result = Dict{Symbol, Vector{Int}}()

    haskey(sys.models, vertex_model) ||
        error("Vertex model :$vertex_model is not in the system (group names are not vertex models)")

    for (model_name, model_vec) in pairs(sys.models)
        isempty(model_vec) && continue
        # Terminals that attach to vertex_model, directly or through a group
        matching_terminals = _vertex_terminals(get_base_model_type(model_vec), vertex_model)
        isempty(matching_terminals) && continue

        terminal_vertices = [terminal_vertex_indices(sys.order, model_name, t, vertex_model) for t in matching_terminals]
        indices = [i for i in 1:length(model_vec) if any(tv -> tv[i] == vertex_idx, terminal_vertices)]

        if !isempty(indices)
            result[model_name] = indices
        end
    end

    return result
end

"""
    build_reverse_terminal_map(sys, vertex_model) -> Dict{Int, Dict{Symbol, Vector{Int}}}

Pre-compute full reverse lookup: vertex_idx => {model_name => [instance indices]}.

Used by GPU for O(1) lookup during batched contingency application.
Builds the map once during setup, then use for fast indexing per scenario.

# Arguments
- `sys::SystemModel`: System with models and index order
- `vertex_model::Symbol`: Vertex model in the system (e.g. `:Bus`); group names error

# Returns
`Dict{Int, Dict{Symbol, Vector{Int}}}` mapping vertex indices to connected devices.

# Example
```julia
bus_map = build_reverse_terminal_map(sys, :Bus)
bus_map[5][:PQ]  # => [1, 3]  PQ instances on bus 5
```
"""
function build_reverse_terminal_map(
    sys::SystemModel,
    vertex_model::Symbol
)
    result = Dict{Int, Dict{Symbol, Vector{Int}}}()

    haskey(sys.models, vertex_model) ||
        error("Vertex model :$vertex_model is not in the system (group names are not vertex models)")

    n_vertices = length(sys.models[vertex_model])

    # Initialize empty dict for each vertex
    for v in 1:n_vertices
        result[v] = Dict{Symbol, Vector{Int}}()
    end

    # Iterate all models and populate the map
    for (model_name, model_vec) in pairs(sys.models)
        isempty(model_vec) && continue
        # Terminals that attach to vertex_model, directly or through a group
        matching_terminals = _vertex_terminals(get_base_model_type(model_vec), vertex_model)
        isempty(matching_terminals) && continue

        terminal_vertices = [terminal_vertex_indices(sys.order, model_name, t, vertex_model) for t in matching_terminals]

        for i in 1:length(model_vec)
            # Add each instance once per distinct vertex it touches (0: attached elsewhere)
            for v_idx in unique(tv[i] for tv in terminal_vertices)
                v_idx == 0 && continue
                push!(get!(result[v_idx], model_name, Int[]), i)
            end
        end
    end

    return result
end

"""
    disable_devices!(sys, devices::Dict{Symbol, Vector{Int}})

Set status=0/false for specified device instances.

Uses `get_status_field` to find the correct status field per model type.
Models without a status field are always in service and are left as they are.
Errors if a model name is not in the system.

# Arguments
- `sys::SystemModel`: System with models to modify
- `devices::Dict{Symbol, Vector{Int}}`: Model names to instance indices

# Example
```julia
devices = find_devices_on_vertex(sys, :Bus, isolated_bus_idx)
disable_devices!(sys, devices)
```
"""
function disable_devices!(sys::SystemModel, devices::Dict{Symbol, Vector{Int}})
    for (model_name, indices) in devices
        model_vec = sys.models[model_name]
        model_type = get_base_model_type(model_vec)
        status_field = get_status_field(model_type)
        isnothing(status_field) && continue

        status_arr = getproperty(model_vec, status_field)
        for i in indices
            if status_arr[i] isa Bool
                status_arr[i] = false
            else
                status_arr[i] = zero(eltype(status_arr))
            end
        end
    end
end

"""
    _find_bridges(incidence::SparseMatrixCSC{Int8}, n_vertices::Int) -> BitVector

Find bridge edges using Tarjan's algorithm on the graph encoded by the incidence matrix.

A bridge is an edge whose removal disconnects the graph. For N-1 contingency
screening, bridges identify outages that create islands.

Uses iterative DFS to avoid stack overflow on large systems.

# Arguments
- `incidence`: L×N incidence matrix from `build_incidence_matrix`
- `n_vertices`: Number of vertices (buses)

# Returns
BitVector of length L where `result[l]` = true if edge l is a bridge.
"""
function _find_bridges(incidence::SparseMatrixCSC{Int8}, n_vertices::Int)
    n_edges = size(incidence, 1)

    if n_edges == 0
        return BitVector()
    end

    from_idx, to_idx = extract_from_to_idx(incidence)

    # Detect parallel edges (same bus pair → never bridges)
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for l in 1:n_edges
        pair = minmax(from_idx[l], to_idx[l])
        edge_count[pair] = get(edge_count, pair, 0) + 1
    end

    has_parallel = falses(n_edges)
    for l in 1:n_edges
        pair = minmax(from_idx[l], to_idx[l])
        if edge_count[pair] > 1
            has_parallel[l] = true
        end
    end

    # Build adjacency list: adj_list[v] = [(neighbor, edge_idx), ...]
    adj_list = [Vector{Tuple{Int, Int}}() for _ in 1:n_vertices]
    for l in 1:n_edges
        i, j = from_idx[l], to_idx[l]
        push!(adj_list[i], (j, l))
        push!(adj_list[j], (i, l))
    end

    # Tarjan's bridge algorithm with iterative DFS
    is_bridge = falses(n_edges)
    visited = falses(n_vertices)
    disc = zeros(Int, n_vertices)
    low = zeros(Int, n_vertices)
    parent_edge = zeros(Int, n_vertices)

    # Stack for iterative DFS: (node, neighbor_index, phase)
    stack = Vector{Tuple{Int, Int, Int}}()
    timer = 0

    for start in 1:n_vertices
        visited[start] && continue

        timer += 1
        visited[start] = true
        disc[start] = timer
        low[start] = timer

        push!(stack, (start, 1, 0))

        while !isempty(stack)
            u, idx, phase = pop!(stack)
            neighbors = adj_list[u]

            if phase == 1 && idx > 1
                # Returning from child - update low and check bridge
                v, edge_idx = neighbors[idx - 1]
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]
                    is_bridge[edge_idx] = true
                end
            end

            # Process remaining neighbors
            while idx <= length(neighbors)
                v, edge_idx = neighbors[idx]

                if !visited[v]
                    # Tree edge - visit child
                    visited[v] = true
                    timer += 1
                    disc[v] = timer
                    low[v] = timer
                    parent_edge[v] = edge_idx

                    # Push current node to resume later, then child
                    push!(stack, (u, idx + 1, 1))
                    push!(stack, (v, 1, 0))
                    break
                elseif edge_idx != parent_edge[u]
                    # Back edge - update low
                    low[u] = min(low[u], disc[v])
                end
                idx += 1
            end
        end
    end

    # Parallel edges are never bridges
    for l in 1:n_edges
        if has_parallel[l]
            is_bridge[l] = false
        end
    end

    return is_bridge
end

#==============================================================================
# Slack/Reference Vertex Identification
==============================================================================#

"""
    _find_slack_per_island(vertices, island_map, n_islands) -> Vector{Int}

Identify one slack/reference vertex per island.

Searches for vertices with type code 3 (slack bus in PSS/E convention).
Issues warnings for:
- Islands with multiple slack vertices
- Islands with no slack vertex

# Arguments
- `vertices`: StructArray of vertex model instances (Bus, Node, etc.)
- `island_map::Vector{Int}`: Maps vertex index → island number
- `n_islands::Int`: Total number of islands

# Returns
Vector of slack vertex internal indices (one per island).
Entry is `0` if island has no slack vertex.

# Example
```julia
slack_buses = _find_slack_per_island(sys.models.Bus, island_map, n_islands)
# slack_buses[1] = 5  # Island 1's slack is vertex 5
# slack_buses[2] = 0  # Island 2 has no slack (warning issued)
```
"""
function _find_slack_per_island(vertices, island_map, n_islands)
    slack_vertices = zeros(Int, n_islands)
    is_slack = Array(vertices.ide .== BUS_SLACK)

    for i in eachindex(vertices)
        # Check for slack bus
        if is_slack[i]
            island_idx = island_map[i]
            if slack_vertices[island_idx] != 0
                @warn "Island $island_idx has multiple slack vertices: " *
                      "$(slack_vertices[island_idx]) and $i"
            else
                slack_vertices[island_idx] = i
            end
        end
    end

    # Warn about islands without slack
    for i in 1:n_islands
        if slack_vertices[i] == 0
            @warn "Island $i has no slack vertex - power flow may not converge"
        end
    end

    return slack_vertices
end
