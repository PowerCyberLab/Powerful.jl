"""
Outage (contingency) actions and their resolution to model instances.

Shared by contingency analysis on CPU status dicts and batched GPU screening, so
every path resolves an outage to the same `(model, index)` and de-energizes the
same islands.
"""

"""
    OutageAction

Abstract type for contingency actions (outages).
"""
abstract type OutageAction end

"""
    BranchOutage(from_bus, to_bus, circuit)

Represents a branch (line) outage.
"""
struct BranchOutage <: OutageAction
    from_bus::Int
    to_bus::Int
    circuit::String
end

"""
    GeneratorOutage(bus, unit_id)

Represents a generator (PV) outage.
"""
struct GeneratorOutage <: OutageAction
    bus::Int
    unit_id::String
end

"""
    Transformer2WOutage(from_bus, to_bus, circuit)

Represents a 2-winding transformer outage.
"""
struct Transformer2WOutage <: OutageAction
    from_bus::Int
    to_bus::Int
    circuit::String
end

"""
    Contingency(name, actions)

A named contingency consisting of one or more outage actions.

# Example
```julia
Contingency("LINE_1_2", [BranchOutage(1, 2, "1")])            # N-1 branch
Contingency("GEN_101", [GeneratorOutage(101, "1")])          # N-1 generator
Contingency("XFMR_1_2", [Transformer2WOutage(1, 2, "1")])     # N-1 transformer
Contingency("N2_TEST", [BranchOutage(1, 2, "1"), GeneratorOutage(101, "1")])
```
"""
struct Contingency
    name::String
    actions::Vector{<:OutageAction}
end

outage_model(::BranchOutage) = :Branch
outage_model(::GeneratorOutage) = :PV
outage_model(::Transformer2WOutage) = :Transformer2W

outage_key(a::BranchOutage) = (a.from_bus, a.to_bus, a.circuit)
outage_key(a::GeneratorOutage) = (a.bus, a.unit_id)
outage_key(a::Transformer2WOutage) = (a.from_bus, a.to_bus, a.circuit)

"""
    resolve_outage(sys::SystemModel, action::OutageAction) -> (model::Symbol, idx::Int)

Model and internal index of the instance taken out by `action`. Errors if the
instance does not exist.
"""
function resolve_outage(sys::SystemModel, action::OutageAction)
    model = outage_model(action)
    return (model, get_index_by_key(sys, model, outage_key(action)))
end

"""
    apply_contingency!(status, sys, cont::Contingency, topo, bus_device_map)
        -> Dict{Symbol, Vector{Int}}

Take the instances removed by the actions of `cont` out of service in `status` (a
status dict of the base case, e.g. from `copy_status_dict`), then de-energize the
islands this leaves without a slack bus (see [`deenergize_islands!`](@ref)).
Returns the instances switched off, keyed by model (`:Bus` for de-energized buses).
The system itself is not modified.

`topo` is `analyze_topology(sys; build_incidence=true)` and `bus_device_map` is
`build_reverse_terminal_map(sys, :Bus)`, both of the base case. The islands are only
searched when the network can split: when a bridge is removed, or two or more
branches or transformers.
"""
function apply_contingency!(status::AbstractDict, sys::SystemModel, cont::Contingency,
                            topo::TopologyInfo, bus_device_map)
    switched = Dict{Symbol, Vector{Int}}()
    n_edges_out = 0
    bridge_out = false
    for action in cont.actions
        model, idx = resolve_outage(sys, action)
        status[model][idx] || continue  # already out of service
        status[model][idx] = false
        push!(get!(switched, model, Int[]), idx)

        l = _find_edge_row(topo, model, idx)
        if l !== nothing
            n_edges_out += 1
            bridge_out |= topo.is_bridge[l]
        end
    end
    (bridge_out || n_edges_out >= 2) && deenergize_islands!(switched, status, topo, bus_device_map)
    return switched
end

"""
    deenergize_islands!(switched, status, topo, bus_device_map) -> switched

Take out of service, in `status`, every bus that the in-service branches and
transformers no longer connect to a slack bus, and every device on such a bus; these
islands have no reference and cannot be solved. The instances switched off are
appended to `switched` (keyed by model, `:Bus` for the buses).
"""
function deenergize_islands!(switched::Dict{Symbol, Vector{Int}}, status::AbstractDict,
                             topo::TopologyInfo, bus_device_map)
    incidence = topo.incidence
    isnothing(incidence) && error("TopologyInfo has no incidence matrix. Use analyze_topology(sys; build_incidence=true)")
    from_idx, to_idx = extract_from_to_idx(incidence)
    n_vertices = size(incidence, 2)

    # Breadth-first search from the slack buses over in-service edges
    energized = falses(n_vertices)
    queue = [b for b in topo.slack_buses if b > 0 && status[:Bus][b]]
    energized[queue] .= true
    while !isempty(queue)
        v = pop!(queue)
        for k in nzrange(incidence, v)
            l = rowvals(incidence)[k]
            model, idx = topo.branch_indices[l]
            status[model][idx] || continue
            u = from_idx[l] == v ? to_idx[l] : from_idx[l]
            if !energized[u]
                energized[u] = true
                push!(queue, u)
            end
        end
    end

    for v in 1:n_vertices
        (energized[v] || !status[:Bus][v]) && continue
        status[:Bus][v] = false
        push!(get!(switched, :Bus, Int[]), v)
        haskey(bus_device_map, v) || continue  # bus hosts no devices
        for (m, idxs) in bus_device_map[v], i in idxs
            status[m][i] || continue
            status[m][i] = false
            push!(get!(switched, m, Int[]), i)
        end
    end
    return switched
end
