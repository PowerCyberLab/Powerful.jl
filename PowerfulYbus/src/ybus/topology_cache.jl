# ============================================================================
# Topology Cache for Pre-computed Y Matrices
# ============================================================================
#
# Pre-computes the Y matrix of each network topology.
#
# ============================================================================

"""
    TopologyKey

Identifier for a specific network topology.

# Fields
- `name::Symbol`: Human-readable name (e.g., :prefault, :fault_on, :postfault)
- `lines_out::Vector{Int}`: Branch instance indices that are out of service
- `faults::Vector{Tuple{Int, ComplexF64}}`: Active faults as (internal bus index, admittance)

# Usage
```julia
# Pre-fault topology
prefault = TopologyKey(:prefault)

# Fault at internal bus index 3
fault_on = TopologyKey(:fault_on; faults=[(3, 1e7im)])

# Post-fault with line 5 removed
postfault = TopologyKey(:postfault; lines_out=[5])
```
"""
struct TopologyKey
    name::Symbol
    lines_out::Vector{Int}
    faults::Vector{Tuple{Int, Complex{Float64}}}  # (internal bus index, y_fault)

    TopologyKey(name::Symbol; lines_out=Int[], faults=Tuple{Int, Complex{Float64}}[]) =
        new(name, lines_out, faults)
end

Base.:(==)(a::TopologyKey, b::TopologyKey) =
    a.name == b.name && a.lines_out == b.lines_out && a.faults == b.faults

Base.hash(k::TopologyKey, h::UInt) = hash(k.name, hash(k.lines_out, hash(k.faults, h)))

# ============================================================================
# Applying a topology to the full Y matrix
# ============================================================================

"""
    _branch_admittance(sys, line_idx, T) -> (from, to, y_series, y_shunt)

Internal from/to bus indices, series admittance 1/(r+jx) and per-end shunt
admittance jb/2 of Branch instance `line_idx`.
"""
function _branch_admittance(sys, line_idx::Int, ::Type{T}) where T
    branch = sys.models[:Branch]
    from = get_foreign_indices(sys, :Branch, :i)[line_idx]
    to = get_foreign_indices(sys, :Branch, :j)[line_idx]
    y_series = one(Complex{T}) / complex(T(branch.r[line_idx]), T(branch.x[line_idx]))
    y_shunt = im * T(branch.b[line_idx]) / 2
    return from, to, y_series, y_shunt
end

"""
    _fault_admittance(sys, fault_idx, T) -> (bus, y_fault)

Internal bus index and shunt admittance gf + j*bf of Fault instance `fault_idx`.
"""
function _fault_admittance(sys, fault_idx::Int, ::Type{T}) where T
    fault = sys.models[:Fault]
    bus = get_foreign_indices(sys, :Fault, :i)[fault_idx]
    gf, bf = fault_admittance(fault.rf[fault_idx], fault.xf[fault_idx])
    return bus, complex(T(gf), T(bf))
end

"""
    _remove_branch_from_ymatrix!(Y, from_int, to_int, y_series, y_shunt)

Remove a branch contribution (trip/open): subtract (y_series + y_shunt) from both
diagonals and add y_series back on the off-diagonals.
"""
@inline function _remove_branch_from_ymatrix!(Y::SparseMatrixCSC{Complex{T}}, from_int::Int, to_int::Int,
                                              y_series::Complex{T}, y_shunt::Complex{T}) where T
    Y[from_int, from_int] -= (y_series + y_shunt)
    Y[to_int, to_int] -= (y_series + y_shunt)
    Y[from_int, to_int] += y_series
    Y[to_int, from_int] += y_series
    nothing
end

"""
    apply_topology(base::FullYMatrix{T}, sys, key::TopologyKey) -> FullYMatrix{T}

Full Y matrix for topology `key`: the base matrix with the branches in
`key.lines_out` removed and the fault admittances in `key.faults` added.
"""
function apply_topology(base::FullYMatrix{T}, sys, key::TopologyKey) where T
    Y = copy(base.Y)
    for line_idx in key.lines_out
        _remove_branch_from_ymatrix!(Y, _branch_admittance(sys, line_idx, T)...)
    end
    for (bus, y_fault) in key.faults
        Y[bus, bus] += y_fault
    end
    return FullYMatrix{T}(Y, base.n_bus, base.gen_buses)
end

# ============================================================================
# Topology cache: one Y matrix per topology
# ============================================================================

"""
    TopologyCache{T, SYS, S}

Cache of network admittance matrices, one per topology.

# Fields
- `base::FullYMatrix{T}`: Pre-fault Y built at construction; every topology,
  including ones added later, is derived from it
- `topologies::Dict{TopologyKey, FullYMatrix{T}}`: Y matrix per topology
- `current::Base.RefValue{TopologyKey}`: Currently active topology
- `sys`: Reference to system model

# Usage
```julia
cache = TopologyCache(sys)
add_topology!(cache, TopologyKey(:fault_on; faults=[(3, 1e7im)]))
switch_topology!(cache, :fault_on)
Y = current_Y(cache).Y
```

Switching topologies (e.g. to apply or clear a fault) is a dictionary lookup rather
than a matrix rebuild, since every topology the simulation can reach is pre-computed
when the cache is built or when the fault or switching event is added to it.
"""
struct TopologyCache{T, SYS}
    base::FullYMatrix{T}
    topologies::Dict{TopologyKey, FullYMatrix{T}}
    current::Base.RefValue{TopologyKey}
    sys::SYS
end

"""
    TopologyCache(sys; events=nothing)

Build a topology cache holding the pre-fault topology, plus the `topology_key` of
each of `events` if given.

# Arguments
- `sys`: SystemModel containing network data
- `events`: Optional switching schedule defining topologies to pre-compute
"""
function TopologyCache(sys; events=nothing)
    base = build_full_ybus(sys)
    prefault_key = TopologyKey(:prefault)
    topologies = Dict{TopologyKey, typeof(base)}(prefault_key => apply_topology(base, sys, prefault_key))
    cache = TopologyCache(base, topologies, Ref(prefault_key), sys)
    if !isnothing(events)
        for event in events
            add_topology!(cache, event.topology_key)
        end
    end
    return cache
end

"""
    _store_topology!(cache, key, base, sys) -> FullYMatrix

Derive and store topology `key` from `base`, unless the cache already has it.
"""
function _store_topology!(cache::TopologyCache, key::TopologyKey, base::FullYMatrix, sys)
    haskey(cache.topologies, key) && return cache.topologies[key]
    return cache.topologies[key] = apply_topology(base, sys, key)
end

"""
    switch_topology!(cache::TopologyCache, key::TopologyKey)
    switch_topology!(cache::TopologyCache, name::Symbol)

Switch to a different pre-computed topology. Errors if it is not in the cache.
"""
function switch_topology!(cache::TopologyCache, key::TopologyKey)
    if !haskey(cache.topologies, key)
        error("Topology $key not found in cache. Available: $(keys(cache.topologies))")
    end
    cache.current[] = key
    return nothing
end

function switch_topology!(cache::TopologyCache, name::Symbol)
    for key in keys(cache.topologies)
        if key.name == name
            cache.current[] = key
            return nothing
        end
    end
    error("Topology with name :$name not found in cache")
end

"""
    current_Y(cache::TopologyCache) -> FullYMatrix

The Y matrix of the currently active topology.
"""
current_Y(cache::TopologyCache) = cache.topologies[cache.current[]]

"""
    add_topology!(cache::TopologyCache, key::TopologyKey) -> FullYMatrix

Add topology `key`, derived from the cache's base Y, and return its Y matrix.
"""
add_topology!(cache::TopologyCache, key::TopologyKey) = _store_topology!(cache, key, cache.base, cache.sys)

export TopologyKey, TopologyCache
export switch_topology!, current_Y, add_topology!
