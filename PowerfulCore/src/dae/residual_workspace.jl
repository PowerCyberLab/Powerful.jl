"""
Residual calculation workspace with zero-allocation design.

This module provides pre-allocated workspaces for efficient residual evaluation
in power system DAE models. Uses NamedTuples for type-stable, dictionary-like
access with zero runtime overhead.
"""

# ============================================================================
# Helper Functions for GPU Detection (Weak Dependency Support)
# ============================================================================

"""
    is_gpu_array_type(::Type) -> Bool

Whether an array type lives on a GPU. `false` here; a GPU package extends it for its
array types (so the GPU library need not be loaded to ask).
"""
is_gpu_array_type(::Type) = false

"""
    device_sparse(array_type, S::SparseMatrixCSC) -> AbstractSparseMatrix

`S` on the device of `array_type`: `S` itself for CPU arrays; a GPU package returns a
device sparse matrix for its array types, so a Jacobian of a GPU system is assembled and
factorized on the GPU.
"""
device_sparse(::Type{<:AbstractArray}, S::SparseMatrixCSC) = S


# ============================================================================
# Helper Functions for Constructor
# ============================================================================

"""
    _build_storage(::Type{T}, items, metadata, am, array_type) where T

Build address and value storage for model items (variables or residuals).

Returns two vectors of pairs: (addr_pairs, val_pairs)
- Addresses are transferred to GPU if array_type <: CuVector
- Values are allocated with the specified array type

# Arguments
- `T`: Numeric type for values
- `items`: Vector of items (from metadata.vars or metadata.residuals)
- `metadata`: ModelMetadata containing model name
- `am`: AddressManager containing address assignments
- `array_type`: Array type constructor (Vector, CuVector, etc.)

# Implementation Note
This unified function handles both variables and residuals, which follow identical
storage patterns. Previously implemented as separate _build_variable_storage and
_build_residual_storage functions, but consolidated to eliminate duplication.
"""
function _build_storage(
    ::Type{T},
    items::Vector,
    metadata::ModelMetadata,
    am::AddressManager,
    array_type::Type{<:AbstractVector}
) where T<:Number
    AT = array_type{T}
    model_name = metadata.name

    addr_pairs = Pair{Symbol, Any}[]
    val_pairs = Pair{Symbol, AT}[]

    for item in items
        item_name = item.name

        addrs_cpu = item_addresses(am, model_name, item)

        # Transfer to GPU if needed (one-time cost)
        # Convert addresses to target array type (Vector or CuVector)
        AddrArrayType = array_type{UInt64}
        addrs = AddrArrayType(addrs_cpu)

        # Allocate value storage (on the target array type directly: zeros on the host
        # converted to it would allocate twice)
        push!(addr_pairs, item_name => addrs)
        push!(val_pairs, item_name => fill!(AT(undef, length(addrs)), zero(T)))
    end

    return addr_pairs, val_pairs
end

"""
    _build_var_gather_tuple(var_addrs, var_vals, metadata) -> Tuple

Build pre-computed tuple for zero-allocation variable gathering in TDS.

Returns a tuple where each element is `(addrs, vals, is_state::Bool)`:
- `addrs`: Address vector for this variable
- `vals`: Value vector for this variable
- `is_state`: true if differential state (gather from x), false if algebraic (gather from y)

This eliminates runtime Symbol lookup via `getfield(ws.var_addrs, var_name)`.

Note: Uses explicit loop instead of `map` to reduce method compilations (~78% fewer
gather/scatter compilations), though JIT time improvement is negligible due to
NamedTuple type specialization being the dominant cost.
"""
function _build_var_gather_tuple(
    var_addrs::NamedTuple,
    var_vals::NamedTuple,
    metadata::ModelMetadata
)
    n = length(metadata.vars)
    entries = Vector{Any}(undef, n)
    @inbounds for i in 1:n
        var_meta = metadata.vars[i]
        var_name = var_meta.name
        addrs = getfield(var_addrs, var_name)
        vals = getfield(var_vals, var_name)
        is_state = is_state_var(var_meta)
        entries[i] = (addrs, vals, is_state)
    end
    return Tuple(entries)
end

"""
    build_res_scatter_tuple(res_addrs, res_vals, metadata; skip=Returns(false)) -> Tuple

Build pre-computed tuple for zero-allocation residual scattering in TDS.

Returns a tuple where each element is `(addrs, vals, is_diff::Bool, mode::Symbol)`:
- `addrs`: Address vector for this residual
- `vals`: Value vector for this residual
- `is_diff`: true if differential residual (scatter to f), false if algebraic (scatter to g)
- `mode`: Scatter mode (:private, :add, or :set)

This eliminates runtime Symbol lookup via `getfield(ws.res_addrs, res_name)`.

# Arguments
- `res_addrs`: NamedTuple of residual address vectors
- `res_vals`: NamedTuple of residual value vectors
- `metadata`: ModelMetadata with residual specifications
- `skip`: Optional predicate function `(res_meta) -> Bool`. Residuals for which
  `skip(res_meta)` returns `true` are excluded from the tuple. Default: include all.

# Example
```julia
# Standard scatter tuple (all residuals)
scatter_tup = build_res_scatter_tuple(res_addrs, res_vals, metadata)

# Filtered scatter tuple (skip bus injections for Y-bus formulation)
internal_tup = build_res_scatter_tuple(res_addrs, res_vals, metadata;
    skip = res -> res.source_model === :bus_group)
```

Note: Uses explicit loop instead of `map` to reduce method compilations (~78% fewer
gather/scatter compilations), though JIT time improvement is negligible due to
NamedTuple type specialization being the dominant cost.
"""
function build_res_scatter_tuple(
    res_addrs::NamedTuple,
    res_vals::NamedTuple,
    metadata::ModelMetadata;
    skip::Function = Returns(false)
)
    entries = Any[]
    @inbounds for i in 1:length(metadata.residuals)
        res_meta = metadata.residuals[i]

        # Apply filter - skip this residual if predicate returns true
        skip(res_meta) && continue

        res_name = res_meta.name
        addrs = getfield(res_addrs, res_name)
        vals = getfield(res_vals, res_name)
        is_diff = is_state_res(res_meta)

        # A residual no instance is delivered through is not scattered; one only some are
        # (the other form of an injection) scatters those
        unrouted = count(==(UNROUTED), addrs)
        unrouted == length(addrs) && !isempty(addrs) && continue
        unrouted > 0 && (addrs = PartialAddresses(collect(UInt, addrs)))

        # Determine scatter mode from access pattern
        mode = if is_private_res(res_meta)
            :private
        elseif is_shared_res(:any)(res_meta)
            res_meta.access.mode  # :add or :set
        else
            error("Unknown residual access pattern: $(res_meta.access)")
        end

        push!(entries, (addrs, vals, is_diff, mode))
    end
    return Tuple(entries)
end

# ============================================================================
# Core Workspace Structure
# ============================================================================

"""
    ResidualBuffers{T, AT, MaskT, VarAddrs, VarVals, ResAddrs, ResVals}

Pre-allocated workspace for a single model type's residual calculation.

Uses NamedTuples for type-stable access to variable and residual arrays.
All memory is allocated during construction - zero allocation during calculation.

# Type Parameters
- `T`: Float type (Float64, Float32, etc.)
- `AT`: Array type for data storage (Vector{T} for CPU, CuVector{T} for GPU)
- `VarAddrs`: NamedTuple type for variable addresses
- `VarVals`: NamedTuple type for variable values
- `ResAddrs`: NamedTuple type for residual addresses
- `ResVals`: NamedTuple type for residual values
- `System`: SystemConfig type

# Fields
- `var_addrs`: NamedTuple mapping variable names to address vectors (always CPU)
- `var_vals`: NamedTuple mapping variable names to value vectors (can be CPU/GPU)
- `res_addrs`: NamedTuple mapping residual names to address vectors (always CPU)
- `res_vals`: NamedTuple mapping residual names to value vectors (can be CPU/GPU)
- `setpoints`: NamedTuple of initialization setpoints (reference, not typed for simplicity)
- `automata`: the `AutomatonState` of the model's automata (modes, discrete quantities,
  memory and the values of its guards), empty for a model without automata; the kernels
  read its discrete quantities and memory
- `system`: SystemConfig of the system the workspace belongs to (`sys.settings`)
- `n_instances`: Number of model instances

# Note on active_mask
The `active_mask` (which instances are in-service) is NOT stored in this struct.
Instead, it is passed as a parameter to `scatter_residuals!` and `calc_residuals!`.
This enables parallel ACCC with thread-local status arrays - each thread can have
its own isolated status without affecting others.

# Setpoints
For dynamic models (e.g., GENCLS), `setpoints` stores initialization setpoints
used in holding equations (e.g., P_m0, E_q0). These are passed by reference from
the TDS initialization result, not copied.
For models without setpoints, this is an empty NamedTuple `()`.

# System settings
`ws.system` holds the system-wide settings of the system (`sys.settings`), type-stable.

# Example
For a PQ load model on CPU (no setpoints):
```julia
ws = ResidualBuffers{Float64, Vector{Float64}, ...}(...)
v_values = ws.var_vals.v  # Vector{Float64}
ws.setpoints  # NamedTuple()
ws.system.frequency_base  # 60.0 Hz
```

For GENCLS on CPU (with setpoints):
```julia
ws = ResidualBuffers{Float64, Vector{Float64}, ...}(...)
ws.setpoints.P_m0  # Vector{Float64} - mechanical power setpoints
ws.setpoints.E_q0  # Vector{Float64} - internal voltage setpoints
omega_base = ws.system.omega_base  # System omega
```
"""
struct ResidualBuffers{T<:Number,
                         AT<:AbstractVector{T},
                         VarAddrs<:NamedTuple,
                         VarVals<:NamedTuple,
                         VarGatherTuple<:Tuple,
                         ResAddrs<:NamedTuple,
                         ResVals<:NamedTuple,
                         ResScatterTuple<:Tuple,
                         Setpoints<:NamedTuple,
                         Automata<:AutomatonState,
                         System<:SystemConfig}
    # --- Variable storage (for gather operations) ---
    var_addrs::VarAddrs       # NamedTuple of address vectors (UInt64 indices into x/y)
    var_vals::VarVals         # NamedTuple of value vectors (workspace for gathered values)
    var_gather_tuple::VarGatherTuple  # Pre-computed tuple: ((addrs, vals, is_state), ...)

    # --- Residual storage (for scatter operations) ---
    res_addrs::ResAddrs       # NamedTuple of address vectors (UInt64 indices into f/g)
    res_vals::ResVals         # NamedTuple of value vectors (workspace for computed residuals)
    res_scatter_tuple::ResScatterTuple  # Pre-computed tuple: ((addrs, vals, is_diff, mode), ...)

    # --- Instance state ---
    n_instances::Int          # Number of model instances

    # --- Configuration ---
    setpoints::Setpoints      # Initialization setpoints (type-stable, empty NamedTuple for models without)
    automata::Automata        # AutomatonState of the model's automata (empty without any)
    system::System            # SystemConfig of the system (sys.settings)
end

# ============================================================================
# Constructor from ModelMetadata
# ============================================================================

"""
    ResidualBuffers(T, metadata, am, model_vec, array_type; system, setpoints, automata)

Construct pre-allocated workspace for a model's residual calculations.

# Arguments
- `T`: Float type (e.g., Float64)
- `metadata`: ModelMetadata with variable/residual specs
- `am`: AddressManager with address assignments
- `model_vec`: Model StructArray (e.g., PQVec)
- `array_type`: Vector (CPU) or CuVector (GPU)
- `system`: SystemConfig of the system (`sys.settings`)
- `setpoints`: initialization setpoints (optional)
- `automata`: automaton state of the Float64 workspace to share (optional): the workspaces
  of every element type built for the same system share the modes and discrete quantities,
  so that Jacobians see those the Float64 evaluation uses; the memory and values are of
  element type `T`, this workspace's own (`AutomatonState(automata, T)`); otherwise a new `AutomatonState` for
  a model that declares automata (which checks them, `check_automata`), and an empty
  one for a model without
"""
function ResidualBuffers(
    ::Type{T},
    metadata::ModelMetadata,
    am::AddressManager,
    model_vec::StructArray,
    array_type::Type{<:AbstractVector};
    system::SystemConfig,
    setpoints::NamedTuple=NamedTuple(),
    automata::Union{Nothing, AutomatonState}=nothing
) where {T<:Number}
    n_instances = length(model_vec)
    AT = array_type{T}

    # Build variable and residual storage using unified function
    var_addr_pairs, var_val_pairs = _build_storage(T, metadata.vars, metadata, am, array_type)
    res_addr_pairs, res_val_pairs = _build_storage(T, metadata.residuals, metadata, am, array_type)

    # Build output variable storage (diagnostic outputs without DAE addresses)
    # These are like P_e, Q_e in GENCLS - computed but not part of DAE residuals
    for output_name in metadata.observable_vars
        output_array = fill!(AT(undef, n_instances), zero(T))
        push!(res_val_pairs, output_name => output_array)
    end

    # Use provided setpoints reference if given, otherwise allocate new arrays
    # Reference semantics: when setpoints kwarg is provided, we store that reference directly
    # Allocation: only when no setpoints provided (e.g., power flow models without holding equations)
    if isempty(setpoints)
        # Allocate new arrays for setpoints from metadata
        # Setpoints are typed as T to support ForwardDiff for input Jacobian computation (Bode plots)
        setpoint_pairs = Pair{Symbol, AT}[]
        for init_param in metadata.initialized_params
            setpoint_array = fill!(AT(undef, n_instances), zero(T))
            push!(setpoint_pairs, init_param.name => setpoint_array)
        end
        setpoints_nt = NamedTuple(setpoint_pairs)
    else
        # Use provided reference directly (no copy)
        setpoints_nt = setpoints
    end

    # Automaton state: shared when given (memory and values in this element type),
    # otherwise new for a model with automata and empty without. Always an
    # AutomatonState (never nothing): a workspace type with a field that may be nothing
    # is inferred as a union of two workspace types, which the setup of every workspace
    # then splits at run time
    automata_state = if !isnothing(automata)
        AutomatonState(automata, T)
    elseif has_automata(metadata)
        AutomatonState(metadata, n_instances, array_type)
    else
        AutomatonState()
    end

    # Construct NamedTuples
    var_addrs_nt = NamedTuple(var_addr_pairs)
    var_vals_nt = NamedTuple(var_val_pairs)
    res_addrs_nt = NamedTuple(res_addr_pairs)
    res_vals_nt = NamedTuple(res_val_pairs)

    # Build pre-computed tuples for zero-allocation TDS gather/scatter
    # These eliminate runtime Symbol lookup in hot path
    var_gather_tup = _build_var_gather_tuple(var_addrs_nt, var_vals_nt, metadata)
    res_scatter_tup = build_res_scatter_tuple(res_addrs_nt, res_vals_nt, metadata)

    # Return workspace with full type parameters
    # setpoints_nt is either the provided reference or newly allocated arrays
    return _new_residual_buffers(T, AT, var_addrs_nt, var_vals_nt, var_gather_tup, res_addrs_nt, res_vals_nt,
                                 res_scatter_tup, n_instances, setpoints_nt, automata_state, system)
end

# The workspace of the given parts, of the type their types make (a function barrier: the
# parts are built with types known at run time only, and the workspace type, with its
# eleven parameters, would be instantiated at run time from them at every construction,
# a type-cache lookup whose cost grows with the parameters; here it is a constant of each
# compiled method)
@noinline function _new_residual_buffers(::Type{T}, ::Type{AT}, var_addrs::VA, var_vals::VV, var_gather::VG,
                                         res_addrs::RA, res_vals::RV, res_scatter::RS, n_instances::Int,
                                         setpoints::SP, automata::AU,
                                         system::SY) where {T, AT, VA, VV, VG, RA, RV, RS, SP, AU, SY}
    return ResidualBuffers{T, AT, VA, VV, VG, RA, RV, RS, SP, AU, SY}(
        var_addrs, var_vals, var_gather,        # variable storage
        res_addrs, res_vals, res_scatter,       # residual storage
        n_instances,                            # instance state
        setpoints, automata, system)            # configuration
end

"""
    ResidualBuffers(ws::ResidualBuffers, ::Type{T}; setpoints_at_values = false) -> ResidualBuffers

A copy of `ws` in element type `T`: every number it holds (variables, residuals,
setpoints, and the discrete quantities, memory and values of its automata) in buffers of
its own of element type `T` at the values of `ws`; its addresses, automaton
modes and configuration shared. With `setpoints_at_values`, the setpoints are those of `ws`
themselves, of its element type and shared (a kernel then reads them as the numbers they
are at that moment, and may branch on them). Built without the address manager and type-stable where
the type of `ws` is known: for evaluating a model's kernels with numbers of another type
on the fly (the trace of the variables its automaton values read, PowerfulTDS
`trace_value_reads`).
"""
function ResidualBuffers(ws::ResidualBuffers{S, AS}, ::Type{T}; setpoints_at_values::Bool = false) where {S, AS, T<:Number}
    AT = typeof(similar(AS(undef, 0), T))
    typed(v) = copyto!(similar(v, T), v)
    var_vals = map(typed, ws.var_vals)
    res_vals = map(typed, ws.res_vals)
    var_gather = map(entry -> (entry[1], _beside(entry[2], ws.var_vals, var_vals), entry[3]), ws.var_gather_tuple)
    res_scatter = map(entry -> (entry[1], _beside(entry[2], ws.res_vals, res_vals), entry[3], entry[4]),
                      ws.res_scatter_tuple)
    a = ws.automata
    automata = AutomatonState(a.modes, map(typed, a.discrete), map(typed, a.memory), map(typed, a.values),
                              _fresh_affine(a.affine, T), a.suspended)
    return _new_residual_buffers(T, AT, ws.var_addrs, var_vals, var_gather, ws.res_addrs, res_vals, res_scatter,
                                 ws.n_instances, setpoints_at_values ? ws.setpoints : map(typed, ws.setpoints),
                                 automata, ws.system)
end

"""
    value_buffers(ws::ResidualBuffers, ::Type{T}) -> ResidualBuffers

A copy of `ws` for evaluating the values of its automata (`calc_automaton_values!`) in element type `T`: its
variables and the values of its automata in buffers of element type `T` of its own (the variables at the values of
`ws`, the values `NaN`), everything else of `ws` itself, shared (residuals, setpoints, the discrete quantities and
memory of its automata: the kernel of the values reads the discrete state as it is, promoted). For the kernel of the
values evaluated in another element type at a state of `ws` (PowerfulTDS, the rounding of the values at a start).
"""
function value_buffers(ws::ResidualBuffers{S, AS}, ::Type{T}) where {S, AS, T<:Number}
    AT = typeof(similar(AS(undef, 0), T))
    var_vals = map(v -> copyto!(similar(v, T), v), ws.var_vals)
    var_gather = map(entry -> (entry[1], _beside(entry[2], ws.var_vals, var_vals), entry[3]), ws.var_gather_tuple)
    a = ws.automata
    automata = AutomatonState(a.modes, a.discrete, a.memory, map(v -> fill!(similar(v, T), T(NaN)), a.values),
                              _fresh_affine(a.affine, T), a.suspended)
    return _new_residual_buffers(T, AT, ws.var_addrs, var_vals, var_gather, ws.res_addrs, ws.res_vals,
                                 ws.res_scatter_tuple, ws.n_instances, ws.setpoints, automata, ws.system)
end

# The buffer of `new` named as the buffer `v` of `old`
_beside(v, old::NamedTuple, new::NamedTuple) = values(new)[findfirst(b -> b === v, values(old))]

# ============================================================================
# Generic Operations (Zero Allocation)
# ============================================================================

"""
    gather_variables!(ws::ResidualBuffers, x::AbstractVector, y::AbstractVector)

Gather variables from DAE vectors (x for differential states, y for algebraic) into workspace buffers.

Unified version that handles both power flow (x empty) and TDS (x contains states).
Uses pre-computed `var_gather_tuple` for ZERO ALLOCATION iteration.

# Arguments
- `ws::ResidualBuffers` - Workspace with variable buffers and pre-computed gather tuple
- `x::AbstractVector` - Differential state vector (delta, omega, ...) - empty for PF
- `y::AbstractVector` - Algebraic variable vector (v, theta, ...) - all variables for PF

# Notes
- For power flow (algebraic-only): x is empty, all variables gathered from y
- For TDS (differential-algebraic): state variables gathered from x, algebraic from y
- The `is_state` flag in var_gather_tuple determines which vector to use
"""
@inline function gather_variables!(
    ws::ResidualBuffers{T},
    x::AbstractVector,
    y::AbstractVector
) where T
    # Use pre-computed tuple for zero-allocation iteration
    # Each entry has (addrs, vals, is_state) - is_state determines x vs y
    _gather_vars!(ws.var_gather_tuple, x, y)
    nothing
end

"""
    gather_variables!(gather::Tuple, x::AbstractVector, y::AbstractVector)

Gather the variables of `gather` from `x` (states) and `y` (algebraic variables) into
their workspace buffers: entries of a workspace's `var_gather_tuple`, all of them or
those of some variables (`variable_gather`). Zero allocation.
"""
@inline function gather_variables!(gather::Tuple, x::AbstractVector, y::AbstractVector)
    _gather_vars!(gather, x, y)
    nothing
end

"""
    variable_gather(ws::ResidualBuffers, metadata::ModelMetadata, names) -> Tuple

The entries of `ws.var_gather_tuple` of the variables `names` of `metadata`, in their
declaration order: what `gather_variables!(gather, x, y)` gathers those variables with,
leaving the others' buffers as they are. Throws for a name that is not a variable of
`metadata`.
"""
function variable_gather(ws::ResidualBuffers, metadata::ModelMetadata, names)
    declared = [var.name for var in metadata.vars]
    unknown = setdiff(names, declared)
    isempty(unknown) || throw(ArgumentError(
        "$(metadata.name) has no variable $(join(unknown, ", ")) (its variables: $(join(declared, ", ")))"))
    return Tuple(ws.var_gather_tuple[k] for k in eachindex(declared) if declared[k] in names)
end

"""
    gather_residuals!(ws::ResidualBuffers{T}, f::AbstractVector, g::AbstractVector) where T

Gather accumulated residual values from global residual vectors into local workspace buffers.

This function is used for models with SharedRes(:set) residuals that need to READ
accumulated residuals from previous models before computing their own residuals.

Example: PV generators in CurrentInjection formulation need to read accumulated
branch currents (I_r, I_i) before computing power balance equations (P = V·I).

Uses pre-computed `res_scatter_tuple` for zero-allocation iteration. Gathers from
`f` (differential) or `g` (algebraic) based on `is_diff` flag in the tuple.

# Arguments
- `ws`: Workspace with pre-allocated residual buffers
- `f`: Global differential residual vector (e.g., `dae.f`)
- `g`: Global algebraic residual vector (e.g., `dae.g`)

# Notes
- Only gathers SharedRes(:set) residuals (ones that overwrite)
- Called BEFORE calc_residuals! for models that need to read accumulated residuals
- Zero allocation via pre-computed tuple iteration
"""
@inline function gather_residuals!(
    ws::ResidualBuffers{T},
    f::AbstractVector,
    g::AbstractVector
) where T
    # Use pre-computed res_scatter_tuple for zero-allocation iteration
    # Only gathers SharedRes(:set) residuals that need to read accumulated values
    _gather_res_set!(ws.res_scatter_tuple, f, g)
    nothing
end

"""
    zero_residuals!(ws::ResidualBuffers{T}) where T

Zero out local residual buffers before calculation.

Sets all values in `ws.res_vals` to zero. Should be called before
`calc_residuals!` to ensure clean slate.

# Arguments
- `ws`: Workspace with pre-allocated residual buffers
"""
@inline function zero_residuals!(ws::ResidualBuffers{T}) where T
    # Iterate over all residuals in the NamedTuple
    for (res_name, vals) in pairs(ws.res_vals)
        fill!(vals, zero(T))
    end

    nothing
end

"""
    scatter_residuals!(ws, f, g, active_mask)

Scatter residual values from workspace buffers to DAE vectors.

Uses pre-computed tuples for type-stable, zero-allocation iteration.
Respects `active_mask` to skip inactive instances.

# Two-Pass Algorithm
1. First pass: Scatter `:private` and `:add` modes (accumulation)
2. Second pass: Scatter `:set` modes (overwrite control equations)

# Arguments
- `ws`: Workspace with residual buffers
- `f`: Differential residual vector (can be empty for algebraic-only problems)
- `g`: Algebraic residual vector
- `active_mask`: Boolean vector indicating which instances are in-service

# See Also
- [`scatter_residuals_accumulate!`](@ref): Cold-path variant that skips `:set` mode
"""
@inline function scatter_residuals!(
    ws::ResidualBuffers{T},
    f::AbstractVector,
    g::AbstractVector,
    active_mask::AbstractVector{Bool}
) where T
    # Two-pass scatter using pre-computed tuples
    # Pass 1: Scatter :private and :add modes
    scatter_res_add!(ws.res_scatter_tuple, f, g, active_mask)
    # Pass 2: Scatter :set modes (overwrite control equations)
    scatter_res_set!(ws.res_scatter_tuple, f, g, active_mask)

    nothing
end

"""
    scatter_residuals_accumulate!(ws, f, g, active_mask)

Cold-path scatter that only does the `:add` pass, skipping `:set` mode entries.

Used by `evaluate_generator_contributions!` to accumulate all model contributions
without overwriting bus balance equations. This allows reading the actual P/Q
injections from generators after power flow converges.

Unlike `scatter_residuals!`, this function:
- Only scatters `:add` and `:private` modes
- Skips `:set` modes entirely (treats them as `:add` implicitly)

This is a cold-path function called once after power flow converges.

# Arguments
- `ws`: Workspace with residual buffers
- `f`: Differential residual vector
- `g`: Algebraic residual vector
- `active_mask`: Boolean vector indicating which instances are in-service
"""
@inline function scatter_residuals_accumulate!(
    ws::ResidualBuffers{T},
    f::AbstractVector,
    g::AbstractVector,
    active_mask::AbstractVector{Bool}
) where T
    # Only the :add pass - skips :set modes
    scatter_res_add!(ws.res_scatter_tuple, f, g, active_mask)
    nothing
end

# ============================================================================
# TDS Gather/Scatter - Handle differential (x,f) and algebraic (y,g) separately
# ============================================================================

# --- Recursive helpers for zero-allocation iteration over pre-computed tuples ---
# NOTE: if/else avoids Union type inference (reduces compilations, not JIT time)

# Base case: empty tuple
@inline _gather_vars!(::Tuple{}, x, y) = nothing

# Recursive case: gather one variable, then recurse
@inline function _gather_vars!(tup::Tuple, x, y)
    addrs, vals, is_state = first(tup)
    if is_state
        gather!(vals, x, addrs)
    else
        gather!(vals, y, addrs)
    end
    _gather_vars!(Base.tail(tup), x, y)
end

# Scatter one residual entry into f (differential) or g (algebraic)
@inline function _scatter_entry!(entry, f, g, active_mask)
    addrs, vals, is_diff, mode = entry
    if is_diff
        scatter!(f, vals, addrs, active_mask, mode)
    else
        scatter!(g, vals, addrs, active_mask, mode)
    end
    nothing
end

# --- Two-pass scatter helpers for correct ordering (add/private first, set last) ---

"""
    scatter_res_add!(res_scatter_tuple::Tuple, f, g, active_mask)

First pass of the residual scatter: write every `:private` and `:add` entry of a
pre-built scatter tuple (see [`build_res_scatter_tuple`](@ref)) into `f`
(differential entries) or `g` (algebraic entries), skipping inactive instances.
`:set` entries are left for [`scatter_res_set!`](@ref).
"""
function scatter_res_add! end

"""
    scatter_res_set!(res_scatter_tuple::Tuple, f, g, active_mask)

Second pass of the residual scatter: overwrite the addresses of every `:set` entry
([`SharedRes`](@ref)`(:set)`) with its values, after all contributions were
accumulated by [`scatter_res_add!`](@ref).
"""
function scatter_res_set! end

# Base case: empty tuple
@inline scatter_res_add!(::Tuple{}, f, g, active_mask) = nothing

# Recursive case: scatter only :add and :private modes, skip :set
@inline function scatter_res_add!(tup::Tuple, f, g, active_mask)
    entry = first(tup)
    entry[4] !== :set && _scatter_entry!(entry, f, g, active_mask)
    scatter_res_add!(Base.tail(tup), f, g, active_mask)
end

# Base case: empty tuple
@inline scatter_res_set!(::Tuple{}, f, g, active_mask) = nothing

# Recursive case: scatter only :set mode, skip :add and :private
@inline function scatter_res_set!(tup::Tuple, f, g, active_mask)
    entry = first(tup)
    entry[4] === :set && _scatter_entry!(entry, f, g, active_mask)
    scatter_res_set!(Base.tail(tup), f, g, active_mask)
end

# ============================================================================
# Zero-Allocation Gather Helpers (for SharedRes(:set) residuals)
# ============================================================================

# Base case: empty tuple
@inline _gather_res_set!(::Tuple{}, f, g) = nothing

# Recursive case: gather only :set mode residuals
@inline function _gather_res_set!(tup::Tuple, f, g)
    addrs, vals, is_diff, mode = first(tup)
    if mode === :set
        if is_diff
            gather!(vals, f, addrs)
        else
            gather!(vals, g, addrs)
        end
    end
    _gather_res_set!(Base.tail(tup), f, g)
end

