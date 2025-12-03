import Base: iterate, length

"""
Helper functions to work with components
"""
function component_types(sys::SystemModel)
    # Get the concrete type of each component
    map(typeof, sys.models)
end

# Helper type with explicit tuple type
struct ComponentIterator{T<:Tuple, S<:SystemModel{T}}
    system::S
    n::Int  # Current iteration number
end

Base.length(sys::SystemModel) = length(sys.models)

# Type-stable iteration
function Base.iterate(iter::ComponentIterator{T, S}) where {T, S}
    # Now compiler knows exact tuple type T
    if iter.n > fieldcount(T)  # More precise than length check
        return nothing
    end
    
    # Compiler can now determine exact type of T.parameters[iter.n]
    return (iter.system.components[iter.n], 
            ComponentIterator{T, S}(iter.system, iter.n + 1))
end

# Make SystemModel iterable using our type-stable iterator
Base.iterate(sys::SystemModel{T}) where {T} = 
    iterate(ComponentIterator{T, typeof(sys)}(sys, 1))

# This is the corrected method
Base.iterate(sys::SystemModel{T}, state::ComponentIterator{T,S}) where {T,S} = 
    iterate(state)

# ============================================================================
# System Setup and Validation
# ============================================================================

"""
    setup!(sys::SystemModel)

Automatically assign addresses for all models in the system.

This function performs three phases:
1. **Phase 1**: Assign internal addresses for all model variables/residuals
2. **Phase 2**: Retrieve external addresses (mapping between models)
3. **Phase 3**: Basic validation

After calling setup!, the system is ready for creating DAE and workspace objects.

# Arguments
- `sys::SystemModel`: System to set up

# Returns
- `sys`: The same system (for chaining)

# Example
```julia
# Automatic setup (recommended)
sys = SystemModel(am, models; auto_setup=true)  # setup! called automatically

# Manual setup
sys = SystemModel(am, models; auto_setup=false)
setup!(sys)

# Check if setup
if !is_setup(sys)
    setup!(sys)
end
```

# Errors
- If a model has 0 instances, it will be skipped with a warning
- If Bus model is missing or empty, validation warnings will be issued
- If external address retrieval fails (e.g., missing source model), an error is thrown

# See Also
- [`is_setup`](@ref): Check if system has been set up
- [`assign_address!`](@ref): Low-level function for manual address assignment
- [`retrieve_addresses`](@ref): Low-level function for external address retrieval
"""
function setup!(sys::SystemModel)
    am = sys.address_manager

    # Phase 1: Assign internal addresses
    for (model_name, model_vec) in pairs(sys.models)
        n_instances = length(model_vec)

        # Skip models with no instances
        if n_instances == 0
            @warn "Model $model_name has 0 instances, skipping address assignment"
            continue
        end

        base_type = get_base_model_type(model_vec)
        metadata = model_metadata(base_type)

        assign_address!(
            am,
            model_name,
            metadata,
            n_instances
        )
    end

    # Phase 2: Retrieve external addresses
    for (_, model_vec) in pairs(sys.models)
        # Skip if no instances
        if length(model_vec) == 0
            continue
        end

        base_type = get_base_model_type(model_vec)
        metadata = model_metadata(base_type)

        retrieve_addresses(am, metadata, sys.models)
    end

    # Phase 3: Basic validation
    _validate_system_basic(sys)

    return sys
end

"""
    is_setup(sys::SystemModel) -> Bool

Check if a system has been set up (addresses assigned).

A system is considered "set up" if all non-empty models have addresses assigned
in the address manager.

# Arguments
- `sys::SystemModel`: System to check

# Returns
- `true` if system is set up, `false` otherwise

# Example
```julia
sys = SystemModel(am, models; auto_setup=false)
@assert !is_setup(sys)

setup!(sys)
@assert is_setup(sys)

# Conditional setup
if !is_setup(sys)
    @info "Setting up system..."
    setup!(sys)
end
```
"""
function is_setup(sys::SystemModel)
    # Check if all non-empty models have addresses assigned
    for (model_name, model_vec) in pairs(sys.models)
        # Skip empty models
        if length(model_vec) == 0
            continue
        end

        # Check if this model has been assigned
        if !haskey(sys.address_manager.assignments, model_name)
            return false
        end
    end

    return true
end

"""
    _validate_system_basic(sys::SystemModel)

Perform basic validation checks on the system.

Internal function called by setup!. Performs sanity checks and issues warnings
for common problems.

# Checks
- Bus model exists and is non-empty
- (Future: Check for orphaned devices, topology validation, etc.)

# Arguments
- `sys::SystemModel`: System to validate
"""
function _validate_system_basic(sys::SystemModel)
    # Check 1: Bus model must exist
    if !haskey(sys.models, :Bus)
        @warn "System does not contain a Bus model. Most models require buses to function."
        return  # Can't do further validation without buses
    end

    # Check 2: At least one bus should exist
    if length(sys.models.Bus) == 0
        error("Bus model exists but contains no buses. Cannot create valid power system.")
    end

    # Future validations:
    # - Check for orphaned devices (devices referencing non-existent buses)
    # - Check for duplicate keys within models
    # - Validate topology connectivity
    # - Check parameter bounds

    nothing
end

export setup!, is_setup

# ============================================================================
# GPU Transfer
# ============================================================================

"""
    to_gpu(sys::SystemModel) -> SystemModel

Transfer a SystemModel to GPU by converting all model data to CuArrays.

Uses `StructArrays.replace_storage` to transfer each model's arrays to GPU.
Returns a new SystemModel with GPU-backed model vectors and `array_type = CuVector`.

# Arguments
- `sys`: CPU-based SystemModel

# Returns
- New SystemModel with GPU arrays

# Example
```julia
# Create system on CPU
sys = SystemModel(am, (Bus = bus, PQ = pq))

# Transfer to GPU
sys_gpu = to_gpu(sys)

# Create GPU DAE and workspace
dae_gpu = DAE(Float64, sys_gpu)
ws_gpu = SystemResidualWorkspace(Float64, sys_gpu)
```

# Notes
- Original system is not modified
- All model vectors are transferred to GPU
- Address manager is shared (addresses are just UInt64)
"""
function to_gpu(sys::SystemModel)
    # Transfer each model to GPU using StructArrays.replace_storage
    gpu_models = NamedTuple(
        name => StructArrays.replace_storage(CUDA.CuArray, model)
        for (name, model) in pairs(sys.models)
    )

    # Create new SystemModel with GPU arrays
    # Don't auto-setup since addresses are already assigned
    return SystemModel(sys.address_manager, gpu_models;
                      auto_setup=false,
                      array_type=CUDA.CuVector)
end

export to_gpu

# ============================================================================
# System Residual Workspace
# ============================================================================

# ============================================================================
# System-Level Workspace
# ============================================================================

"""
    SystemResidualWorkspace{T, WS, Models, Metadata}

System-level workspace containing pre-allocated workspaces for all model types.

Uses NamedTuples for type-stable access to model-specific workspaces. All memory
is allocated during construction - the hot path `evaluate_system_residuals!`
performs zero allocation.

# Type Parameters
- `T`: Float type (Float64, Float32, etc.)
- `WS`: NamedTuple type for model workspaces
- `Models`: NamedTuple type for model vectors
- `Metadata`: NamedTuple type for model metadata

# Fields
- `model_workspaces`: NamedTuple mapping model names to ResidualWorkspace
- `model_vectors`: NamedTuple mapping model names to model vectors (cached)
- `model_metadata`: NamedTuple mapping model names to ModelMetadata (cached)

# Example
```julia
# Construction (once, allocates all memory)
sys = SystemModel(am, (PQ = pq_vec, Bus = bus_vec, Branch = branch_vec))
sys_ws = SystemResidualWorkspace(Float64, sys)

# Hot path (zero allocation)
for iter in 1:max_iter
    evaluate_system_residuals!(sys_ws, dae)
    # ... check convergence, solve Jacobian ...
end
```
"""
struct SystemResidualWorkspace{T<:Number, WS<:NamedTuple, Models<:NamedTuple, Metadata<:NamedTuple}
    model_workspaces::WS
    model_vectors::Models
    model_metadata::Metadata
end

# ============================================================================
# Constructor
# ============================================================================

"""
    SystemResidualWorkspace(::Type{T}, sys::SystemModel) where T

Construct system workspace with all model workspaces pre-allocated.

Iterates through all models in the system, creates a ResidualWorkspace
for each, and caches model vectors and metadata. This function performs ALL
memory allocation - subsequent residual evaluations are zero-allocation.

The array type (CPU vs GPU) is inferred from the model data.
Use `to_gpu(sys)` to transfer the system to GPU before creating workspaces.

# Arguments
- `T`: Float type for calculations
- `sys`: SystemModel containing address manager and model vectors

# Returns
- `SystemResidualWorkspace`: Fully allocated workspace ready for calculations

# Example
```julia
# CPU version
sys = SystemModel(am, (PQ = pq_vec, Bus = bus_vec))
sys_ws = SystemResidualWorkspace(Float64, sys)

# GPU version: transfer system first
sys_gpu = to_gpu(sys)
sys_ws_gpu = SystemResidualWorkspace(Float64, sys_gpu)
```
"""
function SystemResidualWorkspace(
    ::Type{T},
    sys::SystemModel
) where T<:Number

    # Verify system is set up
    if !is_setup(sys)
        error("""
        SystemModel must be set up before creating SystemResidualWorkspace.

        The system has not been initialized with address assignments.

        Solution 1 (Recommended): Create SystemModel with automatic setup:
            sys = SystemModel(am, models)  # auto_setup=true by default

        Solution 2: Call setup! manually:
            sys = SystemModel(am, models; auto_setup=false)
            setup!(sys)
        """)
    end

    am = sys.address_manager

    # NOTE: retrieve_addresses is now called in setup!()
    # No need to call it here - addresses are already retrieved

    # Build workspace for each model
    workspace_pairs = Pair{Symbol, ResidualWorkspace}[]
    model_pairs = Pair{Symbol, Any}[]  # Any to support both AbstractModel and StructArray
    metadata_pairs = Pair{Symbol, ModelMetadata}[]

    for (model_name, model_vec) in pairs(sys.models)
        try
            # Use helper functions (cleaner and eliminates code duplication!)
            base_type = get_base_model_type(model_vec)
            metadata = model_metadata(base_type)

            # Create workspace for this model with system's array type
            ws = ResidualWorkspace(T, metadata, am, model_vec, sys.array_type)

            # Store in pairs
            push!(workspace_pairs, model_name => ws)
            push!(model_pairs, model_name => model_vec)
            push!(metadata_pairs, model_name => metadata)

        catch e
            @warn "Failed to create workspace for model $model_name" exception=e
            rethrow(e)
        end
    end

    # Convert to NamedTuples for type stability
    workspaces_nt = NamedTuple(workspace_pairs)
    models_nt = NamedTuple(model_pairs)
    metadata_nt = NamedTuple(metadata_pairs)

    return SystemResidualWorkspace{T, typeof(workspaces_nt), typeof(models_nt), typeof(metadata_nt)}(
        workspaces_nt,
        models_nt,
        metadata_nt
    )
end

# ============================================================================
# System Residual Evaluation (Hot Path - Zero Allocation)
# ============================================================================

"""
    evaluate_system_residuals!(sys_ws::SystemResidualWorkspace{T}, dae::DAE{T}) where T

Evaluate residuals for all models in the system - ZERO ALLOCATION.

This is the hot path called in every Newton-Raphson iteration. Orchestrates
the gather-compute-scatter pattern across all model types:

1. Zero out global residual vector
2. For each model:
   - Gather variables from DAE into local workspace
   - Zero local residual buffers
   - Calculate model-specific residuals
   - Scatter results back to global DAE

# Arguments
- `sys_ws`: System workspace with pre-allocated buffers
- `dae`: DAE system with global variable and residual vectors

# Performance Notes
- Zero allocation after initial workspace construction
- Type-stable dispatch on each model's calc_residuals!
- Models are evaluated sequentially (can be parallelized if SharedRes uses atomics)
- Compiler inlines and optimizes the entire call chain

# Example
```julia
# In Newton-Raphson loop
for iter in 1:max_iter
    evaluate_system_residuals!(sys_ws, dae)  # ← Zero allocation here

    if norm(dae.g, Inf) < tol
        break
    end

    calculate_jacobian!(J, sys_ws, dae)
    Δy = J \\ dae.g
    dae.y .-= Δy
end
```
"""
function evaluate_system_residuals!(
    sys_ws::SystemResidualWorkspace{T},
    dae::DAE{AT}
) where {T, AT}
    # Zero out global residual vector
    fill!(dae.g, zero(eltype(dae.g)))

    # Evaluate each model sequentially
    # Note: Can be parallelized with Threads.@threads if SharedRes uses atomic operations
    for (model_name, workspace) in pairs(sys_ws.model_workspaces)
        model_vec = getfield(sys_ws.model_vectors, model_name)
        metadata = getfield(sys_ws.model_metadata, model_name)

        # Generic gather-compute-scatter pattern
        evaluate_model_residuals!(workspace, model_vec, dae, metadata)
    end

    # Synchronize GPU kernels if needed (only once at the end)
    if AT <: CUDA.CuArray
        CUDA.synchronize()
    end

    nothing
end

"""
    evaluate_model_residuals!(ws::ResidualWorkspace{T}, model_vec::AbstractModel,
                              dae::DAE{T}, metadata::ModelMetadata) where T

Evaluate residuals for a single model type - zero allocation.

Generic function that works for any model by dispatching to model-specific
calc_residuals! implementations.

# Steps
1. Gather variables from DAE.y into workspace
2. Zero local residual buffers
3. Calculate residuals (model-specific dispatch)
4. Scatter residuals to DAE.g

# Arguments
- `ws`: Pre-allocated workspace for this model
- `model_vec`: Model vector with data
- `dae`: DAE system
- `metadata`: Model metadata for scatter operation
"""
@inline function evaluate_model_residuals!(
    ws::ResidualWorkspace{T},
    model_vec::StructArray{M},
    dae::DAE,
    metadata::ModelMetadata
) where {T, M<:AbstractModel}
    # Step 1: Gather variables from global DAE (zero alloc)
    gather_variables!(ws, dae)

    # Step 2: Zero local residual buffers - SKIPPED
    # Not needed! Kernels use branchless masking and write to ALL elements:
    #   - Active elements: mask * computed_value
    #   - Inactive elements: 0 * computed_value = 0
    # zero_residuals!(ws)

    # Step 3: Calculate residuals (model-specific, zero alloc)
    # Dispatches to model-specific calc_residuals! based on model_vec type
    calc_residuals!(ws, model_vec)

    # Step 4: Scatter results back to global DAE (zero alloc)
    scatter_residuals!(ws, dae, metadata)

    nothing
end

# ============================================================================
# Parallel Evaluation (Future Extension)
# ============================================================================

"""
    evaluate_system_residuals_parallel!(sys_ws::SystemResidualWorkspace{T}, dae::DAE{T}) where T

Parallel version of residual evaluation (requires atomic scatter for SharedRes).

NOT YET IMPLEMENTED - requires atomic operations for SharedRes accumulation.
"""
function evaluate_system_residuals_parallel!(
    ::SystemResidualWorkspace{T},
    ::DAE{T}
) where T
    error("Parallel residual evaluation not yet implemented. " *
          "Requires atomic operations for SharedRes scatter.")

    # Future implementation:
    # fill!(dae.g, zero(T))
    #
    # model_names = collect(keys(sys_ws.model_workspaces))
    # Threads.@threads for model_name in model_names
    #     ws = sys_ws.model_workspaces[model_name]
    #     model_vec = sys_ws.model_vectors[model_name]
    #     metadata = sys_ws.model_metadata[model_name]
    #
    #     gather_variables!(ws, dae)  # Thread-safe (read-only)
    #     zero_residuals!(ws)
    #     calc_residuals!(ws, model_vec)
    #     scatter_residuals_atomic!(ws, dae, metadata)  # Needs atomic +=
    # end
    #
    # nothing
end

# ============================================================================
# Exports
# ============================================================================

export setup!, is_setup
export SystemResidualWorkspace
export evaluate_system_residuals!
