"""
Differential-Algebraic Equations

Represents the system of differential-algebraic equations for power flow analysis:
- Algebraic equations: g(y) = 0
- Future: Differential equations: f(x, y) = dx/dt

# Type Parameters
- `AT`: Array type (Vector{T} for CPU, CuVector{T} for GPU)

# Fields
- `m`: Number of algebraic variables and residuals
- `g`: Algebraic residual vector (size m)
- `y`: Algebraic variable vector (size m)

# Example
```julia
# Create DAE system from SystemModel (CPU)
sys = SystemModel(am, models)
dae = DAE(Float64, sys)

# GPU version: transfer system first, then create DAE
sys_gpu = to_gpu(sys)
dae_gpu = DAE(Float64, sys_gpu)

# Solve power flow
for iter in 1:max_iter
    evaluate_system_residuals!(sys_ws, dae)
    if norm(dae.g, Inf) < tol
        break
    end
    # ... Jacobian and Newton update
end
```
"""
struct DAE{AT<:AbstractVector}
    m::Int  # Number of algeb variables & algeb residuals
    g::AT   # Algeb residuals (can be CPU or GPU)
    y::AT   # Algeb variables (can be CPU or GPU)
end

"""
    DAE(::Type{T}, sys::SystemModel) where T

Construct DAE system from SystemModel with proper sizing.

Determines the size of the algebraic system from the address manager,
allocates y and g vectors, and initializes variables from model data.

The array type (CPU vs GPU) is determined by inspecting the first model's
array type. Use `to_gpu(sys)` before calling this constructor to create
GPU-backed DAE.

# Arguments
- `T`: Float type for calculations
- `sys`: SystemModel with address manager and model data

# Returns
- `DAE{AT}`: Initialized DAE system ready for power flow analysis

# Example
```julia
# CPU version
sys = SystemModel(am, (Bus = bus_vec, PQ = pq_vec))
dae = DAE(Float64, sys)

# GPU version: transfer system first
sys_gpu = to_gpu(sys)
dae_gpu = DAE(Float64, sys_gpu)
```
"""
function DAE(::Type{T}, sys::SystemModel) where T<:Number
    am = sys.address_manager

    # Determine size from address manager
    m = 0
    for (key, addrs) in am.addresses
        addr_type, _, _ = key
        if addr_type isa AlgebVar || addr_type isa AlgebRes
            if !isempty(addrs)
                m = max(m, maximum(addrs))
            end
        end
    end

    # Get array type from system (explicitly stored)
    AT = sys.array_type{T}

    # Allocate vectors with system's array type
    y = AT(undef, m)
    g = AT(undef, m)
    fill!(y, zero(T))
    fill!(g, zero(T))

    # Initialize variables from model data
    if haskey(sys.models, :Bus)
        bus = sys.models.Bus

        # Get voltage magnitude and angle addresses
        if haskey(am.assignments[:Bus].by_type[AlgebVar()], :v)
            v_addrs = am.assignments[:Bus].by_type[AlgebVar()][:v]
            # Direct array assignment (works for both CPU and GPU)
            y[v_addrs] .= bus.vm
        end

        if haskey(am.assignments[:Bus].by_type[AlgebVar()], :theta)
            theta_addrs = am.assignments[:Bus].by_type[AlgebVar()][:theta]
            # Convert degrees to radians and assign
            y[theta_addrs] .= deg2rad.(bus.va)
        end
    end

    return DAE{typeof(y)}(m, g, y)
end

