"""
    DAE{AT}

Paired vectors for DAE/algebraic systems: variables and their residuals.

In **power flow** (purely algebraic): This holds the current operating point during
Newton iteration. The goal is to find y such that g(y) = 0. Differential fields (n, x, f)
are unused (n=0, empty vectors).

In **time-domain simulation** (differential-algebraic): This holds both differential
states (x, f) and algebraic variables (y, g). The system is M dx/dt = f(x,y), 0 = g(x,y),
where `M` is the mass matrix.

# Type Parameters
- `AT`: Array type (Vector{T} for CPU, CuVector{T} for GPU)

# Fields
- `m`: Number of algebraic variables and residuals
- `y`: Algebraic variable vector (unknowns to solve for)
- `g`: Algebraic residual vector (equation evaluations, should be zero at solution)
- `n`: Number of differential states and residuals (0 for power flow)
- `x`: Differential state vector (empty for power flow)
- `f`: Differential residual vector (dx/dt evaluations, empty for power flow)

# Example
```julia
# Power flow (algebraic only, CPU default)
dae = DAE(am)  # Float64, Vector by default
set_algebraic_guess!(dae.y, am, sys, PowerInjection())  # optional initial guess

# Power flow (GPU)
dae = DAE(am; array_type=CuArray)

# Time-domain simulation (differential-algebraic)
dae = DAE(am; allocate_differential=true)
```
"""
struct DAE{AT<:AbstractVector}
    # Algebraic system (always present)
    m::Int  # Number of algebraic variables & algebraic residuals
    y::AT   # Algebraic variables (can be CPU or GPU)
    g::AT   # Algebraic residuals (can be CPU or GPU)

    # Differential system (optional - for TDS)
    n::Int  # Number of differential states & differential residuals
    x::AT   # Differential states (can be CPU or GPU)
    f::AT   # Differential residuals (dx/dt, can be CPU or GPU)
end

"""
    DAE(am; eltype=Float64, array_type=Vector, allocate_differential=false)

Construct DAE with sizing from AddressManager. Initialization is caller's responsibility.

# Arguments
- `am`: AddressManager with assigned addresses
- `eltype`: Element type (default: Float64)
- `array_type`: Vector (default) or CuArray
- `allocate_differential`: Allocate x, f vectors for TDS

# Example
```julia
dae = DAE(am)                              # CPU, Float64
dae = DAE(am; array_type=CuArray)          # GPU
dae = DAE(am; allocate_differential=true)  # TDS
```
"""
function DAE(am::AddressManager;
             eltype::Type{T}=Float64,
             array_type::Type{<:AbstractVector}=Vector,
             allocate_differential::Bool=false) where T<:Number
    # Determine algebraic size from address manager
    m = max_address_in_space(am, AlgebVar())

    # Construct array type with element type T
    AT = array_type{T}

    # Allocate algebraic vectors
    y = AT(undef, m)
    g = AT(undef, m)
    fill!(y, zero(T))
    fill!(g, zero(T))

    # Allocate differential vectors (for TDS)
    if allocate_differential
        n = max_address_in_space(am, StateVar())
        x = AT(undef, n)
        f = AT(undef, n)
        fill!(x, zero(T))
        fill!(f, zero(T))
        # Note: Initialization of x is handled by TDS separately
    else
        # Power flow: no differential states
        n = 0
        x = AT(undef, 0)
        f = AT(undef, 0)
    end

    return DAE{typeof(y)}(m, y, g, n, x, f)
end

"""
    max_address_in_space(am::AddressManager, addr_type) -> Int

Compute the maximum address value for a given address type.

Helper function for determining vector sizes from AddressManager.

# Arguments
- `am`: AddressManager with assigned addresses
- `addr_type`: Address type to query (StateVar, AlgebVar, StateRes, AlgebRes)

# Returns
- Maximum address in that space, or 0 if no addresses assigned

# Example
```julia
m = max_address_in_space(am, AlgebVar())  # Size of algebraic space
n = max_address_in_space(am, StateVar())   # Size of differential space
```
"""
function max_address_in_space(am::AddressManager, addr_type)
    max_addr = 0
    for (key, addrs) in am.addresses
        key_type, _, _ = key
        if key_type == addr_type && !isempty(addrs)
            max_addr = max(max_addr, maximum(addrs))
        end
    end
    return max_addr
end
