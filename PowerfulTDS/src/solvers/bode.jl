"""
Linearization of TDS Problems for Frequency Response Analysis

Computes the linearized state-space model of a DAE system at steady state.

For a DAE system:
    M·ẋ = f(x, y, u)
    0 = g(x, y, u)

The linearized state-space form is:
    ẋ = A·x + B·u
    z = C·x + D·u

where:
    A = M⁻¹(f_x - f_y·g_y⁻¹·g_x)  (state matrix, from eigenvalue.jl)
    B = M⁻¹(f_u - f_y·g_y⁻¹·g_u)  (input matrix)
    C = output selection matrix
    D = direct feedthrough (usually 0)

This module provides:
- `linearize`: Compute state-space model (returns ControlSystemsBase.StateSpace)
- `compute_input_jacobian`: Compute ∂f/∂u and ∂g/∂u
- `compute_B_matrix`: Compute reduced input matrix B
"""

using LinearAlgebra: diag
using DifferentiationInterface
using ADTypes: AutoForwardDiff
using ControlSystemsBase: ss
using PowerfulCore: create_adcache, get_float_ws
using PowerfulCore: AddressManager, StateVar, get_addresses, get_var_address_type
using PowerfulCore: SystemBuffers, ForwardDiffCache, get_typed_workspace!
using PowerfulCore: ParameterSlice, gather_parameters, scatter_parameters!

# ============================================================================
# Setpoint Inputs
# ============================================================================

"""
    setpoint_slices(adcache, inputs) -> Vector{ParameterSlice}

The setpoint inputs described by `inputs`, a vector of `(model_name, setpoint,
instances)` tuples, as consecutive slices of the input vector `u`: one input per
instance, in the order given. Errors if a model has no workspace in `adcache`, has no
such setpoint, or has no such instance.
"""
function setpoint_slices(adcache, inputs::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{Int}}})
    slices = ParameterSlice[]
    n = 0
    for (model_name, setpoint, instances) in inputs
        slice = ParameterSlice(adcache, model_name, setpoint, instances, n + 1; source = :setpoint)
        push!(slices, slice)
        n = last(slice.range)
    end
    return slices
end

# ============================================================================
# Input Jacobian Computation
# ============================================================================

"""
    InputJacobianFunctor

Functor for computing ∂[f;g]/∂u where u are setpoint parameters.

This functor is called by ForwardDiff to compute Jacobians. It:
1. Writes the input values (potentially Dual) into the workspace of their type
   (`scatter_parameters!`)
2. Evaluates residuals with fixed (x, y) states
3. Packs results into output vector

# Fields
- `adapter_ref`: Ref{Any} wrapping TDSAdapter (type barrier for compilation)
- `adcache_ref`: Ref{Any} wrapping ADCache
- `inputs`: Setpoint inputs (see `setpoint_slices`)
- `n_diff`: Number of differential states
- `n_alg`: Number of algebraic states
- `x_fixed`: Fixed differential state values at equilibrium
- `y_fixed`: Fixed algebraic state values at equilibrium
"""
struct InputJacobianFunctor
    adapter_ref::Ref{Any}
    adcache_ref::Ref{Any}
    inputs::Vector{ParameterSlice}
    n_diff::Int
    n_alg::Int
    x_fixed::Vector{Float64}
    y_fixed::Vector{Float64}
    t::Float64                # the time of the fixed state
end

"""
    (functor::InputJacobianFunctor)(du::AbstractVector{T}, u::AbstractVector{T}) where T

Evaluate residuals [f; g] at fixed states (x, y) with input parameters u.

Writes the inputs into the workspace of their element type, then evaluates residuals.
Works uniformly for Float64 and Dual types.
"""
function (functor::InputJacobianFunctor)(du::AbstractVector{T}, u::AbstractVector{T}) where T
    n_diff = functor.n_diff
    n_alg = functor.n_alg

    # Unwrap from type barrier
    adapter = functor.adapter_ref[]::TDSAdapter
    adcache = functor.adcache_ref[]

    # Write the input values where evaluations of type T read them
    scatter_parameters!(adcache, functor.inputs, u)
    ws = get_typed_workspace!(adcache, T)

    if T === Float64
        # Float64 path: use adapter's preallocated buffers
        f_buf = adapter.f_buf
        g_buf = adapter.g_buf

        evaluate_system_residuals!(ws, f_buf, g_buf,
                                    functor.x_fixed, functor.y_fixed, adapter.status, functor.t)
    else
        # Dual path: use typed buffers from adcache
        buffers = get_tds_dual_buffers(adcache, T, n_diff, n_alg)
        x_buf = buffers.x_buf::Vector{T}
        y_buf = buffers.y_buf::Vector{T}
        f_buf = buffers.f_buf::Vector{T}
        g_buf = buffers.g_buf::Vector{T}

        # Copy fixed states as Dual with zero partials (they don't depend on u)
        @inbounds for i in 1:n_diff
            x_buf[i] = T(functor.x_fixed[i])
        end
        @inbounds for i in 1:n_alg
            y_buf[i] = T(functor.y_fixed[i])
        end

        # Evaluate with Dual workspace
        evaluate_system_residuals!(ws, f_buf, g_buf, x_buf, y_buf, adapter.status, functor.t)
    end

    # Pack: du = [f; g]
    @inbounds for i in 1:n_diff
        du[i] = f_buf[i]
    end
    @inbounds for i in 1:n_alg
        du[n_diff + i] = g_buf[i]
    end

    nothing
end

"""
    compute_input_jacobian(prob::TDSProblem, inputs; verbose=false) -> (f_u, g_u, adapter)

Compute the input Jacobians ∂f/∂u and ∂g/∂u at steady state using ForwardDiff.

Uses automatic differentiation to compute exact derivatives of residuals with
respect to the setpoint inputs.

# Arguments
- `prob::TDSProblem`: TDS problem at equilibrium
- `inputs`: Setpoint inputs as `(model_name, setpoint, instances)` tuples (see `linearize`)

# Keyword Arguments
- `verbose::Bool`: Print progress (default: false)

# Returns
- `f_u::Matrix{Float64}`: ∂f/∂u (n_diff × n_inputs)
- `g_u::Matrix{Float64}`: ∂g/∂u (n_alg × n_inputs)
- `adapter::TDSAdapter`: TDS adapter (for mass matrix access)
"""
function compute_input_jacobian(
    prob::TDSProblem,
    inputs::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{Int}}};
    verbose::Bool = false
)
    @vprint "Computing input Jacobian ∂[f;g]/∂u via ForwardDiff..."

    # Create ADCache and adapter
    verbose && print("  Creating ADCache and TDSAdapter... ")
    t_setup = @elapsed begin
        adcache = create_adcache(prob)
        adapter = TDSAdapter(adcache, prob)  # populates setpoints
    end
    @vprint "done ($(round(t_setup*1000, digits=1)) ms)"

    n_diff = adapter.n_diff
    n_alg = adapter.n_alg
    n_total = n_diff + n_alg
    setpoints = setpoint_slices(adcache, inputs)
    p0 = gather_parameters(adcache, setpoints)
    n_inputs = length(p0)

    @vprint "  States: $n_diff differential, $n_alg algebraic"
    @vprint "  Inputs: $n_inputs parameters"

    # Get steady-state values: the adapter's initial state, the state at the start of the problem
    x_fixed = copy(adapter.x_buf)
    y_fixed = copy(adapter.y_buf)

    # Create functor for ForwardDiff
    functor = InputJacobianFunctor(
        Ref{Any}(adapter),
        Ref{Any}(adcache),
        setpoints,
        n_diff, n_alg,
        x_fixed, y_fixed,
        Float64(prob.tspan[1])
    )

    # Compute Jacobian via ForwardDiff
    verbose && print("  Computing input Jacobian ($n_inputs inputs)... ")
    t_jac = @elapsed begin
        du_temp = zeros(n_total)
        backend = AutoForwardDiff()
        jac_prep = DifferentiationInterface.prepare_jacobian(functor, du_temp, backend, p0)

        J_input = zeros(n_total, n_inputs)
        DifferentiationInterface.jacobian!(functor, du_temp, J_input, jac_prep, backend, p0)
    end
    @vprint "done ($(round(t_jac*1000, digits=1)) ms)"

    # Extract f_u and g_u
    f_u = J_input[1:n_diff, :]
    g_u = J_input[n_diff+1:end, :]

    return f_u, g_u, adapter
end

# ============================================================================
# B Matrix Computation
# ============================================================================

"""
    compute_B_matrix(prob::TDSProblem, inputs; verbose=false) -> ZeroMassReductionResult

Reduce the DAE linearized at the operating point of `prob`, with the setpoint inputs
`inputs` (see `compute_input_jacobian`): the state matrix `A`, the input matrix `B = M⁻¹(f_u - f_y·g_y⁻¹·g_u)` and
the sensitivities of the algebraic variables and zero-mass states to the retained
states and the inputs (see `reduce_zero_mass_states`).

# Arguments
- `prob::TDSProblem`: TDS problem at equilibrium
- `inputs`: Setpoint inputs as `(model_name, setpoint, instances)` tuples
"""
function compute_B_matrix(
    prob::TDSProblem,
    inputs::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{Int}}};
    verbose::Bool = false
)
    @vprint "="^60
    @vprint "Computing B Matrix (Input Jacobian)"
    @vprint "="^60

    # Step 1: Compute input Jacobian
    f_u, g_u, adapter = compute_input_jacobian(prob, inputs; verbose=verbose)

    n_diff = adapter.n_diff
    am = prob.address_manager

    # Step 2: Compute state Jacobian
    verbose && print("  Computing state Jacobian... ")

    t_jac = @elapsed J_full, _ = dae_jacobian(adapter, prob.tspan[1])
    @vprint "done ($(round(t_jac*1000, digits=1)) ms)"

    # Extract sub-matrices
    jac_sub = extract_jacobian_submatrices(J_full, n_diff)

    # Get mass matrix
    M = adapter.mass_matrix
    M_diag = diag(M)[1:n_diff]

    # Collect state names
    state_names = collect_state_names(am, n_diff)

    # Step 3: Compute reduced matrices with zero-mass handling
    verbose && print("  Computing A and B matrices... ")
    t_ab = @elapsed begin
        reduction = reduce_zero_mass_states(
            jac_sub.fx, jac_sub.fy, jac_sub.gx, jac_sub.gy, M_diag, state_names;
            f_u = f_u, g_u = g_u
        )
    end
    @vprint "done ($(round(t_ab*1000, digits=1)) ms, $(length(reduction.zero_mass_indices)) zero-mass states)"

    @vprint "  Final dimensions: A=$(size(reduction.A)), B=$(size(reduction.B))"
    @vprint "="^60

    return reduction
end

# ============================================================================
# Linearization (Main Entry Point)
# ============================================================================

"""
    linearize(prob::TDSProblem, inputs, outputs; verbose=false) -> StateSpace

Linearize a TDS problem about its steady-state operating point, returning a
continuous-time state-space model for use with ControlSystemsBase.jl: frequency
response (Bode, Nyquist plots), stability analysis, or controller design.

# Arguments
- `prob::TDSProblem`: a TDS problem already at equilibrium (e.g. `prob.tspan[1]`, before
  any event).
- `inputs`, `outputs`: each a `Vector` of `(model_name, field_name, instances)` tuples.
  - `model_name`, `field_name`: the model type and field, as `Symbol`s — for example
    `(:TGOV1, :pref0, ...)` for a governor's power reference, or `(:GENCLS, :omega, ...)`
    for a classical machine's rotor speed. An input must be a setpoint parameter of the
    model (not a state variable); an output may be any variable: a differential state, a
    controller state with a zero time constant, or an algebraic variable, including a
    network quantity the model reads such as its terminal voltage `:v`.
  - `instances`: which model instances to use, as a range or vector of integer indices
    into that model's array (1-based). `1:1` selects the first instance; `[1, 3, 5]`
    selects three specific ones. To find the index for a particular bus or generator
    (rather than "the first one"), use `PowerfulCore.get_index_by_key(prob.sys, :GENCLS,
    (bus_number, id))`, which returns the array index for that external key.

# Keyword Arguments
- `verbose::Bool`: print progress (default: false)

# Example
```julia
using ControlSystemsBase: bode, bodeplot

prob = TDSProblem(init; tspan=(0.0, 10.0))

sys = linearize(prob,
    [(:TGOV1, :pref0, 1:1)],    # input: first governor's power reference
    [(:GENCLS, :omega, 1:1)]    # output: first machine's rotor speed
)

frequencies = 10.0 .^ range(-2, 2, length=100)
mag, phase, w = bode(sys, frequencies)
```

# Extended help
An output that is not one of the retained differential states is expressed through
those states and the inputs (`w = -(W⁻¹W_x x_r + W⁻¹W_u u)`, see
`reduce_zero_mass_states`), so such an output carries a direct feedthrough term `D`.
"""
function linearize(
    prob::TDSProblem,
    inputs::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{Int}}},
    outputs::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{Int}}};
    verbose::Bool = false
)
    @vprint "="^60
    @vprint "Linearizing TDS System at Steady State"
    @vprint "="^60

    # Step 1: Compute A and B matrices for the setpoint inputs
    reduction = compute_B_matrix(prob, inputs; verbose=verbose)
    @vprint "  $(size(reduction.B, 2)) input(s)"

    # Step 2: Output matrices C and D
    verbose && print("  Setting up outputs... ")
    C, D = _build_output_matrices(outputs, prob.address_manager, reduction)
    @vprint "$(size(C, 1)) output(s)"

    @vprint "  State-space: A=$(size(reduction.A)), B=$(size(reduction.B)), C=$(size(C)), D=$(size(D))"
    @vprint "="^60

    # Return ControlSystemsBase StateSpace object
    return ss(reduction.A, reduction.B, C, D)
end

"""
    _build_output_matrices(outputs, am, reduction) -> (C, D)

Output matrices of the reduced system `ẋ_r = A x_r + B u`, one row per output
instance. A retained state is read directly; an algebraic variable or a zero-mass
state is an eliminated quantity `w_k = -(W⁻¹W_x)[k, :] x_r - (W⁻¹W_u)[k, :] u`.
"""
function _build_output_matrices(
    outputs::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{Int}}},
    am::AddressManager,
    reduction::ZeroMassReductionResult
)
    n_states = length(reduction.retained_indices)
    n_alg = size(reduction.W_inv_Wx, 1) - length(reduction.zero_mass_indices)
    C_rows = Vector{Float64}[]
    D_rows = Vector{Float64}[]

    for (model_name, field_name, instances) in outputs
        addr_type = get_var_address_type(am, model_name, field_name)
        addrs = get_addresses(am, addr_type, model_name, field_name)
        for inst in instances
            checkbounds(Bool, addrs, inst) || throw(ArgumentError(
                "Output ($model_name, $field_name): instance $inst out of range 1:$(length(addrs))"))
            addr = Int(addrs[inst])

            state_pos = addr_type isa StateVar ? findfirst(==(addr), reduction.retained_indices) : nothing
            if state_pos !== nothing
                c_row = zeros(n_states)
                c_row[state_pos] = 1.0
                push!(C_rows, c_row)
                push!(D_rows, zeros(size(reduction.B, 2)))
            else
                # Eliminated: algebraic variables first, then zero-mass states
                k = addr_type isa StateVar ?
                    n_alg + findfirst(==(addr), reduction.zero_mass_indices) : addr
                push!(C_rows, -reduction.W_inv_Wx[k, :])
                push!(D_rows, -reduction.W_inv_Wu[k, :])
            end
        end
    end

    isempty(C_rows) && throw(ArgumentError("linearize: no outputs given"))
    return permutedims(reduce(hcat, C_rows)), permutedims(reduce(hcat, D_rows))
end

# ============================================================================
# Exports
# ============================================================================

export linearize
