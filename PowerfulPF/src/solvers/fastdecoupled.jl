"""
Fast Decoupled Power Flow solver - Allocation-free implementation.

Exploits weak P-θ and Q-V coupling for 2-3x speedup on well-conditioned systems.

# Key Features
- Constant B' and B" matrices (factored once, reused across iterations)
- Decoupled P-θ and Q-V subsystems (smaller linear solves)
- Zero allocations in iteration loop
- KLU sparse direct solver with cached factorizations

# References
- Stott & Alsac (1974). "Fast Decoupled Load Flow". IEEE Trans. PAS-93(3).
"""

using LinearAlgebra
using SparseArrays
using KLU
using PowerfulCore: @vprint

# ============================================================================
# Workspace Structure
# ============================================================================

"""
    FastDecoupledWorkspace{PROB, AT, BT}

Pre-allocated workspace for allocation-free Fast Decoupled Power Flow.

All buffers allocated during setup, iterations are allocation-free.

# Type Parameters
- `PROB`: ACPowerFlowProblem type
- `AT`: Array type (Vector{Float64} or CuArray)
- `BT`: B-matrix type (SparseMatrixCSC)

# Fields
- `prob::PROB`: Power flow problem (sys, address_manager, formulation)
- `adcache::ForwardDiffCache`: Residual evaluator for mismatch computation
- `iter_dae::DAE{AT}`: Current iteration state (θ, V in y vector)
- `B_prime::BT`: Constant P-θ coupling matrix
- `B_double_prime::BT`: Constant Q-V coupling matrix
- `Bp_solver::KLULinearSolver`: KLU factorization cache of B' (with the setup's ordering)
- `Bpp_solver::KLULinearSolver`: KLU factorization cache of B"
- `mismatch_P::AT`: Active power mismatch buffer (ΔP)
- `mismatch_Q::AT`: Reactive power mismatch buffer (ΔQ)
- `rhs_P::AT`: Right-hand side for P-θ solve (ΔP/V)
- `rhs_Q::AT`: Right-hand side for Q-V solve (ΔQ/V)
- `delta_theta::AT`: Angle correction (output buffer)
- `delta_V::AT`: Voltage magnitude correction (output buffer)
- `variant::Symbol`: :XB or :BX
- `pq_bus_indices::Vector{Int}`: Internal indices of PQ buses
- `pv_bus_indices::Vector{Int}`: Internal indices of PV buses
- `pqpv_bus_indices::Vector{Int}`: PQ and PV buses in ascending order (rows of B')
- `slack_bus_index::Int`: Internal index of slack bus
- `theta_addrs::Vector{UInt64}`: DAE addresses for θ variables
- `v_addrs::Vector{UInt64}`: DAE addresses for V variables
- `P_res_addrs::Vector{UInt64}`: DAE addresses for P residuals
- `Q_res_addrs::Vector{UInt64}`: DAE addresses for Q residuals
"""
struct FastDecoupledWorkspace{PROB, AT<:AbstractVector, BT<:AbstractMatrix}
    # Problem and evaluators
    prob::PROB
    adcache::ForwardDiffCache
    iter_dae::DAE{AT}

    # Constant matrices (built once, factored once)
    B_prime::BT
    B_double_prime::BT
    Bp_solver::KLULinearSolver  # KLU factorization cache (factored on first solve)
    Bpp_solver::KLULinearSolver

    # Pre-allocated iteration buffers (allocation-free hot loop)
    mismatch_P::AT  # Full P mismatch vector
    mismatch_Q::AT  # Full Q mismatch vector
    rhs_P::AT       # P-θ RHS: ΔP/V (reduced: excludes slack)
    rhs_Q::AT       # Q-V RHS: ΔQ/V (reduced: PQ buses only)
    delta_theta::AT # P-θ solution (reduced size)
    delta_V::AT     # Q-V solution (reduced size)

    # Configuration
    variant::Symbol

    # Bus type information (internal indices)
    pq_bus_indices::Vector{Int}
    pv_bus_indices::Vector{Int}
    pqpv_bus_indices::Vector{Int}  # PQ and PV buses in ascending order: the rows of B'
    slack_bus_index::Int

    # Address mappings (DAE positions)
    theta_addrs::Vector{UInt64}  # All buses
    v_addrs::Vector{UInt64}      # All buses
    P_res_addrs::Vector{UInt64}  # All buses
    Q_res_addrs::Vector{UInt64}  # All buses

    # Model status for residual evaluation
    status::Dict{Symbol, Vector{Bool}}
end

# ============================================================================
# Setup Function
# ============================================================================

"""
    setup_powerflow_fastdecoupled(prob::ACPowerFlowProblem; kwargs...)
        -> FastDecoupledWorkspace

One-time setup for Fast Decoupled Power Flow with allocation-free iterations.

Constructs and factors B' and B" matrices once. All subsequent iterations
reuse these factorizations with zero allocations.

# Arguments
- `prob::ACPowerFlowProblem`: Power flow problem (must use PowerInjection formulation)

# Keyword Arguments
- `variant::Symbol=:XB`: Variant (:XB neglects R, :BX includes R)
- `ordering::Symbol=:AMD`: KLU fill-reducing ordering (`:AMD`, `:COLAMD`, `:natural`)
- `verbose::Bool=false`: Print setup progress

# Returns
- `FastDecoupledWorkspace`: Pre-allocated workspace for allocation-free iterations

# Example
```julia
prob = ACPowerFlowProblem(sys, PowerInjection())
workspace = setup_powerflow_fastdecoupled(prob; variant=:XB, verbose=true)

y0 = copy(workspace.iter_dae.y)
sol = solve_powerflow_fastdecoupled!(workspace, y0; abstol=1e-6)
```

# Performance
Setup cost is amortized over iterations:
- Setup: ~5-10ms (one-time)
- Iteration: ~1-2ms (allocation-free)

For repeated solves, workspace can be reused with different initial conditions.
"""
function setup_powerflow_fastdecoupled(
    prob::ACPowerFlowProblem;
    variant::Symbol = :XB,
    ordering::Symbol = :AMD,
    verbose::Bool = false
)
    # Validate inputs
    if !(prob.formulation isa PowerInjection)
        error("Fast Decoupled Power Flow requires PowerInjection formulation, got $(typeof(prob.formulation))")
    end

    if !(variant in [:XB, :BX])
        error("Fast Decoupled variant must be :XB or :BX, got :$variant")
    end

    sys = prob.sys
    am = prob.address_manager

    @vprint "="^60
    @vprint "Fast Decoupled Power Flow Setup"
    @vprint "Variant: $variant"
    @vprint "="^60

    # Step 1: Create the ForwardDiffCache and DAE for residual evaluation
    verbose && print("\nStep 1/6: Creating residual workspace... ")
    t_cache = @elapsed begin
        adcache = create_adcache(am, sys, prob.formulation, ForwardDiffBackend())
        iter_dae = DAE(am; array_type=sys.array_type)
        set_algebraic_guess!(iter_dae.y, am, sys, prob.formulation)
        status = build_status_dict(sys)
    end
    @vprint "done ($(round(t_cache*1000, digits=1)) ms)"

    # Step 2: Identify bus types from sys.models.Bus
    verbose && print("Step 2/6: Identifying bus types... ")
    t_bustypes = @elapsed begin
        pq_indices = get_pq_buses(sys)
        pv_indices = get_pv_buses(sys)
        slack_index = get_slack_bus(sys)
    end
    n_pq = length(pq_indices)
    n_pv = length(pv_indices)
    n_bus = length(sys.models.Bus)
    @vprint "done ($(round(t_bustypes*1000, digits=1)) ms)"
    @vprint "  Total buses: $n_bus"
    @vprint "  PQ buses: $n_pq"
    @vprint "  PV buses: $n_pv"
    @vprint "  Slack bus: 1"

    # Step 3: Get address mappings for Bus variables and residuals
    verbose && print("Step 3/6: Extracting address mappings... ")
    t_addrs = @elapsed begin
        theta_addrs, v_addrs, P_res_addrs, Q_res_addrs = _get_bus_addresses(am, sys)
    end
    @vprint "done ($(round(t_addrs*1000, digits=1)) ms)"

    # Step 4: Build B' matrix (P-θ subsystem, excludes slack)
    verbose && print("Step 4/6: Building B' matrix (P-θ subsystem)... ")
    t_Bp = @elapsed begin
        B_prime = _build_B_prime_matrix(sys, variant, slack_index)
    end
    n_P_eqs = n_pq + n_pv  # All non-slack buses
    @vprint "done ($(round(t_Bp*1000, digits=1)) ms)"
    @vprint "  Size: $(size(B_prime, 1)) x $(size(B_prime, 2)) (reduced from $n_bus)"
    @vprint "  Nonzeros: $(nnz(B_prime))"
    @vprint "  Density: $(round(100.0 * nnz(B_prime) / (size(B_prime,1)^2), digits=2))%"

    # Step 5: Build B" matrix (Q-V subsystem, PQ buses only)
    # For FDXB: B' uses :XB (neglect r), B" uses :BX (include r)
    variant_Bpp = :BX
    verbose && print("Step 5/6: Building B\" matrix (Q-V subsystem, variant=$variant_Bpp)... ")
    t_Bpp = @elapsed begin
        B_double_prime = _build_B_double_prime_matrix(sys, variant_Bpp, pq_indices)
    end
    @vprint "done ($(round(t_Bpp*1000, digits=1)) ms)"
    @vprint "  Size: $(size(B_double_prime, 1)) x $(size(B_double_prime, 2)) (PQ only)"
    @vprint "  Nonzeros: $(nnz(B_double_prime))"
    @vprint "  Density: $(round(100.0 * nnz(B_double_prime) / (size(B_double_prime,1)^2), digits=2))%"

    # Step 6: Allocate iteration buffers (allocation-free hot loop)
    verbose && print("Step 6/6: Allocating iteration buffers... ")
    t_buf = @elapsed begin
        # Full mismatch vectors (all buses)
        mismatch_P = zeros(Float64, n_bus)
        mismatch_Q = zeros(Float64, n_bus)

        # Reduced RHS and solution vectors
        rhs_P = zeros(Float64, n_P_eqs)       # ΔP/V (PQ + PV buses)
        rhs_Q = zeros(Float64, n_pq)          # ΔQ/V (PQ buses only)
        delta_theta = zeros(Float64, n_P_eqs) # Δθ (PQ + PV buses)
        delta_V = zeros(Float64, n_pq)        # ΔV (PQ buses only)
    end
    @vprint "done ($(round(t_buf*1000, digits=1)) ms)"
    @vprint "  Allocated: 2x$n_bus + 2x$n_P_eqs + 2x$n_pq = $(2*n_bus + 2*n_P_eqs + 2*n_pq) elements"

    @vprint "\n" * "="^60
    @vprint "Setup complete! Ready for allocation-free iterations."
    @vprint "="^60
    @vprint ""

    return FastDecoupledWorkspace(
        prob,
        adcache,
        iter_dae,
        B_prime,
        B_double_prime,
        KLULinearSolver(ordering),  # Bp_solver (factored on first solve)
        KLULinearSolver(ordering),  # Bpp_solver
        mismatch_P,
        mismatch_Q,
        rhs_P,
        rhs_Q,
        delta_theta,
        delta_V,
        variant,
        pq_indices,
        pv_indices,
        sort!(vcat(pq_indices, pv_indices)),
        slack_index,
        theta_addrs,
        v_addrs,
        P_res_addrs,
        Q_res_addrs,
        status
    )
end

# ============================================================================
# Allocation-Free Iteration Function
# ============================================================================

"""
    solve_powerflow_fastdecoupled!(workspace::FastDecoupledWorkspace, initial_guess; kwargs...)
        -> PowerFlowSolution

Solve power flow using Fast Decoupled method with **zero allocations** per iteration.

All buffers pre-allocated in workspace, factorizations reused across iterations.

# Arguments
- `workspace`: Pre-allocated workspace from `setup_powerflow_fastdecoupled`
- `initial_guess`: Initial guess for variables (θ, V)

# Keyword Arguments
- `abstol::Real=1e-4`: Absolute tolerance on max(|ΔP|, |ΔQ|)
- `maxiters::Int=100`: Maximum iterations
- `verbose::Bool=false`: Print iteration progress

# Returns
- `PowerFlowSolution`: Solution with convergence status

# Example
```julia
workspace = setup_powerflow_fastdecoupled(prob; variant=:XB)
y0 = copy(workspace.iter_dae.y)

# Allocation-free solve
sol = solve_powerflow_fastdecoupled!(workspace, y0; abstol=1e-6)

# Benchmark iterations only (no setup overhead)
using BenchmarkTools
@btime solve_powerflow_fastdecoupled!(\$workspace, \$y0; abstol=1e-6, verbose=false)
```
"""
function solve_powerflow_fastdecoupled!(
    workspace::FastDecoupledWorkspace,
    initial_guess::AbstractVector;
    abstol::Real = 1e-4,
    maxiters::Int = 100,
    verbose::Bool = false
)
    # Unpack workspace (no allocations)
    (; prob, adcache, iter_dae, B_prime, B_double_prime,
       Bp_solver, Bpp_solver, mismatch_P, mismatch_Q, rhs_P, rhs_Q,
       delta_theta, delta_V, variant,
       pq_bus_indices, pv_bus_indices, pqpv_bus_indices,
       theta_addrs, v_addrs, P_res_addrs, Q_res_addrs, status) = workspace

    # Initialize iteration state (in-place copy)
    copyto!(iter_dae.y, initial_guess)

    # CRITICAL: FDPF requires PV bus voltages to be fixed at their setpoints
    # Overwrite PV bus voltages with their generator setpoints (not flat-start values)
    sys = prob.sys
    if haskey(sys.models, :PV)
        pv = sys.models.PV
        pv_bus_idx = get_foreign_indices(sys, :PV, :i)  # PV.i → internal Bus indices
        active = in_service(pv)
        for i in eachindex(pv)
            if active[i]  # Only for active PV generators
                iter_dae.y[v_addrs[pv_bus_idx[i]]] = pv.vs[i]  # Use PV setpoint
            end
        end
    end

    @vprint "\n" * "="^60
    @vprint "Fast Decoupled Power Flow Iterations"
    @vprint "Variant: $variant, abstol: $abstol"
    @vprint "="^60
    @vprint ""

    # Factor matrices on first solve (amortized cost), with the setup's KLU ordering
    if !is_factorized(Bp_solver)
        verbose && print("Factoring B' matrix... ")
        t_fact = @elapsed begin
            factorize!(Bp_solver, B_prime)
        end
        @vprint "done ($(round(t_fact*1000, digits=1)) ms, one-time cost)"
    end

    if !is_factorized(Bpp_solver)
        verbose && print("Factoring B\" matrix... ")
        t_fact = @elapsed begin
            factorize!(Bpp_solver, B_double_prime)
        end
        @vprint "done ($(round(t_fact*1000, digits=1)) ms, one-time cost)"
        @vprint ""
    end

    local max_mismatch = Inf
    local converged = false
    local final_iter = 0

    for iter in 1:maxiters
        final_iter = iter
        @vprint "Iteration $iter:"

        # 1. Compute power mismatches using the ForwardDiffCache (in-place)
        verbose && print("  [1/5] Computing mismatches... ")
        t_mis = @elapsed begin
            # Evaluate g(y) in-place, at t = 0: the operating point, whose residuals do not
            # depend on time
            adcache(iter_dae.g, iter_dae.y, status, 0.0)

            # Extract P and Q mismatches from residuals (in-place)
            _extract_mismatches!(mismatch_P, mismatch_Q, iter_dae.g,
                                P_res_addrs, Q_res_addrs)

            # - Check P mismatches at PQ + PV buses (all non-slack)
            # - Check Q mismatches at PQ buses ONLY (PV buses maintain voltage setpoint)
            # (isolated buses are neither PQ nor PV; an empty set has zero mismatch)
            max_P = max(maximum(i -> abs(mismatch_P[i]), pq_bus_indices; init=zero(eltype(mismatch_P))),
                        maximum(i -> abs(mismatch_P[i]), pv_bus_indices; init=zero(eltype(mismatch_P))))
            max_Q = maximum(i -> abs(mismatch_Q[i]), pq_bus_indices; init=zero(eltype(mismatch_Q)))
            max_mismatch = max(max_P, max_Q)
        end
        @vprint "done ($(round(t_mis*1000, digits=2)) ms)"
        @vprint "        max|ΔP|, max|ΔQ| = $(max_mismatch)"

        # Check convergence
        if max_mismatch < abstol
            @vprint ""
            @vprint "="^60
            @vprint "✓ CONVERGED in $iter iterations!"
            @vprint "  Final max mismatch: $(max_mismatch)"
            @vprint "="^60
            converged = true
            break
        end

        # 2. Solve P-θ subsystem: B' Δθ = ΔP/V (in-place)
        verbose && print("  [2/5] Solving P-θ subsystem... ")
        t_Ptheta = @elapsed begin
            _solve_P_theta_subsystem!(delta_theta, mismatch_P, rhs_P,
                                     iter_dae.y, v_addrs,
                                     Bp_solver, pqpv_bus_indices)
        end
        @vprint "done ($(round(t_Ptheta*1000, digits=2)) ms)"

        # 3. Update angles (in-place)
        verbose && print("  [3/5] Updating θ... ")
        t_update_theta = @elapsed begin
            _update_angles!(iter_dae.y, delta_theta, theta_addrs, pqpv_bus_indices)
        end
        @vprint "done ($(round(t_update_theta*1000, digits=2)) ms)"

        # 4. Solve Q-V subsystem: B" ΔV = ΔQ/V (in-place)
        verbose && print("  [4/5] Solving Q-V subsystem... ")
        t_QV = @elapsed begin
            _solve_Q_V_subsystem!(delta_V, mismatch_Q, rhs_Q,
                                 iter_dae.y, v_addrs,
                                 Bpp_solver, pq_bus_indices)
        end
        @vprint "done ($(round(t_QV*1000, digits=2)) ms)"

        # 5. Update voltage magnitudes (in-place, PQ buses only)
        verbose && print("  [5/5] Updating V... ")
        t_update_V = @elapsed begin
            _update_voltages!(iter_dae.y, delta_V, v_addrs, pq_bus_indices)
        end
        @vprint "done ($(round(t_update_V*1000, digits=2)) ms)"

        @vprint "  Total iteration time: $(round((t_mis + t_Ptheta + t_update_theta + t_QV + t_update_V)*1000, digits=2)) ms"
        @vprint ""
    end

    if !converged
        @vprint "="^60
        @vprint "✗ Did NOT converge in $maxiters iterations"
        @vprint "  Final max mismatch: $(max_mismatch)"
        @vprint "="^60
    end

    return PowerFlowSolution(
        converged,
        max_mismatch,
        final_iter,
        converged ? :Success : :MaxIters,
        nothing,
        prob,
        deepcopy(iter_dae),  # Deep copy to preserve solution state
        Ref{Union{Nothing, Vector{eltype(iter_dae.y)}}}(nothing)  # Lazy cache for generator residuals
    )
end

# ============================================================================
# Convenience Wrapper
# ============================================================================

"""
    solve_powerflow_fastdecoupled(prob::ACPowerFlowProblem; kwargs...) -> PowerFlowSolution
    solve_powerflow_fastdecoupled(sys::SystemModel; initial_guess=nothing, kwargs...) -> PowerFlowSolution

Convenience wrapper: setup + solve in one call. The `SystemModel` form builds a
PowerInjection problem; the problem form keeps the problem's addressing and
`initial_guess`.

Keywords: `variant=:XB`, `abstol=1e-4`, `maxiters=100`, `verbose=false`, `ordering=:AMD`.

For repeated solves or benchmarking, use explicit workspace API instead.

# Example
```julia
sys = SystemModel(case)
sol = solve_powerflow_fastdecoupled(sys; variant=:XB, abstol=1e-6, verbose=true)
```
"""
function solve_powerflow_fastdecoupled(
    prob::ACPowerFlowProblem;
    variant::Symbol = :XB,
    abstol::Real = 1e-4,
    maxiters::Int = 100,
    verbose::Bool = false,
    ordering::Symbol = :AMD
)
    workspace = setup_powerflow_fastdecoupled(prob; variant=variant,
                                             ordering=ordering,
                                             verbose=verbose)

    return solve_powerflow_fastdecoupled!(workspace, _initial_guess(prob, workspace);
                                         abstol=abstol,
                                         maxiters=maxiters,
                                         verbose=false)  # Don't double-print
end

function solve_powerflow_fastdecoupled(
    sys::SystemModel;
    initial_guess::Union{Nothing, Vector} = nothing,
    kwargs...
)
    prob = ACPowerFlowProblem(sys, PowerInjection(); initial_guess=initial_guess)
    return solve_powerflow_fastdecoupled(prob; kwargs...)
end

# ============================================================================
# Helper Functions (Allocation-Free Implementations)
# ============================================================================

"""
    build_subsystem_rhs!(rhs, mismatch, voltage, v_addrs, pqpv_indices)

Build right-hand side for P-θ subsystem: `rhs[k] = ΔP/V` of bus `pqpv_indices[k]`.

Maps from full system (all buses) to reduced subsystem (PQ + PV buses in ascending
order, the row order of B').
"""
@inline function build_subsystem_rhs!(
    rhs::AbstractVector,
    mismatch::AbstractVector,
    voltage::AbstractVector,
    v_addrs::Vector{UInt64},
    pqpv_indices::Vector{Int}
)
    @inbounds for (k, bus_idx) in enumerate(pqpv_indices)
        rhs[k] = mismatch[bus_idx] / voltage[v_addrs[bus_idx]]
    end
    nothing
end

"""
    scatter_subsystem_update!(y, delta, addrs, pqpv_indices)

Scatter subsystem solution back to full state vector: `y[addrs[bus]] -= delta[k]` for
bus `pqpv_indices[k]` (PQ + PV buses in ascending order).
"""
@inline function scatter_subsystem_update!(
    y::AbstractVector,
    delta::AbstractVector,
    addrs::Vector{UInt64},
    pqpv_indices::Vector{Int}
)
    @inbounds for (k, bus_idx) in enumerate(pqpv_indices)
        y[addrs[bus_idx]] -= delta[k]
    end
    nothing
end

"""
    _get_bus_addresses(am::AddressManager, sys::SystemModel)
        -> (theta_addrs, v_addrs, P_res_addrs, Q_res_addrs)

Extract DAE addresses for Bus variables and residuals.

Returns vectors of addresses indexed by internal bus index.
"""
function _get_bus_addresses(am::AddressManager, ::SystemModel)
    # Get address assignments for Bus model
    bus_assignments = am.assignments[:Bus]

    # Variable addresses (AlgebVar)
    theta_addrs = collect(bus_assignments.by_type[AlgebVar()][:theta])
    v_addrs = collect(bus_assignments.by_type[AlgebVar()][:v])

    # Residual addresses (AlgebRes)
    P_res_addrs = collect(bus_assignments.by_type[AlgebRes()][:p])
    Q_res_addrs = collect(bus_assignments.by_type[AlgebRes()][:q])

    return theta_addrs, v_addrs, P_res_addrs, Q_res_addrs
end

"""
    _series_susceptance(r, x, variant) -> Float64

Series susceptance magnitude used by FDPF: `1/x` for XB (resistance neglected),
`x/(r²+x²)` for BX.
"""
@inline function _series_susceptance(r, x, variant::Symbol)
    variant === :XB && return 1.0 / x
    variant === :BX && return x / (r^2 + x^2)
    error("Unknown variant: $variant")
end

"""
    _series_B(sys, variant) -> SparseMatrixCSC

Series-susceptance matrix of in-service Branch and Transformer2W devices in FULL
bus coordinates (including slack), built from the incidence matrix.

B approximates J_Pθ: negative diagonals, positive off-diagonals.
"""
function _series_B(sys::SystemModel, variant::Symbol)
    A, branch_indices = build_incidence_matrix(sys)
    b = [-_series_susceptance(series_impedance(sys, model_name, k)..., variant)
         for (model_name, k) in branch_indices]
    return build_B_from_incidence(A, b)
end

"""
    _build_B_prime_matrix(sys, variant, slack_index) -> SparseMatrixCSC

Build B' matrix for P-θ subsystem (excludes slack bus).

The FULL n×n matrix is built first so that non-slack buses connected to the slack
keep those branches' series susceptance in their diagonal; then the slack
row/column is removed. Shunts are added to B" only.

Matrix size: (n_pq + n_pv) × (n_pq + n_pv)
"""
function _build_B_prime_matrix(sys::SystemModel, variant::Symbol, slack_index)
    non_slack_indices = [i for i in 1:length(sys.models.Bus) if i != slack_index]
    return _series_B(sys, variant)[non_slack_indices, non_slack_indices]
end

"""
    _build_B_double_prime_matrix(sys, variant, pq_indices) -> SparseMatrixCSC

Build B" matrix for Q-V subsystem (PQ buses only).

Built for ALL buses first (a PQ bus's diagonal includes every connected branch,
even if the other endpoint is PV or Slack), then the PQ×PQ submatrix is extracted.
Adds branch charging/end shunts and the susceptance of the fixed shunt models
(`fixed_shunt_models`: Shunt, SwitchedShunt) to the diagonal.

Matrix size: n_pq × n_pq
"""
function _build_B_double_prime_matrix(sys::SystemModel, variant::Symbol, pq_indices)
    shunt_diag = zeros(length(sys.models.Bus))

    branches = sys.models.Branch
    if !isempty(branches)
        from_idx = get_foreign_indices(sys, :Branch, :i)
        to_idx = get_foreign_indices(sys, :Branch, :j)
        active = in_service(branches)
        for k in eachindex(branches)
            active[k] || continue
            shunt_diag[from_idx[k]] -= branches.b[k] / 2.0 + branches.bi[k]
            shunt_diag[to_idx[k]] -= branches.b[k] / 2.0 + branches.bj[k]
        end
    end

    for name in fixed_shunt_models(sys)
        shunts = sys.models[name]
        params = shunt_admittance_params(shunts)
        bus_idx = get_foreign_indices(sys, name, :i)
        active = in_service(shunts)
        for k in eachindex(shunts)
            active[k] || continue
            shunt_diag[bus_idx[k]] -= shunt_admittance(params, k)[2]
        end
    end

    B_full = _series_B(sys, variant) + spdiagm(0 => shunt_diag)
    return B_full[pq_indices, pq_indices]
end

"""
    _extract_mismatches!(mismatch_P, mismatch_Q, g, P_res_addrs, Q_res_addrs)

Extract P and Q mismatches from full residual vector (in-place, allocation-free).

The residual g represents the power balance equation: g = P_calc - P_spec.
The mismatch for FDPF is ΔP = P_spec - P_calc = -g (to drive g to zero).
"""
@inline function _extract_mismatches!(mismatch_P::AbstractVector, mismatch_Q::AbstractVector,
                                     g::AbstractVector,
                                     P_res_addrs::Vector{UInt64}, Q_res_addrs::Vector{UInt64})
    n_bus = length(P_res_addrs)

    @inbounds for i in 1:n_bus
        mismatch_P[i] = g[P_res_addrs[i]]
        mismatch_Q[i] = g[Q_res_addrs[i]]
    end

    nothing
end

"""
    _solve_P_theta_subsystem!(delta_theta, mismatch_P, rhs_P, y, v_addrs, Bp_solver, pqpv_indices)

Solve P-θ subsystem: B' Δθ = ΔP/V (in-place, allocation-free).

Steps:
1. Build RHS: rhs_P = ΔP/V for PQ and PV buses (excludes slack)
2. Solve: delta_theta = B' \\ rhs_P using cached factorization
"""
@inline function _solve_P_theta_subsystem!(delta_theta::AbstractVector,
                                          mismatch_P::AbstractVector,
                                          rhs_P::AbstractVector,
                                          y::AbstractVector,
                                          v_addrs::Vector{UInt64},
                                          Bp_solver::KLULinearSolver,
                                          pqpv_indices::Vector{Int})
    # Build RHS: rhs_P[k] = mismatch_P[bus_idx] / V[bus_idx]
    # CRITICAL: Order must match B' matrix row order (PQ + PV buses, ascending)
    build_subsystem_rhs!(rhs_P, mismatch_P, y, v_addrs, pqpv_indices)

    linear_solve!(Bp_solver, delta_theta, rhs_P)

    nothing
end

"""
    _update_angles!(y, delta_theta, theta_addrs, pqpv_indices)

Update voltage angles (in-place, allocation-free).

θ[bus] += Δθ for PQ and PV buses (slack angle is reference, doesn't change)
"""
@inline function _update_angles!(y::AbstractVector,
                                delta_theta::AbstractVector,
                                theta_addrs::Vector{UInt64},
                                pqpv_indices::Vector{Int})
    # delta_theta follows the B' row order (PQ + PV buses, ascending); slack unchanged
    scatter_subsystem_update!(y, delta_theta, theta_addrs, pqpv_indices)

    nothing
end

"""
    _solve_Q_V_subsystem!(delta_V, mismatch_Q, rhs_Q, y, v_addrs, Bpp_solver, pq_indices)

Solve Q-V subsystem: B" ΔV = ΔQ/V (in-place, allocation-free).

The solution delta_V represents the absolute voltage change ΔV.
Only PQ buses (PV voltages are fixed by control).
"""
@inline function _solve_Q_V_subsystem!(delta_V::AbstractVector,
                                      mismatch_Q::AbstractVector,
                                      rhs_Q::AbstractVector,
                                      y::AbstractVector,
                                      v_addrs::Vector{UInt64},
                                      Bpp_solver::KLULinearSolver,
                                      pq_indices::Vector{Int})
    # Build RHS: rhs_Q[k] = mismatch_Q[bus_idx] / V[bus_idx] for PQ buses
    @inbounds for i in eachindex(pq_indices)
        bus_idx = pq_indices[i]
        V = y[v_addrs[bus_idx]]
        rhs_Q[i] = mismatch_Q[bus_idx] / V
    end

    linear_solve!(Bpp_solver, delta_V, rhs_Q)

    nothing
end

"""
    _update_voltages!(y, delta_V, v_addrs, pq_indices)

Update voltage magnitudes for PQ buses (in-place, allocation-free).

delta_V contains the voltage magnitude correction: V_new = V_old + ΔV
"""
@inline function _update_voltages!(y::AbstractVector,
                                  delta_V::AbstractVector,
                                  v_addrs::Vector{UInt64},
                                  pq_indices::Vector{Int})
    @inbounds for i in eachindex(pq_indices)
        bus_idx = pq_indices[i]
        addr = v_addrs[bus_idx]
        y[addr] -= delta_V[i]
    end

    nothing
end

# ============================================================================
# Exports
# ============================================================================

export FastDecoupledWorkspace,
       setup_powerflow_fastdecoupled,
       solve_powerflow_fastdecoupled!,
       solve_powerflow_fastdecoupled
