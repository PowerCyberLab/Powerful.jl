"""
PowerfulPF NonlinearSolve.jl Extension

Provides NonlinearSolve.jl-based power flow solver backend (RobustNR).
Loaded automatically when `using NonlinearSolve` is called.
"""
module PowerfulPFNonlinearSolveExt

using PowerfulPF
using NonlinearSolve: NonlinearSolve, NonlinearFunction, NonlinearProblem, ReturnCode, remake

using LinearAlgebra: norm
using SparseArrays: nnz
using DifferentiationInterface: DifferentiationInterface

# PowerfulPF internal types and functions
using PowerfulPF: ACPowerFlowProblem, PowerFlowSolution, MultiIslandSolution
using PowerfulPF: NewtonMethod, NonlinearSolveBackend, AlgorithmConfig
using PowerfulPF: ResidualFunctor, _solve_multi_island, _initial_guess

# PowerfulCore types and functions
using PowerfulCore: DAE, ForwardDiffCache, AbstractADCache, SystemModel
using PowerfulCore: create_adcache, set_algebraic_guess!, build_status_dict, set_hold!
using PowerfulCore: jacobian_sparsity_pattern, create_di_backend
using PowerfulCore: is_gpu_array_type
using PowerfulCore: @vprint

# Import functions to extend
import PowerfulPF: setup_powerflow, solve_powerflow!, _solve

# =============================================================================
# NonlinearSolve Workspace - Explicit setup/iteration separation
# =============================================================================

"""
    NonlinearSolveWorkspace{PROB, JT, BT, PT, CT}

Pre-allocated workspace for NonlinearSolve.jl-based power flow solver.

Separates one-time setup (sparsity detection, coloring, problem construction) from
allocation-free solve calls. This enables precise benchmarking and efficient repeated solves.

# Type Parameters
- `PROB`: ACPowerFlowProblem type
- `JT`: Jacobian type (SparseMatrixCSC or Matrix)
- `BT`: Backend type for differentiation (AutoSparse or Nothing)
- `PT`: Prep type for Jacobian (PreparedJacobian or Nothing)
- `CT`: residual evaluator (AbstractADCache) type

# Fields
- `pf_prob::PROB`: Power flow problem (contains sys, address_manager, formulation)
- `adcache::CT`: Residual evaluator (stateless)
- `iter_dae::DAE{Vector{Float64}}`: Workspace owns iteration state (for initial guess access)
- `nl_prob::NonlinearProblem`: Pre-constructed NonlinearSolve problem
- `J::JT`: Pre-allocated sparse Jacobian matrix
- `jac_prep::PT`: Prepared Jacobian (coloring, etc.)
- `jac_backend::BT`: Differentiation backend
- `g_temp::Vector{Float64}`: Temporary buffer for Jacobian computation
- `solver`: Default solver algorithm

# Example
```julia
using NonlinearSolve  # Triggers extension loading

# One-time setup (expensive: sparsity, coloring, problem construction)
workspace = setup_powerflow(prob; verbose=true)

# Solve multiple times with different initial conditions (fast!)
y0_1 = copy(workspace.iter_dae.y)
sol1 = solve_powerflow!(workspace, y0_1; abstol=1e-8)

y0_2 = y0_1 .* 1.05  # Different initial guess
sol2 = solve_powerflow!(workspace, y0_2; abstol=1e-10)

# Benchmark iterations only (no setup overhead)
@btime solve_powerflow!(\$workspace, \$y0_1; abstol=1e-8)
```

# See Also
- [`setup_powerflow`](@ref): One-time setup function
- [`solve_powerflow!`](@ref): Solve using pre-allocated workspace
- [`NewtonWorkspace`](@ref): Hand-baked Newton solver workspace
"""
struct NonlinearSolveWorkspace{PROB, JT, BT, PT, CT<:AbstractADCache}
    pf_prob::PROB  # ACPowerFlowProblem (contains sys, address_manager, formulation)
    adcache::CT  # Residual evaluator (stateless)
    iter_dae::DAE{Vector{Float64}}  # Workspace owns iteration state (for initial guess access)
    nl_prob::NonlinearProblem  # NonlinearSolve problem
    J::JT
    jac_prep::PT
    jac_backend::BT
    g_temp::Vector{Float64}
    solver::Any
end

"""
    setup_powerflow(prob::ACPowerFlowProblem; kwargs...) -> NonlinearSolveWorkspace

Perform all one-time setup for NonlinearSolve.jl-based power flow solver.

This function separates setup cost from iteration cost, enabling:
- Precise benchmarking of iterations only
- Reusing workspace for multiple solves with different initial conditions
- Clear visibility into setup vs iteration allocation

# Arguments
- `prob::ACPowerFlowProblem`: Power flow problem with system, formulation, and address manager

# Keyword Arguments
- `solver=nothing`: NonlinearSolve algorithm (default: TrustRegion)
- `verbose::Bool=false`: Print detailed setup progress

# Returns
- `NonlinearSolveWorkspace`: Pre-allocated workspace ready for repeated solves

# Example
```julia
using NonlinearSolve  # Required for this backend

prob = ACPowerFlowProblem(sys)
workspace = setup_powerflow(prob; verbose=true)

y0 = copy(workspace.iter_dae.y)
sol1 = solve_powerflow!(workspace, y0; abstol=1e-8)
sol2 = solve_powerflow!(workspace, y0; abstol=1e-10)

@btime solve_powerflow!(\$workspace, \$y0; abstol=1e-8)
```

# See Also
- [`solve_powerflow!`](@ref): Solve using pre-allocated workspace
- [`setup_powerflow_newton`](@ref): Hand-baked Newton solver setup
"""
function setup_powerflow(
    prob::ACPowerFlowProblem;
    solver = nothing,
    verbose::Bool = false
)
    sys = prob.sys
    am = prob.address_manager

    if is_gpu_array_type(sys.array_type)
        error("""
        GPU power flow not yet supported. Use CPU version:

            sys_cpu = SystemModel(models)  # CPU version
            prob = ACPowerFlowProblem(sys_cpu, PowerInjection())
            solution = solve(prob, RobustNR())
        """)
    end

    @vprint "Setting up NonlinearSolve Power Flow Solver"

    verbose && print("\nStep 1/4: Creating workspace... ")
    t_cache = @elapsed begin
        adcache = create_adcache(am, sys, prob.formulation, prob.backend)
        iter_dae = DAE(am; array_type=sys.array_type)
        set_algebraic_guess!(iter_dae.y, am, sys, prob.formulation)
        # Create status dict for residual evaluation
        base_status = build_status_dict(sys)
    end
    @vprint "done ($(round(t_cache, digits=3))s)"
    @vprint "  Problem size: $(iter_dae.m) variables"

    # Type barrier: wrap adcache in Ref{Any} to prevent closure type explosion
    # This breaks the compile-time type chain for NonlinearFunction while preserving runtime performance.
    # Cost: one dynamic dispatch per call (~nanoseconds) vs potential compilation time savings.
    adcache_ref = Ref{Any}(adcache)

    function residual_function!(out, u, _p)
        adcache_unwrapped = adcache_ref[]::ForwardDiffCache
        adcache_unwrapped(out, u, base_status, 0.0)   # the operating point, at t = 0
        nothing
    end

    verbose && print("\nStep 2/4: Computing sparsity pattern... ")
    t_sparsity = @elapsed begin
        sparsity = jacobian_sparsity_pattern(am, sys, prob.formulation)
        J = Float64.(sparsity)
    end
    @vprint "done ($(round(t_sparsity, digits=3))s)"
    @vprint "  Nonzeros: $(nnz(J)) / $(iter_dae.m^2)"
    @vprint "  Density: $(round(100.0 * nnz(J) / (iter_dae.m^2), digits=2))%"

    verbose && print("\nStep 3/4: Preparing Jacobian (coloring)... ")
    t_prep = @elapsed begin
        jac_backend = create_di_backend(prob.backend, sparsity)

        g_temp = similar(iter_dae.g)
        # Use ResidualFunctor for consistent type with Enzyme
        residual_functor = ResidualFunctor(adcache, base_status)
        jac_prep = DifferentiationInterface.prepare_jacobian(
            residual_functor, zero(g_temp), jac_backend, zero(iter_dae.y); strict=Val(false)
        )
    end
    @vprint "done ($(round(t_prep, digits=3))s)"

    # Type barrier: wrap jac_prep in Ref{Any} to prevent closure type explosion
    # jac_prep is SMCPushforwardSparseJacobianPrep{...} which embeds the full cache type
    jac_prep_ref = Ref{Any}(jac_prep)
    jac_backend_ref = Ref{Any}(jac_backend)
    g_temp_ref = Ref{Any}(g_temp)
    residual_functor_ref = Ref{Any}(residual_functor)

    function sparse_jacobian!(J_out, u, _p)
        jac_prep_unwrapped = jac_prep_ref[]
        jac_backend_unwrapped = jac_backend_ref[]
        g_temp_unwrapped = g_temp_ref[]
        functor_unwrapped = residual_functor_ref[]
        DifferentiationInterface.jacobian!(
            functor_unwrapped, g_temp_unwrapped, J_out, jac_prep_unwrapped, jac_backend_unwrapped, u
        )
        nothing
    end

    verbose && print("\nStep 4/4: Constructing NonlinearProblem... ")
    t_prob = @elapsed begin
        nonlin_func = NonlinearFunction(
            residual_function!;
            jac = sparse_jacobian!,
            jac_prototype = J
        )
        nlprob = NonlinearProblem(nonlin_func, iter_dae.y, nothing)
    end
    @vprint "done ($(round(t_prob, digits=3))s)"

    # Select default solver
    solver_alg = isnothing(solver) ? NonlinearSolve.TrustRegion() : solver

    @vprint "Setup complete (Sparse Jacobian + TrustRegion)"

    return NonlinearSolveWorkspace(
        prob,  # Store ACPowerFlowProblem
        adcache,
        iter_dae,
        nlprob,  # NonlinearProblem
        J,
        jac_prep,
        jac_backend,
        g_temp,
        solver_alg
    )
end

"""
    solve_powerflow!(workspace::NonlinearSolveWorkspace, initial_guess::Vector; kwargs...)
        -> PowerFlowSolution

Solve power flow using pre-allocated NonlinearSolve workspace.

This function performs ONLY the solve, assuming all setup has been done via
`setup_powerflow`. This enables precise benchmarking and efficient repeated solves.

# Arguments
- `workspace::NonlinearSolveWorkspace`: Pre-allocated workspace from `setup_powerflow`
- `initial_guess::Vector`: Initial guess for algebraic variables

# Keyword Arguments
- `abstol::Real=1e-8`: Absolute tolerance for residual norm
- `reltol::Real=1e-8`: Relative tolerance
- `maxiters::Int=50`: Maximum iterations
- `solver=nothing`: Override default solver (from workspace)

# Returns
- `PowerFlowSolution`: Solution object with convergence status and results
"""
function solve_powerflow!(
    workspace::NonlinearSolveWorkspace,
    initial_guess::Vector;
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
    maxiters::Int = 50,
    solver = nothing
)
    # Use provided solver or workspace default
    solver_alg = isnothing(solver) ? workspace.solver : solver

    # Update problem with new initial guess; out-of-service instances hold their own
    # variables at it
    set_hold!(workspace.adcache, initial_guess)
    nl_prob = remake(workspace.nl_prob; u0 = initial_guess)

    # Solve
    sol = NonlinearSolve.solve(
        nl_prob,
        solver_alg;
        abstol = abstol,
        reltol = reltol,
        maxiters = maxiters
    )

    # Compute final residuals
    final_residuals = similar(sol.u)
    final_status = build_status_dict(workspace.pf_prob.sys)
    workspace.adcache(final_residuals, sol.u, final_status, 0.0)   # the operating point, at t = 0
    final_norm = norm(final_residuals, Inf)

    converged = sol.retcode == ReturnCode.Success
    iterations = sol.stats.nsteps

    # Create DAE with final solution (NonlinearSolve.jl returns sol.u, not workspace.iter_dae)
    # Copy solution into a new DAE structure
    final_dae = deepcopy(workspace.iter_dae)
    copyto!(final_dae.y, sol.u)
    copyto!(final_dae.g, final_residuals)

    return PowerFlowSolution(
        converged,
        final_norm,
        iterations,
        sol.retcode,
        sol.stats,
        workspace.pf_prob,
        final_dae,  # Store final solution DAE
        Ref{Union{Nothing, Vector{eltype(sol.u)}}}(nothing)  # Lazy cache for generator residuals
    )
end

# =============================================================================
# Solve dispatch for NonlinearSolveBackend
# =============================================================================

function _solve(
    prob::ACPowerFlowProblem,
    alg::AlgorithmConfig{<:NewtonMethod, NonlinearSolveBackend};
    solver=nothing,
    verbose::Bool=false,
    kwargs...
)
    # Check for multi-island first
    if prob.sys isa Vector{SystemModel}
        return _solve_multi_island(prob, alg; solver=solver, verbose=verbose, kwargs...)
    end

    # Single-island: setup workspace (includes validation)
    workspace = setup_powerflow(prob; solver=solver, verbose=verbose)
    y0 = _initial_guess(prob, workspace)

    # Solve using NonlinearSolve.jl (remaining kwargs passed through)
    return solve_powerflow!(
        workspace,
        y0;
        abstol = alg.method.abstol,
        reltol = alg.method.reltol,
        maxiters = alg.method.maxiters,
        kwargs...
    )
end

end # module
