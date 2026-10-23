"""
Power flow solvers with automatic differentiation.

# Two-Level API

## Level 1: Modern Composable API (Recommended)
Use the DiffEq-style API from `power_flow_api.jl` for maximum flexibility:
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, RobustNR(); abstol=1e-8)  # NonlinearSolve.jl backend (robust)
sol = solve(prob, Newton(); abstol=1e-8) # Native Newton backend (max performance)
```

## Level 2: Direct Performance API (Advanced Users)
Direct access to high-performance Newton solver without Problem objects.
Use for benchmarking or when you need maximum control.

### Workspace API (Benchmarking, Repeated Solves)
```julia
# Separate setup from iterations for precise timing
prob = ACPowerFlowProblem(sys)
workspace = setup_powerflow_newton(prob)

# Copy initial guess (required - workspace vectors are mutated)
y0 = copy(workspace.iter_dae.y)

# Fast repeated solves (zero setup overhead)
sol = solve_powerflow_newton!(workspace, y0; abstol=1e-8)

# Benchmark iterations only (excludes one-time setup cost)
@btime solve_powerflow_newton!(\$workspace, \$y0; abstol=1e-8)
```

**Performance Note**: For most users, use `solve(prob, Newton())` from the high-level API.
The workspace API is for benchmarking or when you need to reuse setup across multiple solves.
"""

using LinearAlgebra
using KLU
using SparseArrays
using DifferentiationInterface
using DifferentiationInterface: AutoSparse
using ADTypes: AutoForwardDiff
using SparseMatrixColorings: GreedyColoringAlgorithm
import PowerfulCore: compute_observations
using PowerfulCore: @vprint, build_status_dict, set_hold!

# Note: PowerFlowSolution struct is defined in solution_types.jl (included early for circular dependency resolution)

# ============================================================================
# ResidualFunctor - Consistent type for Enzyme compatibility
# ============================================================================

"""
    ResidualFunctor

Functor wrapping ADCache and status for consistent types with Enzyme.

Anonymous closures like `(g, y) -> adcache(g, y, status, 0.0)` get different
compiler-generated types each time they're created. This breaks Enzyme's
`BatchDuplicated` which requires homogeneous types.

This struct provides a stable type that can be reused between
`prepare_jacobian` and `jacobian!` calls.

It evaluates at `t = 0`: a power flow's state is the operating point a simulation starts
from, at the start of the time axis, and its residuals do not depend on time.
"""
struct ResidualFunctor{C, S}
    adcache::C
    status::S
end

function (f::ResidualFunctor)(g::AbstractVector, y::AbstractVector)
    f.adcache(g, y, f.status, 0.0)   # the operating point, at t = 0 (see the docstring)
    return nothing
end

"""
    get_generator_residuals(sol::PowerFlowSolution{T}) where T

Lazy accessor for generator residual contributions after power flow converges.

Returns residual vector containing actual power injections from generators (PV and Slack)
that were computed by running residual evaluation in "observation mode" (force_accumulate=true).

# Physical Interpretation
For converged power flow, this function extracts:
- **PV generators**: Actual Q injection at bus to maintain voltage V = Vs
- **Slack generators**: Actual P and Q injections to maintain voltage and angle

The residuals are evaluated with generators accumulating (not replacing) bus power
balance equations, allowing extraction of their actual contributions.

# Returns
- `Vector{T}`: Residual vector of same length as solution, with generator contributions
  at bus injection positions. Non-generator positions contain load/branch contributions.

# Caching Behavior
- First call: Computes generator residuals and caches result
- Subsequent calls: Returns cached result (O(1) operation)

# Example
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())

# Extract generator contributions (computed and cached on first call)
g_gen = get_generator_residuals(sol)

# Access PV reactive power at specific bus
q_addrs = get_addresses(sol.prob.address_manager, :Bus, :q_inj)
q_gen_bus5 = g_gen[q_addrs[5]]
```

# See Also
- `evaluate_generator_contributions!`: Direct evaluation without caching
"""
function get_generator_residuals(sol::PowerFlowSolution{T}) where T
    # Check if already cached
    if sol._generator_residuals[] === nothing
        # Compute generator contributions
        prob = sol.prob
        sys_ws = SystemBuffers(T, prob.address_manager, prob.sys, prob.formulation)
        status = build_status_dict(prob.sys)
        g_gen = evaluate_generator_contributions!(sys_ws, sol.dae.y, status)

        # Cache the result
        sol._generator_residuals[] = g_gen
    end

    return sol._generator_residuals[]
end

# =============================================================================
# NonlinearSolve Stubs (real implementations in PowerfulPFNonlinearSolveExt)
# =============================================================================

"""
    setup_powerflow(prob::ACPowerFlowProblem; kwargs...)

Setup NonlinearSolve.jl-based power flow workspace. Requires `using NonlinearSolve`.

See also: [`setup_powerflow_newton`](@ref) for the native Newton solver (no extra dependency).
"""
function setup_powerflow end

"""
    solve_powerflow!(workspace, initial_guess::Vector; kwargs...)

Solve power flow using NonlinearSolve workspace. Requires `using NonlinearSolve`.

See also: [`solve_powerflow_newton!`](@ref) for the native Newton solver.
"""
function solve_powerflow! end

# =============================================================================
# Helper functions for solve_powerflow_newton
# =============================================================================

"""
    _validate_newton_setup(sys, use_gpu, verbose)

Validate preconditions for Newton solver.
"""
function _validate_newton_setup(sys::SystemModel, use_gpu::Bool, verbose::Bool)
    # Note: Setup validation removed - SystemModel setup happens at problem construction
    # via ACPowerFlowProblem which calls setup_with_formulation!()

    if use_gpu && isempty(methods(gpu_linear_solver))
        error("GPU linear solver requested, but `gpu_linear_solver` has no method: " *
              "no GPU package providing it is loaded.")
    end

    if use_gpu && !is_gpu_array_type(sys.array_type)
        @vprint "Hybrid mode: CPU Jacobian computation → GPU linear solve → CPU solution"
    end
end

"""
    _setup_newton_workspace(am, sys, formulation, backend, status, verbose)

Create stateless autodiff cache and setup DifferentiationInterface for Jacobian computation.
The cache provides efficient residual evaluation with reusable computational buffers.

Returns: (adcache, jac_prep, jac_backend, J, residual_functor)
"""
function _setup_newton_workspace(am::AddressManager, sys::SystemModel, formulation::FormulationTrait, backend::ADBackend,
                                 status::Dict{Symbol, <:AbstractVector{Bool}}, verbose::Bool)
    verbose && print("Step 1/5: Creating workspace... ")
    t_setup = @elapsed begin
        adcache = create_adcache(am, sys, formulation, backend)
    end

    # Determine problem size from address manager
    n_vars = Int(max_address_in_space(am, AlgebVar()))

    @vprint "done ($(round(t_setup, digits=3))s)"
    @vprint "  Problem size: $n_vars variables"

    verbose && print("Step 2/5: Setting up DifferentiationInterface (sparsity + coloring)... ")
    t_di_setup = @elapsed begin
        # The functor shares the workspace's status dict, so outages applied to it
        # reach both residual and Jacobian. Sparsity comes from the model structure,
        # a superset of the non-zeros under any status (conservative for contingencies).
        # Use ResidualFunctor for consistent type with Enzyme
        residual_functor = ResidualFunctor(adcache, status)

        # Prototype vectors match the system array type (CPU or GPU)
        ArrayType = sys.array_type{Float64}
        y_proto = ArrayType(undef, n_vars)

        # Use strict=Val(false) because we may use different status dicts at jacobian! time
        # (for parallel ACCC). The underlying adcache is the same, so sparsity is a superset.
        (; J, backend, prep) = prepare_sparse_jacobian(
            residual_functor, zero(y_proto), zero(y_proto), adcache; strict=Val(false)
        )
        jac_backend, jac_prep = backend, prep
    end
    @vprint "done ($(round(t_di_setup, digits=3))s)"
    @vprint "  Nonzeros: $(nnz(J)) / $(n_vars^2)"
    @vprint "  Density: $(round(100.0 * nnz(J) / (n_vars^2), digits=2))%"

    @vprint "Step 3/5: Matrix coloring completed by DifferentiationInterface"
    @vprint "  ✓ Automatic chunking handled internally"

    return (adcache, jac_prep, jac_backend, J, residual_functor)
end

# =============================================================================
# Newton Workspace - Explicit setup/iteration separation
# =============================================================================

"""
    NewtonWorkspace{SYS, JT, BT, PT}

Pre-allocated workspace for Newton-Raphson power flow solver.

Separates one-time setup (sparsity detection, coloring, allocation) from
allocation-free Newton iterations. Uses plain vectors for iteration state,
eliminating aliasing issues and ensuring zero allocations in the hot loop.

# Type Parameters
- `SYS`: SystemModel type (matches CachedAutodiffWorkspace)
- `JT`: Jacobian type (SparseMatrixCSC)
- `BT`: Backend type for differentiation (AutoSparse)
- `PT`: Prep type for Jacobian (PreparedJacobian)

# Fields
- `cache::CachedAutodiffWorkspace{SYS}`: Residual evaluator
- `iter_dae::DAE{Vector{Float64}}`: Current iteration state (y, g)
- `ls_dae::DAE{Vector{Float64}}`: Line search trial state (y_ls, g_ls)
- `dy::Vector{Float64}`: Newton update direction (RHS buffer)
- `J::SparseMatrixCSC`: Sparse Jacobian matrix
- `jac_prep`: Prepared Jacobian (coloring, etc.)
- `jac_backend`: Differentiation backend
- `linear_solver`: Factorization cache (`KLULinearSolver`, or the GPU solver of `gpu_linear_solver`)

# Example
```julia
# One-time setup (expensive: sparsity, coloring, allocation)
prob = ACPowerFlowProblem(sys)
workspace = setup_powerflow_newton(prob; verbose=true)

# Solve multiple times with different initial conditions (allocation-free!)
y0_1 = copy(workspace.iter_dae.y)
sol1 = solve_powerflow_newton!(workspace, y0_1; abstol=1e-8)

y0_2 = y0_1 .* 1.1  # Different initial guess
sol2 = solve_powerflow_newton!(workspace, y0_2; abstol=1e-10)

# Benchmark iterations only (no setup overhead)
@btime solve_powerflow_newton!(\$workspace, \$y0_1; abstol=1e-8, verbose=false)
```

# Design Notes
Power flow solves algebraic equations only. When state variables (x, f) are added
for dynamic simulation, they will be added as explicit fields to this workspace
or a separate DynamicSimWorkspace will be created.
"""
struct NewtonWorkspace{PROB, AT<:AbstractVector, JT, PT, BT, CT<:AbstractADCache, FT, ST<:AbstractVector{Bool}, LS<:AbstractNewtonLinearSolver}
    # Power flow problem (contains sys, address_manager, formulation)
    prob::PROB

    # Residual evaluator (stateless, provides callable interface)
    adcache::CT

    # Iteration state (owned by workspace, single source of truth)
    iter_dae::DAE{AT}  # Current iterate: iter_dae.y and iter_dae.g

    # Line search trial state (bundled for scalability)
    ls_dae::DAE{AT}  # Line search trial point (ls_dae.y, ls_dae.g)
    dy::AT             # Newton direction (RHS buffer)

    # Jacobian components (pre-allocated, pre-colored)
    J::JT            # Sparse Jacobian matrix
    jac_prep::PT     # Prepared Jacobian (coloring info)
    jac_backend::BT  # Differentiation backend

    # Residual functor for Jacobian computation (consistent type for Enzyme)
    residual_functor::FT

    # Model status (active/inactive) for contingency analysis (CPU Vector or GPU CuArray)
    status::Dict{Symbol, ST}

    # Linear solver with its factorization cache (KLU, or GPU from gpu_linear_solver)
    linear_solver::LS
end

"""
    Base.show(io::IO, ws::NewtonWorkspace)

Display problem size and configuration.
"""
function Base.show(io::IO, ws::NewtonWorkspace)
    if get(io, :compact, false)
        print(io, "NewtonWorkspace(...)")
    else
        n_vars = length(ws.iter_dae.y)
        backend = nameof(typeof(ws.linear_solver))
        formulation_name = nameof(typeof(ws.prob.formulation))

        println(io, "NewtonWorkspace:")
        println(io, "  Variables: $n_vars")
        println(io, "  Jacobian: $(n_vars) x $(n_vars) (sparse)")
        println(io, "  Backend: $backend")
        println(io, "  Formulation: $formulation_name")
    end
end

"""
    setup_powerflow_newton(prob::ACPowerFlowProblem; kwargs...) -> NewtonWorkspace

Perform all one-time setup for Newton-Raphson power flow solver.

This function separates setup cost from iteration cost, enabling:
- Precise benchmarking of iterations only
- Reusing workspace for multiple solves with different initial conditions
- Clear visibility into setup vs iteration allocation

# Arguments
- `prob::ACPowerFlowProblem`: Power flow problem with system, formulation, and address manager

# Keyword Arguments
- `use_gpu::Bool=false`: Use GPU linear solver
- `ordering::Symbol=:AMD`: KLU fill-reducing ordering (`:AMD`, `:COLAMD`, `:natural`)
- `status`: Model status dict used by both the residual and the Jacobian, owned by the
  workspace (default: `build_status_dict(prob.sys)`, a copy of the system's status);
  change it to apply outages without modifying the system
- `verbose::Bool=false`: Print detailed setup progress

# Returns
- `NewtonWorkspace`: Pre-allocated workspace ready for allocation-free iterations

# Example
```julia
# Create power flow problem with formulation
prob = ACPowerFlowProblem(sys, PowerInjection())

# One-time setup with timing details
workspace = setup_powerflow_newton(prob; verbose=true)

# Multiple solves (fast, no setup overhead!)
sol1 = solve_powerflow_newton!(workspace, initial_guess_1; abstol=1e-8)
sol2 = solve_powerflow_newton!(workspace, initial_guess_2; abstol=1e-10)

# Benchmark iterations precisely (excludes setup)
using BenchmarkTools
y0 = copy(workspace.iter_dae.y)
@btime solve_powerflow_newton!(\$workspace, \$y0; abstol=1e-8, verbose=false)
```

# See Also
- [`solve_powerflow_newton!`](@ref): Allocation-free iteration function
- [`solve`](@ref): High-level API with algorithm presets
- [`NewtonWorkspace`](@ref): Workspace structure documentation
"""
function setup_powerflow_newton(
    prob::ACPowerFlowProblem;
    use_gpu::Bool = false,
    ordering::Symbol = :AMD,
    status::Dict{Symbol, <:AbstractVector{Bool}} = build_status_dict(prob.sys),
    verbose::Bool = false
)
    sys = prob.sys
    am = prob.address_manager

    _validate_newton_setup(sys, use_gpu, verbose)

    @vprint "Setting up Newton-Raphson Power Flow Solver"

    # Model active/inactive status is shared by the residual and the Jacobian functor
    adcache, jac_prep, jac_backend, J, residual_functor = _setup_newton_workspace(am, sys, prob.formulation, prob.backend, status, verbose)

    verbose && print("Step 4/5: Allocating iteration and line search DAE buffers... ")
    t_buf = @elapsed begin
        # Create iteration DAE (workspace owns iteration state)
        iter_dae = DAE(am; array_type=sys.array_type)
        set_algebraic_guess!(iter_dae.y, am, sys, prob.formulation)

        # Create line search DAE (trial state for line search)
        ls_dae = DAE(am; array_type=sys.array_type)
        set_algebraic_guess!(ls_dae.y, am, sys, prob.formulation)

        n = iter_dae.m
        # Create dy with same array type as DAE (Vector for CPU, CuArray for GPU)
        dy = similar(iter_dae.y)
    end
    @vprint "done ($(round(t_buf, digits=3))s)"
    @vprint "  Created: iter_dae, ls_dae ($n variables each)"
    @vprint "  Allocated: dy"
    @vprint "  Iteration state owned by workspace, line search bundled in ls_dae"

    linear_solver = use_gpu ? gpu_linear_solver() : KLULinearSolver(ordering)
    @vprint "Setup complete (Sparse Jacobian + $(nameof(typeof(linear_solver))))"

    return NewtonWorkspace(
        prob,  # Store ACPowerFlowProblem
        adcache,
        iter_dae,
        ls_dae,
        dy,
        J,
        jac_prep,
        jac_backend,
        residual_functor,
        status,
        linear_solver
    )
end

# =============================================================================
# Allocation-free Newton iterations
# =============================================================================

"""
    solve_powerflow_newton!(workspace::NewtonWorkspace, initial_guess::Vector; kwargs...)
        -> PowerFlowSolution

Solve power flow using pre-allocated workspace (allocation-free iterations).

This function performs ONLY Newton-Raphson iterations, assuming all setup
has been done via `setup_powerflow_newton`. All buffers are pre-allocated,
making this suitable for:
- Precise iteration benchmarks (excludes setup overhead)
- Multiple solves with different initial conditions
- Repeated solves in optimization loops

# Arguments
- `workspace::NewtonWorkspace`: Pre-allocated workspace from `setup_powerflow_newton`
- `initial_guess::Vector`: Initial guess for all variables. **IMPORTANT**: Must be a copy of
  `workspace.iter_dae.y` or your custom guess. Do NOT pass the workspace vector directly as it
  will be mutated during iterations

# Keyword Arguments
- `abstol::Real=1e-8`: Absolute tolerance for residual norm
- `maxiters::Int=50`: Maximum Newton iterations
- `damping::Real=1.0`: Step size damping factor (0.5 for cautious, 1.0 for full Newton)
- `jacobian_every::Int=1`: Refresh the Jacobian every N iterations (N>1: dishonest Newton)
- `verbose::Bool=false`: Print iteration details with timing

# Returns
- `PowerFlowSolution`: Solution object with convergence status and results

# Example
```julia
# One-time setup
prob = ACPowerFlowProblem(sys)
workspace = setup_powerflow_newton(prob; verbose=true)

# IMPORTANT: Copy the initial guess (workspace vectors are mutated during solve)
y0 = copy(workspace.iter_dae.y)

# Solve with different tolerances (fast, no setup overhead!)
sol_loose = solve_powerflow_newton!(workspace, y0; abstol=1e-6)
sol_tight = solve_powerflow_newton!(workspace, y0; abstol=1e-10)

# Try different initial conditions
y0_perturbed = y0 .* 1.05  # 5% perturbation
sol_perturbed = solve_powerflow_newton!(workspace, y0_perturbed; abstol=1e-8)

# Benchmark iterations only (excludes setup)
using BenchmarkTools
@btime solve_powerflow_newton!(\$workspace, \$y0; abstol=1e-8, verbose=false)
```

# Notes
- Uses DAE accessor functions to be future-proof when state variables are added
- All iteration buffers pre-allocated in workspace (allocation-free after warmup)
- Linear solver factorization cached and reused across iterations

# See Also
- [`setup_powerflow_newton`](@ref): One-time setup function
- [`solve`](@ref): High-level API (recommended for most users)
- [`NewtonWorkspace`](@ref): Workspace structure documentation
"""
function solve_powerflow_newton!(
    workspace::NewtonWorkspace,
    initial_guess::AbstractVector;
    abstol::Real = 1e-8,
    maxiters::Int = 50,
    damping::Real = 1,
    jacobian_every::Int = 1,
    verbose::Bool = false
)
    (; prob, adcache, iter_dae, ls_dae, dy, J, jac_prep, jac_backend,
       residual_functor, status, linear_solver) = workspace

    # Initialize iteration state from initial guess
    # workspace.iter_dae is the single source of truth for iteration state
    copyto!(iter_dae.y, initial_guess)
    # Out-of-service instances hold their own variables at the initial guess
    set_hold!(adcache, initial_guess)

    # Fresh factorization (symbolic + numeric) for each solve
    reset_factorization!(linear_solver)
    @vprint "Solving with $(nameof(typeof(linear_solver)))..."

    residual! = (g, y) -> adcache(g, y, status, 0.0)   # the operating point, at t = 0 (`ResidualFunctor`)
    # ls_dae.g doubles as the Jacobian output scratch
    jacobian! = (Jm, y) -> DifferentiationInterface.jacobian!(
        residual_functor, ls_dae.g, Jm, jac_prep, jac_backend, y)
    bufs = NewtonBuffers(iter_dae.g, ls_dae.g, ls_dae.y, dy, J)
    policy = NewtonPolicy(; damping = Float64(damping), jacobian_every)

    local result
    t_solve = @elapsed result = newton_solve!(residual!, jacobian!, iter_dae.y, bufs, linear_solver, policy;
                                               abstol, maxiters, verbose)
    res_norm = result.residual_norm

    if result.converged
        @vprint "Converged in $(result.iterations) iterations ($(round(t_solve*1000, digits=1)) ms)"
        return PowerFlowSolution(
            true,
            res_norm,
            result.iterations,
            :Success,
            nothing,
            prob,
            deepcopy(iter_dae),  # Deep copy to preserve solution state
            Ref{Union{Nothing, Vector{eltype(iter_dae.y)}}}(nothing)  # Lazy cache
        )
    end

    retcode = result.status === :diverged ? :Diverged : :MaxIters
    @vprint "Did NOT converge in $maxiters iterations ($(round(t_solve*1000, digits=1)) ms), residual=$(res_norm)"

    return PowerFlowSolution(
        false,
        res_norm,
        result.iterations,
        retcode,
        nothing,
        prob,
        deepcopy(iter_dae),  # Deep copy to preserve solution state
        Ref{Union{Nothing, Vector{eltype(iter_dae.y)}}}(nothing)  # Lazy cache
    )
end

# =============================================================================
# Exports
# =============================================================================
       setup_powerflow_newton, solve_powerflow_newton!,
       NewtonWorkspace, get_generator_residuals

# =============================================================================
# Solution Convenience Methods
# =============================================================================

"""
    compute_observations(sol::PowerFlowSolution, ::Type{ModelT}, var::Symbol) -> Vector{Float64}

Generic fallback implementation for simple observations.

Automatically handles:
1. **State variables**: Extract from solution vector (e.g., :theta, :v for Bus)
2. **Parameters**: Return from model struct (e.g., :pg for PV)

Models only need custom implementations for computed observations (e.g., power flows).
"""
function compute_observations(sol::PowerFlowSolution, ::Type{ModelT}, var::Symbol) where {ModelT <: AbstractModel}
    model_name = nameof(ModelT)

    # Check if model exists in system
    if !haskey(sol.prob.sys.models, model_name)
        error("Model $model_name not found in system")
    end

    model = sol.prob.sys.models[model_name]
    meta = model_metadata(ModelT)
    am = sol.prob.address_manager

    # Try to find observation in metadata
    obs_idx = findfirst(o -> o.name == var, meta.observations)
    if obs_idx === nothing
        error("Unknown observation :$var for model $model_name")
    end

    # Check if it's a state variable we can extract
    for model_var in meta.vars
        if model_var.name == var
            # Extract from solution vector
            addrs = get_addresses(am, model_name, var)
            return sol.dae.y[addrs]
        end
    end

    # Check if it's a model parameter we can return directly
    if hasfield(ModelT, var)
        return getproperty(model, var)
    end

    # If we reach here, model needs custom implementation
    error("Model $model_name must implement compute_observations for derived observation :$var")
end

#==============================================================================
# Ergonomic Accessors for PowerFlowSolution
==============================================================================#

"""
    get_var(sol::PowerFlowSolution, model::Symbol, var::Symbol, ext_key)

Ergonomic wrapper for querying solution variables.

# Example
```julia
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, RobustNR())
v = get_var(sol, :Bus, :v, (4001,))  # Clean and simple!
```
"""
function get_var(sol::PowerFlowSolution, model::Symbol, var::Symbol, ext_key)
    return get_var(sol.prob, sol.dae.y, model, var, ext_key)
end

"""
    get_var_by_index(sol::PowerFlowSolution, model::Symbol, var::Symbol, idx::Int)

Get variable from solution by internal index.
"""
function get_var_by_index(sol::PowerFlowSolution, model::Symbol, var::Symbol, idx::Int)
    return get_var_by_index(sol.prob, sol.dae.y, model, var, idx)
end

"""
    get_all_vars(sol::PowerFlowSolution, model::Symbol, var::Symbol)

Get all variable values from solution.
"""
function get_all_vars(sol::PowerFlowSolution, model::Symbol, var::Symbol)
    return get_all_vars(sol.prob, sol.dae.y, model, var)
end

"""
    get_vars_dict(sol::PowerFlowSolution, model::Symbol, var::Symbol)

Get all variable values from solution as dictionary.
"""
function get_vars_dict(sol::PowerFlowSolution, model::Symbol, var::Symbol)
    return get_vars_dict(sol.prob, sol.dae.y, model, var)
end

