# solve is imported from CommonSolve in the parent module (PowerfulPF.jl)
# This extends CommonSolve.solve for our problem types

"""
Power Flow API - DifferentialEquations.jl-style interface

This module provides a modern, composable API for power flow analysis that matches
Julia ecosystem conventions (DifferentialEquations.jl, NonlinearSolve.jl, etc.).

# Core Pattern

```julia
# 1. Define problem
prob = ACPowerFlowProblem(sys)

# 2. Choose algorithm
alg = Newton()

# 3. Solve
sol = solve(prob, alg)
```

# Design Philosophy

The API separates concerns into composable layers:

1. **Problem Formulation**: What to solve (AC/DC power flow)
2. **Solution Method**: How to solve (Newton, Fast Decoupled)
3. **Implementation Backend**: Which package (Native, NonlinearSolve.jl)
4. **Linear Solver**: Factorization method (KLU, CUDSS)

This design provides:
- **Extensibility**: Easy to add new algorithms
- **Consistency**: Matches Julia ecosystem (DiffEq, Optimization, NonlinearSolve)
- **Performance**: Concrete types enable compiler optimizations
"""

# ============================================================================
# Address Reordering (AMD/CustomPermutation)
# ============================================================================

"""
    apply_addressing_strategy!(am::AddressManager, sys::SystemModel, formulation::FormulationTrait)

Apply AMD or custom permutation reordering to addresses.

This function handles address reordering strategies that require sparsity analysis,
which is only available after `setup_with_formulation!` has assigned basic addresses.

Must be called after `setup_with_formulation!`.

# Supported Strategies
- `NaturalAddressing`: No-op (addresses remain in natural order)
- `AMDAddressing`: Computes sparsity pattern, AMD permutation, and remaps addresses
- `CustomPermutation`: Applies user-provided permutation

# Example
```julia
am = AddressManager(strategy=AMDAddressing())
setup_with_formulation!(am, sys, formulation)      # Basic addresses
apply_addressing_strategy!(am, sys, formulation)   # AMD reordering
```
"""
function apply_addressing_strategy!(am::AddressManager, sys::SystemModel, formulation::FormulationTrait)
    perm = addressing_permutation(am.strategy, am, sys, formulation)
    isnothing(perm) && return am  # NaturalAddressing: no remapping

    remap_addresses!(am, perm)
    am.applied_permutation = perm
    build_reverse_map!(am)  # addresses changed
    return am
end

"""
    addressing_permutation(strategy, am, sys, formulation) -> Union{Nothing, Vector{Int}}

Permutation of the algebraic addresses that `strategy` prescribes for the natural
addresses in `am`, or `nothing` to keep the natural order.
"""
addressing_permutation(::NaturalAddressing, am, sys, formulation) = nothing
addressing_permutation(::AMDAddressing, am, sys, formulation) =
    compute_amd_ordering(jacobian_sparsity_pattern(am, sys, formulation))
addressing_permutation(s::CustomPermutation, am, sys, formulation) = collect(Int, s.permutation)

export apply_addressing_strategy!

# ============================================================================
# Level 1: Problem Types
# ============================================================================

"""
    PowerFlowProblem

Abstract base type for all power flow problem formulations.

Concrete subtypes:
- `ACPowerFlowProblem`: AC power flow (nonlinear); with the `DCPowerFlow()` formulation,
  the linearized DC power flow
"""
abstract type PowerFlowProblem end

"""
    ACPowerFlowProblem <: PowerFlowProblem

AC power flow problem with full nonlinear equations.

Solves the nonlinear system g(y) = 0 where:
- y contains bus voltage magnitudes and angles
- g contains power balance equations at each bus

Automatically detects and handles multi-island systems when auto_split=true (default).

# Fields
- `sys::Union{SystemModel, Vector{SystemModel}}`: Power system (single or multi-island)
- `initial_guess::Union{Nothing, Vector{Float64}}`: Starting point for iteration

# Constructors

```julia
# Use system's current state as initial guess (default)
# Automatically detects islands
prob = ACPowerFlowProblem(sys)

# Provide custom initial guess (single-island systems only)
prob = ACPowerFlowProblem(sys; initial_guess=y0)

# Disable auto-split for manual island handling
prob = ACPowerFlowProblem(sys; auto_split=false)
```

# Example

```julia
using PowerfulCore, PowerfulPF, PowerfulCases

# Load and set up system
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# Create problem (auto-detects islands)
prob = ACPowerFlowProblem(sys)

# Solve with default algorithm
sol = solve(prob)  # Works for single or multi-island!

# Or specify algorithm
sol = solve(prob, Newton(abstol=1e-10))
```
"""
struct ACPowerFlowProblem{F<:FormulationTrait, S, B<:ADBackend} <: PowerFlowProblem
    sys::S  # SystemModel or Vector{SystemModel}
    formulation::F
    address_manager::AddressManager
    initial_guess::Union{Nothing, Vector{Float64}}
    backend::B  # Autodiff backend (ForwardDiffBackend or EnzymeBackend)

    # Constructor with formulation support and auto-split
    function ACPowerFlowProblem(
        sys::SystemModel,
        formulation::FormulationTrait=PowerInjection();
        initial_guess=nothing,
        addressing::AddressingStrategy=NaturalAddressing(),
        auto_split::Bool=true,
        backend::ADBackend=ForwardDiffBackend()
    )
        check_formulation_support(sys, formulation)

        # Auto-detect islands if requested
        if auto_split
            result = split_islands(sys)
            if result isa Vector{SystemModel}
                # Each island gets its own address layout at solve time; there is no
                # whole-system layout a single guess vector could follow.
                isnothing(initial_guess) || error(
                    "initial_guess is not supported for a system with $(length(result)) islands: " *
                    "each island is addressed separately. Build an ACPowerFlowProblem per island " *
                    "(split_islands(sys)) and pass each its own guess.")
                # Multi-island: create problem with vector of systems
                # AddressManager is unused for multi-island (created per-island during solve)
                am_placeholder = AddressManager(strategy=addressing)
                return new{typeof(formulation), Vector{SystemModel}, typeof(backend)}(
                    result, formulation, am_placeholder, initial_guess, backend
                )
            end
        end

        # Single island - existing logic
        am = AddressManager(strategy=addressing)
        setup_with_formulation!(am, sys, formulation)
        apply_addressing_strategy!(am, sys, formulation)  # AMD/CustomPermutation reordering
        return new{typeof(formulation), SystemModel, typeof(backend)}(sys, formulation, am, initial_guess, backend)
    end

    # Inner constructor for GPU transfer (bypasses auto_split and address manager setup)
    # Used by `to_gpu` methods to create a GPU problem from a CPU problem
    function ACPowerFlowProblem(
        sys::S,
        formulation::F,
        address_manager::AddressManager,
        initial_guess,
        backend::B,
        ::Val{:direct}  # Marker to select this constructor
    ) where {F<:FormulationTrait, S, B<:ADBackend}
        return new{F, S, B}(sys, formulation, address_manager, initial_guess, backend)
    end
end

"""
    Base.show(io::IO, prob::ACPowerFlowProblem)

Display formulation and problem size.
"""
function Base.show(io::IO, prob::ACPowerFlowProblem{F, S, B}) where {F, S, B}
    if get(io, :compact, false)
        print(io, "ACPowerFlowProblem{$(nameof(F))}(...)")
    else
        formulation_name = nameof(F)
        is_multi = S <: Vector{SystemModel}

        println(io, "ACPowerFlowProblem:")
        println(io, "  Formulation: $formulation_name")
        if is_multi
            println(io, "  Islands: $(length(prob.sys))")
        else
            n_alg_vars = haskey(prob.address_manager.next_idx, AlgebVar()) ?
                         prob.address_manager.next_idx[AlgebVar()] - 1 : 0
            n_diff_vars = haskey(prob.address_manager.next_idx, StateVar()) ?
                         prob.address_manager.next_idx[StateVar()] - 1 : 0

            println(io, "  Variables: $n_alg_vars algebraic, $n_diff_vars differential")
        end
    end
end

# ============================================================================
# Multi-Island Support
# ============================================================================
#
# Multi-island logic for power flow problems with disconnected components.
#
# **Architecture**:
# - Problem construction: ACPowerFlowProblem auto-detects islands via split_islands()
# - Solution type: MultiIslandSolution aggregates results from all islands (see solution_types.jl)
# - Accessors: Transparent variable access across islands (automatic lookup)
# - Solving: solve() dispatches to solve each island independently
#
# **Related code**:
# - Island detection: src/topology/topology.jl::split_islands()
# - Report merging: src/io/reports.jl::collect_model_table(::MultiIslandSolution)
# - Solution types: src/solvers/solution_types.jl (PowerFlowSolution, MultiIslandSolution)
#
# **Example workflow**:
# ```julia
# sys = SystemModel(case)              # May have disconnected islands
# prob = ACPowerFlowProblem(sys)       # Auto-detects islands
# sol = solve(prob, RobustNR())        # Returns MultiIslandSolution if needed
# v = get_var(sol, :Bus, :v, (4001,))     # Finds correct island automatically
# ```
# ============================================================================

# ============================================================================
# Level 2: Solution Method Types
# ============================================================================

"""
    PowerFlowAlgorithm

Abstract base type for all power flow solution algorithms.
"""
abstract type PowerFlowAlgorithm end

"""
    ACMethod <: PowerFlowAlgorithm

Abstract type for AC power flow solution methods.
"""
abstract type ACMethod <: PowerFlowAlgorithm end

"""
    NewtonMethod <: ACMethod

Newton-Raphson method for AC power flow (internal type).

**Note**: Users should call `Newton()` preset constructor, not this type directly.

The workhorse algorithm for power flow analysis. Uses full Jacobian matrix
with quadratic convergence near the solution.

# Fields
- `maxiters::Int`: Maximum iterations (default: 50)
- `abstol::Float64`: Absolute tolerance on residual norm (default: 1e-8)
- `reltol::Float64`: Relative tolerance (default: 1e-8); used by the NonlinearSolve backend (`RobustNR`), the native solver checks `abstol` only

# Algorithm Details
At each iteration k:
1. Evaluate residuals: g(y_k)
2. Compute Jacobian: J(y_k) = ∂g/∂y
3. Solve linear system: J·Δy = -g
4. Update: y_{k+1} = y_k + Δy
5. Check convergence: ||g(y_{k+1})||∞ < abstol
"""
Base.@kwdef struct NewtonMethod <: ACMethod
    maxiters::Int = 50
    abstol::Float64 = 1e-8
    reltol::Float64 = 1e-8
end

"""
    FastDecoupled <: ACMethod

Fast Decoupled Power Flow method with constant Jacobian matrices.

Exploits the weak coupling between P-θ and Q-V in power systems to use
constant B' and B" matrices, factored once and reused across all iterations.

# Variants
- `:XB` - Traditional fast decoupled (neglects R, uses X only)
- `:BX` - Alternate formulation (includes R effects)

# Performance
2-3x faster than Newton-Raphson for well-conditioned transmission networks.
Typical convergence: 10-20 iterations vs 3-5 for Newton (but smaller subsystems).

# When to Use
✓ Large transmission networks (>1000 buses)
✓ Well-conditioned systems (high X/R ratio)
✓ Repeated solves (factorizations reused)
✗ Stressed systems near voltage collapse
✗ Distribution networks (low X/R ratio)
"""
struct FastDecoupled <: ACMethod
    variant::Symbol  # :XB or :BX
    maxiters::Int
    abstol::Float64

    function FastDecoupled(; variant=:XB, maxiters=30, abstol=1e-6)
        if !(variant in [:XB, :BX])
            error("FastDecoupled variant must be :XB or :BX")
        end
        new(variant, maxiters, abstol)
    end
end

"""
    DishonestNewton(; update_freq=3, maxiters=50, abstol=1e-8) <: ACMethod

Dishonest (modified) Newton method: the Jacobian is refreshed and refactorized
only every `update_freq` iterations, reusing the factorization in between. Trades
convergence speed for reduced factorization cost; useful for very large systems
where Jacobian factorization dominates. Uses KLU on CPU, single-island systems.
"""
struct DishonestNewton <: ACMethod
    update_freq::Int  # Update Jacobian every N iterations
    maxiters::Int
    abstol::Float64

    function DishonestNewton(; update_freq=3, maxiters=50, abstol=1e-8)
        update_freq >= 1 || error("update_freq must be ≥ 1, got $update_freq")
        new(update_freq, maxiters, abstol)
    end
end

# ============================================================================
# Level 3: Implementation Backend Types
# ============================================================================

"""
    ImplementationBackend

Abstract type for choosing which package implements the solver.
"""
abstract type ImplementationBackend end

"""
    NativeBackend <: ImplementationBackend

Use Powerful.jl's native Newton-Raphson implementation.

Optimized for power system problems with:
- Sparse Jacobian with matrix coloring for ForwardDiff
- KLU sparse factorization with fill-reducing ordering
- Zero-allocation residual evaluation
- Optional GPU acceleration

This is the recommended backend for most use cases.
"""
struct NativeBackend <: ImplementationBackend end

"""
    NonlinearSolveBackend <: ImplementationBackend

Use NonlinearSolve.jl for robust solving with automatic algorithm selection.

Provides access to NonlinearSolve.jl's suite of algorithms:
- Trust region methods
- Levenberg-Marquardt
- Pseudo-transient continuation
- Automatic method selection

Recommended when native Newton fails or for research/comparison purposes.
"""
struct NonlinearSolveBackend <: ImplementationBackend end

# ============================================================================
# Level 5: Linear Solver Types
# ============================================================================

"""
    LinearSolver

Abstract type for linear system solvers (for J·Δy = -g in Newton method).
"""
abstract type LinearSolver end

"""
    KLUSolver <: LinearSolver

Sparse LU factorization via KLU.jl with fill-reducing ordering.

KLU is optimized for circuit simulation and works excellently for power
system Jacobians.

# Ordering Options
- `:AMD` - Approximate Minimum Degree (default, best for most cases)
- `:COLAMD` - Column AMD (alternative)
- `:natural` - No reordering (fastest factorization, but larger fill)

# Fields
- `ordering::Symbol`: Fill-reducing ordering strategy (default: :AMD)
"""
Base.@kwdef struct KLUSolver <: LinearSolver
    ordering::Symbol = :AMD

    function KLUSolver(ordering::Symbol)
        if !(ordering in [:AMD, :COLAMD, :natural])
            error("KLUSolver ordering must be :AMD, :COLAMD, or :natural")
        end
        new(ordering)
    end
end

"""
    CUDSSSolver <: LinearSolver

GPU sparse direct solver via CUDSS (CUDA Direct Sparse Solver).

Requires a CUDA-capable GPU and a GPU package that provides `gpu_linear_solver`.
Works in hybrid mode: CPU computes Jacobian, GPU performs factorization and
triangular solves.

# Usage
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton(linear_solver=CUDSSSolver()))
```

# Performance
Best for large systems (>10k buses) where GPU factorization speedup
outweighs PCIe transfer overhead. Data transfers happen automatically
each Newton iteration.
"""
struct CUDSSSolver <: LinearSolver end

# ============================================================================
# Algorithm Configuration (Composition)
# ============================================================================

"""
    AlgorithmConfig{M, B, L}

Full algorithm configuration composing method, backend, and linear solver.

Most users should use preset constructors like `Newton()` or
`RobustNR()` instead of constructing this directly.

# Type Parameters
- `M`: Solution method type
- `B`: Implementation backend type
- `L`: Linear solver type

# Example (advanced users only)

```julia
# Full custom configuration
alg = AlgorithmConfig(
    NewtonMethod(abstol=1e-10),
    NativeBackend(),
    KLUSolver(ordering=:COLAMD)
)
sol = solve(prob, alg)
```
"""
struct AlgorithmConfig{M<:PowerFlowAlgorithm,
                       B<:ImplementationBackend,
                       L<:LinearSolver}
    method::M
    backend::B
    linear_solver::L
end

# ============================================================================
# Preset Algorithm Constructors (90% Use Case)
# ============================================================================

# Default algorithm for each problem type
default_algorithm(::ACPowerFlowProblem) = Newton()

"""
    Newton(; kwargs...) -> AlgorithmConfig

Fast native Newton-Raphson with KLU sparse solver (recommended default).

# Keyword Arguments
- `maxiters::Int=50`: Maximum iterations
- `abstol::Float64=1e-8`: Convergence tolerance on ||g||∞
- `reltol::Float64=1e-8`: Relative tolerance (ignored by the native solver, which checks `abstol`)
- `ordering::Symbol=:AMD`: KLU fill-reducing ordering (:AMD, :COLAMD, :natural)
- `linear_solver::Union{Nothing, LinearSolver}=nothing`: Override linear solver (e.g., CUDSSSolver())

# Algorithm Stack
- Method: Newton-Raphson (full Jacobian, quadratic convergence)
- Globalization: None (pure Newton steps)
- Backend: Native (PowerfulPF's optimized implementation)
- Linear solver: KLU with AMD ordering (default), or user-specified

# Example

```julia
prob = ACPowerFlowProblem(sys)

# Use defaults (abstol=1e-8, maxiters=50, AMD ordering)
sol = solve(prob, Newton())

# Tighter tolerance for sensitivity analysis
sol = solve(prob, Newton(abstol=1e-12))

# More iterations for difficult cases
sol = solve(prob, Newton(maxiters=100))

# Try different ordering
sol = solve(prob, Newton(ordering=:COLAMD))

# Use GPU linear solver (requires CUDA + CUDSS)
sol = solve(prob, Newton(linear_solver=CUDSSSolver()))
```
"""
function Newton(;
    maxiters::Int=50,
    abstol::Float64=1e-8,
    reltol::Float64=1e-8,
    ordering::Symbol=:AMD,
    linear_solver::Union{Nothing, LinearSolver}=nothing
)
    solver = if isnothing(linear_solver)
        KLUSolver(; ordering)
    else
        linear_solver
    end

    return AlgorithmConfig(
        NewtonMethod(; maxiters, abstol, reltol),
        NativeBackend(),
        solver
    )
end

"""
    RobustNR(; kwargs...) -> AlgorithmConfig

Robust Newton via NonlinearSolve.jl with automatic algorithm selection.

Uses NonlinearSolve.jl's default algorithm which automatically selects
trust region, Levenberg-Marquardt, or other robust methods based on
problem characteristics.

# Keyword Arguments
- `maxiters::Int=50`: Maximum iterations
- `abstol::Float64=1e-8`: Convergence tolerance
- `reltol::Float64=1e-8`: Relative tolerance, passed to NonlinearSolve

# Algorithm Stack
- Method: Newton-Raphson (interface to NonlinearSolve)
- Globalization: Automatic (trust region, line search, etc.)
- Backend: NonlinearSolve.jl
- Linear solver: Automatic (typically KLU)

# When to Use
- Native Newton fails to converge
- Poor initial guess (e.g., high-load stressed systems)
- Research/comparison purposes
- Exploring NonlinearSolve.jl's features

# Example

```julia
prob = ACPowerFlowProblem(sys)

# Robust solver (slower but more reliable)
sol = solve(prob, RobustNR())

# Compare with native Newton
sol_native = solve(prob, Newton())
sol_robust = solve(prob, RobustNR())
```
"""
function RobustNR(;
    maxiters::Int=50,
    abstol::Float64=1e-8,
    reltol::Float64=1e-8
)
    return AlgorithmConfig(
        NewtonMethod(; maxiters, abstol, reltol),
        NonlinearSolveBackend(),
        KLUSolver()
    )
end

# Future presets (placeholders)

"""
    FastDecoupledXB(; kwargs...) -> AlgorithmConfig

Fast Decoupled Power Flow with XB formulation (neglects resistance).

XB variant uses susceptance-only approximation: B'[i,j] = -b_ij = -Im(Y_ij).
Best for transmission networks with high X/R ratio.

# Keyword Arguments
- `maxiters::Int=100`: Maximum iterations (FDPF typically needs 2-4x more than Newton)
- `abstol::Float64=1e-6`: Convergence tolerance on max(|ΔP|, |ΔQ|)
- `ordering::Symbol=:AMD`: KLU fill-reducing ordering (:AMD, :COLAMD, :natural)

# Example
```julia
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, FastDecoupledXB())  # 2-3x faster than Newton for large systems

# Tighter tolerance
sol = solve(prob, FastDecoupledXB(abstol=1e-8, maxiters=50))
```
"""
function FastDecoupledXB(; maxiters::Int=100, abstol::Float64=1e-6, ordering::Symbol=:AMD)
    return AlgorithmConfig(
        FastDecoupled(variant=:XB, maxiters=maxiters, abstol=abstol),
        NativeBackend(),
        KLUSolver(ordering=ordering)
    )
end

"""
    FastDecoupledBX(; kwargs...) -> AlgorithmConfig

Fast Decoupled Power Flow with BX formulation (includes resistance).

BX variant includes resistance effects: better for moderate X/R ratios.

# Keyword Arguments
- `maxiters::Int=100`: Maximum iterations
- `abstol::Float64=1e-6`: Convergence tolerance
- `ordering::Symbol=:AMD`: KLU ordering

# Example
```julia
prob = ACPowerFlowProblem(sys, PowerInjection())
sol = solve(prob, FastDecoupledBX())
```
"""
function FastDecoupledBX(; maxiters::Int=100, abstol::Float64=1e-6, ordering::Symbol=:AMD)
    return AlgorithmConfig(
        FastDecoupled(variant=:BX, maxiters=maxiters, abstol=abstol),
        NativeBackend(),
        KLUSolver(ordering=ordering)
    )
end

# ============================================================================
# Unified Solve Interface
# ============================================================================

"""
    solve(prob::PowerFlowProblem, alg=default_algorithm(prob); kwargs...)

Solve a power flow problem with the specified algorithm.

This is the main entry point for the modern API. Matches the interface used by
DifferentialEquations.jl, NonlinearSolve.jl, and Optimization.jl.

# Arguments
- `prob::PowerFlowProblem`: Problem to solve (e.g., ACPowerFlowProblem)
- `alg`: Algorithm to use (defaults to problem-specific default)
- `kwargs...`: Additional solver-specific options

# Returns
Solution object with fields:
- `converged::Bool`: Whether solver converged
- `iterations::Int`: Number of iterations taken
- `residual_norm::Float64`: Final residual norm
- Additional fields depending on solver

# Example

```julia
using PowerfulCore, PowerfulPF, PowerfulCases

# Load system
case = load("ieee14")
sys = SystemModel(parse_network(case.raw))

# Create problem
prob = ACPowerFlowProblem(sys)

# Solve with default algorithm
sol = solve(prob)

# Solve with specific algorithm
sol = solve(prob, Newton(abstol=1e-10))

# Solve with robust method
sol = solve(prob, RobustNR())

# Check convergence
if sol.converged
    println("Converged in ", sol.iterations, " iterations")
    println("Final residual: ", sol.residual_norm)
else
    @warn "Did not converge"
end
```

See also: [`ACPowerFlowProblem`](@ref), [`Newton`](@ref), [`RobustNR`](@ref)
"""
function solve(prob::PowerFlowProblem, alg::PowerFlowAlgorithm=default_algorithm(prob); kwargs...)
    _validate_problem_algorithm(prob, alg)
    return _solve(prob, alg; kwargs...)
end

function solve(prob::PowerFlowProblem, alg::AlgorithmConfig; kwargs...)
    _validate_problem_algorithm(prob, alg)
    return _solve(prob, alg; kwargs...)
end

# ============================================================================
# GPU Support Trait
# ============================================================================

"""
    supports_gpu(alg) -> Bool

Whether `alg` can solve a problem whose system lives on the GPU. Only native
Newton with the CUDSS linear solver (`Newton(linear_solver=CUDSSSolver())`) can;
the NonlinearSolve backend (`RobustNR`), fast decoupled and KLU run on the CPU. A
bare method (e.g. `NewtonMethod()`, `DishonestNewton()`) solves with KLU.
"""
supports_gpu(alg::AlgorithmConfig) = supports_gpu(alg.method, alg.backend, alg.linear_solver)
supports_gpu(::NewtonMethod, ::NativeBackend, ::CUDSSSolver) = true
supports_gpu(::PowerFlowAlgorithm, ::ImplementationBackend, ::LinearSolver) = false
supports_gpu(::PowerFlowAlgorithm) = false

# ============================================================================
# Internal Validation
# ============================================================================

"""
    _validate_problem_algorithm(prob, alg)

Validate that algorithm is compatible with problem type and execution backend.
"""
function _validate_problem_algorithm(prob::ACPowerFlowProblem, alg::AlgorithmConfig)
    if !(alg.method isa ACMethod)
        error("ACPowerFlowProblem requires an ACMethod algorithm, got $(typeof(alg.method))")
    end

    # Check GPU compatibility
    if prob.sys isa SystemModel && prob.sys.array_type != Vector
        if !supports_gpu(alg)
            error("$(nameof(typeof(alg.method))) with $(nameof(typeof(alg.backend))) and " *
                  "$(nameof(typeof(alg.linear_solver))) does not support GPU execution. " *
                  "Use Newton(linear_solver=CUDSSSolver()) for GPU problems.")
        end
    end

    return nothing
end

function _validate_problem_algorithm(prob::ACPowerFlowProblem, alg::ACMethod)
    # Check GPU compatibility for bare method (no config)
    if prob.sys isa SystemModel && prob.sys.array_type != Vector
        if !supports_gpu(alg)
            error("Algorithm $(typeof(alg)) does not support GPU execution. " *
                  "Use Newton(linear_solver=CUDSSSolver()) for GPU problems.")
        end
    end
    return nothing
end

# ============================================================================
# Internal Dispatch to Implementations
# ============================================================================

"""
    _solve(prob, alg; kwargs...)

Internal dispatch to actual solver implementations based on algorithm configuration.
"""

# Native Newton dispatch (KLU on CPU, CUDSS on GPU)
function _solve(
    prob::ACPowerFlowProblem,
    alg::AlgorithmConfig{<:NewtonMethod, NativeBackend, <:Union{KLUSolver, CUDSSSolver}};
    verbose::Bool = false,
    kwargs...
)
    if prob.sys isa Vector{SystemModel}
        alg.linear_solver isa CUDSSSolver &&
            error("Multi-island GPU solving not yet supported. Split islands manually or use KLUSolver.")
        return _solve_multi_island(prob, alg; verbose=verbose, kwargs...)
    end

    # Single-island: preserve problem (including address_manager with its strategy)
    workspace = setup_powerflow_newton(
        prob;
        _newton_setup_kwargs(alg.linear_solver)...,
        verbose = verbose
    )

    y0 = _initial_guess(prob, workspace)

    return solve_powerflow_newton!(
        workspace, y0;
        abstol = alg.method.abstol,
        maxiters = alg.method.maxiters,
        verbose = verbose
    )
end

# setup_powerflow_newton keywords selecting the linear solver
_newton_setup_kwargs(ls::KLUSolver) = (ordering = ls.ordering,)
_newton_setup_kwargs(::CUDSSSolver) = (use_gpu = true,)

# Dishonest Newton: native Newton with Jacobian refresh every `update_freq` iterations
function _solve(prob::ACPowerFlowProblem, alg::DishonestNewton; verbose::Bool = false, kwargs...)
    prob.sys isa Vector{SystemModel} &&
        error("DishonestNewton does not support multi-island systems yet. Use Newton().")

    workspace = setup_powerflow_newton(prob; verbose = verbose)
    return solve_powerflow_newton!(
        workspace, _initial_guess(prob, workspace);
        abstol = alg.abstol,
        maxiters = alg.maxiters,
        jacobian_every = alg.update_freq,
        verbose = verbose
    )
end

"""
    _initial_guess(prob, workspace) -> Vector

Initial iterate for a solve: the problem's `initial_guess` if given
(length-checked), otherwise the workspace's DAE initialization.
"""
function _initial_guess(prob::ACPowerFlowProblem, workspace)
    isnothing(prob.initial_guess) && return copy(workspace.iter_dae.y)
    length(prob.initial_guess) == length(workspace.iter_dae.y) ||
        error("Initial guess length ($(length(prob.initial_guess))) must match " *
              "number of variables ($(length(workspace.iter_dae.y)))")
    return copy(prob.initial_guess)
end

# Note: _solve for NonlinearSolveBackend is provided by PowerfulPFNonlinearSolveExt.
# If not loaded, the MethodError hint registered in __init__() provides guidance.

# Fast Decoupled dispatch
function _solve(
    prob::ACPowerFlowProblem,
    alg::AlgorithmConfig{<:FastDecoupled, NativeBackend, KLUSolver};
    verbose::Bool=false,
    kwargs...
)
    # Fast Decoupled requires PowerInjection formulation (checked in setup)
    # Dispatch to fastdecoupled.jl implementation
    solve_powerflow_fastdecoupled(
        prob;
        variant = alg.method.variant,
        abstol = alg.method.abstol,
        maxiters = alg.method.maxiters,
        verbose = verbose,
        ordering = alg.linear_solver.ordering,
        kwargs...
    )
end

# Convenience: Allow bare method without full AlgorithmConfig
function _solve(prob::ACPowerFlowProblem, alg::NewtonMethod; kwargs...)
    # Convert to full config with defaults
    full_alg = Newton(
        maxiters = alg.maxiters,
        abstol = alg.abstol,
        ordering = :AMD
    )
    return _solve(prob, full_alg; kwargs...)
end

"""
    _solve_multi_island(prob, alg; kwargs...) -> MultiIslandSolution

Solve each island of a multi-island problem with the single-island dispatch of
`alg` (native or NonlinearSolve backend) and aggregate the results.
"""
function _solve_multi_island(prob::ACPowerFlowProblem, alg; kwargs...)
    islands = prob.sys  # Vector{SystemModel}

    # Create island problems preserving addressing strategy, formulation, and backend
    island_probs = [
        ACPowerFlowProblem(
            island,
            prob.formulation;
            addressing = prob.address_manager.strategy,
            auto_split = false,  # Don't recursively split islands
            backend = prob.backend
        )
        for island in islands
    ]

    # Solve all islands using single-island dispatch (preserves all configuration)
    solutions = map(island_prob -> _solve(island_prob, alg; kwargs...), island_probs)

    return _aggregate_island_solutions(islands, solutions)
end

# Helper: Aggregate island solutions into unified MultiIslandSolution
"""
    _aggregate_island_solutions(islands, solutions) -> MultiIslandSolution

Aggregate multiple island solutions into a unified MultiIslandSolution.

Computes overall convergence status (ALL islands must converge), maximum iterations,
and maximum residual norm across all islands.
"""
function _aggregate_island_solutions(islands::Vector{SystemModel}, solutions::Vector)
    all_converged = all(sol.converged for sol in solutions)
    max_iterations = maximum(sol.iterations for sol in solutions)
    max_residual = maximum(sol.residual_norm for sol in solutions)
    retcode = all_converged ? :Success : :SomeIslandsFailed

    return MultiIslandSolution(
        islands,
        solutions,
        all_converged,
        max_iterations,
        max_residual,
        retcode,
        nothing  # stats
    )
end

# ============================================================================
# Exports
# ============================================================================

# Problem types
export PowerFlowProblem, ACPowerFlowProblem

# Algorithm types (user-facing)
export PowerFlowAlgorithm, ACMethod
export Newton, FastDecoupled, DishonestNewton


# Backend types (for advanced users)
export ImplementationBackend, NativeBackend, NonlinearSolveBackend

# Linear solver types (for advanced users)
export LinearSolver, KLUSolver

# Configuration (for advanced users)
export AlgorithmConfig

# Preset constructors (main user interface)
export RobustNR
export FastDecoupledXB, FastDecoupledBX  # Future

# Main interface
export solve
