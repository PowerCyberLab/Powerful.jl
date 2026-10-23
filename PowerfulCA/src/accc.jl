# AC Contingency Calculation (ACCC) - Sequential and Parallel Implementations
#
# Architecture:
# - Both solve_accc (sequential) and solve_accc_parallel share _accc_setup and
#   _solve_contingency!, which run each contingency on a power-flow NewtonWorkspace
# - Status dicts are used instead of modifying the system directly (thread-safe)
#
# Each worker owns a full NewtonWorkspace (adcache, jac_prep, status dict, DAE
# vectors, J matrix, factorization): adcache and jac_prep hold mutable internal
# buffers (SystemBuffers, DifferentiationInterface caches) and are NOT thread-safe.

using SparseArrays: nnz

"""
    solve_accc(sys::SystemModel, contingencies; kwargs...) -> NamedTuple

Sequential AC Contingency Calculation with infrastructure reuse.

Returns a named tuple with:
- `base`: Base case PowerFlowSolution
- `contingencies`: Vector of per-contingency results (name, converged, iterations,
  residual_norm, y, out_of_service, singular); `out_of_service` holds the instances
  the contingency switched off by model, including the buses (`:Bus`) and devices of
  islands it leaves without a slack bus

# Keyword Arguments
- `abstol::Real=1e-8`: Absolute tolerance for convergence
- `maxiters::Int=20`: Maximum Newton iterations per contingency
- `damping::Real=1.0`: Step damping factor
- `ordering::Symbol=:AMD`: KLU fill-reducing ordering (`:AMD`, `:COLAMD`, `:natural`)
- `verbose::Bool=false`: Print progress

# Example
```julia
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))

contingencies = [
    Contingency("LINE_1_2", [BranchOutage(1, 2, "1")]),
    Contingency("GEN_3", [GeneratorOutage(3, "1")]),
    Contingency("N2_TEST", [BranchOutage(1, 2, "1"), GeneratorOutage(3, "1")]),
]

result = solve_accc(sys, contingencies; verbose=true)

for r in result.contingencies
    println("\$(r.name): \$(r.converged ? "OK" : "FAILED")")
end
```

# See Also
- [`solve_accc_parallel`](@ref): Multi-threaded version for large contingency lists
"""
function solve_accc(
    sys::SystemModel,
    contingencies;
    abstol::Real = 1e-8,
    maxiters::Int = 20,
    damping::Real = 1.0,
    ordering::Symbol = :AMD,
    verbose::Bool = false
)
    n_contingencies = length(contingencies)
    verbose && println("Sequential ACCC: $n_contingencies contingencies")

    ctx = _accc_setup(sys, 1; ordering, verbose)
    worker_ws = only(ctx.workers)

    verbose && print("  [3/3] Solving $n_contingencies contingencies... ")
    results = Vector{NamedTuple}(undef, n_contingencies)
    t_cont = @elapsed for (i, cont) in enumerate(contingencies)
        results[i] = _solve_contingency!(worker_ws, cont, ctx; abstol, maxiters, damping)
    end
    verbose && println("done ($(round(t_cont, digits=3))s)")

    verbose && _print_accc_summary(results, ctx.t_setup + t_cont)

    return (base = ctx.base_sol, contingencies = results)
end

"""
    _accc_setup(sys, n_workers; ordering, verbose) -> NamedTuple

Shared ACCC setup: solve the base case (warm start for every contingency), build
the base status, topology and bus→device map for island de-energization, and one
NewtonWorkspace per worker.
"""
function _accc_setup(sys::SystemModel, n_workers::Int; ordering::Symbol, verbose::Bool)
    verbose && print("  [1/3] Solving base case... ")
    t_base = @elapsed begin
        base_prob = ACPowerFlowProblem(sys)
        base_sol = solve(base_prob, Newton())
    end
    verbose && println("done ($(round(t_base, digits=3))s)")

    base_sol.converged || error("Base case did not converge")

    verbose && print("  [2/3] Setting up $n_workers worker workspace(s)... ")
    t_workers = @elapsed begin
        base_status = build_status_dict(sys)
        topo = analyze_topology(sys; build_incidence=true)
        bus_device_map = build_reverse_terminal_map(sys, :Bus)
        # Each worker owns its status dict, copied from the shared base
        workers = [setup_powerflow_newton(base_prob; ordering, status=copy_status_dict(base_status))
                   for _ in 1:n_workers]
    end
    verbose && println("done ($(round(t_workers, digits=3))s)")
    n_vars = length(base_sol.dae.y)
    verbose && println("        Jacobian: $(n_vars)x$(n_vars), $(nnz(first(workers).J)) nonzeros")
    verbose && println("        Bridges: $(count(topo.is_bridge)) / $(length(topo.is_bridge))")

    return (; base_sol, y_warm = copy(base_sol.dae.y), base_status, topo, bus_device_map,
            workers, t_setup = t_base + t_workers)
end

"""
    _solve_contingency!(worker_ws, cont, ctx; abstol, maxiters, damping) -> NamedTuple

Solve one contingency on a worker's NewtonWorkspace: reset its status to the base
case, apply the outages and de-energize the islands they leave without a slack bus
(`apply_contingency!`), and run Newton from the base-case solution.

Returns `(name, converged, iterations, residual_norm, y, out_of_service, singular)`:
`out_of_service` holds the instances the contingency switched off, keyed by model
(`:Bus` for de-energized buses). A singular Jacobian is a contingency without a
solution: `converged = false`, `singular = true` and `y = nothing`.
"""
function _solve_contingency!(worker_ws, cont::Contingency, ctx; abstol::Real, maxiters::Int, damping::Real)
    status = worker_ws.status
    for (k, v) in ctx.base_status
        copyto!(status[k], v)
    end

    out_of_service = apply_contingency!(status, worker_ws.prob.sys, cont, ctx.topo, ctx.bus_device_map)

    sol = try
        solve_powerflow_newton!(worker_ws, ctx.y_warm; abstol, maxiters, damping)
    catch err
        err isa SingularException || rethrow()
        nothing
    end

    isnothing(sol) && return (name = cont.name, converged = false, iterations = 0,
        residual_norm = NaN, y = nothing, out_of_service, singular = true)
    return (
        name = cont.name,
        converged = sol.converged,
        iterations = sol.iterations,
        residual_norm = sol.residual_norm,
        y = sol.dae.y,
        out_of_service,
        singular = false,
    )
end

function _print_accc_summary(results, t_total)
    n_contingencies = length(results)
    n_converged = count(r -> r.converged, results)
    n_failed = n_contingencies - n_converged
    avg_iters = sum(r -> r.iterations, results) / n_contingencies

    println("  Summary:")
    println("    Converged: $n_converged / $n_contingencies")
    if n_failed > 0
        println("    Failed: $n_failed")
    end
    println("    Avg iterations: $(round(avg_iters, digits=1))")
    println("    Total time: $(round(t_total, digits=3))s")
end

"""
    solve_accc_parallel(sys::SystemModel, contingencies; kwargs...) -> NamedTuple

Parallel AC Contingency Calculation using multi-threading.

Shares the base-case solution, topology and bus→device map across all threads,
while each worker task owns a NewtonWorkspace with its own status dictionary,
AD cache, Jacobian preparation, DAE vectors and factorization.

# Arguments
- `sys::SystemModel`: System model
- `contingencies`: Vector of `Contingency` objects

# Keyword Arguments
- `abstol::Real=1e-8`: Absolute tolerance for convergence
- `maxiters::Int=20`: Maximum Newton iterations per contingency
- `damping::Real=1.0`: Step damping factor
- `ordering::Symbol=:AMD`: KLU fill-reducing ordering (`:AMD`, `:COLAMD`, `:natural`)
- `verbose::Bool=false`: Print progress

# Returns
Named tuple with:
- `base`: Base case PowerFlowSolution
- `contingencies`: Vector of per-contingency results

# Example
```julia
using Powerful
using PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))

contingencies = [
    Contingency("LINE_1_2", [BranchOutage(1, 2, "1")]),
    Contingency("LINE_1_5", [BranchOutage(1, 5, "1")]),
    Contingency("LINE_2_3", [BranchOutage(2, 3, "1")]),
]

# Parallel solve using all available threads
result = solve_accc_parallel(sys, contingencies; verbose=true)

# Results
println("Processed \$(length(contingencies)) contingencies using \$(Threads.nthreads()) threads")
for r in result.contingencies
    status = r.converged ? "OK" : "FAILED"
    println("  \$(r.name): \$status (\$(r.iterations) iters)")
end
```

# Performance Notes
- Set `JULIA_NUM_THREADS` environment variable before starting Julia
- Typical speedup is near-linear up to 8-16 threads for large contingency lists
- For small contingency counts (<10), sequential `solve_accc` may be faster due to threading overhead

# See Also
- [`solve_accc`](@ref): Sequential version
- [`Contingency`](@ref): Contingency definition
- [`BranchOutage`](@ref), [`GeneratorOutage`](@ref): Outage types
"""
function solve_accc_parallel(
    sys::SystemModel,
    contingencies;
    abstol::Real = 1e-8,
    maxiters::Int = 20,
    damping::Real = 1.0,
    ordering::Symbol = :AMD,
    verbose::Bool = false
)
    n_contingencies = length(contingencies)
    n_threads = Threads.nthreads()

    verbose && println("Parallel ACCC: $n_contingencies contingencies, $n_threads threads")

    ctx = _accc_setup(sys, n_threads; ordering, verbose)

    verbose && print("  [3/3] Solving contingencies in parallel... ")
    results = Vector{NamedTuple}(undef, n_contingencies)

    errors = Vector{Union{Nothing, Exception}}(nothing, n_contingencies)

    # One task per workspace, each pulling contingency indices from a shared
    # atomic counter. Workspaces are owned by tasks, not indexed by threadid(),
    # which is neither stable under task migration nor bounded by nthreads().
    next_idx = Threads.Atomic{Int}(1)
    t_parallel = @elapsed @sync for worker_ws in ctx.workers
        Threads.@spawn while true
            i = Threads.atomic_add!(next_idx, 1)
            i > n_contingencies && break
            try
                results[i] = _solve_contingency!(worker_ws, contingencies[i], ctx; abstol, maxiters, damping)
            catch e
                errors[i] = e
                results[i] = (
                    name = contingencies[i].name,
                    converged = false,
                    iterations = 0,
                    residual_norm = Inf,
                    y = Float64[],
                )
            end
        end
    end

    # Check for errors and report
    first_error_idx = findfirst(!isnothing, errors)
    if first_error_idx !== nothing
        error("Error in contingency $(contingencies[first_error_idx].name): $(errors[first_error_idx])")
    end
    verbose && println("done ($(round(t_parallel, digits=3))s)")

    if verbose
        _print_accc_summary(results, ctx.t_setup + t_parallel)
        println("    Throughput: $(round(n_contingencies / t_parallel, digits=1)) contingencies/s on $n_threads threads")
    end

    return (base = ctx.base_sol, contingencies = results)
end
