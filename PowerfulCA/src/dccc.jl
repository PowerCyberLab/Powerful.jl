"""
    dccc.jl

DC Contingency Calculation (Screening) using LODF.

Fast contingency screening using DC power flow and Line Outage Distribution Factors.
Much faster than full AC contingency calculation, suitable for N-1 screening of large
contingency lists to identify critical contingencies for detailed AC analysis.

Key equations:
    P_base[l] = b[l] * (θ_from - θ_to)             # DC power flow
    P_pre     = P_base - Σ_g pg[g] * PTDF[:, bus(g)] # generator outages, slack absorbs
    P_post    = P_pre + H[:,K] * (I - H[K,K])⁻¹ * P_pre[K]  # outage of branch set K

with `H = PTDF * Aᵀ`. The last line is the multiple-outage compensation form of the
LODF (Güler, Gross and Liu, "Generalized Line Outage Distribution Factors",
IEEE Trans. Power Systems 22(2), 2007); for a single outage it reduces to
`P_post[m] = P_pre[m] + LODF[m,l] * P_pre[l]`.
"""

using PowerfulPF: get_slack_bus, dc_bus_injections

# ============================================================================
# DC Screening Result Types
# ============================================================================

"""
    BranchViolation

A branch thermal limit violation detected during contingency screening.

# Fields
- `branch_idx::Int`: Internal branch index (in incidence ordering)
- `model_name::Symbol`: Model type (:Branch or :Transformer2W)
- `instance_idx::Int`: Instance index within model
- `flow::Float64`: Post-outage power flow (MW or pu)
- `limit::Float64`: Thermal limit (same units as flow)
- `loading::Float64`: Percentage loading (flow/limit * 100)
"""
struct BranchViolation
    branch_idx::Int
    model_name::Symbol
    instance_idx::Int
    flow::Float64
    limit::Float64
    loading::Float64
end

"""
    DCScreeningResult

Result for a single contingency from DC screening.

# Fields
- `name::String`: Contingency name
- `outaged_branches::Vector{Int}`: Incidence rows of the in-service branches the
  contingency takes out (empty for a generator-only contingency)
- `outaged_generators::Vector{Int}`: PV instance indices of the in-service
  generators the contingency takes out
- `violations::Vector{BranchViolation}`: Detected thermal violations
- `max_loading::Float64`: Maximum branch loading (%) after outage
- `is_radial::Bool`: True if the branch outages disconnect the network (flows undefined)
"""
struct DCScreeningResult
    name::String
    outaged_branches::Vector{Int}
    outaged_generators::Vector{Int}
    violations::Vector{BranchViolation}
    max_loading::Float64
    is_radial::Bool
end

"""
    DCContingencyScreenResult

Complete result from DC contingency screening.

# Fields
- `base_flows::Vector{Float64}`: Base case DC power flows (same ordering as incidence)
- `branch_indices::Vector{Tuple{Symbol, Int}}`: Branch model mapping
- `branch_limits::Vector{Float64}`: Thermal limits for each branch
- `results::Vector{DCScreeningResult}`: Per-contingency screening results
- `critical_contingencies::Vector{Int}`: Indices of contingencies with violations
"""
struct DCContingencyScreenResult{T<:Real}
    base_flows::Vector{T}
    branch_indices::Vector{Tuple{Symbol, Int}}
    branch_limits::Vector{T}
    results::Vector{DCScreeningResult}
    critical_contingencies::Vector{Int}
end

# ============================================================================
# DC Power Flow for Base Case
# ============================================================================

"""
    solve_dc_powerflow(sys::SystemModel) -> (theta, P_flows, branch_indices)

Solve DC power flow to get bus angles and branch flows.

# Returns
- `theta::Vector{Float64}`: Bus voltage angles (radians)
- `P_flows::Vector{Float64}`: Branch active power flows (pu)
- `branch_indices::Vector{Tuple{Symbol, Int}}`: Branch model mapping
"""
function solve_dc_powerflow(sys::SystemModel)
    slack_idx = get_slack_bus(sys)

    (; A, branch_indices, b, from_indices, to_indices, non_slack, F) = dc_network(sys, slack_idx)
    L = length(b)
    N = size(A, 2)

    # Net power injections at each bus (shared DC model, V = 1 pu)
    P_inj_reduced = dc_bus_injections(sys)[non_slack]

    # Solve for angles: B * theta = P
    theta = _expand_reduced!(zeros(Float64, N), F \ P_inj_reduced, non_slack)

    # Compute branch flows: P[l] = b[l] * (theta_from - theta_to)
    P_flows = _dc_branch_flows!(zeros(Float64, L), b, from_indices, to_indices, theta)

    return (theta = theta, flows = P_flows, branch_indices = branch_indices, incidence = A, susceptances = b, from_indices = from_indices, to_indices = to_indices, factorization = F, non_slack = non_slack, slack_idx = slack_idx)
end

# ============================================================================
# Branch Limit Gathering
# ============================================================================

"""
    gather_branch_limits(sys::SystemModel, branch_indices::Vector{Tuple{Symbol, Int}}) -> Vector{Float64}

Gather thermal limits (Rate A) for each branch in incidence order, converted to per-unit.

# Arguments
- `sys::SystemModel`: System model
- `branch_indices`: Mapping from incidence row to (model_name, instance_idx)

# Returns
Vector of thermal limits in per-unit (same as DC power flow output).
Limits from PSS/E (MVA) are divided by system MVA base.
"""
function gather_branch_limits(sys::SystemModel, branch_indices::Vector{Tuple{Symbol, Int}})
    n_branches = length(branch_indices)
    limits = Vector{Float64}(undef, n_branches)

    # Get system MVA base for per-unit conversion
    sbase = sys.parsed_data.system_base

    for (l, (model_name, instance_idx)) in enumerate(branch_indices)
        model = sys.models[model_name]

        if model_name == :Branch
            # rate_a is in MVA, convert to per-unit
            limits[l] = model.rate_a[instance_idx] / sbase
        elseif model_name == :Transformer2W
            # rata1 is in MVA, convert to per-unit
            limits[l] = model.rata1[instance_idx] / sbase
        else
            error("Unknown branch model type: $model_name")
        end
    end

    return limits
end

# ============================================================================
# Main DC Screening Function
# ============================================================================

"""
    screen_contingencies_dc(sys::SystemModel, contingencies; kwargs...) -> DCContingencyScreenResult

Fast DC contingency screening using LODF.

Estimates post-outage branch flows and identifies thermal limit violations
without solving full power flow for each contingency.

# Arguments
- `sys::SystemModel`: System model
- `contingencies`: Vector of `Contingency` objects

# Keyword Arguments
- `loading_threshold::Real=100.0`: Report violations above this loading (%)
- `use_lazy_ptdf::Bool=false`: Use LazyPTDF for memory efficiency (large systems)
- `verbose::Bool=false`: Print progress

# Returns
`DCContingencyScreenResult` with base flows, per-contingency results, and critical contingencies.

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

result = screen_contingencies_dc(sys, contingencies; verbose=true)

# Show critical contingencies
for idx in result.critical_contingencies
    r = result.results[idx]
    println("\$(r.name): \$(length(r.violations)) violations, max loading \$(round(r.max_loading, digits=1))%")
end
```

# Algorithm
1. Solve DC power flow to get base case flows
2. Compute the PTDF matrix (or use LazyPTDF)
3. For each contingency:
   - Shift each lost generator's output to the slack through its PTDF column
   - Remove the outaged branches with the multiple-outage LODF (see file header)
   - Check for thermal violations
4. Return results sorted by severity

# Limitations
- The slack bus picks up all lost generation
- Uses DC power flow approximation (assumes flat voltage profile)
- Thermal limits use Rate A (normal operation)

# See Also
- [`solve_accc`](@ref): Full AC contingency calculation
- [`estimate_post_outage_flows`](@ref): Single contingency flow estimation
"""
function screen_contingencies_dc(
    sys::SystemModel,
    contingencies;
    loading_threshold::Real = 100.0,
    use_lazy_ptdf::Bool = false,
    verbose::Bool = false
)
    n_contingencies = length(contingencies)
    verbose && println("DC Contingency Screening: $n_contingencies contingencies")

    # Step 1: Solve DC power flow
    verbose && print("  [1/3] Solving DC power flow... ")
    t_dc = @elapsed begin
        dc_result = solve_dc_powerflow(sys)
        base_flows = dc_result.flows
        branch_indices = dc_result.branch_indices
        A = dc_result.incidence
    end
    verbose && println("done ($(round(t_dc, digits=3))s)")

    # Step 2: PTDF columns (full matrix, or on demand from the DC factorization)
    verbose && print("  [2/3] Computing PTDF... ")
    t_ptdf = @elapsed begin
        ptdf = if use_lazy_ptdf
            LazyPTDF(
                dc_result.factorization,
                A, dc_result.susceptances, branch_indices,
                dc_result.from_indices, dc_result.to_indices,
                dc_result.slack_idx, dc_result.non_slack
            )
        else
            compute_ptdf(sys; slack_bus = dc_result.slack_idx)
        end
    end
    verbose && println("done ($(round(t_ptdf, digits=3))s)")

    # Step 3: Gather branch limits (per-unit, same as the flows)
    limits = gather_branch_limits(sys, branch_indices)
    L = length(base_flows)

    # Step 4: Screen contingencies
    verbose && print("  [3/3] Screening contingencies... ")
    results = Vector{DCScreeningResult}(undef, n_contingencies)
    critical_indices = Int[]

    t_screen = @elapsed for (i, cont) in enumerate(contingencies)
        outaged_branches = Int[]
        outaged_generators = Int[]
        for action in cont.actions
            _collect_outage!(outaged_branches, outaged_generators, sys, branch_indices, action)
        end
        unique!(outaged_branches)
        unique!(outaged_generators)

        # Generator outages: the slack picks up the lost injection
        pre_flows = _flows_after_generator_outages(base_flows, ptdf, sys, outaged_generators)
        post_flows = _flows_after_branch_outages(pre_flows, ptdf, dc_result.from_indices, dc_result.to_indices, outaged_branches)
        if isnothing(post_flows)
            results[i] = DCScreeningResult(
                cont.name, outaged_branches, outaged_generators, BranchViolation[], Inf, true
            )
            push!(critical_indices, i)
            continue
        end

        violations = BranchViolation[]
        max_loading = 0.0

        for m in 1:L
            m in outaged_branches && continue  # Outaged branches carry no flow

            post_flow = post_flows[m]
            limit_m = limits[m]

            # Compute loading percentage
            loading = if limit_m > 0
                abs(post_flow) / limit_m * 100.0
            else
                0.0  # No limit specified
            end

            max_loading = max(max_loading, loading)

            # Check for violation
            if loading > loading_threshold
                model_name, instance_idx = branch_indices[m]
                push!(violations, BranchViolation(
                    m, model_name, instance_idx,
                    post_flow, limit_m, loading
                ))
            end
        end

        results[i] = DCScreeningResult(
            cont.name, outaged_branches, outaged_generators, violations, max_loading, false
        )

        if !isempty(violations)
            push!(critical_indices, i)
        end
    end
    verbose && println("done ($(round(t_screen, digits=3))s)")

    # Summary
    if verbose
        n_critical = length(critical_indices)
        n_radial = count(r -> r.is_radial, results)
        total_violations = sum(r -> length(r.violations), results)

        println("  Summary:")
        println("    Critical contingencies: $n_critical / $n_contingencies")
        if n_radial > 0
            println("    Radial outages: $n_radial (network disconnection)")
        end
        println("    Total violations: $total_violations")
        println("    Total time: $(round(t_dc + t_ptdf + t_screen, digits=3))s")
    end

    return DCContingencyScreenResult(
        base_flows, branch_indices, limits, results, critical_indices
    )
end

"""
    _collect_outage!(branches, generators, sys, branch_indices, action)

Record the device taken out by `action`: its incidence row for a Branch or
Transformer2W, its PV index for a generator. Devices already out of service are
not part of the DC case and are skipped. Errors if the device does not exist.
"""
function _collect_outage!(branches, generators, sys::SystemModel, branch_indices,
                          action::Union{BranchOutage, Transformer2WOutage})
    row = findfirst(==(resolve_outage(sys, action)), branch_indices)
    isnothing(row) || push!(branches, row)
    return nothing
end

function _collect_outage!(branches, generators, sys::SystemModel, branch_indices,
                          action::GeneratorOutage)
    model, idx = resolve_outage(sys, action)
    in_service(sys.models[model])[idx] && push!(generators, idx)
    return nothing
end

_ptdf_column(ptdf::PTDF, bus::Int) = view(ptdf.matrix, :, bus)
_ptdf_column(ptdf::LazyPTDF, bus::Int) = get_ptdf_column(ptdf, bus)

"""
    _flows_after_generator_outages(flows, ptdf, sys, generators) -> Vector

Flows after taking the PV generators `generators` offline, with the slack bus
picking up their output: `flows - Σ pg[g] * PTDF[:, bus(g)]`.
"""
function _flows_after_generator_outages(flows::Vector{Float64}, ptdf, sys::SystemModel, generators::Vector{Int})
    isempty(generators) && return flows
    pv = sys.models[:PV]
    pv_bus = get_foreign_indices(sys, :PV, :i)
    shifted = copy(flows)
    for g in generators
        shifted .-= pv.pg[g] .* _ptdf_column(ptdf, pv_bus[g])
    end
    return shifted
end

"""
    _flows_after_branch_outages(pre_flows, ptdf, from_indices, to_indices, K) -> flows or nothing

Flows after removing the branches in rows `K`, from the intact-network flows
`pre_flows`: `pre + H[:,K] (I - H[K,K])⁻¹ pre[K]` with `H[:,l] = PTDF[:,from(l)] -
PTDF[:,to(l)]`, and zero on the removed branches. Returns `nothing` when
`I - H[K,K]` is singular, i.e. the outages split the network.
"""
function _flows_after_branch_outages(pre_flows::Vector{Float64}, ptdf, from_indices, to_indices, K::Vector{Int})
    isempty(K) && return pre_flows
    H_K = stack(_ptdf_column(ptdf, from_indices[l]) .- _ptdf_column(ptdf, to_indices[l]) for l in K)
    M = I - H_K[K, :]
    minimum(svdvals(M)) < 1e-10 && return nothing
    post = pre_flows .+ H_K * (M \ pre_flows[K])
    post[K] .= 0.0
    return post
end

# ============================================================================
# Reporting
# ============================================================================

"""
    print_dc_screening_report(result::DCContingencyScreenResult; max_violations::Int=10)

Print a summary report of DC contingency screening results.

# Arguments
- `result`: DCContingencyScreenResult from screen_contingencies_dc
- `max_violations::Int=10`: Maximum violations to show per contingency
"""
function print_dc_screening_report(result::DCContingencyScreenResult; max_violations::Int = 10)
    println("=" ^ 70)
    println("DC Contingency Screening Report")
    println("=" ^ 70)
    println()

    n_total = length(result.results)
    n_critical = length(result.critical_contingencies)

    println("Overview:")
    println("  Total contingencies: $n_total")
    println("  Critical (violations): $n_critical")
    println()

    if isempty(result.critical_contingencies)
        println("No thermal violations detected.")
        return
    end

    # Sort critical contingencies by max loading
    sorted_critical = sort(
        result.critical_contingencies,
        by = i -> result.results[i].max_loading,
        rev = true
    )

    println("Critical Contingencies (sorted by max loading):")
    println("-" ^ 70)

    for (rank, idx) in enumerate(sorted_critical)
        r = result.results[idx]

        if r.is_radial
            println("$rank. $(r.name) [RADIAL - Network Disconnection]")
        else
            println("$rank. $(r.name)")
            println("   Max loading: $(round(r.max_loading, digits=1))%")
            println("   Violations: $(length(r.violations))")

            # Show top violations
            sorted_violations = sort(r.violations, by = v -> v.loading, rev = true)
            for (j, v) in enumerate(sorted_violations)
                if j > max_violations
                    println("   ... and $(length(r.violations) - max_violations) more")
                    break
                end
                println("   - Branch $(v.branch_idx) ($(v.model_name)): $(round(v.loading, digits=1))% ($(round(v.flow, digits=3)) pu / $(round(v.limit, digits=3)) pu)")
            end
        end
        println()
    end
end
