"""
Eigenvalue Analysis for Small-Signal Stability

This module provides eigenvalue analysis (small-signal stability) for TDS problems.
It computes the state matrix, eigenvalues, and participation factors to identify
potential instability modes.

The state matrix is computed as:
    As = T^{-1} * (fx - fy * gy^{-1} * gx)

where:
- T = diagonal mass matrix for differential states
- fx, fy = ∂f/∂x, ∂f/∂y (differential equation Jacobians)
- gx, gy = ∂g/∂x, ∂g/∂y (algebraic equation Jacobians)

States with zero mass coefficients (M=0) are treated as pseudo-algebraic
and are reduced following the ANDES approach.
"""

using LinearAlgebra: eigen, I, diag
using PrettyTables
using Printf: @sprintf
using DifferentiationInterface
using PowerfulCore: create_adcache

# ============================================================================
# Result Struct
# ============================================================================

"""
    EigenvalueResult{T}

Result of eigenvalue analysis containing eigenvalues, eigenvectors, participation
factors, and stability statistics.

# Fields
- `eigenvalues::Vector{Complex{T}}` - Complex eigenvalues (μ)
- `right_eigenvectors::Matrix{Complex{T}}` - Right eigenvectors (N), each column = one mode
- `left_eigenvectors::Matrix{Complex{T}}` - Left eigenvectors (W)
- `participation_factors::Matrix{T}` - Normalized participation factors |W| * |N|
- `state_names::Vector{Symbol}` - Variable names for each state

# Statistics
- `n_positive::Int` - Count of eigenvalues with positive real part (unstable)
- `n_zeros::Int` - Count of eigenvalues with near-zero real part
- `n_negative::Int` - Count of eigenvalues with negative real part (stable)

# Zero Mass Handling
- `zero_mass_indices::Vector{Int}` - Indices of states with zero mass coefficient (reduced)
- `n_states::Int` - Number of states after reduction
"""
struct EigenvalueResult{T, E<:Number, MR<:AbstractMatrix{E}, ML<:AbstractMatrix}
    eigenvalues::Vector{E}
    right_eigenvectors::MR
    left_eigenvectors::ML  # Can be Matrix or Transpose
    participation_factors::Matrix{T}
    state_names::Vector{Symbol}

    # Statistics
    n_positive::Int
    n_zeros::Int
    n_negative::Int

    # Zero mass handling
    zero_mass_indices::Vector{Int}
    n_states::Int
end

# ============================================================================
# Main Function
# ============================================================================

"""
    compute_eigenvalues(prob::TDSProblem; tol=1e-6, verbose=false) -> EigenvalueResult

Compute eigenvalues and participation factors for small-signal stability analysis.

The state matrix is computed from the DAE Jacobian at the initial steady-state point
using the same AD infrastructure as TDS solve. States with zero mass coefficients
are treated as pseudo-algebraic and reduced.

# Arguments
- `prob::TDSProblem` - Initialized TDS problem (from power flow solution)
- `tol::Float64` - Tolerance for classifying eigenvalues as zero (default: 1e-6)
- `verbose::Bool` - Print progress information (default: false)

# Returns
`EigenvalueResult` containing eigenvalues, eigenvectors, participation factors, and statistics

# Example
```julia
# Create TDS problem from power flow
prob = TDSProblem(setup_tds(sol_pf, dyr); tspan=(0.0, 10.0))

# Compute eigenvalues
result = compute_eigenvalues(prob)

# Check stability
if result.n_positive > 0
    @warn "System has \$(result.n_positive) unstable modes"
end

# Print report
print_eigenvalue_report(result)
```
"""
function compute_eigenvalues(
    prob::TDSProblem;
    tol::Float64 = 1e-6,
    verbose::Bool = false
)
    @vprint "="^60
    @vprint "Eigenvalue Analysis"
    @vprint "="^60

    # 1. Create ADCache and TDSAdapter (same as solve())
    verbose && print("  Creating ADCache and TDSAdapter... ")
    t_setup = @elapsed begin
        adcache = create_adcache(prob)
        adapter = TDSAdapter(adcache, prob)
    end
    @vprint "done ($(round(t_setup*1000, digits=1)) ms)"

    # 2. Get dimensions and initial state
    n_diff = adapter.n_diff
    n_alg = adapter.n_alg
    am = prob.address_manager

    @vprint "  States: $n_diff differential, $n_alg algebraic"

    # 3-4. Full DAE Jacobian at the operating point, the adapter's initial state, which is the
    # state at the start of the problem
    verbose && print("  Computing Jacobian... ")
    t_jac = @elapsed J_full, _ = dae_jacobian(adapter, prob.tspan[1])
    @vprint "done ($(round(t_jac*1000, digits=1)) ms, nnz=$(nnz(J_full)))"

    # 5-11. Reduction to the retained states, eigenvalues and participation factors
    M_diag = diag(adapter.mass_matrix)[1:n_diff]
    result, _ = _eigen_analysis(J_full, M_diag, am, n_diff; tol, verbose)
    return result
end

"""
    _eigen_analysis(J_full, M_diag, am, n_diff; tol=1e-6, verbose=false)
        -> (EigenvalueResult, ZeroMassReductionResult)

Eigenvalues, eigenvectors and participation factors of the DAE with Jacobian `J_full`
(`∂[f; g]/∂[x; y]`, the first `n_diff` rows and columns differential) and mass
coefficients `M_diag` of the differential states, reduced to the states with nonzero
mass (`reduce_zero_mass_states`).
"""
function _eigen_analysis(J_full::AbstractMatrix, M_diag::AbstractVector, am, n_diff::Int;
                         tol::Float64 = 1e-6, verbose::Bool = false)
    # 5. Extract sub-matrices
    jac_sub = extract_jacobian_submatrices(J_full, n_diff)

    # 7. Collect state names
    all_state_names = collect_state_names(am, n_diff)

    # 8. Compute state matrix with zero mass handling
    verbose && print("  Computing state matrix... ")
    t_as = @elapsed reduction = reduce_zero_mass_states(
        jac_sub.fx, jac_sub.fy, jac_sub.gx, jac_sub.gy, M_diag, all_state_names
    )
    As = reduction.A
    zero_mass_indices = reduction.zero_mass_indices
    reduced_state_names = reduction.reduced_state_names
    n_states = size(As, 1)
    @vprint "done ($(round(t_as*1000, digits=1)) ms, $(length(zero_mass_indices)) zero-mass states reduced)"

    if n_states == 0
        error("No non-zero mass states found. Cannot compute eigenvalues.")
    end

    # 9. Eigenvalue decomposition
    verbose && print("  Computing eigenvalues... ")
    t_eig = @elapsed begin
        eig_result = eigen(As)
        mu = eig_result.values
        N = eig_result.vectors
    end
    @vprint "done ($(round(t_eig*1000, digits=1)) ms, $n_states eigenvalues)"

    # 10. Compute participation factors
    verbose && print("  Computing participation factors... ")
    t_pf = @elapsed pfactors, W = _compute_participation_factors(mu, N)
    @vprint "done ($(round(t_pf*1000, digits=1)) ms)"

    # 11. Count statistics
    n_positive = count(real.(mu) .> tol)
    n_zeros = count(abs.(real.(mu)) .<= tol)
    n_negative = count(real.(mu) .< -tol)

    @vprint "="^60
    @vprint "Results: $n_positive unstable, $n_zeros marginal, $n_negative stable"
    @vprint "="^60

    result = EigenvalueResult(
        mu, N, W, pfactors, reduced_state_names,
        n_positive, n_zeros, n_negative,
        zero_mass_indices, n_states
    )
    return result, reduction
end

"""
    _compute_participation_factors(mu, N) -> (pfactors, W)

Compute participation factors from eigenvalues and right eigenvectors.

Uses the relationship: W = inv(N)^T (from orthogonality W^T * N = I)
Participation factor: pf_{ij} = |W_{ij}| * |N_{ij}|, normalized by column.
"""
function _compute_participation_factors(mu, N)
    n_states = length(mu)
    ET = eltype(N)  # Could be Float64 or ComplexF64

    # Left eigenvector: W = inv(N)^T (from orthogonality W^T * N = I)
    # N \ I gives inv(N), then transpose
    WT = N \ Matrix{ET}(I, n_states, n_states)
    W = transpose(WT)

    # Participation factors: |W| * |N| (element-wise product of magnitudes)
    pfactor = abs.(W) .* abs.(N)

    # Normalize columns (each mode sums to 1)
    for j in 1:n_states
        col_sum = sum(pfactor[:, j])
        if col_sum > 0
            pfactor[:, j] ./= col_sum
        end
    end

    return round.(real.(pfactor), digits=5), W
end

# ============================================================================
# Reporting Functions
# ============================================================================

"""
    print_eigenvalue_report(result::EigenvalueResult; io::IO=stdout, max_modes=50)

Print formatted eigenvalue analysis report.

# Arguments
- `result::EigenvalueResult` - Result from `compute_eigenvalues`
- `io::IO` - Output stream (default: stdout)
- `max_modes::Int` - Maximum number of modes to show in participation factor table (default: 50)
"""
function print_eigenvalue_report(result::EigenvalueResult; io::IO = stdout, max_modes::Int = 50)
    mu = result.eigenvalues
    n_states = result.n_states

    # Header
    println(io, "\n" * "="^70)
    println(io, "Eigenvalue Analysis Report")
    println(io, "="^70)

    # Statistics
    println(io, "\nSTATISTICS")
    println(io, "-"^70)
    stats_data = [
        result.n_positive result.n_zeros result.n_negative n_states
    ]
    stats_labels = ["Positive (unstable)", "Zero (marginal)", "Negative (stable)", "Total"]
    pretty_table(io, stats_data; column_labels=stats_labels, display_size=(-1, -1))

    # Warning for unstable modes
    if result.n_positive > 0
        println(io, "\n⚠ WARNING: $(result.n_positive) eigenvalue(s) with positive real part detected!")
        println(io, "  The system may be unstable.")
    end

    # Zero mass states (if any)
    if !isempty(result.zero_mass_indices)
        println(io, "\nZERO MASS STATES (treated as pseudo-algebraic)")
        println(io, "-"^70)
        println(io, "  $(length(result.zero_mass_indices)) state(s) with M=0 were reduced.")
    end

    # Eigenvalue data table
    println(io, "\nEIGENVALUE DATA")
    println(io, "-"^70)

    # Compute derived quantities
    mu_real = real.(mu)
    mu_imag = imag.(mu)
    freq = zeros(n_states)      # Damped frequency (Hz)
    ufreq = zeros(n_states)     # Undamped frequency (Hz)
    damping = zeros(n_states)   # Damping ratio (%)

    for i in 1:n_states
        if mu_imag[i] != 0
            ufreq[i] = abs(mu[i]) / (2π)
            freq[i] = abs(mu_imag[i]) / (2π)
            damping[i] = -mu_real[i] / abs(mu[i]) * 100
        end
    end

    # Find most associated state for each mode
    assoc_state = Symbol[]
    for i in 1:n_states
        pf_col = result.participation_factors[:, i]
        max_idx = argmax(pf_col)
        push!(assoc_state, result.state_names[max_idx])
    end

    # Create mode labels with stability markers
    mode_labels = String[]
    for i in 1:n_states
        marker = ""
        if mu_real[i] > 1e-6
            marker = " **"  # Unstable
        elseif abs(mu_real[i]) <= 1e-6
            marker = " *"   # Marginal
        end
        push!(mode_labels, "#$i$marker")
    end

    # Build eigenvalue table
    eig_data = hcat(
        mode_labels,
        assoc_state,
        [@sprintf("%.6e", r) for r in mu_real],
        [@sprintf("%.6e", i) for i in mu_imag],
        [@sprintf("%.4f", f) for f in freq],
        [@sprintf("%.4f", u) for u in ufreq],
        [@sprintf("%.2f", d) for d in damping]
    )
    eig_headers = ["Mode", "Most Associated", "Real", "Imag", "Freq (Hz)", "Nat.Freq (Hz)", "Damping (%)"]
    pretty_table(io, eig_data; column_labels=eig_headers, display_size=(-1, -1))

    println(io, "\n  * = marginal (|Re| ≤ tol), ** = unstable (Re > tol)")

    # Participation factors (if small enough)
    if n_states <= max_modes
        println(io, "\nPARTICIPATION FACTORS")
        println(io, "-"^70)
        println(io, "  Rows = states, Columns = modes")
        println(io, "  Values show how strongly each state participates in each mode")
        println(io)

        # Build participation factor table
        pf_data = result.participation_factors
        pf_headers = vcat([:State], [Symbol("#$i") for i in 1:n_states])

        # Combine state names with participation factors
        full_data = hcat(result.state_names, pf_data)
        pretty_table(io, full_data; column_labels=pf_headers, display_size=(-1, -1))
    else
        println(io, "\nNote: Participation factors table omitted ($n_states modes > $max_modes limit)")
    end

    println(io, "\n" * "="^70)
end

"""
    write_eigenvalue_report(result::EigenvalueResult, filename::String; kwargs...)

Write eigenvalue analysis report to file.

# Arguments
- `result::EigenvalueResult` - Result from `compute_eigenvalues`
- `filename::String` - Output file path
- `kwargs...` - Passed to `print_eigenvalue_report`
"""
function write_eigenvalue_report(result::EigenvalueResult, filename::String; kwargs...)
    open(filename, "w") do io
        print_eigenvalue_report(result; io=io, kwargs...)
    end
    println("Eigenvalue report written to $filename")
end

# ============================================================================
# Exports
# ============================================================================

export compute_eigenvalues, EigenvalueResult
export print_eigenvalue_report, write_eigenvalue_report
