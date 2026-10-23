"""
Eigenvalue Sensitivities to Model Parameters

Derivatives of the eigenvalues of the linearized TDS system with respect to parameters
of its dynamic models (controller gains, time constants, inertias, ...), for tuning a
controller to a required small-signal performance by gradient-based optimization.

A parameter changes the eigenvalues directly, through the Jacobian, and through the
operating point (`parameter_seeding.jl`). Both paths are differentiated with
forward-mode AD: the initialization runs with the parameters as Dual numbers, and the
Jacobian of the DAE is evaluated with the Dual operating point, giving the Jacobian and
its derivative along each parameter. The Dual numbers have one of two fixed widths
(`SENSITIVITY_WIDTHS`), and more parameters than the wider one takes go through in
passes, so two compilations serve any number of parameters. The derivative of the state
matrix follows from the reduction to the states with mass (`reduced_state_matrix_derivative`), and that of each eigenvalue
from its eigenvectors:

    dλ_i = (N⁻¹ dA N)_ii

with `N` the right eigenvectors, which requires distinct eigenvalues.
"""

using PowerfulCore: register_parameters!, scatter_parameters!, gather_parameters
using PowerfulCore: create_di_backend, ForwardDiffBackend, jacobian_sparsity_pattern

"""
    EigenvalueSensitivity

Eigenvalues of the linearized TDS system and their derivatives with respect to model
parameters.

# Fields
- `eigen::EigenvalueResult`: Eigenvalues, eigenvectors and participation factors at the
  operating point (as `compute_eigenvalues`)
- `parameters::Vector{Tuple{Symbol, Symbol, Int}}`: `(model, field, instance)` of each
  parameter, in the order of the columns
- `dlambda::Matrix{ComplexF64}`: `dlambda[i, j]` is the derivative of eigenvalue `i`
  (`eigen.eigenvalues[i]`) with respect to parameter `j`
"""
struct EigenvalueSensitivity{E<:EigenvalueResult}
    eigen::E
    parameters::Vector{Tuple{Symbol, Symbol, Int}}
    dlambda::Matrix{ComplexF64}
end

"""
    damping_ratio_sensitivity(s::EigenvalueSensitivity) -> Matrix{Float64}

Derivatives of the damping ratios `ζ = -σ/|λ|` of the eigenvalues `λ = σ + jω` with
respect to the parameters, `dζ = ω (σ dω - ω dσ) / |λ|³`: zero for a real eigenvalue,
undefined (`NaN`) for a zero eigenvalue.
"""
function damping_ratio_sensitivity(s::EigenvalueSensitivity)
    λ = s.eigen.eigenvalues
    σ, ω = real.(λ), imag.(λ)
    dσ, dω = real.(s.dlambda), imag.(s.dlambda)
    return ω .* (σ .* dω .- ω .* dσ) ./ abs.(λ) .^ 3
end

"""
    eigenvalue_sensitivity(sol_pf, dyn_models, params; config=SimConfig(),
                           formulation=PowerInjection(), addressing=NaturalAddressing(),
                           tol=1e-6, verbose=false) -> EigenvalueSensitivity

Eigenvalues of the TDS system initialized from the power flow `sol_pf` with the dynamic
models `dyn_models` (as `setup_tds`), and their derivatives with respect to the
parameters `params`: `(model, field, instances)` tuples of fields of the dynamic models,
e.g. `[(:SEXS, :K, 1:5), (:TGOV1, :R, [2])]`, one parameter per instance.

The power flow is held; the states and setpoints of the dynamic models follow the
parameters through the initialization. Errors if a parameter is not a real field of one
of `dyn_models`, is given twice, or enters the equations of another model through its
external parameters.

# Example
```julia
dyr = parse_dynamics(file(load("ieee14"), :dyr, variant = "genrou"), sys)
s = eigenvalue_sensitivity(sol_pf, dyr, [(:SEXS, :K, 1:5)])

# Least-damped oscillatory mode and the gradient of its damping ratio
λ = s.eigen.eigenvalues
osc = findall(>(0.1), imag.(λ))
i = osc[argmin(-real.(λ[osc]) ./ abs.(λ[osc]))]
dζ_dK = damping_ratio_sensitivity(s)[i, :]
```
"""
function eigenvalue_sensitivity(
    sol_pf,
    dyn_models::NamedTuple,
    params::Vector{<:Tuple{Symbol, Symbol, <:AbstractVector{<:Integer}}};
    config::SimConfig = SimConfig(),
    formulation::FormulationTrait = PowerInjection(),
    addressing::AddressingStrategy = NaturalAddressing(),
    tol::Float64 = 1e-6,
    verbose::Bool = false
)
    param_list = _parameter_list(dyn_models, params)
    D = _sensitivity_dual_type(length(param_list))
    return _eigenvalue_sensitivity(sol_pf, dyn_models, param_list, D;
                                   config, formulation, addressing, tol, verbose)
end

# `eigenvalue_sensitivity` with the Dual type `D` of the chosen width (a function barrier)
function _eigenvalue_sensitivity(sol_pf, dyn_models, param_list, ::Type{D};
                                 config, formulation, addressing, tol, verbose) where D
    np = length(param_list)

    # Operating point
    prob = TDSProblem(setup_tds(sol_pf, dyn_models; config); formulation, addressing)
    cache = create_adcache(prob)
    adapter = TDSAdapter(cache, prob)
    n_diff = adapter.n_diff
    sparsity = something(tds_jacobian_sparsity(adapter),
                         jacobian_sparsity_pattern(prob.address_manager, prob.sys, prob.formulation;
                                                   include_differential = true))
    backend = create_di_backend(ForwardDiffBackend(), sparsity)

    # The initialization with the parameters of each pass seeded, and the model fields
    # and setpoints that move with them. The fields of every pass are registered before
    # the first Jacobian, so the Dual workspaces are built once.
    passes = collect(Iterators.partition(1:np, ForwardDiff.npartials(D)))
    verbose && print("  Initializing with $np seeded parameter(s) in $(length(passes)) pass(es)... ")
    t_init = @elapsed seeded = [_seeded_operating_point(sol_pf, dyn_models, param_list[chunk], D,
                                                        prob, adapter, cache; config)
                                for chunk in passes]
    all_slices = reduce(vcat, [p.slices for p in seeded])
    register_parameters!(cache, all_slices)
    values_float = gather_parameters(cache, all_slices)
    @vprint "done ($(round(t_init*1000, digits=1)) ms)"

    local eig, reduction, blocks, W_fact, N, N_inv, M_diag
    dlambda = Matrix{ComplexF64}(undef, 0, np)
    for (pass, chunk) in enumerate(passes)
        verbose && print("  Jacobian along parameters $(first(chunk))-$(last(chunk))... ")
        # Fields moved by another pass carry that pass's derivatives in the Dual workspaces:
        # reset every registered field to its value, with zero derivatives
        scatter_parameters!(cache, all_slices, values_float)
        t_jac = @elapsed J_d = _jacobian_along_parameters(seeded[pass], cache, adapter, sparsity, backend, prob.tspan[1])
        @vprint "done ($(round(t_jac*1000, digits=1)) ms)"
        M_d = seeded[pass].M_d

        # Eigenvalues at the operating point (the values of every pass)
        if pass == 1
            J = ForwardDiff.value.(J_d)
            M_diag = ForwardDiff.value.(M_d)
            eig, reduction = _eigen_analysis(J, M_diag, prob.address_manager, n_diff; tol, verbose)
            jac_sub = extract_jacobian_submatrices(J, n_diff)
            blocks = _reduction_blocks(jac_sub.fx, jac_sub.fy, jac_sub.gx, jac_sub.gy,
                                       reduction.retained_indices, reduction.zero_mass_indices)
            W_fact = lu(blocks[3])
            N = eig.right_eigenvectors
            N_inv = transpose(eig.left_eigenvectors)   # left_eigenvectors = inv(N)ᵀ
            dlambda = Matrix{ComplexF64}(undef, length(eig.eigenvalues), np)
        end

        # Derivatives of the eigenvalues along each parameter of the pass
        for (j, k) in enumerate(chunk)
            dJ = extract_jacobian_submatrices(ForwardDiff.partials.(J_d, j), n_diff)
            dA = reduced_state_matrix_derivative(reduction, blocks, W_fact, M_diag,
                                                 (dJ.fx, dJ.fy, dJ.gx, dJ.gy),
                                                 ForwardDiff.partials.(M_d, j))
            dlambda[:, k] = vec(sum(transpose(N_inv) .* (dA * N); dims = 1))
        end
    end

    return EigenvalueSensitivity(eig, param_list, dlambda)
end

"""
    _jacobian_along_parameters(seeded, cache, adapter, sparsity, backend, t) -> SparseMatrixCSC

The DAE Jacobian at the seeded operating point `seeded` (`_seeded_operating_point`) and its
time `t` (the start of the problem, where the initialization puts the operating point), as
Dual numbers whose partials are its derivatives along the seeded parameters: the moving
model fields and setpoints are written, as Dual values, into the workspace of the
Jacobian's element type (registered beforehand).
"""
function _jacobian_along_parameters(seeded, cache, adapter, sparsity, backend, t::Real)
    (; u_d, slices, values) = seeded
    D = eltype(u_d)
    closure = let cache = cache, adapter = adapter, slices = slices, values = values, t = Float64(t)
        (du, u) -> begin
            scatter_parameters!(cache, slices, convert.(eltype(du), values))
            _evaluate_packed!(du, adapter, u, t, nothing)
        end
    end
    du_d = zeros(D, length(u_d))
    prep = DifferentiationInterface.prepare_jacobian(closure, du_d, backend, u_d)
    J_d = SparseMatrixCSC{D, Int}(sparsity)
    DifferentiationInterface.jacobian!(closure, du_d, J_d, prep, backend, u_d)
    return J_d
end

export EigenvalueSensitivity, eigenvalue_sensitivity, damping_ratio_sensitivity
