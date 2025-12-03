# Address Reordering for Optimal Sparse Factorization

## Motivation

The current address assignment is sequential: all `theta[1:n]`, then all `v[1:n]`. This ordering is simple but not optimal for sparse matrix factorization. AMD (Approximate Minimum Degree) ordering can significantly reduce fill-in during factorization, improving both memory usage and computational efficiency.

## Key Insight

For power flow problems, the Jacobian sparsity pattern can be determined from **connectivity alone**:
- Bus variables: `theta[i]`, `v[i]` for each bus `i`
- Bus residuals: `p[i]`, `q[i]` for each bus `i`
- Self-coupling: Each bus's residuals depend on its own variables (diagonal blocks)
- Network coupling: Branches/transformers create off-diagonal dependencies between connected buses

This allows us to:
1. Build the sparsity pattern **before** address assignment (from connectivity data)
2. Run AMD to get optimal ordering
3. Assign addresses according to this ordering
4. Provide permutation to KLU to skip internal AMD

## Architecture

### Three-Stage Process

```julia
# Stage 1: Connectivity-based sparsity (BEFORE addressing)
pattern = build_topology_sparsity(models)  # Only needs bus/branch connectivity

# Stage 2: AMD ordering
perm = amd_ordering(pattern)  # Get optimal permutation

# Stage 3: Readdressing
assign_address!(am, models, perm)  # Assign addresses using permutation
```

### Connectivity-Based Sparsity Pattern

```julia
function build_topology_sparsity(models::NamedTuple) -> SparseMatrixCSC{Bool}
    """
    Build Jacobian sparsity pattern from network topology alone.

    Does NOT require addresses to be assigned yet!

    Algorithm:
    1. Extract bus numbers from Bus model
    2. For each bus i: Add diagonal blocks (theta[i], v[i] affect p[i], q[i])
    3. For each branch (i,j): Add off-diagonal blocks (coupling between buses)
    4. Return sparsity pattern as Boolean sparse matrix
    """
end
```

### Permuted Address Assignment

```julia
function assign_address_permuted!(
    am::AddressManager,
    model_name::Symbol,
    metadata::ModelMetadata,
    count::Int,
    permutation::Vector{Int}  # <-- NEW: permutation from AMD
)
    """
    Assign addresses according to AMD permutation.

    Instead of: [1, 2, 3, ..., n]
    Use:        [perm[1], perm[2], ..., perm[n]]

    This reorders variables/residuals to match optimal AMD ordering.
    """
end
```

### KLU Integration

```julia
# Store permutation in workspace
workspace.permutation = perm

# First factorization: provide permutation to skip AMD
if isnothing(solver_fact)
    solver_fact = klu(J)
    # Customize: Use our permutation instead of internal AMD
    K = solver_fact
    K.common.ordering = 3  # User-provided ordering
    klu_analyze!(K, P, Q)  # Provide permutation
else
    klu!(solver_fact, J)  # Refactor with same pattern
end
```

## Implementation Challenges

### Challenge 1: Variable Numbering in Pattern Matrix

The sparsity pattern is built BEFORE addresses are assigned. We need a temporary numbering scheme:

**Solution**: Use a "logical" numbering based on bus indices and variable types:
```julia
# For each bus i (assuming n buses):
logical_idx(bus_i, :theta) = 2*(i-1) + 1
logical_idx(bus_i, :v)     = 2*(i-1) + 2
# Pattern matrix is 2n × 2n (p and q residuals map to same indices)
```

Then AMD permutation operates on this logical space. When assigning actual addresses, we apply the permutation.

### Challenge 2: Models with Different Key Schemes

Different models use different keys:
- Bus: `i` (bus number)
- PQ: `(i, id)` (bus + load ID)
- Branch: `(i, j, ckt)` (from-bus, to-bus, circuit)

**Solution**: Readdressing operates at the **bus variable level** only:
```julia
# Only reorder bus variables (theta, v) and bus residuals (p, q)
# Other models (PQ, PV, Branch) reference bus variables via indexer
# Their addressing is derived from bus addressing
```

### Challenge 3: Permutation Storage and Use

The permutation must be:
1. Computed once during system setup
2. Stored in the address manager
3. Provided to KLU during first factorization

**Solution**: Add permutation field to AddressManager:
```julia
mutable struct AddressManager
    addresses::Dict{...}
    next_idx::Dict{...}
    assignments::Dict{...}
    retrieved::Dict{...}
    permutation::Union{Nothing, Vector{Int}}  # <-- NEW
end
```

## Sparsity Pattern Construction

### Diagonal Blocks (Self-Coupling)

Each bus `i` has a 2×2 diagonal block in the Jacobian:
```
∂p[i]/∂theta[i]   ∂p[i]/∂v[i]
∂q[i]/∂theta[i]   ∂q[i]/∂v[i]
```

### Off-Diagonal Blocks (Network Coupling)

Each branch connecting buses `i` and `j` creates four 2×2 blocks:
```
# Bus i residuals depend on bus j variables
∂p[i]/∂theta[j]   ∂p[i]/∂v[j]
∂q[i]/∂theta[j]   ∂q[i]/∂v[j]

# Bus j residuals depend on bus i variables
∂p[j]/∂theta[i]   ∂p[j]/∂v[i]
∂q[j]/∂theta[i]   ∂q[j]/∂v[i]
```

### Implementation

```julia
function build_topology_sparsity(models::NamedTuple)
    buses = models[:Bus]
    n_buses = length(buses)
    bus_numbers = buses.i  # Bus numbering from data

    # Build mapping: bus_number → index in array
    bus_to_idx = Dict(buses.i[k] => k for k in 1:n_buses)

    I = Int[]  # Row indices
    J = Int[]  # Column indices

    # Diagonal blocks: Each bus coupled to itself
    for idx in 1:n_buses
        # theta[idx] and v[idx] affect p[idx] and q[idx]
        for res_var in 1:2, col_var in 1:2
            push!(I, 2*(idx-1) + res_var)  # p or q residual
            push!(J, 2*(idx-1) + col_var)  # theta or v variable
        end
    end

    # Off-diagonal blocks: Branch/transformer coupling
    for branch_model_name in [:Branch, :Transformer]
        if haskey(models, branch_model_name)
            branches = models[branch_model_name]
            for k in 1:length(branches)
                i_bus = bus_to_idx[branches.i[k]]
                j_bus = bus_to_idx[branches.j[k]]

                # Bus i depends on bus j
                for res_var in 1:2, col_var in 1:2
                    push!(I, 2*(i_bus-1) + res_var)
                    push!(J, 2*(j_bus-1) + col_var)
                end

                # Bus j depends on bus i
                for res_var in 1:2, col_var in 1:2
                    push!(I, 2*(j_bus-1) + res_var)
                    push!(J, 2*(i_bus-1) + col_var)
                end
            end
        end
    end

    # Build sparse pattern
    m = 2 * n_buses
    return sparse(I, J, ones(Bool, length(I)), m, m)
end
```

## AMD Integration

```julia
using AMD  # AMD.jl package

function compute_amd_permutation(sparsity::SparseMatrixCSC{Bool})
    # AMD expects symmetric matrix for best results
    # Make symmetric if not already
    pattern = sparsity .| sparsity'

    # Compute AMD ordering
    perm = amd(pattern)

    return perm
end
```

## Expected Benefits

1. **Reduced fill-in**: AMD minimizes new nonzeros during factorization
2. **Faster factorization**: Less fill-in → faster numerical factorization
3. **Lower memory**: Smaller factorized matrix
4. **Skip AMD in KLU**: Provide pre-computed ordering, save ~10-20ms per factorization

### Estimated Performance Impact

For large systems (70k bus):
- **Pattern construction**: ~1ms (lightweight connectivity analysis)
- **AMD ordering**: ~50-100ms (one-time cost during setup)
- **Saved per iteration**: ~10-20ms (skip KLU's internal AMD)
- **Net benefit**: Significant for multi-iteration solves

For small systems (14 bus):
- Negligible overhead, small benefit
- Main value is in algorithm correctness and large-system scalability

## Usage API

```julia
# Option 1: Automatic (default behavior)
sys = SystemModel(am, models; reorder=true)  # Uses AMD automatically

# Option 2: Manual control
pattern = build_topology_sparsity(models)
perm = compute_amd_permutation(pattern)
sys = SystemModel(am, models; permutation=perm)

# Option 3: Custom permutation
perm_custom = custom_ordering(models)  # User-defined
sys = SystemModel(am, models; permutation=perm_custom)
```

## Testing Strategy

1. **Correctness**: Verify reordered system produces same power flow solution
2. **Sparsity accuracy**: Compare topology-based pattern vs. current address-based pattern
3. **Performance**: Benchmark factorization time with/without reordering
4. **KLU integration**: Verify permutation is correctly provided to KLU

## Future Extensions

1. **Island-aware ordering**: Respect electrical islands in permutation
2. **Geographical ordering**: Use physical bus locations
3. **Hierarchical ordering**: Partition by voltage levels (transmission/distribution)
4. **Custom constraints**: Allow users to fix certain variables at specific addresses
