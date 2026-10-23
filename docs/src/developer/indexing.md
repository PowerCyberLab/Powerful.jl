# The Three-Layer Addressing System

Maps between three coordinate systems so that computation is efficient while users keep
working with the IDs from their data files.

## Why Three Layers?

**Problem 1: External IDs are sparse**
- Bus numbers: 4001, 4502, 4211 (not 1, 2, 3)
- Can't be used as array indices
- Need dense indices: 1, 2, 3...

**Problem 2: Array access needs contiguous indices**
- Model data is stored in StructArrays (one array per field)
- Solution: Internal indices (1, 2, 3...) into those arrays

**Problem 3: Solvers need flat vectors**
- All angles in positions 1-14
- All voltages in positions 15-28
- Or a reordered layout for better sparse factorization
- Solution: DAE addresses specify the solver memory layout

**Each layer solves exactly one problem.**

## Layer 1: External Keys → Internal Indices

**Mapper:** `IndexOrder`, built by `ext2int` (`PowerfulCore/src/dae/ext2int.jl`)

**Keys are tuples** built from the fields named by each model's `get_key_name()` trait:
- `Bus`: `(i,)` - e.g. `(101,)`
- `PQ`: `(i, id)` - e.g. `(101, "1")`
- `Branch`: `(i, j, ckt)` - e.g. `(101, 205, "1")`

**Example:**
- Bus numbers [101, 205, 103] from file
- Map to internal [1, 2, 3]
- Lookup: `order.key_indices[:Bus][(101,)] → 1`

**When created:** SystemModel construction (once)

**Storage:**
```julia
sys.order.key_indices[:Bus]   # Dict: key tuple → internal index
sys.order.key_indices[:PQ]    # Same for each model type
```

The stored keys use the field types of the model (`Int32` bus numbers, `String3` IDs).
To look up a key given as plain Julia values, use `get_index_by_key`, which converts it
with `normalize_key` and errors if the instance does not exist:

```julia
using PowerfulCore: get_index_by_key, get_field_by_key

idx = get_index_by_key(sys, :PQ, (101, "1"))      # internal index of load (101, "1")
pl = get_field_by_key(sys, :PQ, (101, "1"), :pl)  # its pl field (per-unit)
```

Models store external IDs; only the index maps know internal positions:
```julia
sys.models.PQ.i[3]                          # external bus number of load 3
bus_idx = get_index_by_key(sys, :Bus, sys.models.PQ.i[3])
sys.models.Bus.basekv[bus_idx]              # base kV of that bus
```

**Handles:** Sparse, non-contiguous external IDs → dense, contiguous internal indices

## Layer 2: Internal Indices → DAE Addresses

**Mapper:** `AddressManager`

**Example (14-bus system with PowerInjection formulation):**
- 28 algebraic variables total (14 θ + 14 V)
- Natural addressing: θ at 1-14, V at 15-28
- Bus 1 (internal index 1): θ at address 1, V at address 15
- Bus 2 (internal index 2): θ at address 2, V at address 16

**When created:** Problem construction (`ACPowerFlowProblem`, `TDSProblem`), per formulation

**Storage:**
```julia
using PowerfulCore: AlgebVar, get_addresses

am = prob.address_manager
am.assignments[:Bus].by_type[AlgebVar()][:v]   # addresses of all bus voltages
get_addresses(am, :Bus, :v)                     # same, via the accessor
get_addresses(am, AlgebVar(), :Bus, :v)         # typed lookup
```

Element i of the address vector is the DAE address of instance i's variable.

**Why formulation-specific:**
- PowerInjection: θ, V variables
- CurrentInjection: V_r, V_i variables
- Same system, different variable sets → different addresses

**Addressing strategies** (the `addressing=` keyword of `ACPowerFlowProblem` and `TDSProblem`):
- `NaturalAddressing()` - Group by variable (all θ, then all V); the default
- `AMDAddressing()` - Approximate minimum degree reordering to reduce sparse fill-in
- `CustomPermutation(perm)` - User-provided ordering

## Layer 3: Foreign References (Precomputed Links)

**Problem:** Models reference other models repeatedly

**Example:**
- PV generator at bus 4001
- Needs the bus's data and variables during setup and evaluation
- Full path: External 4001 → Internal index → Addresses

**Solution:** Precompute the external → internal hop for every reference field.

**Storage:** `sys.order.foreign_refs[model][ref_key]` is a `GroupLinkRef` with two vectors:
```julia
ref = sys.order.foreign_refs[:PQ][:i]
ref.target_type   # [:Bus, :Bus, ...]  model referenced by each PQ instance
ref.target_idx    # internal index in that model, one entry per PQ instance
```

A `GroupLinkRef` can be heterogeneous: a governor that controls "the generator at
(BUS, ID)" may point to a `GENCLS` for one instance and a `GENROU` for another.
`is_heterogeneous(ref)` tells the two cases apart. Convenience accessors:

```julia
using PowerfulCore: get_foreign_indices, get_foreign_link_ref

bus_idx = get_foreign_indices(sys, :PV, :i)      # internal Bus indices of all PV generators
from_bus = get_foreign_indices(sys, :Branch, :i)
to_bus = get_foreign_indices(sys, :Branch, :j)
ref = get_foreign_link_ref(sys, :TGOV1, (:BUS, :ID))  # full GroupLinkRef
```

**Created:** SystemModel construction, from the models' `Terminal` declarations
(`build_terminal_refs`, called by `ext2int`)

**Benefit:** Eliminates dictionary lookups in hot paths

### Controller Links (TDS Only)

Controllers (TGOV1, SEXS) use the same `foreign_refs` map through a terminal that
references the controlled group, e.g. `Terminal(:generators, [:BUS, :ID])`. The
`GroupLinkRef` stored under the composite key `(:BUS, :ID)` can point to different
generator models.

**Group key names vs. arity:** Group interfaces declare `key_fields` (e.g. `[:BUS, :ID]`
for `:generators`). If a model implementing the group uses different key field names,
a warning is emitted and keys are matched by position. Prefer matching names when
adding new models.

### Pipeline: Where `foreign_refs` Are Built and Used

- **ext2int:** Builds `foreign_refs` from terminal declarations; entries are `GroupLinkRef` keyed by the terminal's reference field.
- **Derived-system builders:** Rebuild the index for the derived system, so dynamic model terminals and controller terminals get their own entries.
- **External parameter retrieval:** `retrieve_external_params` uses `get_group_ref(order, model, key)` to load `ExtParam` data.
- **Address retrieval:** `setup_with_formulation!` passes `foreign_refs` to `retrieve_addresses` to resolve foreign variables and residuals.
- **Initialization:** `InitInput`s of kind `:from_group_init` (e.g. TGOV1's `tm0` from the generator's `P_m`) are resolved through the controller's group link.

### Dynamic Model Linking Flow

Dynamic models use different linking mechanisms:

| Phase | Models | Traits Used | Output |
|-------|--------|-------------|--------|
| TDS system preparation | GENCLS, GENROU, ... | `get_replacement_map`, `get_data_inheritance` | `ModelLinks` |
| Index rebuild (ext2int) | TGOV1, SEXS | `get_controlled_group`, `get_terminals` | `GroupLinkRef` |

**Replacement models** (generators) replace static PV/Slack devices. `ModelLinks`
records which static instance each one replaces, and `get_data_inheritance` declares
which data to copy (e.g. GENCLS takes `xd_p` from `PV.zx`).

**Controllers** are not linked by replacement. Their terminal on the controlled group
is resolved by composite key matching in `build_terminal_refs`.

## API Boundaries

### Layer 1 API (IndexOrder)

**Input:** Model collections with external IDs in struct fields

**Output:** `order.key_indices[model][key_tuple] → internal_index`

**Used by:** System for array access, users for specifying devices

### Layer 2 API (AddressManager)

**Input:** SystemModel + formulation trait

**Output:** `am.assignments[model].by_type[addr_type][var_name] → Vector{UInt}`
- Element i = DAE address for instance i's variable

**Used by:** Gather/scatter operations during residual evaluation

### Layer 3 API (Foreign References)

**Input:** Model `Terminal` declarations

**Output:** `order.foreign_refs[model][ref_key] → GroupLinkRef`
- Target model and internal index for each referencing instance

**Used by:** External parameter retrieval, address retrieval and initialization

## Concrete Example: Variable Access

`get_var(sol, :Bus, :v, (14,))` performs these steps on the IEEE 14-bus case:

```@example indexing
using Powerful, PowerfulCases
using PowerfulCore: AlgebVar, get_index_by_key

sys = SystemModel(parse_network(load("ieee14").raw))
sol = solve(ACPowerFlowProblem(sys), Newton())

# Step 1: External key → internal index
idx = get_index_by_key(sys, :Bus, (14,))

# Step 2: Internal index → DAE address
v_addresses = sol.prob.address_manager.assignments[:Bus].by_type[AlgebVar()][:v]
addr = v_addresses[idx]

# Step 3: Read from the DAE vector (per-unit)
value = sol.dae.y[addr]
(idx, Int(addr), value, get_var(sol, :Bus, :v, (14,)))
```

**Step 4 (optional): Convert to physical units.** Bus voltage is a variable, so its
physical value uses the bus base voltage: `value * sys.models.Bus.basekv[idx]` kV.
Reports do this conversion through the `ReportConfig` passed to
`collect_model_table`. Stored per-unit coefficients (`sys.pu_coeffs`) apply to model
parameters and satisfy `pu = physical × coefficient`; see
[Per-Unit & Multi-Island](miscellaneous.md).

## Changing Address Layout

**Flexibility:** Only the AddressManager changes; the external API stays the same

**Example strategies:**
- **Natural:** All θ, then all V
- **AMD:** Reduce sparse fill-in
- **Custom:** Any permutation, e.g. θ,V for bus 1, then θ,V for bus 2

**Implementation:**
```julia
prob1 = ACPowerFlowProblem(sys, PowerInjection(); addressing = NaturalAddressing())
prob2 = ACPowerFlowProblem(sys, PowerInjection(); addressing = AMDAddressing())

# Same system, same external IDs
# Different DAE address layouts
# No model code changes
```

**Why it works:**
- The external ID mapping is unchanged
- Models don't see DAE addresses directly
- `ResidualBuffers` handle gather/scatter

**This is why the AddressManager lives in the problem, not in SystemModel.**

## Extension Points

### Custom Addressing Strategies

`CustomPermutation(perm)` accepts any permutation, where `perm[i]` is the new address
for old address `i`. Strategies are applied in `apply_addressing_strategy!`
(`PowerfulPF/src/solvers/power_flow_api.jl`), which computes the permutation and calls
`remap_addresses!`; a strategy that computes its own ordering needs an
`AddressingStrategy` subtype and a branch there.

### New Models with Cross-References

Declare references through terminals and use them in external parameters:
```julia
const MY_MODEL_TERMINALS = (Terminal(:Bus, :i),)
get_terminals(::Type{<:MyModel}) = MY_MODEL_TERMINALS

const MY_MODEL_EXT_PARAMS = [
    ExtParam(:basekv, MY_MODEL_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
]
```

The system creates `foreign_refs[:MyModel][:i]` during construction.

### Alternative Index Mapping

`ext2int` assigns internal indices in the order instances appear in the model arrays.
Another mapping (sorted by bus number, grouped by area) must still produce an
`IndexOrder` with:
- `key_indices`: key tuple → internal index for every model
- `foreign_refs`: `GroupLinkRef`s for every terminal

## Key Principle

**Maintain clear separation:**
- **External IDs** - User-facing, stable across solutions
- **Internal indices** - Implementation-efficient array positions
- **DAE addresses** - Solver-specific equation vector positions

Each layer has one job and provides a stable interface to the layer above.
