# The Three-Layer Addressing System

Powerful.jl maps between three coordinate systems: the IDs of the data file, the array positions of the model collections, and the positions in the solver vectors. Users work with the first, models are stored in the second, and solvers operate on the third.

## Overview of the Layers

External IDs are sparse. Bus numbers such as 4001, 4502 and 4211 cannot serve as array indices, so every model collection is addressed through dense internal indices 1, 2, 3, ... into its StructArray fields. Solvers in turn need flat vectors, for example all angles in positions 1 to 14 followed by all voltages in positions 15 to 28, or a reordered layout that reduces sparse fill-in. The DAE addresses define this layout.

Each layer has one mapper:

| Layer | Mapper | Created | Maps |
|-------|--------|---------|------|
| 1 | `IndexOrder` (built by `ext2int`) | `SystemModel` construction | External key to internal index |
| 2 | `AddressManager` | Problem construction, per formulation | Internal index to DAE address |
| 3 | Foreign references (`foreign_refs`) | `SystemModel` construction | Referencing instance to the internal index of the referenced instance |

## Layer 1: External Keys to Internal Indices

`IndexOrder` is built by `ext2int` (`PowerfulCore/src/dae/ext2int.jl`) once, when the `SystemModel` is constructed. Keys are tuples built from the fields named by the `get_key_name()` trait of each model:

- `Bus`: `(i,)`, for example `(101,)`
- `PQ`: `(i, id)`, for example `(101, "1")`
- `Branch`: `(i, j, ckt)`, for example `(101, 205, "1")`

For bus numbers `[101, 205, 103]` in the file, the internal indices are `[1, 2, 3]`, and `order.key_indices[:Bus][(101,)]` returns `1`. The mapping of every model type is stored in `sys.order`:

```julia
sys.order.key_indices[:Bus]   # Dict: key tuple → internal index
sys.order.key_indices[:PQ]    # Same for each model type
```

The stored keys use the field types of the model (`Int32` bus numbers, `String3` IDs). The function `get_index_by_key` converts a key given as plain Julia values with `normalize_key` and raises an error if the instance does not exist:

```julia
using PowerfulCore: get_index_by_key, get_field_by_key

idx = get_index_by_key(sys, :PQ, (101, "1"))      # internal index of load (101, "1")
pl = get_field_by_key(sys, :PQ, (101, "1"), :pl)  # its pl field (per-unit)
```

Models store external IDs in their fields, and the index maps provide the internal positions:

```julia
sys.models.PQ.i[3]                          # external bus number of load 3
bus_idx = get_index_by_key(sys, :Bus, sys.models.PQ.i[3])
sys.models.Bus.basekv[bus_idx]              # base kV of that bus
```

## Layer 2: Internal Indices to DAE Addresses

The `AddressManager` is created with the problem (`ACPowerFlowProblem`, `TDSProblem`) and assigns a DAE address to every variable and residual of every instance. For the 14-bus system with the `PowerInjection` formulation, there are 28 algebraic variables (14 θ and 14 V). Natural addressing places θ at addresses 1 to 14 and V at 15 to 28, so bus 1 (internal index 1) has θ at address 1 and V at address 15, and bus 2 has θ at address 2 and V at address 16.

```julia
using PowerfulCore: AlgebVar, get_addresses

am = prob.address_manager
am.assignments[:Bus].by_type[AlgebVar()][:v]   # addresses of all bus voltages
get_addresses(am, :Bus, :v)                     # same, via the accessor
get_addresses(am, AlgebVar(), :Bus, :v)         # typed lookup
```

Element `i` of an address vector is the DAE address of the variable of instance `i`.

The `AddressManager` belongs to the problem, because each formulation has its own variables: `PowerInjection` uses θ and V, and `CurrentInjection` uses V_r and V_i. The same `SystemModel` therefore yields different addresses for different problems.

The `addressing=` keyword of `ACPowerFlowProblem` and `TDSProblem` selects the layout:

- `NaturalAddressing()` groups addresses by variable (all θ, then all V) and is the default.
- `AMDAddressing()` applies an approximate minimum degree reordering to reduce sparse fill-in.
- `CustomPermutation(perm)` applies a user-provided ordering, where `perm[i]` is the new address of old address `i`. An example is θ and V of bus 1, followed by θ and V of bus 2.

```julia
prob1 = ACPowerFlowProblem(sys, PowerInjection(); addressing = NaturalAddressing())
prob2 = ACPowerFlowProblem(sys, PowerInjection(); addressing = AMDAddressing())
```

Both problems share the system and the external IDs and differ only in the DAE address layout. Model code is unchanged, because models access DAE addresses through `ResidualBuffers`, which gather and scatter by address.

## Layer 3: Foreign References

A model that refers to another model, such as a PV generator at bus 4001, needs the data and variables of the referenced instance during setup and evaluation. The path from the external ID 4001 to the internal index to the DAE addresses is resolved once at construction. The result is stored in `sys.order.foreign_refs[model][ref_key]`, a `GroupLinkRef` with two vectors, which avoids dictionary lookups in hot paths:

```julia
ref = sys.order.foreign_refs[:PQ][:i]
ref.target_type   # [:Bus, :Bus, ...]  model referenced by each PQ instance
ref.target_idx    # internal index in that model, one entry per PQ instance
```

The references are built from the `Terminal` declarations of the models by `build_terminal_refs`, which `ext2int` calls. A `GroupLinkRef` may refer to several model types: a governor that controls "the generator at (BUS, ID)" may point to a `GENCLS` for one instance and a `GENROU` for another. The function `is_heterogeneous(ref)` distinguishes the two cases. Convenience accessors return the indices directly:

```julia
using PowerfulCore: get_foreign_indices, get_foreign_link_ref

bus_idx = get_foreign_indices(sys, :PV, :i)      # internal Bus indices of all PV generators
from_bus = get_foreign_indices(sys, :Branch, :i)
to_bus = get_foreign_indices(sys, :Branch, :j)
ref = get_foreign_link_ref(sys, :TGOV1, (:BUS, :ID))  # full GroupLinkRef
```

### Controller links

Controllers (TGOV1, SEXS) use the same `foreign_refs` map through a terminal that references the controlled group, for example `Terminal(:generators, [:BUS, :ID])`. The `GroupLinkRef` stored under the composite key `(:BUS, :ID)` can point to different generator models. Group interfaces declare `key_fields` (for example `[:BUS, :ID]` for `:generators`). If a model implementing the group uses different key field names, a warning is emitted and the keys are matched by position. New models should use the declared names.

### Where the references are used

- `ext2int` builds `foreign_refs` from the terminal declarations. Each entry is a `GroupLinkRef` keyed by the reference field of the terminal.
- The builders of derived systems rebuild the index, so that dynamic model terminals and controller terminals receive their own entries.
- `retrieve_external_params` loads the data of `ExtParam` declarations through `get_group_ref(order, model, key)`.
- `setup_with_formulation!` passes `foreign_refs` to `retrieve_addresses`, which resolves foreign variables and residuals.
- `InitInput`s of kind `:from_group_init` (for example the `tm0` of TGOV1, taken from the `P_m` of the generator) are resolved through the group link of the controller.

### Linking of dynamic models

Dynamic generators and controllers are linked by two mechanisms:

| Phase | Models | Traits used | Output |
|-------|--------|-------------|--------|
| TDS system preparation | GENCLS, GENROU, ... | `get_replacement_map`, `get_data_inheritance` | `ModelLinks` |
| Index rebuild (`ext2int`) | TGOV1, SEXS | `get_controlled_group`, `get_terminals` | `GroupLinkRef` |

Generator models replace static PV and Slack devices. `ModelLinks` records which static instance each generator replaces, and `get_data_inheritance` declares the data to copy (for example, GENCLS takes `xd_p` from `PV.zx`). Controllers are linked by the composite-key match of their terminal on the controlled group, performed in `build_terminal_refs`.

## Example: Variable Access

The call `get_var(sol, :Bus, :v, (14,))` performs the following steps on the IEEE 14-bus case:

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

A fourth, optional step converts the value to physical units. Bus voltage is a variable, so its physical value is `value * sys.models.Bus.basekv[idx]` kV. Reports apply this conversion through the `ReportConfig` passed to `collect_model_table`. The stored per-unit coefficients (`sys.pu_coeffs`) apply to model parameters and satisfy `pu = physical × coefficient`; see [Per-Unit & Multi-Island](miscellaneous.md).

## Extension Points

### Custom addressing strategies

`CustomPermutation(perm)` accepts any permutation. Strategies are applied in `apply_addressing_strategy!` (`PowerfulPF/src/solvers/power_flow_api.jl`), which computes the permutation and calls `remap_addresses!`. A strategy that computes its own ordering requires an `AddressingStrategy` subtype and a branch in that function.

### New models with cross-references

A model declares references through terminals and uses them in external parameters:

```julia
const MY_MODEL_TERMINALS = (Terminal(:Bus, :i),)
get_terminals(::Type{<:MyModel}) = MY_MODEL_TERMINALS

const MY_MODEL_EXT_PARAMS = [
    ExtParam(:basekv, MY_MODEL_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
]
```

The system creates `foreign_refs[:MyModel][:i]` during construction.

### Alternative index mappings

`ext2int` assigns internal indices in the order in which instances appear in the model arrays. An alternative mapping, such as one sorted by bus number or grouped by area, must produce an `IndexOrder` that contains `key_indices` (key tuple to internal index for every model) and `foreign_refs` (a `GroupLinkRef` for every terminal).
