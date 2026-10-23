# Per-Unit Systems and Multi-Island Handling

Two topics every real-world study meets: unit conversion and network topology.

## Per-Unit Conversion System

**Problem:** Solvers use dimensionless per-unit values. Users think in physical units (MW, kV).

**Solution:** Convert on input, store the coefficients, convert back on output.

### Why Per-Unit?

**1. Better conditioned equations**
- Physical units: 500 MW next to a 1.05 pu voltage (orders of magnitude apart)
- Per-unit: 5.0 pu power next to 1.05 pu voltage (similar scales)

**2. Portable parameters**
- A per-unit impedance means the same thing at any voltage level
- Transformer impedances are comparable across voltage levels

**3. User-friendly reports**
- Compute in per-unit
- Report in physical units when requested

### Three-Part Storage

**1. `sys.models` - Per-unit values (for computation)**
```julia
sys.models.PQ.pl[3]  # 0.478 pu (on the 100 MVA system base)
```
- Used by residual calculations
- The solvers never see physical units

**2. `sys.parsed_data` - Original physical values**
```julia
sys.parsed_data.models.PQ.pl[3]  # 47.8 MW (from the file)
```
- Data provenance and debugging
- Kept by the `SystemModel(case)` constructor

**3. `sys.pu_coeffs` - Conversion coefficients**
```julia
sys.pu_coeffs.PQ[:pl][3]  # 0.01 = 1 / (100 MVA)
```
- A `NamedTuple` of `Dict{Symbol, Vector{Float64}}`, one per-instance vector per converted parameter
- Enables exact inverse conversion
- **pu = physical × coefficient**, so **physical = pu / coefficient**

On the IEEE 14-bus case:

```@example perunit
using Powerful, PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))
(pu = sys.models.PQ.pl[3],
 physical = sys.parsed_data.models.PQ.pl[3],
 coeff = sys.pu_coeffs.PQ[:pl][3])
```

### Metadata-Driven Conversion

**Models declare which fields need conversion:**
```julia
const PQ_PARAMS = [
    ModelParam(:pl, PowerNominal, description="Constant power P", units="MW"),
    ModelParam(:ql, PowerNominal, description="Constant power Q", units="MVAR"),
    ModelParam(:yp, PowerNominal, description="Constant admittance P", units="MW"),
    # ...
]
```

**The PUType determines the coefficient** (`calc_pu_coefficients` in
`PowerfulCore/src/system/per_unit.jl`; Sb = system base MVA, Vb = bus base kV,
Sn/Vn = device ratings):

| PUType | Input | Coefficient |
|--------|-------|-------------|
| `PowerNominal` | MW, MVAR | 1 / Sb |
| `VoltageNominal` | kV | 1 / Vb |
| `CurrentNominal` | A | Vb / Sb |
| `ImpedanceNominal` | Ω | Sb / Vb² (= 1/Zb) |
| `AdmittanceNominal` | S | Vb² / Sb (= Zb) |
| `AngleDegrees` | degrees | π / 180 |
| `PowerPU` | pu on device base | Sn / Sb |
| `VoltagePU` | pu on device base | Vn / Vb |
| `ImpedancePU` | pu on device base | Zn / Zb |
| `AdmittancePU` | pu on device base | Zb / Zn |
| `NoPU` | already per-unit | 1 |

**Example:**
```julia
# System base: 100 MVA
# Physical: 47.8 MW
# Coefficient: 1/100 = 0.01
# Per-unit: 47.8 × 0.01 = 0.478 pu
```

### External Base Quantities

**Some conversions need data from other models:**

```julia
const PQ_EXT_PARAMS = [
    ExtParam(:basekv, PQ_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
]
```

An `ExtParam` with reference type `BaseVoltage()` supplies Vb; a `ModelRating` supplies
a device rating:
```julia
const PV_RATINGS = [
    ModelRating(:mbase, NominalPower, description="Generator MVA base"),
]
```

A model whose parameters are per unit on the base of a device it replaces (a converter on
the `mbase` of its generator record, which its own record does not carry) takes the rating
from that device with `from = LinkedDevice()`. The rating is the one the replaced device
declares in the same field:
```julia
ModelRating(:mbase, NominalPower; from = LinkedDevice(), description="MVA base of the replaced generator")
```
`derive_system(sys, keep, new_models; links)` resolves it through the replacement links
(`get_replacement_map`, `create_model_links`), converts the model's parameters to the system
base, and stores the coefficients in `sys.pu_coeffs`. A model with no linked device, or a
device that declares no such rating, is refused, naming the model and instance.

### Machines and converters on their own base

**DYR data is per unit on the machine base.** A PSS/E dynamic record gives a machine's
parameters (`H`, `D`, the reactances, the armature resistance, the time constants), and the
exciter and governor data, per unit on the machine's own MVA base, `MBASE` of the generator
record it belongs to. The network is on the system base `Sb`.

A model keeps what belongs to the device on the device base `Sn` and converts where it meets
the network. It declares the rating and the ratio:
```julia
ModelMetadata(
    ...,
    ratings = [ModelRating(:mbase, NominalPower; from = LinkedDevice())],   # or OwnField(): ModelRating(:mbase, NominalPower)
    derived_params = [RatingRatio(:Mratio, NominalPower)],
)
```
`Mratio` is `Sn/Sb`, for a rating from the model's own field or from the linked device alike.
It is not a field and has no DYR entry. It is kept with the per-unit coefficients
(`sys.pu_coeffs[model][:Mratio]`, the coefficient of a `PowerPU` parameter) and read as
`external.Mratio[i]` in `calc_residuals!` and `calc_automaton_values!` (the external
parameters of the model), and as `gather_init_data(ctx, M).Mratio` in
`initialize_dynamic_states`. A derived ratio without a `NominalPower` rating is an error.

A rating is read when the system is built. `check_events` therefore refuses a
`FieldModification` of a field that a model's `ModelRating` declares, the model's own (a
generator's `mbase`) or the linked device's (the field a machine or converter names), because
`Mratio` and the per-unit coefficients would be stale.

Rule: a value multiplied by `Mratio` is on the device base (a state, or a parameter declared
`NoPU`), never a parameter Powerful converted to the system base. Power and current into the
network are `Mratio` times the device-base values; the system-base powers of the power flow
are divided by it at initialization.

**Machines** (GENCLS, GENROU, GENSAL, GENTPJ1, with SEXS and TGOV1):

| Quantity | Base | Who converts |
|---|---|---|
| parameters (`H`, `D`, reactances, `ra`, `Kis`; time constants, saturation) | machine | nothing (declared on the machine base; `ra` is `ZSORCE`, machine base; GENSAL has no `Xq2`: X''q = X''d) |
| states, `Id`, `Iq`, `P_m`, `E_fd` | machine (`E_fd` in the machine's own field units) | nothing |
| initialization inputs `p_gen`, `q_gen` (power flow) | system | the model: `/ Mratio` |
| `p_inj`, `q_inj`, and the observables `Pe`, `Qe` | system | the model: `Mratio *` the machine-base value |
| SEXS gains and times, `vref`, `E_fd` | dimensionless; field units of the machine | nothing |
| TGOV1 `R`, `Vmax`, `Vmin`, `Dt`, `pref`, `pout` | machine (they are power of `P_m`) | nothing |

**Converters** (a model that declares a linked rating and the ratio, as REGCA1 and REGFM_A1 do):

| Quantity | Base | Who converts |
|---|---|---|
| states, `Pe`, `Qe`, currents `ip`, `iq`, and parameters the model multiplies by `Mratio` | device | nothing (`Pe` stays on the device base, for a controller that reads it) |
| `p_inj`, `q_inj`, `ip_inj`, `iq_inj` (the `bus_group` injections) | system | the model: `Mratio *` the device-base value |
| initialization inputs from the power flow | system | the model: `/ Mratio` |
| parameters declared with a `PUType` (`PowerPU`, `ImpedancePU`, ...) | system after conversion | Powerful (`calc_pu_coefficients`) |

A parameter converted to the system base by its `PUType` is never multiplied by `Mratio`.

The same physical machine on the system base has `H_s = H Sn/Sb`, `D_s = D Sn/Sb`,
`Z_s = Z Sb/Sn`, and a governor with `R_s = R Sb/Sn`, `P_s = P Sn/Sb`
(`Powerful/test/tds/test_machine_base.jl` shows that the two give the same trajectories).

**Retrieval phases:**
- **PrePerUnit** - Retrieved before conversion (dimensionless values, bases)
- **PostPerUnit** - Retrieved after conversion (already in per-unit)

**Example:** Impedance conversion needs the bus base voltage
```julia
# System base: 100 MVA, bus base: 13.8 kV
# Zb = 13.8² / 100 = 1.9044 Ω
# pu = Ω × (1 / Zb)
```

### Conversion Workflow

`SystemModel(case)` (`PowerfulCore/src/system/system_builder.jl`) runs these steps:

```julia
models = load_models(case, all_models())

# 1. External keys → internal indices
models, order = ext2int(models)

# 2. External parameters needed before conversion (base voltages, bus types)
external_pre = retrieve_external_params(models, models, order, PrePerUnit)

# 3. Coefficients and conversion (per model: calc_pu_coefficients, convert_to_per_unit)
parsed = ParsedSystemData(models, system_base)
result = convert_system_to_per_unit(parsed.models, parsed.system_base, external_pre)

# 4. External parameters that must be in per-unit
external_post = retrieve_external_params(result.models, result.models, order, PostPerUnit)

# 5. Assemble
sys = SystemModel(result.models, merge_external(external_pre, external_post), order;
                  settings = SystemConfig(), parsed_data = parsed, pu_coeffs = result.coeffs)
```

Structural fields (IDs, status) are not converted.

**Result:**
- `sys.models` - Per-unit (computation)
- `sys.pu_coeffs` - Coefficients (reports)
- `sys.parsed_data` - Original values (provenance)

### Inverse Conversion for Reports

```julia
bus_table = collect_model_table(sol, :Bus)                       # default ReportConfig
bus_table = collect_model_table(sol, :Bus; config = AllPhysical)  # kV, MW, MVAR
```

`collect_model_table` returns a Tables.jl-compatible NamedTuple of columns. The units
come from the `ReportConfig` (`UnitPreferences`): by default power is reported in
MW/MVAR, voltage in pu and angles in degrees.

**Process** (`_convert_to_report_units` in `PowerfulPF/src/io/reports.jl`):
- **Parameters** use the stored coefficients: physical = pu / coefficient
- **Observations** (computed values such as bus voltage or load consumption) have no
  stored coefficients; they are converted by their `PUType`, e.g. `VoltageNominal`
  with the bus base voltage:

```julia
# Per-unit: V = 1.03
# Base: Vb = 138 kV
# Report: 1.03 × 138 = 142.14 kV
```

**Per-instance coefficients:** Generators have different MVA ratings and buses have
different base voltages, so coefficients are stored per instance.

## Multi-Island Topology

**Problem:** Networks can separate into islands (outages, planned separations). Each
island is solved independently.

### Island Detection

`analyze_topology(sys)` (`PowerfulCore/src/system/topology.jl`) uses breadth-first
search on the network graph:

**1. Build the graph from the Terminal API**
```julia
# Each model declares terminals through get_terminals
# Buses = vertices; multi-terminal devices (Branch, Transformer) = edges
# Only in-service devices are included
```

**2. Find connected components**

**3. Identify the slack vertex per island**
```julia
# Each island needs a slack bus
# A warning is issued for islands with none or several
```

**Result stored in `TopologyInfo`:**
- `n_islands` - Number of islands
- `island_map` - Island of each bus
- `slack_buses` - Slack bus indices per island
- `adjacency` - Adjacency matrix
- `components` - Bus lists per island
- `incidence`, `is_bridge`, `branch_indices` - Branch-bus incidence and bridge flags

`analyze_topology!(sys)` caches the result in `sys.topology`; `analyze_topology(sys)`
always recomputes it.

### Island Splitting

**`split_islands(sys; auto_analyze=true)`** extracts disconnected islands
(`PowerfulCore/src/system/multi_island.jl`). It uses `sys.topology` when cached;
otherwise it runs `analyze_topology(sys)` (with `auto_analyze=true`).

**Single island:**
```julia
result = split_islands(sys)
# Returns the original SystemModel unchanged
```

**Multiple islands:**
```julia
result = split_islands(sys)
# Returns Vector{SystemModel} (one per island)
```

**Each island is an independent SystemModel:**
- Own `IndexOrder`
- Own address space (when an `ACPowerFlowProblem` is created)
- Own model collections
- Same `settings` and array type as the parent

**Splitting process:**

**1. Identify the buses of each island** from `island_map`

**2. Filter models** using their terminals
```julia
# Load at bus 7 → the island containing bus 7
# Branch 3-4 → kept only if both buses are in the island
```

**3. Rebuild index mappings**
```julia
# Island 2 has buses [5, 7, 11] from the original
# External IDs: 5, 7, 11 (unchanged)
# Internal indices: 1, 2, 3 (new, per island)
```

**4. Re-retrieve external parameters**
```julia
# Each island runs ext2int (new foreign_refs), then retrieve_external_params;
# the parent's external parameters can't be reused because internal indices changed
```

**5. Filter per-unit coefficients and parsed data**
```julia
# The parent has entries for all devices; the island only for its devices
# Coefficient values are unchanged (the arrays are sliced)
```

### Solving Multi-Island Systems

**Automatic detection and handling:**
```julia
prob = ACPowerFlowProblem(sys)     # Detects islands (auto_split=true)
sol = solve(prob, Newton())        # Solves each island, returns a MultiIslandSolution
```

**Unified variable access:**
```julia
v = get_var(sol, :Bus, :v, (4001,))
# MultiIslandSolution routes to the correct island
```

**Manual splitting:**
```julia
result = split_islands(sys)
if result isa Vector
    solutions = [solve(ACPowerFlowProblem(island), Newton()) for island in result]
else
    solution = solve(ACPowerFlowProblem(result), Newton())
end
```

**The result type encodes the split:**
- `SystemModel` - Single connected island
- `Vector{SystemModel}` - Multiple islands

### Outages and Topology

```julia
# Take a branch out of service
sys.models.Branch.status[42] = false

# Recompute connectivity with the outage in effect
topo = analyze_topology(sys)
if topo.n_islands > 1
    @warn "Outage creates $(topo.n_islands) islands"
end
```

`split_islands(sys)` also sees the outage when no topology is cached in
`sys.topology`. For systematic contingency studies use `PowerfulCA`, which applies
outages through the solver's status dict without modifying the system.

### Island and Per-Unit Interaction

**Workflow:**
1. Parse the data file (physical units)
2. Create the `SystemModel` (conversion to per-unit on the system base)
3. Split islands if needed

**Each island inherits the parent's per-unit values and coefficients:**
```julia
# Parent: load 3 is 0.478 pu with coefficient 0.01
# Island: load 3 is 0.478 pu with coefficient 0.01
# Physical: 47.8 MW (same before and after)
```

**Island-specific coefficient storage:**
- A parent with 50 loads has 50 coefficient entries per load parameter
- An island with 12 loads has 12
- Indexing is local to the island

Each island's `collect_model_table` uses its own `pu_coeffs`, so reports convert correctly.

## Architecture Implications

### 1. Separation of Data from Computation
- `SystemModel` = data (per-unit values, coefficients, topology)
- `ACPowerFlowProblem` = computation specification
- Island splitting creates new SystemModels (data layer)

### 2. Metadata-Driven Design
- Per-unit conversion uses model metadata (`params`, `ratings`, `ext_params`)
- Island splitting uses terminal declarations (connectivity)
- No hardcoded model-specific knowledge

### 3. Composition Over Special Cases
- Multi-island = `Vector{SystemModel}` (not a special system type)
- Each island is a regular SystemModel
- Composes with existing problem/solver code

### 4. Derived External Parameters
- External parameters are derived from the models and index maps
- Splitting does not copy them; each island retrieves its own

## Implementation Notes

### Adding Models with Per-Unit Conversion

**1. Declare params metadata**
```julia
const MY_MODEL_PARAMS = [
    ModelParam(:power, PowerNominal, description="...", units="MW"),
    ModelParam(:voltage, VoltageNominal, description="...", units="kV"),
    ModelParam(:angle, AngleDegrees, description="...", units="deg"),
]
```

**2. Declare device ratings (if needed)**
```julia
const MY_MODEL_RATINGS = [
    ModelRating(:mbase, NominalPower, description="Device MVA base"),
]
```

**3. Declare external bases (if needed)**
```julia
const MY_MODEL_EXT_PARAMS = [
    ExtParam(:basekv, MY_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
]
```

A conversion that needs a bus base voltage raises an error if the model has no
`BaseVoltage` external parameter.

**4. Test the round trip**
```julia
# pu / coefficient must reproduce the original value in sys.parsed_data
```

### Working with Multi-Island Systems

**1. Use the Terminal API** - declare terminals through `get_terminals` so the model
works with topology analysis and island filtering.

**2. Don't store cross-model references** - store external IDs and let `IndexOrder`
handle the mapping, so filtering does not break references.

**3. Test both connected and disconnected systems** - the model should behave the same;
only instance counts differ.

**4. Handle missing slack buses** - the topology analysis warns about islands without a
slack bus; report such cases clearly instead of letting the solver fail to converge.

## Quick Reference

### Per-Unit Conversion
```julia
using Powerful, PowerfulCases

case = parse_network(load("ieee14").raw)
sys = SystemModel(case)           # Converts to pu, stores coefficients

sys.models.PQ.pl[3]               # 0.478 pu
sys.pu_coeffs.PQ[:pl][3]          # 0.01
sys.parsed_data.models.PQ.pl[3]   # 47.8 MW

sol = solve(ACPowerFlowProblem(sys), Newton())
bus_table = collect_model_table(sol, :Bus; config = AllPhysical)
```

### Multi-Island Handling
```julia
# Automatic detection
sol = solve(ACPowerFlowProblem(sys), Newton())

# Manual splitting
result = split_islands(sys)
if result isa Vector
    for island in result
        sol = solve(ACPowerFlowProblem(island), Newton())
    end
else
    sol = solve(ACPowerFlowProblem(result), Newton())
end

# Connectivity after an outage
sys.models.Branch.status[42] = false
topo = analyze_topology(sys)
```

### Common PUTypes
```julia
PowerNominal        # MW, MVAR → × 1/Sb
VoltageNominal      # kV → × 1/Vb
ImpedanceNominal    # Ω → × Sb/Vb²
AdmittanceNominal   # S → × Vb²/Sb
AngleDegrees        # degrees → radians
NoPU                # No conversion
```
