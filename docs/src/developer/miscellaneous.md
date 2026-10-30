# Per-Unit Systems and Multi-Island Handling

This page describes two parts of system construction: the conversion of input data to per-unit values, and the detection and splitting of electrical islands.

## Per-Unit Conversion

Solvers operate on dimensionless per-unit values, while users specify data in physical units such as MW and kV. Powerful.jl converts the data on input, stores the conversion coefficients, and converts back on output.

Per-unit values give better conditioned equations: a 500 MW load becomes 5.0 pu and is of the same order as a 1.05 pu voltage. A per-unit impedance has the same meaning at every voltage level, so transformer impedances are comparable across the network. Reports can then be produced in physical units on request.

### Storage

Three fields of a `SystemModel` hold the data of the conversion:

- `sys.models` holds the per-unit values used by the residual calculations. The solvers operate on these values only.
- `sys.parsed_data` holds the original physical values as read from the file. It is kept by the `SystemModel(case)` constructor for provenance and debugging.
- `sys.pu_coeffs` holds the conversion coefficients: a `NamedTuple` of `Dict{Symbol, Vector{Float64}}` with one per-instance vector for each converted parameter. The relation is `pu = physical × coefficient`, so `physical = pu / coefficient` recovers the original value exactly.

On the IEEE 14-bus case, load 3 has 0.478 pu in `sys.models.PQ.pl[3]`, 47.8 MW in `sys.parsed_data.models.PQ.pl[3]` and the coefficient 0.01 (= 1 / 100 MVA) in `sys.pu_coeffs.PQ[:pl][3]`:

```@example perunit
using Powerful, PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))
(pu = sys.models.PQ.pl[3],
 physical = sys.parsed_data.models.PQ.pl[3],
 coeff = sys.pu_coeffs.PQ[:pl][3])
```

Coefficients are stored per instance, because generators have different MVA ratings and buses have different base voltages. Structural fields (IDs, status) are not converted.

### Declaring the Conversion

A model declares the fields to convert in its metadata, with the `PUType` of each field:

```julia
const PQ_PARAMS = [
    ModelParam(:pl, PowerNominal, description="Constant power P", units="MW"),
    ModelParam(:ql, PowerNominal, description="Constant power Q", units="MVAR"),
    ModelParam(:yp, PowerNominal, description="Constant admittance P", units="MW"),
    # ...
]
```

The `PUType` determines the coefficient (`calc_pu_coefficients` in `PowerfulCore/src/system/per_unit.jl`). In the table, Sb is the system base MVA, Vb the bus base kV, and Sn, Vn and Zn the rating of the device.

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

For a 47.8 MW load on a 100 MVA system base, the coefficient is 1/100 = 0.01 and the per-unit value is 0.478 pu.

### External Base Quantities

Some conversions need data from other models. A bus base voltage is declared as an `ExtParam` with the reference type `BaseVoltage()`:

```julia
const PQ_EXT_PARAMS = [
    ExtParam(:basekv, PQ_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
]
```

External parameters are retrieved in one of two phases. `PrePerUnit` parameters are retrieved before the conversion (dimensionless values and bases), and `PostPerUnit` parameters after it, when the referenced values are already in per-unit. A model that converts a quantity with a bus base voltage must declare a `BaseVoltage` external parameter; otherwise the conversion raises an error. For example, with a 100 MVA system base and a 13.8 kV bus, Zb = 13.8² / 100 = 1.9044 Ω, and an impedance in ohms is multiplied by 1 / Zb.

A device rating is declared with `ModelRating`:

```julia
const PV_RATINGS = [
    ModelRating(:mbase, NominalPower, description="Generator MVA base"),
]
```

A model whose parameters are per unit on the base of a device it replaces takes the rating from that device with `from = LinkedDevice()`. An example is a converter on the `mbase` of its generator record, which its own record does not carry. The rating is the one that the replaced device declares in the same field:

```julia
ModelRating(:mbase, NominalPower; from = LinkedDevice(), description="MVA base of the replaced generator")
```

`derive_system(sys, keep, new_models; links)` resolves the rating through the replacement links (`get_replacement_map`, `create_model_links`), converts the model parameters to the system base, and stores the coefficients in `sys.pu_coeffs`. A model with no linked device, or a device that declares no such rating, is refused with an error that names the model and the instance.

### Machines and Converters on Their Own Base

The parameters of a dynamic record are per unit on the MVA base of the machine, `MBASE` of the generator record that the record belongs to. This applies to the machine parameters (`H`, `D`, the reactances, the armature resistance, the time constants) and to the exciter and governor data. The network is on the system base `Sb`.

A model keeps the quantities of the device on the device base `Sn` and converts them where they meet the network. It declares the rating and the ratio in its metadata:

```julia
ModelMetadata(
    ...,
    ratings = [ModelRating(:mbase, NominalPower; from = LinkedDevice())],   # or OwnField(): ModelRating(:mbase, NominalPower)
    derived_params = [RatingRatio(:Mratio, NominalPower)],
)
```

`Mratio` is `Sn/Sb`, for a rating taken from the model's own field or from the linked device. It has no DYR entry and is kept with the per-unit coefficients (`sys.pu_coeffs[model][:Mratio]`, the coefficient of a `PowerPU` parameter). It is read as `external.Mratio[i]` in `calc_residuals!` and `calc_automaton_values!`, and as `gather_init_data(ctx, M).Mratio` in `initialize_dynamic_states`. A derived ratio without a `NominalPower` rating is an error.

A rating is read when the system is built. `check_events` therefore refuses a `FieldModification` of a field declared by a `ModelRating`, either the model's own field (the `mbase` of a generator) or the field of the linked device that a machine or converter names. A modification would leave `Mratio` and the per-unit coefficients out of date.

`Mratio` multiplies values on the device base: states and parameters declared `NoPU`. Power and current into the network are `Mratio` times the device-base values, and the system-base powers of the power flow are divided by `Mratio` at initialization. Parameters that Powerful converted to the system base by their `PUType` are used as converted.

The machines GENCLS, GENROU, GENSAL and GENTPJ1, with SEXS and TGOV1, follow this table:

| Quantity | Base | Who converts |
|---|---|---|
| parameters (`H`, `D`, reactances, `ra`, `Kis`; time constants, saturation) | machine | no conversion (declared on the machine base; `ra` is `ZSORCE`, machine base; GENSAL has no `Xq2`: X''q = X''d) |
| states, `Id`, `Iq`, `P_m`, `E_fd` | machine (`E_fd` in the machine's own field units) | no conversion |
| initialization inputs `p_gen`, `q_gen` (power flow) | system | the model: `/ Mratio` |
| `p_inj`, `q_inj`, and the observables `Pe`, `Qe` | system | the model: `Mratio *` the machine-base value |
| SEXS gains and times, `vref`, `E_fd` | dimensionless; field units of the machine | no conversion |
| TGOV1 `R`, `Vmax`, `Vmin`, `Dt`, `pref`, `pout` | machine (they are power of `P_m`) | no conversion |

The converters REGCA1 and REGFM_A1 declare a linked rating and the ratio and follow this table:

| Quantity | Base | Who converts |
|---|---|---|
| states, `Pe`, `Qe`, currents `ip`, `iq`, and parameters the model multiplies by `Mratio` | device | no conversion (`Pe` stays on the device base, for a controller that reads it) |
| `p_inj`, `q_inj`, `ip_inj`, `iq_inj` (the `bus_group` injections) | system | the model: `Mratio *` the device-base value |
| initialization inputs from the power flow | system | the model: `/ Mratio` |
| parameters declared with a `PUType` (`PowerPU`, `ImpedancePU`, ...) | system after conversion | Powerful (`calc_pu_coefficients`) |

The same physical machine on the system base has `H_s = H Sn/Sb`, `D_s = D Sn/Sb` and `Z_s = Z Sb/Sn`, and a governor with `R_s = R Sb/Sn` and `P_s = P Sn/Sb`. The test `Powerful/test/tds/test_machine_base.jl` shows that the two descriptions give the same trajectories.

### Conversion Workflow

`SystemModel(case)` (`PowerfulCore/src/system/system_builder.jl`) performs the following steps:

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

### Inverse Conversion for Reports

`collect_model_table` returns a Tables.jl-compatible NamedTuple of columns. The units come from the `ReportConfig` (`UnitPreferences`): by default, power is reported in MW and MVAR, voltage in pu and angles in degrees.

```julia
bus_table = collect_model_table(sol, :Bus)                       # default ReportConfig
bus_table = collect_model_table(sol, :Bus; config = AllPhysical)  # kV, MW, MVAR
```

The conversion is implemented in `PowerfulPF/src/io/reports.jl`. Parameters use the stored coefficients, so that physical = pu / coefficient. Observations, which are computed values such as the bus voltage or the load consumption, have no stored coefficients and are converted by their `PUType`. A bus voltage of 1.03 pu on a 138 kV base, with `VoltageNominal`, is reported as 1.03 × 138 = 142.14 kV.

### Adding a Model with Per-Unit Conversion

A model declares its `params` with a `PUType` for each field to convert, its `ratings` if it has a device base, and an `ExtParam` for each external base that the `PUType`s require (typically `BaseVoltage()`). The conversion can be tested by checking that `pu / coefficient` reproduces the original value in `sys.parsed_data`.

## Multi-Island Topology

A network can separate into islands through outages or planned separations. Each island is solved independently.

### Island Detection

`analyze_topology(sys)` (`PowerfulCore/src/system/topology.jl`) builds a graph from the Terminal API and finds its connected components with a breadth-first search. The buses are the vertices. Multi-terminal devices (Branch, Transformer) are the edges, and only in-service devices are included. Each model declares its terminals through `get_terminals`. The function then identifies the slack vertex of every island. An island with no slack bus, or with several, produces a warning.

The result is a `TopologyInfo` with these fields:

- `n_islands`: the number of islands.
- `island_map`: the island of each bus.
- `slack_buses`: the slack bus indices per island.
- `adjacency`: the adjacency matrix.
- `components`: the bus lists per island.
- `incidence`, `is_bridge`, `branch_indices`: the branch-bus incidence and the bridge flags.

`analyze_topology!(sys)` caches the result in `sys.topology`, and `analyze_topology(sys)` always recomputes it.

### Island Splitting

`split_islands(sys; auto_analyze=true)` (`PowerfulCore/src/system/multi_island.jl`) extracts the islands. It uses `sys.topology` when cached and otherwise calls `analyze_topology(sys)` (with `auto_analyze=true`). A single island returns the original `SystemModel`, and several islands return a `Vector{SystemModel}`, one per island:

```julia
result = split_islands(sys)
```

Each island is an independent `SystemModel` with its own `IndexOrder`, its own model collections and its own address space (assigned when an `ACPowerFlowProblem` is created). It has the same `settings` and array type as the parent. An island is built in five steps:

1. The buses of each island are taken from `island_map`.
2. The models are filtered by their terminals. A load at bus 7 belongs to the island that contains bus 7, and branch 3-4 is kept only if both buses are in the island.
3. The index mappings are rebuilt. If island 2 has the buses 5, 7 and 11 of the original system, the external IDs remain 5, 7 and 11, and the internal indices are 1, 2 and 3.
4. The external parameters are retrieved again. Each island runs `ext2int`, which creates new `foreign_refs`, followed by `retrieve_external_params`. The external parameters of the parent cannot be reused, because the internal indices have changed.
5. The per-unit coefficients and the parsed data are sliced to the devices of the island. A parent with 50 loads has 50 coefficient entries per load parameter, and an island with 12 loads has 12. Coefficient values are unchanged, so load 3 with 0.478 pu and coefficient 0.01 has the same 47.8 MW in the parent and in the island.

Each island therefore converts its reports correctly with `collect_model_table`, which uses the `pu_coeffs` of the island.

### Solving Multi-Island Systems

`ACPowerFlowProblem` detects islands automatically (`auto_split=true`), and `solve` solves each island and returns a `MultiIslandSolution`. The function `get_var` routes a request to the island that contains the device:

```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())
v = get_var(sol, :Bus, :v, (4001,))
```

Islands can also be solved manually:

```julia
result = split_islands(sys)
if result isa Vector
    solutions = [solve(ACPowerFlowProblem(island), Newton()) for island in result]
else
    solution = solve(ACPowerFlowProblem(result), Newton())
end
```

### Outages and Topology

An outage is applied by setting the status of a device, after which the topology is recomputed:

```julia
# Take a branch out of service
sys.models.Branch.status[42] = false

# Recompute connectivity with the outage in effect
topo = analyze_topology(sys)
if topo.n_islands > 1
    @warn "Outage creates $(topo.n_islands) islands"
end
```

`split_islands(sys)` also sees the outage when no topology is cached in `sys.topology`. Systematic contingency studies use `PowerfulCA`, which applies outages through the status dict of the solver and leaves the system unchanged.

### Writing Models for Island Splitting

A model works with topology analysis and island filtering when it follows these requirements:

1. It declares its terminals through `get_terminals`.
2. It stores external IDs and leaves the mapping to `IndexOrder`, so that filtering preserves its references.
3. It behaves the same on connected and disconnected systems, where only the instance counts differ.
4. Cases with a missing slack bus are reported clearly. The topology analysis warns about islands without a slack bus, which avoids a solver failure with no explanation.
