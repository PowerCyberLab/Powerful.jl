# SystemModel convenience constructors
# NOTE: These extend Foundation.SystemModel with convenience signatures

using PowerFlowData
using ..Foundation

"""
    SystemModel(case::PowerFlowData.Network, model_types=nothing;
                system_base=100.0,
                frequency_base=case.caseid.basfrq, constant_power_vmin=0.7)

Convenience constructor to create a SystemModel from a PowerFlowData case.

# Arguments
- `case`: PowerFlowData Network object (from `parse_network`)
- `model_types`: Optional tuple/vector of model types. If `nothing` (default), uses `all_models()`
- `system_base::Float64=100.0`: System base power in MVA (the models are converted to per-unit on this base)
- `frequency_base::Float64`: System frequency in Hz (default: the case's base frequency)
- `constant_power_vmin::Float64=0.7`: Voltage (pu) below which constant-power load is
  reduced smoothly to zero (see `constant_power_factor`)

The system-wide settings (`sys.settings`, a `SystemConfig`) hold the frequency base,
the MVA base (`system_base`) and the constant-power breakpoint; every analysis of the
system uses them.

DAE addresses are assigned per problem, not here: pass `addressing` (e.g.
`AMDAddressing()`) to `ACPowerFlowProblem` or `TDSProblem`.

# Workflows

## Default Workflow (with per-unit conversion)
```julia
using Powerful
using PowerfulCases

case = parse_network(load("ieee14").raw)
sys = SystemModel(case)  # Auto-converts to per-unit on 100 MVA base

# Access original data (preserved automatically)
sys.parsed_data.models.PQ.pl[1]  # → 21.7 MW (original)
sys.models.PQ.pl[1]              # → 0.217 pu (computational)

# Use different system base
sys = SystemModel(case; system_base=500.0)  # 500 MVA base

# Use AMD addressing for better sparse factorization
prob = ACPowerFlowProblem(sys; addressing=AMDAddressing())
```

## Explicit Workflow (Maximum Control)
```julia
# For advanced users who need full control over each step
case = parse_network(load("ieee14").raw)
models, order = ext2int(load_models(case, all_models()))
external_pre = retrieve_external_params(models, models, order, PrePerUnit)
parsed = ParsedSystemData(models, 100.0)
result = convert_system_to_per_unit(parsed.models, parsed.system_base, external_pre)
external_post = retrieve_external_params(result.models, result.models, order, PostPerUnit; pu_coeffs = result.coeffs)
sys = SystemModel(result.models, merge_external(external_pre, external_post), order;
                  settings=SystemConfig(), parsed_data=parsed, pu_coeffs=result.coeffs)
```

# See Also
- [`ParsedSystemData`](@ref): Container for original data
- [`convert_system_to_per_unit`](@ref): Per-unit conversion
- [`load_models`](@ref): Parse models from case
"""
function Foundation.SystemModel(
    case::PowerFlowData.Network,
    model_types=nothing;
    system_base::Float64=100.0,
    frequency_base::Float64=Float64(case.caseid.basfrq),
    constant_power_vmin::Float64=0.7
)
    # Default to all registered models if not specified
    types = isnothing(model_types) ? all_models() : model_types
    settings = SystemConfig(; frequency_base, mva_base = system_base, constant_power_vmin)

    # Load models from case (3W transformer conversion handled in Transformer2W parser)
    models = load_models(case, types)

    # 1. Build external-to-internal index mappings (MATPOWER ext2int pattern)
    #    Maps external bus numbers to internal array indices
    models_indexed, index_order = ext2int(models)

    # 2. Retrieve external parameters BEFORE per-unit conversion
    #    (dimensionless params like bus.ide, indices, etc.)
    external_pre = retrieve_external_params(models_indexed, models_indexed, index_order, PrePerUnit)

    # 3. Convert to per-unit (using pre-retrieved external params), preserving the original data
    parsed = ParsedSystemData(models_indexed, system_base)
    result = convert_system_to_per_unit(parsed.models, parsed.system_base, external_pre)
    models_pu = result.models
    pu_coeffs = result.coeffs

    # 4. Retrieve external parameters AFTER per-unit conversion
    #    (physical quantities that need to be in per-unit)
    external_post = retrieve_external_params(models_pu, models_pu, index_order, PostPerUnit; pu_coeffs)

    # 5. Merge external parameters from both phases
    external = merge_external(external_pre, external_post)

    # 6. Create SystemModel with external parameters, index mappings, and coefficients
    # Note: No AddressManager created here - that happens at problem construction
    return SystemModel(models_pu, external, index_order;
                       settings, parsed_data=parsed, pu_coeffs=pu_coeffs)
end

"""
    SystemModel(components::NamedTuple)

Convenience constructor for creating a SystemModel from pre-parsed model components.

This constructor automatically:
1. Applies external-to-internal index mapping (ext2int)
2. Retrieves external parameters from model metadata
3. Creates the SystemModel with proper index ordering

# Arguments
- `components`: NamedTuple of model arrays (e.g., from `load_models`)
- `settings`: System-wide settings (default `SystemConfig()`: 60 Hz, 100 MVA)

# Example
```julia
models = load_models(case, all_models())
sys = SystemModel(models)  # Auto-populates external and order
```

# Note
For full control over per-unit conversion and data preservation,
use `SystemModel(case::PowerFlowData.Network)` instead.
"""
function Foundation._build_system_from_namedtuple(components::NamedTuple;
                                                  settings::SystemConfig{Float64} = SystemConfig())
    # Build external-to-internal index mappings
    components_indexed, index_order = ext2int(components)

    external = retrieve_external_params(components_indexed, components_indexed, index_order)

    return SystemModel(components_indexed, external, index_order; settings)
end
