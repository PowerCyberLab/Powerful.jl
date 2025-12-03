"""
Convenience constructors for SystemModel.

These constructors require access to both PowerCore and Models modules,
so they are defined in the main Powerful module rather than in PowerCore.
"""

using .PowerCore: SystemModel, AddressManager
using .Models: all_models
using PowerFlowData

# load_system is defined in formats/psse_raw.jl and is not in a module
# It's already available in the Powerful namespace

"""
    SystemModel(case::PowerFlowData.Network, model_types=nothing; auto_setup::Bool=true)

Convenience constructor to create a SystemModel directly from a PowerFlowData case.

This is the recommended way to create a power system model. It handles all the
setup automatically: parsing the case file, loading models, creating the address
manager, and setting up address assignments.

# Arguments
- `case`: PowerFlowData Network object (from `PowerFlowData.parse_network`)
- `model_types`: Optional tuple/vector of model types to load. If `nothing` (default),
  uses all registered models via `all_models()`.
- `auto_setup`: Automatically assign addresses (default: `true`). If `false`, you must
  call `setup!(sys)` manually before using the system.

# Returns
- Fully constructed SystemModel, ready for creating DAE and workspace objects

# Examples
```julia
using PowerFlowData
using Powerful

# Load case file
case = PowerFlowData.parse_network("cases/ieee14.raw")

# Quick start: Load and setup in one line (uses all registered models)
sys = SystemModel(case)

# Now ready to use
dae = DAE(Float64, sys)
workspace = SystemResidualWorkspace(Float64, sys)

# Use specific models only
sys = SystemModel(case, (Bus, PQ, Branch))

# Manual setup control
sys = SystemModel(case; auto_setup=false)
# ... custom configuration ...
setup!(sys)
```

# See Also
- [`setup!`](@ref): Manual setup function
- [`is_setup`](@ref): Check if system is set up
- [`all_models`](@ref): Get all registered models
- [`load_system`](@ref): Lower-level function to load models from case data
"""
function PowerCore.SystemModel(
    case::PowerFlowData.Network,
    model_types=nothing;
    auto_setup::Bool=true
)
    # Default to all registered models if not specified
    types = isnothing(model_types) ? all_models() : model_types

    # Load models from case
    models = load_system(case, types)

    # Create address manager
    am = AddressManager()

    # Create system with auto-setup
    return SystemModel(am, models; auto_setup=auto_setup)
end
