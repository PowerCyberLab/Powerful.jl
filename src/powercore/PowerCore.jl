"""
Core interface definition for all models
"""

module PowerCore
using DocStringExtensions
using TestItems
using Logging
using StructArrays
using CUDA

include("core_types.jl")


export AbstractModel, AbstractModelVec, AbstractModelNumerical
export ModelMetadata, model_metadata
export ContiguousVariables, ContiguousInstances, LayoutStrategy

include("address_manager.jl")
export assign_address!, get_var_addresses

include("format_traits.jl")

include("model_var.jl")

include("model_utils.jl")

include("dae.jl")
export DAE

include("residual_workspace.jl")
export ResidualWorkspace
export gather_variables!, zero_residuals!, scatter_residuals!

include("model_metadata.jl")

include("system.jl")
# Exports from system.jl (setup, residuals, iteration)
export SystemResidualWorkspace
export evaluate_system_residuals!
export setup!, is_setup

end # module