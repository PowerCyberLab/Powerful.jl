module Models
using DocStringExtensions
using PowerFlowData
using InlineStrings
using TestItems
using StructArrays
using PrecompileTools
using CUDA

using Powerful.PowerCore
import Powerful.PowerCore: get_key_name, supports_format, parse_model, model_metadata, calc_residuals!

# NOTE: Both model_metadata and calc_residuals! are defined in PowerCore
# Models module only provides implementations for specific model types

include("soa_interface.jl")
include("model_traits.jl")
include("model_registry.jl")
include("model_access.jl")

include("topological/topological.jl")

include("bus/bus.jl")
include("pq.jl")
include("pv.jl")
include("branch.jl")
include("transformer.jl")

# Precompilation workload removed - no longer needed with StructArray migration
# Type aliases are defined directly in each model file

end # module