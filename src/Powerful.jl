module Powerful
using DocStringExtensions

include("powercore/PowerCore.jl")
# include("essential/Essential.jl")
include("formats/Formats.jl")

include("models/Models.jl")

# Convenience constructors that need both PowerCore and Models
include("system_constructors.jl")

# Re-export key types
using .PowerCore: DAE, SystemModel
export DAE, SystemModel

end
