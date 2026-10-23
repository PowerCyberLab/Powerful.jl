# ============================================================================
# Non-Conforming Load Type Traits
# ============================================================================
#
# Sub-classification for models with NetworkPartitionRole = NonConformingLoad: how
# each one's current enters the network Newton iteration.
#
# ============================================================================

"""
    NCLoadType

How a non-conforming load's current enters the network Newton iteration, as the model
declares it (`nc_load_type`):

- [`VoltageDependentNCLoad`](@ref): a current the model computes from its bus voltage, with
  its Jacobian (PQ)
- [`CurrentSourceNCLoad`](@ref): a current its states give (REGCA1, REGFM_A1)

# See Also
- [`NonConformingLoad`](@ref): The parent partition role
- [`nc_load_type`](@ref): Query a model's NC load type
"""
abstract type NCLoadType end

"""
    VoltageDependentNCLoad <: NCLoadType

A load whose current is a function of its bus voltage, iterated with its Jacobian by the network
Newton solve: the model implements [`nc_load_current`](@ref) and [`nc_load_current_jacobian`](@ref)
(and, if some instances inject nothing outside Y, [`injects_nonlinear_current`](@ref)). PQ: the
constant-power and constant-current parts of its ZIP load (`calc_load_current`, with the
constant-power breakpoint); its constant-impedance part is in Y.
"""
struct VoltageDependentNCLoad <: NCLoadType end

"""
    CurrentSourceNCLoad <: NCLoadType

A source whose current its states give: `I = (Ip - j Iq) e^{jθ}` at the angle θ of the current
iterate of its terminal voltage, gathered into the `CurrentSourceBuffer`
(`get_source_current_components`). The converters: the REGCA1 family and REGFM_A1 (its current
behind XL).
"""
struct CurrentSourceNCLoad <: NCLoadType end

# ============================================================================
# Trait Query Function
# ============================================================================

"""
    nc_load_type(::Type{<:AbstractModel}) -> NCLoadType

How the current of a model whose role is `NonConformingLoad` enters the network Newton
iteration. Every such model declares it: there is no default. Building the Y-bus buffers
refuses a non-conforming load without one, naming it (`check_ybus_support`).
"""
function nc_load_type end

"""
    declares_nc_load_type(::Type{<:AbstractModel}) -> Bool

Whether a model type declares its non-conforming load type (a method of `nc_load_type`).
"""
declares_nc_load_type(::Type{M}) where {M<:AbstractModel} = hasmethod(nc_load_type, Tuple{Type{M}})

# The types a non-conforming load can declare, named by a refusal of one without a type
const NC_LOAD_TYPES = (VoltageDependentNCLoad, CurrentSourceNCLoad)

# ============================================================================
# Model-specific NC Load Type Implementations
# ============================================================================

# PQ (its ZIP load: the constant-power and constant-current parts iterated, the constant-impedance part in Y)
nc_load_type(::Type{<:PQ}) = VoltageDependentNCLoad()

# The REGCA1 family (renewable generator - pure current source; REGCA1 and its realization at a finite bandwidth)
# Note: AbstractREGCA1 is imported in network_partition_traits.jl
nc_load_type(::Type{<:AbstractREGCA1}) = CurrentSourceNCLoad()

# REGFM_A1 (grid-forming converter: its current behind XL, a current source of the network)
nc_load_type(::Type{<:REGFM_A1}) = CurrentSourceNCLoad()

# Export all trait types and query function
export NCLoadType, VoltageDependentNCLoad, CurrentSourceNCLoad
export nc_load_type, declares_nc_load_type
