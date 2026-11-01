# ============================================================================
# Network Partition Role Traits
# ============================================================================
#
# Classifies models according to their role in the Y-bus network formulation.
# This determines how each model participates in the two-level network solution:
#
# Level 1: I = Y_gg * E_internal (generators with internal voltage)
# Level 2: Non-conforming loads requiring Newton iteration
#
# ============================================================================

"""
    NetworkPartitionRole

Abstract type for classifying models by their role in Y-bus network formulation.

# Subtypes
- [`FixedNetworkElement`](@ref): Elements that contribute constant admittance (Branch, Shunt, Transformer)
- [`MachineWithInternalVoltage`](@ref): Generators with transient/subtransient reactance (GENCLS, GENROU, GENSAL, GENTPJ1)
- [`NonConformingLoad`](@ref): Loads and current sources iterated by the network Newton solve (PQ, converters)
- [`StaticGenerator`](@ref): Power flow generators without a dynamic model (PV, Slack)
- [`NotInYBusFormulation`](@ref): Controllers and other non-network elements (TGOV1, SEXS)
- [`UnsupportedInYBus`](@ref): Models the Y-bus formulation does not represent, with the reason

Every model declares its role (`network_partition_role`): there is no default.

# Network Solution
The network equations are solved at every bus by a sparse Newton iteration
(`solve_network!`):
- Generators contribute Y = 1/(ra + j*xd') to the diagonal and a Norton current
- Non-conforming loads and current sources are nonlinear current injections
- Static generators are constraints at their buses

# See Also
- [`network_partition_role`](@ref): Query a model's partition role
- [`NCLoadType`](@ref): Sub-classification for non-conforming loads
"""
abstract type NetworkPartitionRole end

"""
    FixedNetworkElement <: NetworkPartitionRole

Marker for network elements that contribute constant (or topology-dependent) admittance.

# Characteristics
- Contribute to Y matrix at construction time
- No state-dependent admittance changes during simulation
- Only change Y matrix when topology changes (line trip, transformer tap)

# Models with this role
- Branch (transmission lines)
- Transformer (with fixed or per-topology tap position)
- Shunt (fixed capacitor/reactor banks)
- Fault (contributes large shunt admittance during fault-on topology)

# Example Y Contribution
```
Y_ij -= y_series        # Off-diagonal (negative of series admittance)
Y_ii += y_series + y_shunt_i   # Diagonal at from-bus
Y_jj += y_series + y_shunt_j   # Diagonal at to-bus
```
"""
struct FixedNetworkElement <: NetworkPartitionRole end

"""
    MachineWithInternalVoltage <: NetworkPartitionRole

Marker for synchronous machines modeled with internal voltage behind reactance.

# Characteristics
- Generator contributes admittance Y_gen = 1/(ra + j*xd') or 1/(ra + j*xd") to Y diagonal
- Provides E_internal = |E| * exp(j*delta) for network solution
- Network current computed via: I = Y_gg * E_internal
- Electrical power: P + jQ = E * conj(I)

# Models with this role
- GENCLS: Classical model with xd' (transient reactance)
- GENROU: Round-rotor with xd" (subtransient reactance)
- GENSAL: Salient-pole with xd" (subtransient reactance)
- GENTPJ1: Type-J machine model with xd" (its xq" may differ: the stator current is a
  machine variable, so the Norton source E = V + (ra + j*xd")*I is exact for any xq")

The admittance uses the subtransient reactance `xd"` when the machine model provides
one (GENROU, GENSAL, GENTPJ1), falling back to the transient reactance `xd'` for models
that only have one (GENCLS).

# See Also
- [`get_internal_voltage`](@ref): Returns E_internal from state variables
"""
struct MachineWithInternalVoltage <: NetworkPartitionRole end

"""
    NonConformingLoad <: NetworkPartitionRole

Marker for loads that don't conform to simple constant-admittance representation.

# Characteristics
- Require Newton iteration to solve for bus voltage
- Current depends nonlinearly on voltage: I = I(V)
- A current injection of the network Newton solve (`nc_load_current`), except the
  constant-impedance part of a load, which is in Y

# Models with this role
- PQ: the constant-power and constant-current parts of its ZIP load, I(V) (its constant-impedance part is in Y)
- Converters (REGCA1, REGFM_A1): current sources of their power Pe + j Qe at the terminal voltage

# Sub-classification
Each declares how its current enters the network Newton iteration ([`NCLoadType`](@ref)):
- `VoltageDependentNCLoad`: a current of the bus voltage, with its Jacobian (PQ)
- `CurrentSourceNCLoad`: a current its states give (the converters)

# See Also
- [`NCLoadType`](@ref): Sub-classification for Newton behavior
- [`nc_load_type`](@ref): Query specific NC load type
"""
struct NonConformingLoad <: NetworkPartitionRole end

"""
    StaticGenerator <: NetworkPartitionRole

Marker for power flow generators left static in the TDS system (instances without
a dynamic model that replaces them).

# Characteristics
- Hold their bus as in the phasor TDS formulation: a PV generator keeps its active
  power `pg` and voltage magnitude `vs` (reactive power free); a Slack generator
  keeps the bus voltage `vs∠va`
- Constraints in the network solve (their bus is retained), not admittances or
  current injections

# Models with this role
- PV, Slack
"""
struct StaticGenerator <: NetworkPartitionRole end

"""
    NotInYBusFormulation <: NetworkPartitionRole

Marker for models that don't participate in network formulation.

# Characteristics
- No admittance contribution to Y matrix
- No current injection into network
- Typically controllers or auxiliary models

# Models with this role
- TGOV1: Governor (provides Pm to generator, doesn't inject current)
- SEXS: Exciter (provides Efd to generator, doesn't inject current)
- REECA1: Renewable electrical control (provides the current commands of REGCA1)
- Area, Zone, Owner (`AreaInterchange`, `Zone`, `Owner`): data without equations

# Note
These models still participate in TDS via their residual equations,
but they don't appear in the network Y matrix or current balance.
"""
struct NotInYBusFormulation <: NetworkPartitionRole end

"""
    UnsupportedInYBus <: NetworkPartitionRole

Marker for models the Y-bus formulation does not represent: `YBusTDSProblem` refuses a
system holding one, naming it and `reason` (`check_ybus_support`).

# Models with this role
- VoltagePlayback: it holds its bus voltage to a recording, which the network solve
  does not impose, and its breaks are events of the TDS problem

A three-phase model (`phase_representation`) declares no role: the Y-bus simulation is
positive sequence and refuses it by its phase representation.
"""
struct UnsupportedInYBus <: NetworkPartitionRole
    reason::String
end

# ============================================================================
# Trait Query Function
# ============================================================================

"""
    network_partition_role(::Type{<:AbstractModel}) -> NetworkPartitionRole

Return the network partition role a model type declares.

Every positive-sequence model (`phase_representation`) declares its role: there is no
default. Building the Y-bus buffers refuses a system holding a model without a role, naming
it (`check_ybus_support`).

# Implementation
Models declare their role here:
```julia
network_partition_role(::Type{<:Branch}) = FixedNetworkElement()
network_partition_role(::Type{<:GENCLS}) = MachineWithInternalVoltage()
network_partition_role(::Type{<:PQ}) = NonConformingLoad()
network_partition_role(::Type{<:TGOV1}) = NotInYBusFormulation()
network_partition_role(::Type{<:VoltagePlayback}) = UnsupportedInYBus("why")
```

# Usage
```julia
role = network_partition_role(GENCLS)
if role isa MachineWithInternalVoltage
    # This generator adds its admittance 1/(ra + j*x) to Y (`_add_generator_contributions!`)
    # and provides E_internal (`get_internal_voltage`)
end
```

# See Also
- [`NetworkPartitionRole`](@ref): The abstract type hierarchy
- [`NCLoadType`](@ref): Further classification for `NonConformingLoad`
"""
function network_partition_role end

# The roles a model can declare, named by a refusal of a model without one
const NETWORK_PARTITION_ROLES = (FixedNetworkElement, MachineWithInternalVoltage, NonConformingLoad,
                                 StaticGenerator, NotInYBusFormulation, UnsupportedInYBus)

"""
    declares_network_partition_role(::Type{<:AbstractModel}) -> Bool

Whether a model type declares its role (a method of `network_partition_role`).
"""
declares_network_partition_role(::Type{M}) where {M<:AbstractModel} =
    hasmethod(network_partition_role, Tuple{Type{M}})

"""
    check_ybus_support(sys)

Throw an `ArgumentError` naming each model of `sys` the Y-bus formulation cannot build: a
model that declares no phase representation (`phase_representation`), a three-phase model
(the Y-bus formulation is positive sequence), a positive-sequence model that declares no role
(`network_partition_role`), with the roles it can declare, a model with instances whose
role is `UnsupportedInYBus`, with its reason, and a `NonConformingLoad` that declares no load
type (`nc_load_type`), with the types it can declare. Called first when the Y-bus buffers are built,
before any role is read.
"""
function check_ybus_support(sys)
    refused = String[]
    for (name, model) in pairs(sys.models)
        M = get_base_model_type(model)
        if !declares_phase_representation(M)
            push!(refused, "$name (no phase representation declared: a method of phase_representation " *
                           "returning PositiveSequence or ThreePhase)")
            continue
        elseif phase_representation(M) isa ThreePhase
            push!(refused, "$name: a three-phase model; the Y-bus simulation is positive sequence")
            continue
        elseif !declares_network_partition_role(M)
            push!(refused, "$name (no Y-bus role declared: a method of network_partition_role returning one of " *
                           join((nameof(R) for R in NETWORK_PARTITION_ROLES), ", ") * ")")
            continue
        end
        role = network_partition_role(M)
        role isa UnsupportedInYBus && !isempty(model) && push!(refused, "$name: $(role.reason)")
        role isa NonConformingLoad && !declares_nc_load_type(M) &&
            push!(refused, "$name (a NonConformingLoad without its load type declared: a method of nc_load_type " *
                           "returning one of " * join((nameof(T) for T in NC_LOAD_TYPES), ", ") * ")")
    end
    isempty(refused) || throw(ArgumentError("The Y-bus simulation (PowerfulYbus) does not support " *
                                            join(refused, "; ")))
    return nothing
end

# ============================================================================
# Model Implementations - FixedNetworkElement
# ============================================================================

# Import model types for trait implementations
using PowerfulCore: get_base_model_type, phase_representation, declares_phase_representation, ThreePhase
using PowerfulPF: Branch, Bus, Transformer2W, Shunt, SwitchedShunt, PQ, PV, Slack, AreaInterchange, Zone, Owner
using PowerfulTDS: GENCLS, GENROU, GENSAL, GENTPJ1, TGOV1, SEXS, AbstractREGCA1, REECA1, REGFM_A1, Fault,
    VoltagePlayback

# Branch (transmission lines)
network_partition_role(::Type{<:Branch}) = FixedNetworkElement()

# Transformer (two-winding)
network_partition_role(::Type{<:Transformer2W}) = FixedNetworkElement()

# Shunt (capacitors and reactors)
network_partition_role(::Type{<:Shunt}) = FixedNetworkElement()
network_partition_role(::Type{<:SwitchedShunt}) = FixedNetworkElement()   # held at binit

# Bus - network node (power balance handled by network solve in Y-bus mode)
network_partition_role(::Type{<:Bus}) = FixedNetworkElement()

# Fault - a shunt admittance in Y of the topologies where it is applied (`process_ybus_events`)
network_partition_role(::Type{<:Fault}) = FixedNetworkElement()

# ============================================================================
# Model Implementations - MachineWithInternalVoltage
# ============================================================================

# Classical generator model (xd')
network_partition_role(::Type{<:GENCLS}) = MachineWithInternalVoltage()

# Round-rotor, salient-pole and type-J generator models (xd")
network_partition_role(::Type{<:GENROU}) = MachineWithInternalVoltage()
network_partition_role(::Type{<:GENSAL}) = MachineWithInternalVoltage()
network_partition_role(::Type{<:GENTPJ1}) = MachineWithInternalVoltage()

# ============================================================================
# Model Implementations - NonConformingLoad
# ============================================================================

# PQ load (its ZIP load: see nc_load_type)
network_partition_role(::Type{<:PQ}) = NonConformingLoad()

# REGCA1 (renewable generator - current source)
network_partition_role(::Type{<:AbstractREGCA1}) = NonConformingLoad()

# REGFM_A1 (grid-forming converter: its current behind its coupling reactance, a current source)
network_partition_role(::Type{<:REGFM_A1}) = NonConformingLoad()

# ============================================================================
# Model Implementations - StaticGenerator
# ============================================================================

network_partition_role(::Type{<:PV}) = StaticGenerator()
network_partition_role(::Type{<:Slack}) = StaticGenerator()

# ============================================================================
# Model Implementations - NotInYBusFormulation
# ============================================================================

# Controllers: their outputs are inputs of other models
network_partition_role(::Type{<:TGOV1}) = NotInYBusFormulation()
network_partition_role(::Type{<:SEXS}) = NotInYBusFormulation()
network_partition_role(::Type{<:IEEEX1}) = NotInYBusFormulation()
network_partition_role(::Type{<:REECA1}) = NotInYBusFormulation()

# Area, zone and owner data: no equations
network_partition_role(::Type{<:AreaInterchange}) = NotInYBusFormulation()
network_partition_role(::Type{<:Zone}) = NotInYBusFormulation()
network_partition_role(::Type{<:Owner}) = NotInYBusFormulation()

# ============================================================================
# Model Implementations - UnsupportedInYBus
# ============================================================================

network_partition_role(::Type{<:VoltagePlayback}) = UnsupportedInYBus(
    "the network solve does not hold a bus voltage to a recording, and the playback's breaks are events of " *
    "the TDS problem, which the Y-bus problem does not take")

# Export all trait types and query function
export NetworkPartitionRole
export FixedNetworkElement, MachineWithInternalVoltage, NonConformingLoad, StaticGenerator
export NotInYBusFormulation, UnsupportedInYBus
export network_partition_role, declares_network_partition_role, check_ybus_support
