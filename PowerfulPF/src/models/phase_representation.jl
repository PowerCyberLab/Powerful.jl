# Phase representation of the power flow models (`PowerfulCore.phase_representation`): balanced
# positive-sequence phasors, and the area, zone and owner data of such a network
using PowerfulCore: PositiveSequence

PowerfulCore.phase_representation(::Type{<:Bus}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:PQ}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:PV}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Slack}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Branch}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Transformer2W}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Shunt}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:SwitchedShunt}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:AreaInterchange}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Zone}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Owner}) = PositiveSequence()
