# Phase representation of the dynamic models (`PowerfulCore.phase_representation`): the phasor
# time-domain simulation is positive sequence. Declared here, beside the generated model files
using PowerfulCore: PositiveSequence

PowerfulCore.phase_representation(::Type{<:GENCLS}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:GENROU}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:GENSAL}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:GENTPJ1}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:TGOV1}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:SEXS}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:IEEEX1}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:AbstractREGCA1}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:REGFM_A1}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:REECA1}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:Fault}) = PositiveSequence()
PowerfulCore.phase_representation(::Type{<:VoltagePlayback}) = PositiveSequence()
