# Parameter domains of the hand-maintained dynamic models (`PowerfulCore.check_parameter_domains`):
# none has constraints of its own (their parameters are impedances, time constants and reactances
# that the DSL-generated models (TGOV1, SEXS, REGCA1, REGFM_A1, REECA1) check, not this
# project's). VoltagePlayback's recording is per-instance data, not a bounded parameter.
PowerfulCore.check_parameter_domains(x::GENCLS, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::GENROU, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::GENSAL, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::GENTPJ1, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Fault, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::VoltagePlayback, ::AbstractString) = x
