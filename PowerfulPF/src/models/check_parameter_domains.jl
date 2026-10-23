# Parameter domains of the power flow models (`PowerfulCore.check_parameter_domains`): none has
# constraints of its own (ratings and per-unit coefficients are checked where the system is built,
# not free parameters an event can set to a nonsensical value during a simulation)
PowerfulCore.check_parameter_domains(x::Bus, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::PQ, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::PV, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Slack, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Branch, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Transformer2W, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Shunt, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::SwitchedShunt, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::AreaInterchange, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Zone, ::AbstractString) = x
PowerfulCore.check_parameter_domains(x::Owner, ::AbstractString) = x
