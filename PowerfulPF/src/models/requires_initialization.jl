# Whether the power flow models require initialization (`PowerfulCore.requires_initialization`):
# none does. They have no state of their own and start from the power flow solution: the bus
# voltages of Bus from its snapshot (`init_tds_buffers!`); the other models read them, or have no
# variables (areas, zones, owners)
PowerfulCore.requires_initialization(::Type{<:Bus}) = false
PowerfulCore.requires_initialization(::Type{<:PQ}) = false
PowerfulCore.requires_initialization(::Type{<:PV}) = false
PowerfulCore.requires_initialization(::Type{<:Slack}) = false
PowerfulCore.requires_initialization(::Type{<:Branch}) = false
PowerfulCore.requires_initialization(::Type{<:Transformer2W}) = false
PowerfulCore.requires_initialization(::Type{<:Shunt}) = false
PowerfulCore.requires_initialization(::Type{<:SwitchedShunt}) = false
PowerfulCore.requires_initialization(::Type{<:AreaInterchange}) = false
PowerfulCore.requires_initialization(::Type{<:Zone}) = false
PowerfulCore.requires_initialization(::Type{<:Owner}) = false
