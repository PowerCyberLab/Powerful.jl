# Whether the power flow models are ZIP loads (`PowerfulCore.is_zip_load`), which the time-domain
# simulation converts: PQ, with its ZIP coefficients, is
PowerfulCore.is_zip_load(::Type{<:PQ}) = true
PowerfulCore.is_zip_load(::Type{<:Bus}) = false
PowerfulCore.is_zip_load(::Type{<:PV}) = false
PowerfulCore.is_zip_load(::Type{<:Slack}) = false
PowerfulCore.is_zip_load(::Type{<:Branch}) = false
PowerfulCore.is_zip_load(::Type{<:Transformer2W}) = false
PowerfulCore.is_zip_load(::Type{<:Shunt}) = false
PowerfulCore.is_zip_load(::Type{<:SwitchedShunt}) = false
PowerfulCore.is_zip_load(::Type{<:AreaInterchange}) = false
PowerfulCore.is_zip_load(::Type{<:Zone}) = false
PowerfulCore.is_zip_load(::Type{<:Owner}) = false
