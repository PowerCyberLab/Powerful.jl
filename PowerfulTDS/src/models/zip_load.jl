# Whether the dynamic models are ZIP loads (`PowerfulCore.is_zip_load`): none is.
# Declared here, beside the generated model files
PowerfulCore.is_zip_load(::Type{<:GENCLS}) = false
PowerfulCore.is_zip_load(::Type{<:GENROU}) = false
PowerfulCore.is_zip_load(::Type{<:GENSAL}) = false
PowerfulCore.is_zip_load(::Type{<:GENTPJ1}) = false
PowerfulCore.is_zip_load(::Type{<:TGOV1}) = false
PowerfulCore.is_zip_load(::Type{<:SEXS}) = false
PowerfulCore.is_zip_load(::Type{<:AbstractREGCA1}) = false
PowerfulCore.is_zip_load(::Type{<:REGFM_A1}) = false
PowerfulCore.is_zip_load(::Type{<:REECA1}) = false
PowerfulCore.is_zip_load(::Type{<:Fault}) = false
PowerfulCore.is_zip_load(::Type{<:VoltagePlayback}) = false
