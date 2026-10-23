# GPU support of the dynamic models (`PowerfulCore.gpu_support`): none has a GPU kernel.
# Declared here, beside the generated model files; VoltagePlayback declares its own reason
using PowerfulCore: UnsupportedOnGPU

const NO_GPU_KERNEL = "no GPU kernel: the models of the phasor time-domain simulation run on the CPU only"

PowerfulCore.gpu_support(::Type{<:GENCLS}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:GENROU}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:GENSAL}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:GENTPJ1}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:TGOV1}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:SEXS}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:REGCA1}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:REGFM_A1}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:REECA1}) = UnsupportedOnGPU(NO_GPU_KERNEL)
PowerfulCore.gpu_support(::Type{<:Fault}) = UnsupportedOnGPU(NO_GPU_KERNEL)
