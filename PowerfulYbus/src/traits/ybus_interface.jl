# ============================================================================
# Y-Bus Interface Functions
# ============================================================================
#
# Interface functions that models must implement to participate in Y-bus
# network formulation. Each function corresponds to a specific aspect of
# the network solution.
#
# ============================================================================

using PowerfulCore: ResidualBuffers

# ============================================================================
# Machine Interface (MachineWithInternalVoltage)
# ============================================================================

"""
    get_internal_voltage(ws::ResidualBuffers, model::StructArray, idx::Int) -> Complex

Return the internal voltage E_internal for generator `idx`.

# Required for
Models with `network_partition_role(M) == MachineWithInternalVoltage()`.

# Returns
`Complex`: E_internal = |E| * exp(j*delta) in system reference frame

# Computation
For a synchronous machine:
```
E_internal = E_q * exp(j*delta)
```
where:
- `E_q` is the internal voltage magnitude (from algebraic state)
- `delta` is the rotor angle (from differential state)

# Usage
Internal voltages are gathered into a vector E and used in:
```
I_gen = Y_gg * E_internal
```

# Example Implementation
```julia
function get_internal_voltage(ws::ResidualBuffers, model::StructArray{<:GENCLS}, idx::Int)
    E_q = ws.var_vals.E_q[idx]
    delta = ws.var_vals.delta[idx]
    return E_q * cis(delta)  # cis(x) = cos(x) + i*sin(x)
end
```

# Note
The workspace `ws` contains gathered variable values that are current
for the present timestep after `gather_variables!` has been called.
"""
function get_internal_voltage end

# ============================================================================
# Non-Conforming Load Interface (NonConformingLoad)
# ============================================================================

"""
    nc_load_current(model::StructArray, idx::Int, v_r, v_i, vmin) -> (i_r=..., i_i=...)

Current drawn (consumption convention) by the part of load `idx` that is not a
constant admittance in Y, at the rectangular bus voltage `v_r + j v_i`.

# Required for
Models with `network_partition_role(M) == NonConformingLoad()` whose `nc_load_type`
is `VoltageDependentNCLoad`, together with [`nc_load_current_jacobian`](@ref) and,
if some instances inject nothing outside Y, [`injects_nonlinear_current`](@ref).

# Arguments
- `vmin`: The system's constant-power breakpoint (`sys.settings.constant_power_vmin`)

# Example Implementation (PQ, in network_injections.jl)
```julia
nc_load_current(model::StructArray{<:PQ}, idx, v_r, v_i, vmin) =
    calc_load_current(v_r, v_i, model.pl[idx], model.ql[idx], model.ip[idx],
                      model.iq[idx], zero(v_r), zero(v_r), vmin)  # Z part is in Y
```
"""
function nc_load_current end

"""
    nc_load_current_jacobian(model::StructArray, idx::Int, v_r, v_i, vmin) -> (rr, ri, ir, ii)

Exact rectangular derivative of [`nc_load_current`](@ref):
`(∂i_r/∂v_r, ∂i_r/∂v_i, ∂i_i/∂v_r, ∂i_i/∂v_i)`. A real 2×2 block, since the
current depends on |V| and is not holomorphic in V.
"""
function nc_load_current_jacobian end

# ============================================================================
# Exports
# ============================================================================

export get_internal_voltage
