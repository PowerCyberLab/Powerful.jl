"""
DC power flow physics shared by the `DCPowerFlow` model kernels, the B-matrix
builders (`gather_branch_susceptances`, `dc_network`) and DC contingency screening
(`dc_bus_injections`), so every DC computation uses one model:

- Branch flow `P = b (θ_from - θ_to)` with `b = 1/(x·tap)`; `tap = windv1/windv2` for
  transformers and 1 for lines. Resistance and phase shift are neglected and the
  network is linear (|V| = 1 on the branches), so PTDF/LODF do not depend on the
  operating point.
- Loads and shunts draw their power at the case bus voltage `vm` (`Bus.vm`): a ZIP
  load draws `pl·f(vm) + ip·vm + yp·vm²` (`f = constant_power_factor`, 1 above the
  system's constant-power breakpoint), a shunt `gl·vm²`. This is the consumption at
  the case's operating point, the same as the AC load model at `vm`; for a flat case
  (vm = 1) it equals the |V| = 1 value.
"""

"""
    dc_susceptance(x, tap=1) -> b

DC susceptance `1/(x·tap)` of a branch with series reactance `x` (pu) and turns ratio
`tap` (`windv1/windv2` for a transformer, 1 for a line).
"""
@inline dc_susceptance(x, tap) = inv(x * tap)
@inline dc_susceptance(x) = inv(x)

"""
    dc_load_power(pl, ip, yp, vm, vmin) -> P

Active power drawn by a ZIP load at the case bus voltage `vm` (pu): constant power
`pl·f(vm)` (`f = constant_power_factor(vm, vmin)`, with `vmin` the constant-power
breakpoint), constant current `ip·vm` and constant impedance `yp·vm²`.
"""
@inline dc_load_power(pl, ip, yp, vm, vmin) = pl * constant_power_factor(vm, vmin) + ip * vm + yp * vm * vm

"""
    dc_shunt_power(gl, vm) -> P

Active power drawn by a shunt conductance `gl` at the case bus voltage `vm` (pu).
"""
@inline dc_shunt_power(gl, vm) = gl * vm * vm
