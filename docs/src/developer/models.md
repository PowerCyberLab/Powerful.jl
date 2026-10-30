# Implementing Models

This page describes how to add a device model. A model that follows these patterns
supports power flow, time-domain simulation, automatic differentiation and reporting,
and dispatches on every array type for which a workspace is specialized.

The reference implementations are the hand-written static models in
`PowerfulPF/src/models/`, in particular `pq.jl`, `pv.jl` and `branch.jl`. The code
excerpts below are taken from `pq.jl`.

!!! note "Generated dynamic model files"
    The files `sexs.jl` (`exciters/`), `tgov1.jl` (`governors/`), `reeca1.jl`, `regca1.jl`
    and `regfm_a1.jl` (`renewable/`) in `PowerfulTDS/src/models/` are machine-generated
    and carry a "DO NOT EDIT" header. The other dynamic model files, including `genrou.jl`,
    `gencls.jl`, `gensal.jl` and `gentpj1.jl`, are hand-maintained, as their headers state.
    All of them follow the metadata and kernel pattern of the hand-written models. New
    models are written by hand following this guide.

## The Model Interface

### Struct Definition and Registration

A model is a `@kwdef mutable struct` that subtypes `AbstractModel`, which gives it keyword
constructors and modifiable fields. The following excerpt is abridged from `PQ`:

```julia
@kwdef mutable struct PQ{Tv} <: AbstractModel
    # Structural fields (concrete types)
    i::Int32           # Bus number (external ID)
    id::String3        # Load identifier
    status::Bool       # In-service flag

    # Numerical fields (type parameter)
    pl::Tv             # Constant power P
    ql::Tv             # Constant power Q
    ip::Tv             # Constant current P
    iq::Tv             # Constant current Q
    yp::Tv             # Constant admittance P
    yq::Tv             # Constant admittance Q
    # ... area, zone, owner, scale, intrpt
end
```

Numerical fields use the type parameter `Tv`, so that the same struct works with `Float64`
and with other number types (for example `Float32`, or `ForwardDiff.Dual` when
differentiating with respect to parameters). Structural fields have concrete types
(`Int32`, `String3`, `Bool`) and store external IDs as given in the data file; the index
maps hold the internal positions. A status field marks in-service instances.

The macro `@register_model` makes the model available to `all_models()` for system
discovery, and the type name is exported:

```julia
@register_model PQ
```

Inside a package, `@register_model` collects the types at precompile time. The package
calls `@init_model_registry` before including its model files and
`register_pending_models!(_PENDING_MODELS)` in its `__init__()` (see
`PowerfulPF/src/PowerfulPF.jl`).

### Keys and Connectivity

The key of a model is a composite of its fields that identifies an instance in events,
lookups and linking. The terminals declare which devices the model attaches to. Traits
such as `get_terminals`, `get_key_name`, `model_metadata` and `calc_residuals!` are
defined in `PowerfulCore`; a package extends them with
`import PowerfulCore: model_metadata, calc_residuals!, get_terminals, get_key_name`.

```julia
# Composite key identifying an instance (used by events, lookups, linking)
get_key_name(::Type{T}) where T<:PQ = [:i, :id]

# Single-attachment device
const PQ_TERMINALS = (Terminal(:Bus, :i),)

get_terminals(::Type{<:PQ}) = PQ_TERMINALS
```

A multi-terminal device declares one `Terminal` per connection:

```julia
# Branch: two-terminal device
const BRANCH_TERMINALS = (
    Terminal(:Bus, :i),  # From-bus
    Terminal(:Bus, :j)   # To-bus
)
```

A terminal can also reference a group, with a composite key where needed:

```julia
# Dynamic models attach to whatever implements the bus group (Bus, or another model
# that provides the same interface)
Terminal(:bus_group, :BUS)

# Controllers match their generator by bus number + machine ID
Terminal(:generators, [:BUS, :ID])
```

Composite keys let a controller such as TGOV1 find its target (GENCLS, GENROU, ...)
across heterogeneous model types in a group. Terminals are the single source of truth for
topology analysis, index linking and external parameter retrieval.

Group interfaces define `key_fields` (for example `[:BUS, :ID]` for `:generators`). A model
that implements the group with different key field names is accepted with a warning and
matched by position; new models should use the key names of the group. Controllers expose
their links through terminals on the controlled group, and all links use the same
`foreign_refs`/`GroupLinkRef` path (see [Three-Layer Indexing](indexing.md)).

### Injection Ports

A model delivers its injection to `:bus_group` through residuals of the group. Every bus
model provides `p` and `q` (power). A bus may additionally provide `ip` and `iq` (current in
phase with and in quadrature to the bus voltage, positive into the bus, with `p = v ip` and
`q = v iq`), or `jr` and `ji` (voltage form: a model that is an EMF E behind an impedance Z
declares its short-circuit current J = E/Z in the network frame).

A model declares one form, its natural one, through `InjectionForms` in its metadata, and
computes the others from it. A model of power form needs no declaration because every bus
accepts it. REGCA1 declares the current form with power as a fallback:

```julia
InjectionForms(:current, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj])
```

The machines (GENROU, GENSAL, GENTPJ1, GENCLS) declare the voltage form with power and
current as fallbacks:

```julia
InjectionForms(:voltage, power = [:p_inj, :q_inj], current = [:ip_inj, :iq_inj],
               voltage = [:jr_inj, :ji_inj])
```

A model of voltage form reads its port voltage (`v_port`, `theta_port`), and the bus maps
it to its own voltage, so that the model's current (E − V)/Z realizes the form
algebraically.

Each connection is delivered through exactly one port, which is chosen when the addresses
are retrieved (`injection_routes`). The bus's port in the model's natural form is used if
the bus implements it. A model of voltage form falls back to the current port and then the
power port; a model of power or current form falls back to the other one of the two. A bus
that implements none of them is an error. The residuals of the forms that are not delivered
have the address `UNROUTED` for that instance and are not scattered. REGFM_A1 declares the
same two ports as REGCA1 (current, with power as a fallback), so at a bus that implements
only the power port its injection routes through the power port.

### Metadata Declaration

The function `model_metadata` collects the declarations of the model in a `ModelMetadata`:

```julia
function model_metadata(::Type{PQ})
    metadata = ModelMetadata(
        name = :PQ,
        vars = PQ_VARS,
        observable_vars = PQ_observable_vars,
        residuals = PQ_RESIDUALS,
        params = PQ_PARAMS,
        observations = PQ_OBSERVATIONS,
        ratings = PQ_RATINGS,
        ext_params = PQ_EXT_PARAMS,
        status_field = :status,
    )

    # Formulation overrides
    metadata.vars_by_formulation[CurrentInjection] = PQ_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = PQ_RESIDUALS_CURRENT
    metadata.vars_by_formulation[DCPowerFlow] = PQ_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = PQ_RESIDUALS_DC

    return metadata
end
```

The fields of `ModelMetadata` (`PowerfulCore/src/foundation/core_types.jl`) are:

- `vars`: variables to allocate in (or reference from) the DAE vectors
- `observable_vars`: computed outputs stored in the workspace without DAE addresses
- `residuals`: equations this model contributes to
- `params`: fields that need per-unit conversion
- `observations`: report outputs computed by `compute_observations`
- `ratings`: device ratings used as conversion bases (for example the MVA base)
- `ext_params`: values retrieved from other models
- `initialized_params`, `init_inputs`: setpoints and power flow inputs for dynamic model initialization
- `variable_limits`: bounds of limited variables (`VariableLimit`), checked at initialization
- `automata`, `memory`: discrete behaviour, namely modes, transitions, and values kept between them
- `status_field`: the `Bool` field that indicates in-service
- `evaluation_pass`: evaluation order (default 0; later passes run after earlier ones)

## Variables and Residuals

### Variable Declarations

A variable owned by another model is declared with the name of its source model and the
field that holds the foreign key:

```julia
const PQ_VARS = [
    ModelVar(
        :theta,                # Local name in workspace
        AlgebVar(),            # Address type
        :Bus,                  # Source model that owns it
        :theta,                # Variable name in source model
        :i;                    # Field holding the foreign key
        description="Bus voltage phase angle",
        units="rad",
    ),
    ModelVar(
        :v,
        AlgebVar(),
        :Bus,
        :v,
        :i;
        description="Voltage magnitude",
        units="pu",
    ),
]
```

Dynamic models typically own variables as well:

```julia
ModelVar(:delta, StateVar(); description="Rotor angle", units="rad")

# A differential state with a mass coefficient (M dx/dt = f), from GENCLS
ModelVar(:omega, StateVar(); description="Rotor speed (pu)",
         mass_coefficient=MassCoefficient(m -> 2.0 .* m.H))
```

The address type `AlgebVar()` marks an algebraic variable (stored in `y`) and `StateVar()`
a differential state (stored in `x`).

### Residual Declarations

A residual that contributes to the equations of another model uses the accumulation mode
`SharedRes`:

```julia
const PQ_RESIDUALS = [
    ModelResidual(
        :p,                    # Local residual name
        AlgebRes(),            # Address type
        SharedRes(),           # Accumulation mode (default :add)
        :Bus,                  # Target model
        :p,                    # Residual name in target
        :i,                    # Foreign key to target
        description="Active power injection from load",
    ),
    ModelResidual(
        :q,
        AlgebRes(),
        SharedRes(),           # g[addr] += residual
        :Bus,
        :q,
        :i,
        description="Reactive power injection from load",
    )
]
```

Control equations use the setter mode, as in `PV`:

```julia
ModelResidual(
    :v_ctrl,
    AlgebRes(),
    SharedRes(:set),           # g[addr] = residual (overwrites)
    :Bus,
    :q,                        # Replaces the Q balance at PV buses
    :i,
    description="Voltage magnitude control equation (replaces Q balance)",
)
```

The equations of the model's own variables are private residuals:

```julia
ModelResidual(:ddelta, StateRes(), PrivateRes(); description="Rotor angle")
```

`SharedRes()` (equivalent to `SharedRes(:add)`) accumulates the contributions of all
devices. `SharedRes(:set)` overwrites the accumulated value and is applied after the `:add`
contributions. `PrivateRes()` holds the model's own equation; each owned variable needs one.

!!! warning "Explicit algebraic equation form (required for the Y-bus solver)"
    Controller algebraic equations must use the explicit form `g = f - y`:

    ```julia
    # CORRECT: g_P_m = P_m0 - P_m (explicit form)
    outputs.g_P_m[i] = mask * (params.P_m0[i] - vars.P_m[i])

    # WRONG: g_P_m = P_m - P_m0 (inverted sign)
    outputs.g_P_m[i] = mask * (vars.P_m[i] - params.P_m0[i])
    ```

    The Y-bus solvers (PowerfulYbus) move these algebraic variables by a Newton iteration
    on all the algebraic rows and the network equations together (`y -= J⁻¹ g`, with
    `J = ∂g/∂y` the sparse Jacobian taken by ForwardDiff). For a residual of the form
    `g = (target) - (variable)`, the entry of its own variable is `-1` and the step is
    `y + g = target` where nothing else moves. A row that does not read its own variable is
    also solved, as long as the Jacobian of the whole system is nonsingular (the coupling
    through the network is part of it).

    The residual `:g_foo` pairs with the variable `:foo`. The Y-bus solver discovers these
    pairs from the metadata (`AlgebraicNewton` in `PowerfulYbus/src/integration/`).

    The form applies to generator control inputs (`P_m`, `E_q`), governor outputs (`pref`,
    `pout`) and any controller algebraic variable that holds a setpoint.

### Multi-Formulation Support

A model that supports several formulations declares the variables and residuals of each
formulation:

```julia
# PowerInjection formulation (default)
const PQ_VARS = [ ... polar variables (theta, v) ... ]
const PQ_RESIDUALS = [ ... power residuals (p, q) ... ]

# CurrentInjection formulation override
const PQ_VARS_CURRENT = [
    ModelVar(:v_r, AlgebVar(), :Bus, :v_r, :i; description="Voltage real part", units="pu"),
    ModelVar(:v_i, AlgebVar(), :Bus, :v_i, :i; description="Voltage imaginary part", units="pu"),
]
const PQ_RESIDUALS_CURRENT = [
    ModelResidual(:i_r, AlgebRes(), SharedRes(), :Bus, :i_r, :i,
                  description="Real current injection from load"),
    ModelResidual(:i_i, AlgebRes(), SharedRes(), :Bus, :i_i, :i,
                  description="Imaginary current injection from load"),
]
```

The overrides are registered in `model_metadata` as shown above. See [Multi-Physics Formulations](formulations.md).

## Parameters and External References

Parameters are the fields that need per-unit conversion. Each is declared with its
`PUType` and physical units:

```julia
const PQ_PARAMS = [
    ModelParam(:pl, PowerNominal, description="Constant power P", units="MW"),
    ModelParam(:ql, PowerNominal, description="Constant power Q", units="MVAR"),
    ModelParam(:ip, PowerNominal, description="Constant current P", units="MW"),
    ModelParam(:iq, PowerNominal, description="Constant current Q", units="MVAR"),
    ModelParam(:yp, PowerNominal, description="Constant admittance P", units="MW"),
    ModelParam(:yq, PowerNominal, description="Constant admittance Q", units="MVAR"),
]
```

During `SystemModel` construction the fields are converted from physical units (MW, kV) to
per-unit, the conversion coefficients are stored in `sys.pu_coeffs`, and the original
values are kept in `sys.parsed_data`.

External parameters are values retrieved from other models:

```julia
const PQ_EXT_PARAMS = [
    ExtParam(:basekv, PQ_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage"),
    ExtParam(:ide, PQ_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type (1=PQ, 2=PV, 3=Slack)"),
    ExtParam(:vm, PQ_TERMINALS[1], :vm, Dimensionless(), PostPerUnit(),
             required=true, description="Case bus voltage magnitude (pu), the operating voltage of the DC model"),
]
```

The retrieval phase is `PrePerUnit()` for values needed before the per-unit conversion
(dimensionless values, bases) and `PostPerUnit()` for quantities that must already be in
per-unit. Residual calculations access the values through the `external` argument of
`calc_residuals!`:

```julia
external.vm[i]  # Case voltage of the bus of the i-th load (used by the DC kernel)
```

A value is a struct field when it is intrinsic to the device (load consumption, generator
rating), comes directly from the data file, may be modified independently (parameter
sweeps, events) and is the same across formulations. A value is an external parameter when
it belongs to another model (bus base kV, bus type) or is derived from topology. A PQ load,
for example, stores its power consumption and retrieves the bus base voltage.

## The Residual Kernel Pattern

Every model with residuals consists of three parts: a CPU kernel, a GPU kernel and a data
adapter.

### CPU Kernel

The kernel is a pure function of NamedTuples of arrays (and scalars such as `vmin`):

```julia
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:Vector}  # CPU dispatch
)
    @inbounds @simd for i in eachindex(outputs.p)
        mask = active[i]  # Branchless masking

        # ZIP load: P = pl·f(v) + ip·v + yp·v², f = constant-power factor below vmin
        power = calc_zip_load(vars.v[i], params.pl[i], params.ql[i],
                              params.ip[i], params.iq[i], params.yp[i], params.yq[i], params.vmin)

        # Sign convention: load consumption = negative injection into bus
        outputs.p[i] = mask * (-power.p)
        outputs.q[i] = mask * (-power.q)
    end

    nothing
end
```

Masking is branchless: `mask * value` replaces a branch on `active[i]`, which keeps the
loop vectorizable under `@inbounds @simd`. The argument `::Type{<:Vector}` selects the CPU
method. The physics lives in a separate function (`calc_zip_load`) that the CPU kernel, the
GPU kernel and reporting share. Injections follow the sign convention that positive values
flow into the bus.

### GPU Kernel

A GPU extension adds a method of the same kernel function, dispatched on the GPU array
type, together with the CUDA device function that it launches. The device function indexes
by thread and calls the same physics function as the CPU kernel, and the launcher method
dispatches on `::Type{<:CuArray}`:

```julia
import PowerfulPF: pq_residual_kernel!, calc_zip_load

# Device kernel
function pq_residual_cuda_kernel!(
    outputs_p, outputs_q,
    params_pl, params_ql, params_ip, params_iq, params_yp, params_yq, vmin,
    vars_v, active, n
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= n
        mask = active[i] ? one(eltype(outputs_p)) : zero(eltype(outputs_p))
        power = calc_zip_load(vars_v[i], params_pl[i], params_ql[i],
                              params_ip[i], params_iq[i], params_yp[i], params_yq[i], vmin)
        @inbounds outputs_p[i] = mask * (-power.p)
        @inbounds outputs_q[i] = mask * (-power.q)
    end
    return nothing
end

# Launcher
@inline function pq_residual_kernel!(
    outputs::NamedTuple{(:p, :q)},
    params::NamedTuple{(:pl, :ql, :ip, :iq, :yp, :yq, :vmin)},
    vars::NamedTuple{(:v,)},
    active::AbstractVector{Bool},
    ::Type{<:CUDA.CuArray}  # GPU dispatch
)
    n = length(outputs.p)
    (; threads, blocks) = cuda_kernel_config(n)
    CUDA.@cuda threads=threads blocks=blocks pq_residual_cuda_kernel!(
        outputs.p, outputs.q,
        params.pl, params.ql, params.ip, params.iq, params.yp, params.yq, params.vmin,
        vars.v, active, n
    )
    nothing  # kernels run in order on the default stream
end
```

### Data Adapter

The adapter `calc_residuals!` extracts the arrays from the workspace, the model and the
system settings (`ws.system`) and calls the kernel with NamedTuples. It performs no
computation.

```julia
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pq::StructArray{<:PQ},
    _external::NamedTuple,
    ::PowerInjection,             # Formulation dispatch
    t = zero(T);                  # Simulation time
    active_mask::AbstractVector{Bool}
) where {T, AT}
    pq_residual_kernel!(
        (p = ws.res_vals.p, q = ws.res_vals.q),              # Outputs from workspace
        (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq,
         vmin = ws.system.constant_power_vmin),              # Params from model and system
        (v = ws.var_vals.v,),                                # Variables from workspace
        active_mask,                                         # In-service mask (passed in)
        AT                                                   # Dispatch type parameter
    )
end
```

The keyword `active_mask` is supplied by the system workspace and is not stored in
`ResidualBuffers`. Dynamic models read initialization setpoints from `ws.setpoints`. The
type parameter `AT` selects the CPU or GPU kernel method.

A model with several formulations has one adapter per formulation, dispatched on the
formulation argument (`::PowerInjection`, `::CurrentInjection`). Each adapter calls the
kernel of its formulation (`pq_current_kernel!` for CurrentInjection) with the variables and
residuals declared for that formulation.

## Observations and Reporting

Observations are computed outputs that are not part of the DAE vector. Each is declared
with its quantity type, which selects the conversion to report units:

```julia
const PQ_OBSERVATIONS = [
    ModelObservable(:p_consumed, ActivePower, units="pu",
                    description="Active power consumed", default_output=true),
    ModelObservable(:q_consumed, ReactivePower, units="pu",
                    description="Reactive power consumed", default_output=true),
]
```

The method `compute_observations` evaluates all instances at once. It extracts the
addresses once and reuses the physics function of the residual kernel:

```julia
function compute_observations(sol::PowerFlowSolution, ::Type{<:PQ}, var::Symbol)
    pq = sol.prob.sys.models.PQ
    vmin = sol.prob.sys.settings.constant_power_vmin

    # Voltage magnitudes for all loads, extracted once
    v_addrs = get_addresses(sol.prob.address_manager, :PQ, :v)
    v_vals = sol.dae.y[v_addrs]

    n = length(pq)
    result = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        # Same physics function as the residual kernel
        power = calc_zip_load(v_vals[i], pq.pl[i], pq.ql[i],
                              pq.ip[i], pq.iq[i], pq.yp[i], pq.yq[i], vmin)
        if var == :p_consumed
            result[i] = power.p
        elseif var == :q_consumed
            result[i] = power.q
        else
            error("Unknown PQ output variable: $var")
        end
    end

    return result
end
```

## Testing

Tests use `@testitem` blocks. A model needs a test of the kernel in isolation, a test of
the full residual calculation, and a check that inactive devices are masked.

## Reference Implementations

| Model | Pattern | File |
|-------|---------|------|
| PQ load | Foreign variables (`theta`, `v` of the bus); `SharedRes()` power balance; several formulations; ZIP load physics | `PowerfulPF/src/models/pq.jl` |
| PV generator | Foreign bus variables; `SharedRes()` power injection; `SharedRes(:set)` voltage control that replaces the Q balance; `evaluation_pass = 1`, so that in CurrentInjection it reads the currents accumulated by pass-0 models | `PowerfulPF/src/models/pv.jl` |
| Branch | Foreign variables from both ends (buses `i` and `j`); `SharedRes()` contributions to both buses; the formulations change the expressions | `PowerfulPF/src/models/branch.jl` |
| GENCLS | Owned `StateVar` variables (rotor angle, speed) with private residuals; foreign algebraic variables (terminal voltage) through `:bus_group`; replaces a static PV or Slack generator (`get_replacement_map`, `get_data_inheritance`); initialization from the power flow (`init_inputs`, `initialized_params`) | `PowerfulTDS/src/models/generators/gencls.jl` |
| TGOV1 | Controller of generators; linked through a terminal on the controlled group; reads initialization data from the target's initialization (`InitInput` of kind `:from_group_init`) | `PowerfulTDS/src/models/governors/tgov1.jl` |

A controller declares the group it controls and a terminal with an explicit composite key:

```julia
# Which group this controller targets
get_controlled_group(::Type{<:TGOV1}) = :generators

# Composite key for matching (must be defined explicitly)
get_terminals(::Type{<:TGOV1}) = (Terminal(:generators, [:BUS, :ID]),)
```

## Implementing Anti-Windup Limits

Dynamic models often have state variables that must stay within bounds (for example the
REECA1 power order `PMIN ≤ s5_y ≤ PMAX`). A limited state is an automaton (see
[Implementing Automata](#Implementing-Automata)), and the generated models (SEXS, TGOV1,
REECA1, REGCA1) state their anti-windup that way. The implementation has three parts: the
automaton, the kernel that reads its gate, and the declaration of the bounds.

### The Automaton

Three modes hold a gate that multiplies the rate of the state: `free` (gate 1), and
`upper` and `lower` (gate 0). A free state that reaches a limit with its unlimited rate
pointing beyond it is held there, and a held state is released where its unlimited rate
turns inward:

```julia
reach_upper = Atom(:y_upper, :>=) & Atom(:rate, :>)   # y ≥ upper and rate > 0
reach_lower = Atom(:y_lower, :<=) & Atom(:rate, :<)
automata = [Automaton(:LAG;
    modes = [Mode(:free; holds = [:gate => 1]), Mode(:upper; holds = [:gate => 0]), Mode(:lower; holds = [:gate => 0])],
    transitions = [Transition(:free, :upper; when = reach_upper), Transition(:free, :lower; when = reach_lower),
                   Transition(:upper, :free; when = Atom(:rate, :<)), Transition(:lower, :free; when = Atom(:rate, :>))],
    initial = [InitialMode(:upper; when = reach_upper), InitialMode(:lower; when = reach_lower), InitialMode(:free)])]
```

The transitions are located: the simulation finds where `y - upper` or `lower - y` crosses
zero, and where the unlimited rate turns inward, and restarts the integrator there. At the
start and after every event the transitions apply at once, so a condition that starts on its
other side is no crossing. A held state has a zero rate wherever its unlimited rate points
outward, so the release must be located; a step whose ends both have a zero rate would
pass over the interval where the state should move inward.

### The Kernels

`calc_residuals!` multiplies the rate by the gate (`ws.automata.discrete.gate`), and
`calc_automaton_values!` computes the atoms into `ws.automata.values`: the unlimited rate
and the distances `y - upper` and `y - lower` of the state from its limits. The
expressions must be branch-free and AD-compatible: within an integration step the equations
do not switch, and the Jacobian is that of the mode.

```julia
out.dy[i] = active[i] * p.gate[i] * (u - v.y[i]) / p.T[i]        # residual kernel
vals.rate[i] = (u - v.y[i]) / p.T[i]                              # values kernel
vals.y_upper[i] = v.y[i] - p.upper[i]
vals.y_lower[i] = v.y[i] - p.lower[i]
```

### The Bounds

The limited state is declared with its bounds as a `VariableLimit`
(`variable_limits = [VariableLimit(:LAG_y, :Vmin, :Vmax)]`). The initialization checks the
initial value of the variable against them (`enforce_initial_limits!`). A violated bound
that names a model field is moved to the value with `TDSConfig(adjust_limits = true)`, and
the start is computed again; any other violation is an error. A bound may also be computed
from the parameters (`m -> 2 .* m.Vmax`); such a bound cannot be moved to a violating
initial value.

A limit is tested with a scenario that pushes the state to its limit (the state must stay
at the limit within tolerance) and a scenario that reduces the demand (the state must come
back off the limit). For a lag with a known input the trajectory has a closed form, against
which the switching times and the trajectory are checked
(`Powerful/test/tds/test_automata.jl`).

## Implementing Automata

A model with discrete behaviour (relays, timers, protection, anti-windup) declares automata
and memory in its metadata (`ModelMetadata(; automata, memory)`, see
`PowerfulCore/src/foundation/automata.jl`). An automaton consists of modes that hold values
of the model's discrete quantities, transitions between the modes on conditions, and memory
that the transitions set. A condition is built from atoms (a value compared with zero) or
crossings of guards. Two kernels of the pattern above read the discrete state from the
workspace's `AutomatonState` (`ws.automata`):

- `calc_residuals!` reads the discrete quantities (`ws.automata.discrete`) and memory
  (`ws.automata.memory`) of each instance where its equations depend on them.
- `calc_automaton_values!` computes every value that the automata name into
  `ws.automata.values`, from the same inputs. These are the atoms `a - b`, the crossing
  guards, and the sets and memory starts that are not numbers.

```julia
# A relay on x: on above hi, off below lo; z integrates the time on (r = 1 while on)
automata = [Automaton(:relay; modes = [Mode(:off; holds = [:r => 0]), Mode(:on; holds = [:r => 1])],
    transitions = [Transition(:off, :on; when = Atom(:x_hi, :>)), Transition(:on, :off; when = Atom(:x_lo, :<))],
    initial = [InitialMode(:off)])]
# calc_residuals!:        ws.res_vals.dz .= active_mask .* ws.automata.discrete.r
# calc_automaton_values!: ws.automata.values.x_hi .= ws.var_vals.x .- m.hi
#                         ws.automata.values.x_lo .= ws.var_vals.x .- m.lo
```

The first subsection below states what the kernels must satisfy; the others describe how the
runtime settles the transitions.

### Kernel Requirements

Both kernels compute the quantities of instance `i` from the variables, parameters,
setpoints and external parameters, together with the discrete quantities and memory of
instance `i` itself. A kernel that reads another instance's discrete quantity or memory
breaks this contract (a rate `r[1] * r[2]` is an example). When the contract holds, a
switch of one instance changes only that instance's rates, and instances that switch at
the same instant (identical devices tripping together) do not interact. Two mechanisms rely
on this. Evaluating the automaton with dual numbers at the instant differentiates through
the switches with one switching time per instance, and the settle after a located event
considers only the instances that switched.

Every simulation checks the contract for the values at setup, using their trace (described
below). A model that violates it is accepted with a warning that names the model, the value,
the quantity and the instances; the simulation then settles every instance of that model in
every round, which is correct and slower. An augmentation that integrates sensitivities
through the switches also checks the rates at instants where several instances switch. Each
switch alone from the discrete state before the instant, and each switch left out from the
state after it, must change only its own instance's rates, and a kernel that fails this
check is refused. An interaction that is visible at neither end remains the responsibility
of the kernel author.

The values kernel is traced at setup with the global detector of SparseConnectivityTracer.
The trace records which variables each value reads; the discrete quantities and memory of
each instance enter as separate inputs. The automaton callback then gathers only the
variables that the values read. An augmentation that locates where a guard grazes a limit
reads the entries of the state that each value reads from the same trace, which its
classification of a grazed guard requires. For a kernel that does not branch on parameters
or setpoints, the trace is made once for all solves of a problem (`TDSProblem` `traces`).

The time, variables, discrete quantities and memory enter the trace as tracers, so the
values kernel must avoid `if`, `?:` and `&&` on them. The function `ifelse` is allowed (it
traces both of its branches), as are `min`, `max`, `abs` and arithmetic. Real parameters and
setpoints enter as tracers first, so that the trace holds under whatever events do to them.
A kernel that branches on them cannot be traced this way and is traced with the parameters
and setpoints at their values. That trace holds until an event changes them, so it is made
again after a `FieldModification` or `SetpointChange` of the model (an event of a kind that
does not declare what it changes re-traces every model), and every variable is gathered,
whichever the values read. Evaluating the automaton with dual numbers at an instant uses
what the values read at setup, and throws where a later trace finds that this changed. A
kernel that cannot be traced either way throws at setup, naming the model.

The simulation also traces `calc_residuals!` at setup, for the instances whose discrete
quantities or memory are read by a residual without mass (an algebraic residual, or the
residual of a state whose mass is zero for that instance, as for a zero inertia). Only the
transitions of these instances can change the algebraic equations, and the network must
then be reinitialized (`PowerFlowInit` solves the algebraic variables and the states without
mass together). A settle evaluates nothing where none of those instances switches; an
anti-windup automaton is such a case, because its gate holds a rate with mass. Otherwise
the settle compares the contributions of the switched models to the rows without mass
before and after the transitions. The comparison covers the system's rows without mass
where a model sets a shared residual (`SharedRes(:set)`, which reads what the others
accumulated) or where another model sets a residual to which the model adds.

Parameters enter this second trace at their values, so a residual kernel may branch on them
(for a time constant that may be zero). The trace is made again after every event that
changes the model's fields (`FieldModification`; an event of a kind that does not declare
what it changes re-traces every model). The variables, time and setpoints enter as tracers,
so the kernel must avoid `if`, `?:` and `&&` on them and on the discrete quantities and
memory. A `SetpointChange` needs no new trace.

Where the algebraic equations changed (by a transition or an event), the reinitialization
factorizes its Jacobian afresh with pivoting. Refactorizing with the pivots of an earlier
factorization could meet a zero pivot, because a mode can move the nonzero of a row off its
own variable: the high-voltage clamp of REGCA1 holds `Vt = Volim` in the equation of `dq`.

### Settling an Instant in Rounds

One rule applies to every instant: an event while a model runs, the start, and the coupling
between devices through a network. The instant settles in rounds.

1. In a round, every transition whose condition holds, in every automaton of every
   in-service instance, fires on the point before the round. An automaton fires at most one
   transition per round, and two transitions that hold at once from one mode are refused.
   The sets of the transitions read the values from before the round. The order of
   automata, modes and transitions has no meaning.
2. After the round, the algebraic point is solved again from the modes entered. This point
   is provisional, and the next round reads it. The algebraic check of the settle
   (`AlgebraicRounds`) is made for each round on the instances that took a step. A round
   that leaves the equations without mass as they were leaves the point as it is (it is
   still a solution), so an instant that settles in one round costs the single solve it
   always did. The host performs the solve: the phasor simulation through its
   `initializealg`, the Y-bus DAE through `HostReinit`, and a fixed-step integration
   through its own algebraic solve (`FixedInstant`).
3. The instant is settled when no transition holds at the point. Only the settled point is
   published (the integration saves it), and the sliding check and the sides of the located
   values are made there. A nonfinite value, a right side outside its domain, a rate that is
   not finite, and a solve that fails or is singular each fail at once, in any round (the
   integration terminates with `ConvergenceFailure`, or the error is thrown).
4. The budget (`max_rounds`) counts rounds, once for an instance whatever the number of its
   automata that fire.
5. The start chooses the initial modes and settles by the same rule. The initialization
   does not know which modes the start chooses, so the equations without mass in the initial
   modes may not hold at the point the start was handed (a residual above the integration's
   tolerance). In that case the point is solved from the initial modes first, and again
   between the rounds (`StartResolve`).

   The adapter starts its automata without a network first, because the detection of the
   state-only constraints reads the discrete state and the solve needs the constraints. It
   repeats the start with the solve only where the first start fired a transition or the
   equations of the initial modes do not hold at the point. The callback of the integration
   starts the automata at its first step with the integration's reinitialization. A start
   that fires transitions without changing the equations stays at rest (the values decided as
   zero are not released); a start that solves the point has left the rest
   and releases them.

### Located Values

A value located at its zero (an atom that changed its truth, or a crossing guard) is at zero
only to rounding, and the reinitialization of the network at the instant moves it by the
tolerance of its Newton iterations. The atoms over such a value are therefore decided by the
direction of its motion.

At the instant, every atom over the value holds as the crossing leaves it, whichever
condition reads the atom: `>=` and `>` are true after an upward crossing, `<` and `<=` are
false. A release that reads the complement of its hold's atom (hold on `x - hi >= 0`, release
on `x - hi < 0`) consequently does not fire at the hold.

After the instant, the value is measured from zero, or from where the instant left it where
that is across zero from its side. A later change of sign is a crossing only where the
value's rate points the same way. Otherwise the change comes from what the integration adds
to the value besides its motion. A held integrator (its rate gated to zero) still moves by
the tolerance of the nonlinear solves, and the short steps after an event move a value by
their error, which may be against its rate. Such a value keeps the truth of the crossing
until the next instant measures it again (held still), or until a step starts with the value
back on its side (moving back).

An automaton that changes mode reads its values afresh. A value whose discrete quantities or
memory a transition changes (a timer whose transition records the time it reads) moves by
that change, a jump that its own motion does not produce. It loses its side at once, located
or not, and is decided by its value, using the same trace of the values (the discrete
quantities and memory each reads). A transition located where its value turns back across
the surface at once (a Filippov sliding mode) throws, naming the instance.

### Rates at Located Events

The sliding check and the values held at zero by continuity use the rate of a value,
`ġ = g_u u̇ + g_t`, along the trajectory of the DAE `M u̇ = F(u, t)`. For a value that reads no
variable without mass (an algebraic variable, or a state of zero mass), the rate follows from
the rows with mass alone, `u̇ = F/M`, and needs no solve. A value that reads such a variable
needs the rows without mass completed, `F_ZZ u̇_Z = -(F_ZW u̇_W + ∂F_Z/∂t)`, which factorizes
`F_ZZ` with KLU (the symbolic analysis is kept from one factorization to the next).

The check after the transitions reuses the factorization from the one before them where no
automaton that switched has its discrete state read by a residual without mass, and the
state and time are the same (the equations of the rows without mass are then those of the
event). An `F_ZZ` that KLU finds singular, or whose condition number it estimates at
`1/eps()` or more, has no rates: `SingularRatesError` names its empty rows and columns, so
that no decision rests on noise.

### Values Held at Zero by Continuity

Some values are zero at an instant by the structure of the model, and their computed value
there is rounding of either sign. Machine-generated models declare these values in their
metadata with the following options.

- `Transition(...; continuity = [:surface => [:w, ...]])`: where the transition fires on the
  zero of `surface`, the values `w` vanish (in the mode it leaves) and the transition leaves
  the algebraic state where it is. For example, REGCA1 entering its clamp where `Vt` crosses
  `Volim` holds `dq` at the zero that the clear mode held it at, because the equation of the
  clamp, `Vt = Volim`, already holds there. At the instant the runtime decides the atoms over
  each `w` by the sign of `w`'s rate in the mode entered (the rate the sliding check
  evaluates), and afterwards measures `w` as it measures a located value. A rate that points
  out of the mode fires the next transition at the same instant.
- `Transition(...; crossing_continuity = [:surface => [:w, ...]])`: the same for facts that
  hold only where the transition fires at a crossing located while the simulation was
  integrating its source mode, which is the first round of an instant. The source mode was
  then settled before the instant, so none of its other exits held. In REGCA1's high-voltage
  clamp (`clamp -> clear` when `dq < 0`, `clamp -> floor` when `dq >= 0` and
  `Iq_lag v - dq < Iolim`), `Vt - Volim` is zero at the release of `dq` only from a settled
  clamp: the exit to `floor` not holding there makes the reactive power of `clamp` and `clear`
  agree. A clamp entered earlier in the same instant has no such guarantee. The runtime
  decides the values by their rate as with `continuity`, and applies the marker only in the
  first round of an instant with a located crossing. A transition that fires in a later round,
  at a start, at an event's reinitialization or at another instant does not use it. A
  `continuity` marker holds however the source mode was entered.
- `Transition(...; continuity_writers = [:w => [:automaton, ...]])`: for a value `w` held by
  `continuity` or `crossing_continuity` on this transition, this option lists the automata of
  the same instance whose firing voids the hold, because `w` reads a quantity that one of them
  writes. (REECA1's PI clip holds values that read the limit `Iqmax`, which the automata
  `injection` and `ipmax_hold` can write.) The fact stated by the marker holds for the point
  as it was before such an automaton fired. The runtime holds `w` only where none of the
  listed automata fires in the round of the transition; otherwise it decides `w` by its own
  value, as it does for a value that no marker names. A writer that fired in an earlier round
  leaves the hold intact: the algebraic point is solved again between rounds, and the round of
  the transition reads what the earlier rounds left (see
  [Settling an Instant in Rounds](#Settling-an-Instant-in-Rounds)). A value that is not
  listed has no writers. `check_automata` refuses a value that the transition does not hold,
  an automaton the model does not have, and the transition's own automaton. A model declares
  its writers independently of the order of its automata, since every automaton fires in the
  round in which its transition holds.
- A marker of any kind (`continuity`, `crossing_continuity`) for a value that reads a
  quantity the network solves (a variable without mass such as a bus voltage, or what a value
  computes from one) is skipped in a round in which another transition, of any instance,
  changes a quantity delivered to the network. Such a change is either a discrete quantity or
  memory that a residual without mass reads (a `bus_group` injection, or a row of the network
  that a connection reads) set to another number by the mode entered or by a set, or a state
  governed by the mode entered where such a residual reads one (`AlgebraicCheck` `inputs`).
  The generator derives a marker from the inputs of the model as they are across the
  transition; what other devices do reaches the model through the network, and Powerful
  bounds that effect. For example, REGCA1's LVPL automaton projecting `Ip` in the round its
  high-voltage clamp engages moves `Pe`, and with it the bus voltage that the markers of the
  clamp read. The value is then decided by its own value.
- A rate decides a held value only where the integration resolves its sign: `|ġ|` must
  exceed the first-order bound `Σ_j |∂ġ/∂u_j| tol_j` of its error from the state's
  (`tol_j = abstol_j + reltol_j |u_j|`; see "The sign of a rate the integration does not
  resolve" in the [TDS architecture](tds.md) page). A rate of the size of that bound is the
  error of the integration, and its sign is a coin. For example, a state oscillating at a
  rate of the order of 10² per second has a state error of 10⁻⁶ that moves the rate of a slow
  envelope by about 10⁻⁴. The value then stays at its zero through the instant, kept against
  the release of an instant that changes the equations until the instant ends and the value is
  measured from where it left it. A transition that reads the value strictly therefore does
  not fire on the rounding of the reinitialization, and a crossing that follows is located as
  any other. The release of REGCA1's clamp (`dq` crossing zero, which holds `Volim - v` at
  zero) is such a value. The measured magnitude is a state continuous across it, so its rate
  in `clear` is the clamp's, zero. Deciding this value by the sign of the error entered the
  clamp again at once and left it again, up to the budget of rounds, at tolerances where the
  error was positive. A fixed-step integration has no tolerance and decides by the rate.
- `Automaton(...; at_rest = [:rate, ...])` lists the rates that the initialization makes zero
  (an anti-windup lag's rate in `free`, which is zero at the equilibrium the simulation starts
  from). They are decided as the exact zero they are (`>=` and `<=` hold, `>` and `<` do not)
  at the start and while the system stays at that rest. At rest their rates are zero too, so
  the initial clauses choose, in their order, the mode the model was initialized in: a governor
  initialized on its (adjusted) valve limit starts and stays in `free`. The first instant that
  changes the equations (an event, or a transition anywhere) releases them, and from there on
  they are decided by their values.

### Conditions That Change Between Sampled Points

The located callback detects a change in the truth of a condition where its component changes
sign between two points at which DiffEq evaluates it (the ends of a step and `interp_points`
between them). A condition that holds only between two of these points goes unseen. Two cases
occur: a value crosses its surface and comes back inside the step, or atoms change their
truths one after the other (the state reaches its limit, then its rate turns outward).

At the start of each step (`PowerfulTDS/src/solvers/step_extremum.jl`), every component that
can fire and has one truth at both ends of the step is examined. The examination finds where
the values read by its atoms cross their surfaces inside the step. Where the ends are on
different sides there is one crossing. Where they are on one side there are none or two,
depending on the extremum of the value on the dense output. The extremum lies beyond the
surface when the rates at the two ends point opposite ways, the first toward the surface
(the rates are the derivative of the dense output, carried through one Dual evaluation of
the values). The extremum is located by the Hermite interpolant of the ends and a Newton step on
the rate, and it is shared with the grazing of guards. The truth of the condition between the
crossing times then follows from the truths of its atoms.

Where this truth differs from that at the start, the component is reported folded: as it is
up to a time inside that interval, and as its value there after it. The folded component is a
function that changes sign where the condition first changes its truth. DiffEq locates that
time and moves the integrator back to it, and the event proceeds as any other located
crossing. An extremum exactly at the surface is a tangency and fires nothing.
`TDSConfig(locate_extrema = false)` leaves such conditions unseen.

### Where Automata Run

The phasor simulation (`solve(::TDSProblem)`), the Y-bus simulation (PowerfulYbus) and any
other simulation mode built on `TDSProblem` run the automata of their models with the same
runtime. The Y-bus DAE solve locates the transitions with the automaton callback. The
trapezoidal solver cuts its fixed steps where a condition changes its truth
(`step_automata!`): it locates the atoms and crossings that changed within the bracket and
settles there with the code of the adaptive callback (`settle_instant!`), with the rates of
the values that continuity holds at zero evaluated on the Y-bus DAE.

The Y-bus simulation refuses a machine model with automata, since its machine residuals are
`calc_ybus_residuals!`. It solves the algebraic rows of the models together with the network
equations in one Newton iteration (`AlgebraicNewton`). In REGCA1's clamp, `g_dq` is
`Volim - v` and does not read `dq`, and the Jacobian of the system carries the response of
the network. A model whose rows do not determine their variables throws, naming the rows and
the modes. Events that change a model's data (`FieldModification`) are applied at their
times, and the automata trace again what they read (`retrace_automata!`, the code of the
phasor simulation).

### Governing a State by a Mode

A mode can hold a differential state on a relation that moves, which a constant hold cannot
state. The declaration `Mode(:upper; holds = [:gate => 0], states = [:y => :upper_bound])`
governs `y` while the model is in `upper`. The right side is a number or a value that the
model computes in `calc_automaton_values!`. A bounded state lowers to such modes (`free`,
`upper`, `lower`, `pinned`) with the bound as the right side; the bound may be a parameter, a
region's constant or a state that moves. `Mode.relations`, a branch current's relation, is
refused.

In a mode that governs `y`, the residual kernel must give `y` the rate of its right side,
`Σ_j (∂r/∂x_j) f_j` over the states it reads, gated by what the mode holds (the gate and the
bound's rate as discrete quantities). Powerful derives no rate itself. The multistep methods
preserve a linear invariant, so the state stays on the relation (`y - r = 0` is affine in the
states) to the tolerance of the integration, as a state-only constraint reduction does
(`StateConstraintReduction`, `adapter.jl`).

A governed state is a differential state with mass at every instance. Its right side must be
affine in the differential states of its own model with constant coefficients and must read
nothing else that evolves. Powerful checks this at setup from the trace of the values
(`ValueReads`) and the degree of each value as a polynomial in the variables and the time
(`PolyDegree`), and checks it again after an event that changes the parameters a branching
kernel reads (`retrace_automata!`). The right side must satisfy the following conditions:

- The degree is at most 1 in the variables and the time, so that its coefficients read no
  variable (`s * w` and `t * s` are refused).
- It reads no algebraic variable, variable without mass, quantity of another model, time,
  discrete quantity or memory (`gate * s` has a coefficient that evolves).
- It reads neither the governed state itself nor a state governed by any mode of the model
  (another mode or another automaton), because such a state moves by its relation.
- The parameters it reads are finite where the simulation starts (a coefficient that is not
  finite is refused, naming the right side).

`check_automata` refuses the declarations that need no trace: a name that is not a
differential state of the model, a right side that is not a finite number or a value, and a
state governed by the modes of two automata (two modes active together). Each refusal names
the model, the mode, the state and the reason. `ModelMetadata` has no place for invariants
such as the crossing of moving bounds (`x_bounds_ordered`), and the code generator rejects a
model that declares one before it reaches Powerful.

The transition that enters a mode writes to each state it governs the value of its right
side, at that instant (`project_states!` with the `AlongStates` policy: the governed state
moves and the states the right side reads stay). The function `enforce_state_constraints!`,
which restores the state-only constraints of a reduced DAE after a fault clears, applies the
same routine with the `MinimumNorm` policy over every state of a constraint. A crossing located on
`x - r` lands on the relation only to the located root, so the write is made whether or not it
moves the state. After the write, inside the same settle:

- The values that read the state are decided by their values and lose the sides that a located
  crossing gave them (as for a discrete quantity that a transition changed).
- The values are computed again and the guards decided again, so that an immediate release
  (`upper` to `free`) keeps the written state.
- An algebraic residual that reads the state (traced with the governed states as inputs of
  the residuals, or a variable of another model that is the state) ends the round there, so
  that the network is reinitialized before the next round, as for a mode that changes the
  algebraic equations.

The extremum search of the automaton callback and the continuity of values held at zero take
the state after the write. The values are measured from where the instant leaves them, and
the held ones are decided from their rates in the mode entered, evaluated at the written
state. Only a settle at an instant, reached through a transition, writes. A step the
integrator rejects, a root-finding evaluation and the bracketing steps of the fixed-step
integration only evaluate the components of the callback.

A start refuses a state that is off its relation, in an initial mode or in a mode the first
settle enters, and the error names the state, the relation and the residual. The right side
is computed from the model and cannot be moved to the start, as a computed bound of a
`VariableLimit` cannot be. A bound that is a model field is moved by
`TDSConfig(adjust_limits = true)` before the automata are set up, so the right side is read
after the move.

Two transitions that hold at once from one mode of an automaton are refused: the settle
throws, naming both. The conditions of the transitions from a mode must therefore exclude
one another. An instance takes part in at most `max_rounds` rounds at one instant (100 unless
the runtime's value is set); a round still needed after the budget is spent throws, and a
a count equal to the budget does not.

## How the Core Uses a Model

All integration with the core goes through metadata queries, and the core contains no
code for particular model types. The function `setup_with_formulation!` reads
`model_metadata` for the formulation of the problem to assign addresses. `SystemBuffers`
creates a `ResidualBuffers` per model, with buffer sizes from the instance counts and buffer
fields from the declarations. `evaluate_system_residuals!` calls each model's
`calc_residuals!` in the order given by `evaluation_pass`. Jacobians are computed with
ForwardDiff (Dual numbers) or Enzyme through DifferentiationInterface; with ForwardDiff the
workspace buffers hold `Dual` numbers, so kernels must be generic over the element type. The
`SystemModel` constructor reads the `params` declarations and converts with bases from
ratings and external parameters. `collect_model_table` calls `compute_observations`, and the
`PUType` of each observation selects its conversion to report units.

## Minimal Model Template

```julia
import PowerfulCore: model_metadata, calc_residuals!, get_terminals, get_key_name

# 1. Struct
@kwdef mutable struct MyModel{Tv} <: AbstractModel
    i::Int32
    id::String3
    status::Bool
    param1::Tv
end

# 2. Registration
@register_model MyModel

# 3. Keys and terminals
get_key_name(::Type{<:MyModel}) = [:i, :id]
const MY_TERMINALS = (Terminal(:Bus, :i),)
get_terminals(::Type{<:MyModel}) = MY_TERMINALS

# 4. Variables
const MY_VARS = [
    ModelVar(:v, AlgebVar(), :Bus, :v, :i; description="Bus voltage magnitude", units="pu"),
]

# 5. Residuals
const MY_RESIDUALS = [
    ModelResidual(:p, AlgebRes(), SharedRes(), :Bus, :p, :i;
                  description="Active power injection"),
]

# 6. Metadata
function model_metadata(::Type{MyModel})
    ModelMetadata(
        name = :MyModel,
        vars = MY_VARS,
        observable_vars = Symbol[],
        residuals = MY_RESIDUALS,
        params = MY_PARAMS,
        ext_params = MY_EXT_PARAMS,
        status_field = :status,
    )
end

# 7. Kernel
@inline function my_kernel!(outputs, params, vars, active, ::Type{<:Vector})
    @inbounds @simd for i in eachindex(outputs.p)
        mask = active[i]
        outputs.p[i] = mask * my_physics(params.param1[i], vars.v[i])
    end
    nothing
end

# 8. Adapter
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    model::StructArray{<:MyModel},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    my_kernel!(
        (p = ws.res_vals.p,),
        (param1 = model.param1,),
        (v = ws.var_vals.v,),
        active_mask,
        AT
    )
end
```
