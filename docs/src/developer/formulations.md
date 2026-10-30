# Multi-Physics Architecture

A formulation is a choice of coordinates and balance equations for the network. The same physical system can be solved in several formulations, which share the physics and differ in their numerical properties. The formulation types are defined in `PowerfulCore/src/foundation/formulation_traits.jl`.

## Available Formulations

`PowerInjection` is the default formulation. Its coordinates are polar (θ, V), and its equations are the active and reactive power balance at each bus, P = 0 and Q = 0. The Jacobian contains trigonometric terms.

`CurrentInjection` uses rectangular coordinates (V_r, V_i) and the real and imaginary current balance at each bus, I_r = 0 and I_i = 0. The network current I = Y·V is linear in V, so the Jacobian has no trigonometric terms.

`DCPowerFlow` uses bus angles as its only variables. Its equations are the active power balance with the linearized branch flows P_ij ≈ (θ_i - θ_j) / x. The Jacobian is constant, because the system is linear. This formulation suits fast approximations and screening studies.

## Formulations as Types

Each formulation is a Julia type:

```julia
abstract type FormulationTrait end
struct PowerInjection <: FormulationTrait end
struct CurrentInjection <: FormulationTrait end
struct DCPowerFlow <: FormulationTrait end
```

The compiler selects the method of a function from the formulation argument:

```julia
calc_residuals!(ws, pq, external, PowerInjection(); active_mask)    # P, Q
calc_residuals!(ws, pq, external, CurrentInjection(); active_mask)  # I_r, I_i
```

Dispatch on the formulation type needs no runtime check of the form `if formulation == ...` and gives each path specialized machine code. A new formulation is added by defining a type and implementing the methods for it, which keeps conditionals out of the rest of the code base. Function signatures also show which code is formulation-specific.

## How Models Declare Multi-Formulation Support

### Default Variables and Residuals

The `vars` and `residuals` passed to `ModelMetadata` are the default (PowerInjection)
declarations. From `PowerfulPF/src/models/pq.jl`:

```julia
const PQ_VARS = [
    ModelVar(:theta, AlgebVar(), :Bus, :theta, :i;
             description="Bus voltage phase angle", units="rad"),
    ModelVar(:v, AlgebVar(), :Bus, :v, :i;
             description="Voltage magnitude", units="pu"),
]

const PQ_RESIDUALS = [
    ModelResidual(:p, AlgebRes(), SharedRes(), :Bus, :p, :i,
                  description="Active power injection from load"),
    ModelResidual(:q, AlgebRes(), SharedRes(), :Bus, :q, :i,
                  description="Reactive power injection from load"),
]
```

### Formulation-Specific Overrides

```julia
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

### Registration in Metadata

```julia
function model_metadata(::Type{PQ})
    metadata = ModelMetadata(
        name = :PQ,
        vars = PQ_VARS,              # Default (PowerInjection)
        residuals = PQ_RESIDUALS,
        # ... observable_vars, params, observations, ratings, ext_params, status_field
    )

    # Register overrides, keyed by formulation type
    metadata.vars_by_formulation[CurrentInjection] = PQ_VARS_CURRENT
    metadata.residuals_by_formulation[CurrentInjection] = PQ_RESIDUALS_CURRENT
    metadata.vars_by_formulation[DCPowerFlow] = PQ_VARS_DC
    metadata.residuals_by_formulation[DCPowerFlow] = PQ_RESIDUALS_DC

    return metadata
end

supports_formulation(::Type{PQ}, ::Type{PowerInjection}) = true
supports_formulation(::Type{PQ}, ::Type{CurrentInjection}) = true
supports_formulation(::Type{PQ}, ::Type{DCPowerFlow}) = true
```

`get_metadata_for_formulation(metadata, formulation)` returns the metadata with the
overrides for that formulation applied. A model without an override for a formulation
uses its default variables and residuals, and still needs a `calc_residuals!` method
for that formulation.

## Formulation Affects Variable Layout

The same system has different variables in each formulation. For the IEEE 14-bus case with natural addressing, `PowerInjection` has 28 variables:

```
θ[1..14] at addresses 1-14
V[1..14] at addresses 15-28
```

`CurrentInjection` also has 28 variables, with a different meaning:

```
V_r[1..14] at addresses 1-14
V_i[1..14] at addresses 15-28
```

Because the variable sets differ, the `AddressManager` is formulation-specific. It is created for each problem and stored in the problem, and the `SystemModel` holds no addresses.

## Formulation-Specific Residual Dispatch

The system workspace calls the following method for every model:

```julia
calc_residuals!(ws::ResidualBuffers, model::StructArray, external::NamedTuple,
                formulation::FormulationTrait, t = zero(T); active_mask)
```

- `external` holds the model's external parameters (`sys.external[model_name]`)
- `t` is the simulation time (zero in power flow)
- `active_mask` marks the in-service instances. It is passed with each call, so parallel
  contingency analysis can use thread-local status arrays

Multiple dispatch selects the method for the formulation. From `PowerfulPF/src/models/pq.jl`:

```julia
# PowerInjection implementation
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pq::StructArray{<:PQ},
    _external::NamedTuple,
    ::PowerInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    pq_residual_kernel!(
        (p = ws.res_vals.p, q = ws.res_vals.q),
        (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq,
         vmin = ws.system.constant_power_vmin),
        (v = ws.var_vals.v,),
        active_mask,
        AT
    )
end

# CurrentInjection implementation
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    pq::StructArray{<:PQ},
    _external::NamedTuple,
    ::CurrentInjection,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    pq_current_kernel!(  # Different kernel, different physics
        (i_r = ws.res_vals.i_r, i_i = ws.res_vals.i_i),
        (pl = pq.pl, ql = pq.ql, ip = pq.ip, iq = pq.iq, yp = pq.yp, yq = pq.yq,
         vmin = ws.system.constant_power_vmin),
        (v_r = ws.var_vals.v_r, v_i = ws.var_vals.v_i),
        active_mask,
        AT
    )
end
```

The formulation type is known when the workspace closures are built, so the code path is selected at compile time.

### Evaluation Order

Some models read residuals accumulated by others. In CurrentInjection, the `SharedRes(:set)` residuals of the PV generator convert the current accumulated at its bus by branches and loads into a power balance. PV declares `evaluation_pass = 1`, so it is evaluated after the default pass 0, and its workspace gathers the accumulated residuals (`gather_residuals!`) before `calc_residuals!` runs.

## Adding a New Formulation

### 1. Define the Formulation Trait
```julia
# In PowerfulCore/src/foundation/formulation_traits.jl
struct MyFormulation <: FormulationTrait end
```

### 2. Specify Variables for Models
```julia
# For models needing different variables
const MY_MODEL_VARS = [
    ModelVar(:my_var1, AlgebVar(), :Bus, :my_var1, :i; description="..."),
    ModelVar(:my_var2, AlgebVar(), :Bus, :my_var2, :i; description="..."),
]

# Register in model_metadata
metadata.vars_by_formulation[MyFormulation] = MY_MODEL_VARS
```

The bus model owns the network variables, so its metadata declares the new internal variables (e.g. `ModelVar(:my_var1, AlgebVar(); ...)`), and devices reference them.

### 3. Specify Residuals for Models
```julia
const MY_MODEL_RESIDUALS = [
    ModelResidual(:res1, AlgebRes(), SharedRes(), :Bus, :res1, :i; description="..."),
    ModelResidual(:res2, AlgebRes(), SharedRes(), :Bus, :res2, :i; description="..."),
]

metadata.residuals_by_formulation[MyFormulation] = MY_MODEL_RESIDUALS
```

### 4. Implement Residual Kernels
```julia
@inline function calc_residuals!(
    ws::ResidualBuffers{T, AT},
    model::StructArray{<:MyModel},
    _external::NamedTuple,
    ::MyFormulation,
    t = zero(T);
    active_mask::AbstractVector{Bool}
) where {T, AT}
    my_formulation_kernel!(
        (res1 = ws.res_vals.res1, res2 = ws.res_vals.res2),  # outputs
        (param1 = model.param1,),                             # params
        (my_var1 = ws.var_vals.my_var1, my_var2 = ws.var_vals.my_var2),  # vars
        active_mask,
        AT
    )
end
```

Follow the kernel pattern in [Implementing Models](models.md).

### 5. Test
```julia
prob = ACPowerFlowProblem(sys, MyFormulation())
sol = solve(prob, Newton())

# Compare against PowerInjection where the formulations are mathematically equivalent
# Check the Jacobian structure and convergence
```

The `AddressManager` and the workspaces allocate addresses and buffers for the new formulation from the declared metadata.

## Choosing a Formulation

`PowerInjection` is the standard formulation for transmission power flow. Its residuals are power mismatches, the quantities engineers reason about.

`CurrentInjection` uses rectangular voltages and has no trigonometric functions in the network equations. It is an alternative when `PowerInjection` has convergence difficulty.

`DCPowerFlow` gives fast approximate active-power solutions and suits screening studies with many scenarios. Voltage magnitudes are not variables, and reactive power and losses are ignored.

## Example: Switching Formulations

```@example formulations
using Powerful, PowerfulCases

sys = SystemModel(parse_network(load("ieee14").raw))

sol1 = solve(ACPowerFlowProblem(sys, PowerInjection()), Newton())
sol2 = solve(ACPowerFlowProblem(sys, CurrentInjection()), Newton())

# The formulations are mathematically equivalent. CurrentInjection's bus
# variables are v_r and v_i, so the magnitude is computed from them.
v1 = get_var(sol1, :Bus, :v, (14,))
v2 = hypot(get_var(sol2, :Bus, :v_r, (14,)), get_var(sol2, :Bus, :v_i, (14,)))
(v1, v2)
```

Switching the formulation requires a new problem only; the model code is unchanged.
