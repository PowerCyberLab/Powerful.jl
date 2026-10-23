# ============================================================================
# Interface Documentation - Model Contracts (Documentation Only)
# ============================================================================
#
# This file documents the interface contracts for Powerful.jl models.
# Following Julia's AbstractArray pattern, these are documentation-only
# interfaces with no compile-time verification.
#
# Models should implement the required methods documented here.
# Interface violations will be caught at runtime via dispatch errors.
#
# ============================================================================

"""
    PowerSystemComponent

Base interface for all power system models.
Covers fundamental model identification and metadata.

# Required Methods

Models implementing this interface **must** define:

- `model_metadata(::Type{T}) -> ModelMetadata` where `T <: AbstractModel`

  Return ModelMetadata describing model structure. Every model MUST implement this.
  It declares:
  - Variables (state and algebraic)
  - Residuals (equations)
  - Parameters (for per-unit conversion)
  - Observations (computed outputs)
  - External parameters (data from other models)

- `get_key_name(::Type{T}) -> Vector{Symbol}` where `T <: AbstractModel`

  Return composite key field names for model instances.

  Examples:
  - Bus: `[:i]` (bus number only)
  - PQ: `[:i, :id]` (bus number + load id)
  - Branch: `[:i, :j, :ckt]` (from bus, to bus, circuit id)

# Optional Methods

- `supports_format(::Type{T}, ::Type{F}) -> Union{FormatSupport{F}, NoFormatSupport}`
  where `{T <: AbstractModel, F <: DataFormat}`

  Check if model can be parsed from a data format (PSSE, PSSEDYR, etc.).
  Returns `FormatSupport{F}()` if supported, `NoFormatSupport()` otherwise.
  Used by parser to discover which models to parse from each file type.

- `get_status_field(::Type{T}) -> Union{Symbol, Nothing}` where `T <: AbstractModel`

  Return status field name if model has device status tracking.
  Returns `:status` or `nothing`.
  Status fields indicate if device is in-service (true/1) or out-of-service (false/0).

# See Also

- [`PowerFlowModel`](@ref): Interface for static power flow models
- [`DynamicModel`](@ref): Interface for time-domain simulation models
"""
PowerSystemComponent = nothing  # Marker for documentation

"""
    PowerFlowModel <: PowerSystemComponent

Interface for models participating in power flow analysis.
Includes static loads, generators, network elements.

Inherits all requirements from [`PowerSystemComponent`](@ref).

# Required Methods

In addition to `PowerSystemComponent` requirements, models implementing this interface **must** define:

- `calc_residuals!(ws::ResidualBuffers, model::StructArray, external::NamedTuple, formulation::FormulationTrait)`

  Calculate residual contributions for power flow equations.

  Models implement this for each formulation they support.
  Julia's multiple dispatch selects the right method based on formulation type.

  Common formulations:
  - `PowerInjection`: P-θ, Q-V (polar coordinates)
  - `CurrentInjection`: I_r-V_r, I_i-V_i (rectangular)
  - `DCPowerFlow`: P-θ only (DC approximation)

# Optional Methods

- `supports_formulation(::Type{T}, ::Type{F}) -> Bool` where `{T <: AbstractModel, F <: FormulationTrait}`

  Query if model supports a specific formulation.
  Optional trait for capability introspection.
  If not implemented, caller must rely on Julia's method dispatch.

- `compute_observations(sol::PowerFlowSolution, ::Type{T}, var::Symbol) -> Vector{Float64}` where `T <: AbstractModel`

  Compute derived outputs from solution (vectorized).

  Required only if model declares observations in metadata.
  Examples: power flows in branches, reactive power from PV generators.

  Returns `Vector{Float64}` with values for all model instances.

- `get_terminals(::Type{T}) -> Tuple` where `T <: AbstractModel`

  Declare terminal connections to network vertices (buses/nodes).

  Returns tuple of `Terminal` objects specifying connectivity.
  Used for topology analysis and external parameter retrieval.

  Examples:
  - Load: `(Terminal(:Bus, :i),)` - single attachment
  - Branch: `(Terminal(:Bus, :i), Terminal(:Bus, :j))` - two terminals

# Implementation Requirements

Models implementing `PowerFlowModel` should have:
- Non-empty residuals in their metadata (at least one algebraic equation)
- Implementation of `calc_residuals!` for at least one formulation

# See Also

- [`LoadModel`](@ref): Semantic marker for load models
- [`NetworkElement`](@ref): Semantic marker for network topology elements
"""
PowerFlowModel = nothing  # Marker for documentation

"""
    DynamicModel <: PowerSystemComponent

Interface for models with differential equations.
Includes synchronous machines, controllers, dynamic loads.

Inherits all requirements from [`PowerSystemComponent`](@ref).

# Required Methods

In addition to `PowerSystemComponent` requirements, models implementing this interface **must** define:

- `calc_residuals!(ws::ResidualBuffers, model::StructArray, external::NamedTuple, formulation::FormulationTrait)`

  Calculate DAE residuals (differential + algebraic).

  Computes:
  - State derivatives: f(x, y) for dx/dt = f(x, y)
  - Algebraic equations: g(x, y) = 0

  Note: Uses same 4-parameter signature as `PowerFlowModel` for consistency.
  TDS currently uses `PowerInjection` formulation, but signature is standardized.

- `initialize_dynamic_states(model::StructArray{<:AbstractModel}, links::ModelLinks, source_model::Symbol, sys::SystemModel, pf_sol::PowerFlowSolution) -> NamedTuple`

  Initialize state variables from power flow solution (batch/vectorized).

  Takes power flow solution and model links, computes initial conditions
  for differential states (e.g., rotor angle, speed, flux linkages) for
  all instances at once using vectorized data gathering.

  Returns `NamedTuple` with initial values for all state variables.

# Optional Trait Methods

- `requires_initialization(::Type{T}) -> Bool` where `T <: AbstractModel`

  Trait, declared by every model (no default): `true` if `initialize_dynamic_states`
  should be called. Models without states of their own (e.g., Slack) declare `false`.

- `get_dyr_model_name(::Type{T}) -> Union{String, Nothing}` where `T <: AbstractModel`

  Return DYR model name string for PSS/E DYR file parsing.
  Examples: "GENCLS", "GENROU", "EXST1"
  Returns `nothing` if model is not parseable from DYR files.
  Used by DYR parser to auto-discover available models.

- `get_replacement_map(::Type{T}) -> Dict{Symbol, NamedTuple}` where `T <: AbstractModel`

  Declare which static models this replaces in TDS.

  Returns `Dict{Symbol, NamedTuple}` mapping:
  - Keys: static model types (`:PV`, `:Slack`)
  - Values: `(source_fields=(...), target_fields=(...))`

  Example: `Dict(:PV => (source_fields=(:BUS, :ID), target_fields=(:i, :id)))`
  Used for linking DYR models to power flow generators.

- `get_data_inheritance(::Type{T}) -> Dict{Symbol, Dict{Symbol, Symbol}}` where `T <: AbstractModel`

  Declare data to inherit from linked power flow models.

  Returns `Dict{Symbol, Dict{Symbol, Symbol}}`:
  - Outer keys: source model types (`:PV`, `:Slack`)
  - Inner dicts: destination_field => source_field mappings

  Example: `Dict(:PV => Dict(:xd_p => :zx))  # Inherit PV.zx → GENCLS.xd_p`

- `parse_dynamic_model(::Type{T}, dyr_data, sys_pf::SystemModel, links::Union{ModelLinks, Nothing}) -> StructArray{T}` where `T <: AbstractModel`

  Parse dynamic model from DYR data.

  **For replacement models** (GENCLS, GENROU): Pass `ModelLinks` to inherit data from PV/Slack.
  **For controller models** (TGOV1, EXST1): Pass `nothing` - linking deferred to TDS setup via group reference terminals.

  Returns `StructArray{T}` with all model instances.

# Implementation Requirements

Models implementing `DynamicModel` should have:
- At least one state variable in their metadata
- Implementation of both `calc_residuals!` and `initialize_dynamic_states`
- Appropriate trait methods (e.g., `requires_initialization` returning `true`)

# See Also

- [`GeneratorModel`](@ref): Semantic marker for generator models
"""
DynamicModel = nothing  # Marker for documentation

# ============================================================================
# Semantic Marker Interfaces (Optional)
# ============================================================================

"""
    GeneratorModel <: DynamicModel

Semantic marker for generator models (GENCLS, GENROU, GENSAL, etc.).

Inherits all requirements from [`DynamicModel`](@ref).
No additional requirements - just semantic grouping for dispatch and filtering.
"""
GeneratorModel = nothing  # Marker for documentation

"""
    LoadModel <: PowerFlowModel

Semantic marker for load models (PQ, ZIP, etc.).

Inherits all requirements from [`PowerFlowModel`](@ref).
Can be used for semantic dispatch or filtering.
"""
LoadModel = nothing  # Marker for documentation

"""
    NetworkElement <: PowerFlowModel

Semantic marker for network topology elements (Branch, Transformer, Shunt).

Inherits all requirements from [`PowerFlowModel`](@ref).
Useful for topology analysis and filtering network components.
"""
NetworkElement = nothing  # Marker for documentation

# ============================================================================
# Model Traits Implementation
# ============================================================================

import Base: length

# model_metadata interface is defined in model_metadata.jl
# Individual models implement model_metadata(::Type{ModelType}) in their source files

"""
    length(model::AbstractModel)

Returns the number of devices represented by this model by checking the length of
any of its data fields. All fields in a model must have equal lengths.
"""
function length(model::AbstractModel)
    fields = fieldnames(typeof(model))
    isempty(fields) && return 0
    return length(getfield(model, first(fields)))
end

# ============================================================================
# Time-Domain Simulation (TDS) Traits
# ============================================================================

"""
    requires_initialization(::Type{<:AbstractModel}) -> Bool

Whether the time-domain initialization calls the model's `initialize_dynamic_states`, which
gives its states and setpoints: `true` for the machines, controls and converters; `false` for a
model with no state of its own (the power flow models, whose variables are the bus voltages of
the power flow) or whose states another initializer owns entirely. Every model declares it: there is no default, since a model with
states that is not initialized would start them at zero. The initialization refuses a system
holding a model without a declaration, naming it, before it initializes any model.
"""
function requires_initialization end

"""
    declares_requires_initialization(::Type{<:AbstractModel}) -> Bool

Whether a model type declares whether it requires initialization (a method of
`requires_initialization`).
"""
declares_requires_initialization(::Type{M}) where {M<:AbstractModel} =
    hasmethod(requires_initialization, Tuple{Type{M}})

"""
    check_parameter_domains(x::M, from::AbstractString) where M <: AbstractModel -> M

Check that every parameter of the instance `x` lies in its valid domain (ranges, ordering between
limits, divisors that must not be zero), throwing `ArgumentError` naming `from` (where the values
came from, e.g. "keyword construction" or a `FieldModification`'s description) if not. Returns `x`
unchanged on success.

Every model declares it: there is no default, since a model without a declaration would let an
event set a parameter to an invalid value with no check. A model with no domain constraints of its
own declares `check_parameter_domains(x::M, ::AbstractString) = x` and says so in a comment; one
with constraints (typically generated with the model from its DSL source) checks them.
"""
function check_parameter_domains end

"""
    declares_check_parameter_domains(::Type{M}) where M <: AbstractModel -> Bool

Whether a model type declares its parameter domains (a method of `check_parameter_domains`
taking an instance of it).
"""
declares_check_parameter_domains(::Type{M}) where {M<:AbstractModel} =
    hasmethod(check_parameter_domains, Tuple{M, AbstractString})

"""
    is_zip_load(::Type{<:AbstractModel}) -> Bool

Whether the instances of a model are ZIP loads that the time-domain simulation's
`apply_zip_conversion!` converts: each carries the (bus `i`, `id`) key of its power flow `PQ`
load and the ZIP coefficients `pl`, `ql`, `ip`, `iq`, `yp`, `yq` in their phasor meaning.
Every model declares it: there is no default. The conversion refuses a system holding a model
without a declaration, naming it, before it converts any load.
"""
function is_zip_load end

"""
    declares_zip_load(::Type{<:AbstractModel}) -> Bool

Whether a model type declares whether it is a ZIP load (a method of `is_zip_load`).
"""
declares_zip_load(::Type{M}) where {M<:AbstractModel} = hasmethod(is_zip_load, Tuple{Type{M}})

"""
    PhaseRepresentation

How a model represents the three phases, as it declares it (`phase_representation`):
[`PositiveSequence`](@ref) or [`ThreePhase`](@ref).
"""
abstract type PhaseRepresentation end

"""
    PositiveSequence()

The model is a balanced, positive-sequence phasor model (power flow, phasor TDS, Y-bus).
"""
struct PositiveSequence <: PhaseRepresentation end

"""
    ThreePhase()

The model represents the three phases a, b, c directly, rather than a single-phase
equivalent.
"""
struct ThreePhase <: PhaseRepresentation end

"""
    phase_representation(::Type{<:AbstractModel}) -> PhaseRepresentation

How a model represents the three phases. Every model declares it: there is no default. A
simulation of one representation refuses a model of the other by name (e.g. the Y-bus
simulation, positive sequence, refuses a three-phase model).
"""
function phase_representation end

"""
    declares_phase_representation(::Type{<:AbstractModel}) -> Bool

Whether a model type declares its phase representation (a method of `phase_representation`).
"""
declares_phase_representation(::Type{M}) where {M<:AbstractModel} = hasmethod(phase_representation, Tuple{Type{M}})

"""
    GPUSupport

Whether a model runs on the GPU, as it declares it (`gpu_support`): [`SupportedOnGPU`](@ref)
or [`UnsupportedOnGPU`](@ref).
"""
abstract type GPUSupport end

"""
    SupportedOnGPU()

The model runs on the GPU: a GPU package has its kernels (or the model has no equations, and
its data transfers).
"""
struct SupportedOnGPU <: GPUSupport end

"""
    UnsupportedOnGPU(reason)

The model does not run on the GPU, for `reason`: `to_gpu` refuses a system holding instances
of it, naming it and `reason` (`check_gpu_support`).
"""
struct UnsupportedOnGPU <: GPUSupport
    reason::String
end

"""
    gpu_support(::Type{<:AbstractModel}) -> GPUSupport

Whether a model runs on the GPU. Every model declares it: there is no default. A model is
declared supported by the GPU package that has its kernels, and unsupported, with the
reason, by the package that defines it. `to_gpu` refuses a system holding a model without
a declaration or with instances of an unsupported one, naming it (`check_gpu_support`),
before any transfer.
"""
function gpu_support end

"""
    declares_gpu_support(::Type{<:AbstractModel}) -> Bool

Whether a model type declares its GPU support (a method of `gpu_support`).
"""
declares_gpu_support(::Type{M}) where {M<:AbstractModel} = hasmethod(gpu_support, Tuple{Type{M}})

"""
    get_dyr_model_name(::Type{<:AbstractModel}) -> Union{String, Nothing}

Return the PSS/E DYR model name string (e.g., "GENCLS", "GENROU").

DYR-parseable models override this trait. Default: `nothing`.
Used by the DYR parser to map model names to Julia types.
"""
get_dyr_model_name(::Type{<:AbstractModel}) = nothing

"""
    get_replacement_map(::Type{<:AbstractModel}) -> Dict{Symbol, NamedTuple}

Declare which static models this dynamic model replaces and how to match instances.

Returns `Dict{Symbol, NamedTuple}` where keys are model symbols (`:PV`, `:Slack`)
and values specify field correspondence for linking. Default: empty dict.
"""
function get_replacement_map(::Type{<:AbstractModel})
    return Dict{Symbol, NamedTuple}()
end

"""
    get_data_inheritance(::Type{<:AbstractModel}) -> Dict{Symbol, Dict{Symbol, Symbol}}

Declare which fields to copy from linked power flow models during construction.

Returns `Dict{model_symbol => Dict{dest_field => src_field}}`.
E.g., `Dict(:PV => Dict(:xd_p => :zx))` copies `PV.zx` → `GENCLS.xd_p`. Default: empty.
"""
function get_data_inheritance(::Type{<:AbstractModel})
    return Dict{Symbol, Dict{Symbol, Symbol}}()
end

# Note: get_controlled_group is defined in groups.jl

"""
    find_control_terminal(::Type{M}) where M <: AbstractModel -> Union{Terminal, Nothing}

Convenience helper to locate the terminal on a controller that links to its
`get_controlled_group`. Returns `nothing` when the model is not a controller
or does not declare a matching terminal.
"""
function find_control_terminal(::Type{M}) where M <: AbstractModel
    group = get_controlled_group(M)
    isnothing(group) && return nothing
    terminals = get_terminals(M)
    idx = findfirst(t -> t.vertex_model == group, terminals)
    return isnothing(idx) ? nothing : terminals[idx]
end

"""
    parse_dynamic_model(::Type{M}, dyr_data, sys_pf::SystemModel, links::ModelLinks) -> StructArray{M}

Parse dynamic model from DYR data with linking to power flow models.

- `dyr_data`: Raw DYR parameters (StructVector from PowerDynData)
- `sys_pf`: Power flow system (for data inheritance via `get_data_inheritance`)
- `links`: Instance mappings between DYR and PF models

Dynamic models must implement this trait. Default throws error.
"""
function parse_dynamic_model(
    ::Type{M},
    dyr_data,
    sys_pf::SystemModel,
    links::ModelLinks
) where M
    error("parse_dynamic_model not implemented for $M")
end
