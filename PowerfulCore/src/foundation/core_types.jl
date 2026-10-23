"""
    AbstractModel

Supertype of every model struct (buses, branches, loads, generators, controllers).

A model type describes a single device instance; a system stores all instances of
one model as a `StructArray`. Behavior is attached through traits such as
[`model_metadata`](@ref), [`calc_residuals!`](@ref), [`get_terminals`](@ref) and
[`get_key_name`](@ref), and the type is made known to the system with
[`@register_model`](@ref PowerfulCore.@register_model).
"""
abstract type AbstractModel end

# ============================================================================
# Virtual Bus Numbering
# ============================================================================
"""
Base number for virtual buses created from 3-winding transformer conversions.
Virtual buses are numbered VIRTUAL_BUS_BASE + 1, VIRTUAL_BUS_BASE + 2, ...
This avoids conflicts with real power system bus numbers (typically < 100,000).
"""
const VIRTUAL_BUS_BASE = 1_000_000

# ============================================================================
# Variable Traits
# ============================================================================
"""
    VarProperty

Supertype of the optional properties attached to a [`ModelVar`](@ref):
[`Description`](@ref), [`Units`](@ref) and
[`MassCoefficient`](@ref). Properties are stored in a [`PropertyDict`](@ref) keyed by
their type and read with [`get_property`](@ref).
"""
abstract type VarProperty end

"""
    Description(value::String)

Human-readable description of a model variable ([`VarProperty`](@ref)); read with
[`get_description`](@ref).
"""
struct Description <: VarProperty
    value::String
end

"""
    Units(value::String)

Unit label of a model variable, such as `"pu"` or `"rad"` ([`VarProperty`](@ref));
read with [`get_units`](@ref).
"""
struct Units <: VarProperty
    value::String
end

"""
    MassCoefficient{F}

Mass coefficient property for differential state variables in DAE systems.
Specifies how to obtain the mass coefficient that multiplies the time derivative: M*dx/dt = f(x).

# Type Parameter
- `F`: Either `Symbol` (field lookup) or a `Function` type (computed value)

# Field
- `source::F`: Source of mass coefficient values - either:
  - `Symbol`: Field name in model struct (e.g., `:M`) for direct lookup
  - `Function`: Callable `f(models) -> Vector` that computes mass values from model data

# Constructors
```julia
MassCoefficient(:M)                    # Field lookup: models.M
MassCoefficient(m -> 2.0 .* m.H)       # Function: compute M = 2H inline
```

# Usage in ModelVar
```julia
# Direct field lookup (backward compatible)
ModelVar(:omega, StateVar(); mass_coefficient=:M, ...)

# Function-based (AD-compatible, for derived coefficients)
ModelVar(:omega, StateVar(); mass_coefficient=MassCoefficient(m -> 2.0 .* m.H), ...)
```

# AD Compatibility
Function-based mass coefficients are fully AD-compatible. When model parameters are Dual
numbers (during Jacobian computation), the function receives Dual arrays and returns Duals,
preserving gradient information. Example: `∂M/∂H = 2.0` flows through `m -> 2.0 .* m.H`.

# See Also
- [`eval_mass_coeff`](@ref): Type-stable extraction of mass values
- [`get_mass_coeff`](@ref): Retrieve MassCoefficient from ModelVar
"""
struct MassCoefficient{F} <: VarProperty
    source::F
end

# Convenience constructors
MassCoefficient(field::Symbol) = MassCoefficient{Symbol}(field)
MassCoefficient(func::F) where {F<:Function} = MassCoefficient{F}(func)

"""
    eval_mass_coeff(mc::MassCoefficient, models) -> AbstractVector

Type-stable dispatch for extracting mass values from a MassCoefficient.
Uses compile-time dispatch on `F` type parameter instead of runtime `isa` checks.

# Arguments
- `mc::MassCoefficient{F}`: Mass coefficient specification
- `models`: StructArray of model instances

# Returns
`AbstractVector{T}` of mass coefficient values, one per model instance.
When `models` contains Dual numbers, returns Dual vector (AD-compatible).

# Methods
- `MassCoefficient{Symbol}`: Returns `getproperty(models, mc.source)` (field lookup)
- `MassCoefficient{<:Function}`: Returns `mc.source(models)` (function call)

# Example
```julia
models = StructArray(; H=[5.0, 10.0], D=[1.0, 1.0])

# Symbol-based: direct field access
mc1 = MassCoefficient(:H)
eval_mass_coeff(mc1, models)  # [5.0, 10.0]

# Function-based: computed value
mc2 = MassCoefficient(m -> 2.0 .* m.H)
eval_mass_coeff(mc2, models)  # [10.0, 20.0]
```
"""
@inline eval_mass_coeff(mc::MassCoefficient{Symbol}, models) = getproperty(models, mc.source)
@inline eval_mass_coeff(mc::MassCoefficient{F}, models) where {F<:Function} = mc.source(models)

"""
    PropertyDict

Alias for `Dict{Type{<:VarProperty},VarProperty}`: the properties of a
[`ModelVar`](@ref), keyed by property type so each kind appears at most once.
"""
const PropertyDict = Dict{Type{<:VarProperty},VarProperty}

"""
    VarTrait

Supertype of the ownership traits of a model variable or residual:
[`Internal`](@ref) (allocated by the model itself) or [`External`](@ref)
(owned by another model and referenced through a foreign key).
"""
abstract type VarTrait end

"""
    Internal <: VarTrait

Trait of a variable or residual owned by the model that declares it. The address
manager allocates its addresses when the model is set up.
"""
struct Internal <: VarTrait end

"""
    External <: VarTrait

Trait of a variable or residual owned by another model (for example the bus voltage
read by a load). Its addresses are retrieved from the owning model through the
declared source model, source name and reference key.
"""
struct External <: VarTrait end

"""
    AddressableType

Supertype of the address spaces of the DAE system. Variables live in
[`AddressableVar`](@ref) spaces and equations in [`AddressableRes`](@ref) spaces;
the [`AddressManager`](@ref) numbers each space separately.
"""
abstract type AddressableType end

"""
    AddressableVar <: AddressableType

Supertype of the variable address spaces: [`AlgebVar`](@ref), [`StateVar`](@ref)
and [`ObservableVar`](@ref).
"""
abstract type AddressableVar <: AddressableType end

"""
    AddressableRes <: AddressableType

Supertype of the residual (equation) address spaces: [`AlgebRes`](@ref) and
[`StateRes`](@ref).
"""
abstract type AddressableRes <: AddressableType end

"""
    AlgebVar()

Address space of algebraic variables, stored in the `y` vector of the [`DAE`](@ref PowerfulCore.DAEInfra.DAE).
Bus voltages and generator currents are algebraic variables.
"""
struct AlgebVar <: AddressableVar end

"""
    StateVar()

Address space of differential state variables, stored in the `x` vector of the
[`DAE`](@ref PowerfulCore.DAEInfra.DAE). Rotor angles and controller integrator states are state variables.
"""
struct StateVar <: AddressableVar end

"""
    ObservableVar()

Address space of observables: quantities computed from the variables for reporting
(for example generator electrical power) that take no part in the equations.
"""
struct ObservableVar <: AddressableVar end

"""
    AlgebRes()

Address space of algebraic residuals `g(x, y) = 0`, stored in the `g` vector of the
[`DAE`](@ref PowerfulCore.DAEInfra.DAE).
"""
struct AlgebRes <: AddressableRes end

"""
    StateRes()

Address space of differential residuals, the right-hand sides `f(x, y)` of
`M dx/dt = f(x, y)`, stored in the `f` vector of the [`DAE`](@ref PowerfulCore.DAEInfra.DAE).
"""
struct StateRes <: AddressableRes end

"""
    ModelVar{T<:VarTrait, AT, SM, RK, Props} represents a variable in a power system model

Type parameters:
- T: Variable trait (Internal/External)
- AT: Addressable type (e.g., Algeb, State)
- SM: Source model type or Nothing
- RK: Reference key type or Nothing (fields used to look up foreign reference)
- Props: Type of properties dictionary
"""
struct ModelVar{
    T<:VarTrait,
    AT<:AddressableVar,
    SM<:Union{Symbol,Nothing},
    RK<:Union{Symbol,Vector{Symbol},Nothing},
    Props<:PropertyDict
}
    name::Symbol
    address_type::AT
    source_model::SM
    source_name::Union{Symbol,Nothing}
    ref_key::RK
    properties::Props
end

abstract type ResAccess end

"""
    SharedRes(mode::Symbol = :add)

Residual access pattern for shared address space (e.g., bus power balance).

# Modes
- `:add` - Accumulate with `+=` (default). Used for power balance contributions.
- `:set` - Overwrite with `=`. Used for control equations (e.g., PV voltage control).
         Applied AFTER all `:add` modes to replace accumulated values.

# Examples
```julia
SharedRes()        # Default :add mode for power balance
SharedRes(:add)    # Explicit :add mode
SharedRes(:set)    # Setter mode for control equations
```

# Mutability
Made mutable to support observation mode where SharedRes(:set) is temporarily changed
to SharedRes(:add) to extract generator power contributions.
"""
mutable struct SharedRes <: ResAccess
    mode::Symbol

    # Default constructor with :add mode
    SharedRes(mode::Symbol = :add) = begin
        mode in (:add, :set) || error("SharedRes mode must be :add or :set, got :$mode")
        new(mode)
    end
end

"""
    PrivateRes()

Residual access pattern for an equation owned by a single model instance: each
instance writes its own address and nothing else contributes. Contrast with
[`SharedRes`](@ref).
"""
struct PrivateRes <: ResAccess end

struct ModelResidual{
    T<:VarTrait,
    AT<:AddressableRes,
    RA<:ResAccess,
    SM<:Union{Symbol,Nothing},
    SR<:Union{Symbol,Nothing},
    RK<:Union{Symbol,Vector{Symbol},Nothing}
}
    name::Symbol
    address_type::AT
    access::RA
    source_model::SM
    source_name::SR
    ref_key::RK
    description::String
end

export AddressableType, AddressableVar, AddressableRes
export AlgebVar, StateVar, ObservableVar
export AlgebRes, StateRes
export VarTrait, Internal, External
export VarProperty, Description, Units, MassCoefficient, PropertyDict, eval_mass_coeff
export ModelVar

export SharedRes
export ModelResidual
### ============================= ###

### ========== Address Manager ========== ###

"""
Tracks address assignment status for a specific model type
"""
mutable struct ModelAddressAssignment
    model_name::Symbol
    # maps addressable type => variable name => range
    by_type::Dict{AddressableType,Dict{Symbol,Vector{UInt}}}
    is_complete::Bool  # true if all required variables are allocated
end

# Convenience constructor with default empty dictionaries
function ModelAddressAssignment(model_name::Symbol)
    by_type = Dict{AddressableType, Dict{Symbol,Vector{UInt}}}()
    # Initialize empty dictionaries for each addressable type
    for T in (AlgebVar, AlgebRes, StateVar, StateRes, ObservableVar)
        by_type[T()] = Dict{Symbol,Vector{UInt}}()
    end
    ModelAddressAssignment(model_name, by_type, false)
end

"""
Track retrieved addresses for a specific model instance
"""
struct ModelRetrievedAddresses
    model_name::Symbol
    by_type::Dict{AddressableType,Dict{Symbol,Vector{UInt}}}
end

# Convenience constructor with default empty dictionaries
function ModelRetrievedAddresses(model_name::Symbol)
    by_type = Dict{AddressableType, Dict{Symbol,Vector{UInt}}}()
    for T in (AlgebVar, AlgebRes, StateVar, StateRes, ObservableVar)
        by_type[T()] = Dict{Symbol,Vector{UInt}}()
    end
    ModelRetrievedAddresses(model_name, by_type)
end

"""
Reverse address lookup map: address → (model, variable, instance_index)

Provides efficient reverse lookup from DAE vector addresses to model metadata.
Used for debugging, diagnostics, and validation reporting.

# Fields
- `algeb_vars::Dict{UInt, Tuple{Symbol, Symbol, Int}}`: Algebraic variable addresses
- `state_vars::Dict{UInt, Tuple{Symbol, Symbol, Int}}`: Differential state addresses
- `algeb_res::Dict{UInt, Tuple{Symbol, Symbol, Int}}`: Algebraic residual addresses
- `state_res::Dict{UInt, Tuple{Symbol, Symbol, Int}}`: Differential residual addresses
- `observable_vars::Dict{UInt, Tuple{Symbol, Symbol, Int}}`: Observable variable addresses

Each tuple contains:
- `Symbol`: Model name (e.g., :Bus, :PQ, :GENCLS)
- `Symbol`: Variable/residual name (e.g., :theta, :v, :p, :q)
- `Int`: Instance index in model array (1-based)

# Example
```julia
# After setup_with_formulation!
model, var, inst = am.reverse_map.algeb_res[UInt(42)]
# → (:PQ, :p, 3)  # 3rd PQ load's active power residual
```
"""
struct ReverseAddressMap
    algeb_vars::Dict{UInt, Tuple{Symbol, Symbol, Int}}
    state_vars::Dict{UInt, Tuple{Symbol, Symbol, Int}}
    algeb_res::Dict{UInt, Tuple{Symbol, Symbol, Int}}
    state_res::Dict{UInt, Tuple{Symbol, Symbol, Int}}
    observable_vars::Dict{UInt, Tuple{Symbol, Symbol, Int}}
end

"""
Address reordering strategies for sparse factorization optimization.

System-level address reordering can improve sparse matrix factorization performance
by reducing fill-in and improving cache locality.

# Available Strategies
- `NaturalAddressing`: Sequential addressing in model order (default)
- `AMDAddressing`: Approximate Minimum Degree ordering for minimal fill-in
- `CustomPermutation`: User-provided permutation vector
"""
abstract type AddressingStrategy end

"""
    NaturalAddressing()

Natural sequential addressing in the order models are processed (default).
No reordering is applied.
"""
struct NaturalAddressing <: AddressingStrategy end

"""
    AMDAddressing()

Approximate minimum degree (AMD) reordering of the algebraic addresses, computed
from the Jacobian sparsity pattern when the problem is built. The permutation it
applied is recorded in the problem's `AddressManager` as `applied_permutation`;
to reuse a known ordering, pass it as [`CustomPermutation`](@ref).

# Example
```julia
prob = ACPowerFlowProblem(sys; addressing=AMDAddressing())
perm = prob.address_manager.applied_permutation
```
"""
struct AMDAddressing <: AddressingStrategy end

"""
    CustomPermutation(perm::Vector{Int})

Custom address reordering using a user-provided permutation vector.

# Arguments
- `perm::Vector{Int}`: Permutation of the algebraic address space, where new
  address `k` holds old address `perm[k]` (the convention of `compute_amd_ordering`).
  Its length must equal the number of algebraic addresses; state and observable
  addresses are not permuted (see `remap_addresses!`).
"""
struct CustomPermutation{P<:AbstractVector{Int}} <: AddressingStrategy
    permutation::P
end

"""
Enhanced address manager with explicit address assignment tracking

# Fields
- `addresses`: Maps (AddressableType, model_name, variable_name) → address vector
- `next_idx`: Tracks next available address for each AddressableType
- `assignments`: Tracks assigned addresses per model (internal variables)
- `retrieved`: Tracks retrieved addresses per model (external variables)
- `reverse_map`: Reverse lookup (address → model/var metadata), built during setup
- `strategy`: Address reordering strategy (NaturalAddressing, AMDAddressing, etc.)
- `applied_permutation`: Permutation the strategy applied to the algebraic addresses
  (new address `k` holds natural address `perm[k]`), or `nothing` for natural order
"""
Base.@kwdef mutable struct AddressManager
    addresses::Dict{Tuple{AddressableType,Symbol,Symbol},Vector{UInt}} = Dict()
    next_idx::Dict{AddressableType,Int} = Dict()

    # Address assignment tracking
    assignments::Dict{Symbol,ModelAddressAssignment} = Dict()
    retrieved::Dict{Symbol,ModelRetrievedAddresses} = Dict()

    # Reverse address lookup (built during setup_with_formulation!)
    reverse_map::Union{Nothing, ReverseAddressMap} = nothing

    # Address reordering strategy, and the permutation it applied
    strategy::AddressingStrategy = NaturalAddressing()
    applied_permutation::Union{Nothing, Vector{Int}} = nothing
end

export AddressManager
export AddressingStrategy, NaturalAddressing, AMDAddressing, CustomPermutation
export ModelAddressAssignment, ReverseAddressMap

### ============= End Address Manager ============= ###

### ============= External-to-Internal Indexing (ext2int) ============= ###

"""
    GroupLinkRef

References to targets that may be heterogeneous (different model types).

Used for both:
1. **Control links**: Controllers (TGOV1) → targets in a group (:generators)
2. **Terminal links**: Models → bus-like targets in :bus_group

Each instance can reference a different target model type, enabling polymorphic
connections (e.g., TGOV1[1] → GENCLS, TGOV1[2] → GENROU).

# Fields
- `target_type::Vector{Symbol}`: Target model types for each instance
- `target_idx::Vector{Int}`: Target indices within each model type's array

# Example
```julia
# PQ[1] → Bus[1], PQ[2] → Bus[2], PQ[3] → Bus[3] (homogeneous)
ref = GroupLinkRef([:Bus, :Bus, :Bus], [1, 2, 3])

# TGOV1[1] → GENCLS[5], TGOV1[2] → GENROU[1] (heterogeneous)
ref = GroupLinkRef([:GENCLS, :GENROU], [5, 1])

# Access: get target for instance 2
target_model = ref.target_type[2]  # → :GENROU
target_idx = ref.target_idx[2]     # → 1
```
"""
struct GroupLinkRef
    target_type::Vector{Symbol}
    target_idx::Vector{Int}
end

"""
    GroupRefKey

Alias for `Union{Symbol, Tuple{Vararg{Symbol}}}`, the key of the per-model foreign
reference tables in [`IndexOrder`](@ref): a single reference field such as `:i`, or a
composite key such as `(:BUS, :ID)`.
"""
const GroupRefKey = Union{Symbol, Tuple{Vararg{Symbol}}}

# Helper to create homogeneous GroupLinkRef (all same target type)
function GroupLinkRef(target_model::Symbol, indices::Vector{Int})
    GroupLinkRef(fill(target_model, length(indices)), indices)
end

"""
    is_heterogeneous(ref::GroupLinkRef) -> Bool

Check if a GroupLinkRef references multiple different model types.

Heterogeneous refs come from group references (e.g., :generators → GENCLS, GENROU).
Homogeneous refs come from direct model references (e.g., :Bus → all Bus).

Use this to determine which gather/scatter code path to use after looking up
a foreign reference by ref_key.

# Example
```julia
ref = foreign_refs[:TGOV1][(:BUS, :ID)]
if is_heterogeneous(ref)
    # Gather from multiple model types
else
    # Fast path: all same model type
end
```
"""
is_heterogeneous(ref::GroupLinkRef) = !allequal(ref.target_type)

"""
    IndexOrder

Stores external-to-internal index mappings following MATPOWER's ext2int/int2ext pattern.

This structure is created during the ext2int preprocessing stage and stored in SystemModel.
It enables fast lookups from external bus numbers (and composite keys) to internal array indices.

# Fields
- `key_indices::Dict{Symbol, Dict}`: Model instance lookup tables
  - Key: Model name (:Bus, :PQ, :PV, ...)
  - Value: Dict mapping composite key tuple → internal array index
  - Example: `key_indices[:Bus][(69953,)] → 67900` (bus number → array index)
  - Example: `key_indices[:PQ][(101, "1")] → 5` (bus number, load id → array index)

- `foreign_refs::Dict{Symbol, Dict{GroupRefKey, GroupLinkRef}}`: Pre-computed group references
  - Outer key: Referencing model name (:PQ, :Branch, :TGOV1, ...)
  - Inner key: Field/handle for the reference (e.g., :i, :j, or composite key tuple)
  - Value: GroupLinkRef with target types and indices (supports heterogeneous targets)
  - Example: `foreign_refs[:PQ][:i]` → GroupLinkRef([:Bus, :Bus, :Bus], [1, 2, 3])
  - Example: `foreign_refs[:Branch][:i]` → GroupLinkRef(:Bus, [5, 7, 9])  # from-buses
  - Example: `foreign_refs[:TGOV1][(:BUS, :ID)]` → GroupLinkRef([:GENCLS, :GENROU], [1, 1])

# Purpose
All models store external bus numbers (semantic meaning preserved),
but computation requires internal array indices. This structure provides O(1) lookup via
pre-computed indexer vectors, avoiding Dict lookups in hot paths.

# See Also
- `ext2int`: Builds this structure from parsed models
- `int2ext`: Future function to convert results back to external numbering
"""
struct IndexOrder
    # Model instance key → internal array index
    # key_indices[:ModelName][(external_key_tuple...)] → internal_array_index
    key_indices::Dict{Symbol, Dict}

    # Pre-computed group reference indexers for fast array-based lookup
    # foreign_refs[:ModelName][field_key] = GroupLinkRef(target_types, target_indices)
    # field_key supports Symbols (single-field terminals) and Tuple{Symbol,...} (composite terminals)
    foreign_refs::Dict{Symbol, Dict{GroupRefKey, GroupLinkRef}}
end

IndexOrder() = IndexOrder(
    Dict{Symbol, Dict}(),
    Dict{Symbol, Dict{GroupRefKey, GroupLinkRef}}()
)

export IndexOrder

### ============= End ext2int ============= ###

### ============= Model Metadata and Traits ============= ###

# === Per-Unit Conversion Types === #

"""
    PUType

Supertype of the per-unit conversion kinds declared for model parameters
([`ModelParam`](@ref)) and observables ([`ModelObservable`](@ref)). The kind selects
the coefficient that [`calc_pu_coefficients`](@ref PowerfulCore.System.calc_pu_coefficients) computes from the system base
power `Sb`, the bus base voltage `Vb` and the device ratings `Sn`, `Vn`; the same
coefficient converts results back to physical units in reports.
"""
abstract type PUType end

"""
    NoPU()

Parameter needs no conversion: it is already in per unit on the system base, or unitless.
"""
struct NoPU <: PUType end

# Nominal (physical) units → per-unit conversion
"""
    PowerNominal()

Power given in MW or Mvar; converted to per unit by dividing by the system base `Sb`.
"""
struct PowerNominal <: PUType end

"""
    ActivePower()

Active power (MW in reports). Used to tag observables so that reports convert
them from per unit to MW.
"""
struct ActivePower <: PUType end

"""
    ReactivePower()

Reactive power (Mvar in reports). Used to tag observables so that reports convert
them from per unit to Mvar.
"""
struct ReactivePower <: PUType end

"""
    VoltageNominal()

Voltage given in kV; converted to per unit by dividing by the bus base voltage `Vb`.
"""
struct VoltageNominal <: PUType end

"""
    CurrentNominal()

Current given in physical units; converted to per unit on the base current `Sb/Vb`.
"""
struct CurrentNominal <: PUType end

"""
    ImpedanceNominal()

Impedance given in ohms; converted to per unit on the base impedance `Vb²/Sb`.
"""
struct ImpedanceNominal <: PUType end

"""
    AdmittanceNominal()

Admittance given in siemens; converted to per unit on the base admittance `Sb/Vb²`.
"""
struct AdmittanceNominal <: PUType end

"""
    AngleDegrees()

Angle given in degrees; converted to radians.
"""
struct AngleDegrees <: PUType end

# Already per-unit → base conversion (rare in PSS/E)
"""
    PowerPU()

Power in per unit on the device base `Sn`; converted to the system base (factor `Sn/Sb`).
"""
struct PowerPU <: PUType end

"""
    VoltagePU()

Voltage in per unit on the device voltage rating `Vn`; converted to the bus base
(factor `Vn/Vb`).
"""
struct VoltagePU <: PUType end

"""
    CurrentPU()

Current in per unit on the device base `Sn/Vn`; converted to the system base
`Sb/Vb`.
"""
struct CurrentPU <: PUType end

"""
    ImpedancePU()

Impedance in per unit on the device base `Vn²/Sn`; converted to the system base
`Vb²/Sb`.
"""
struct ImpedancePU <: PUType end

"""
    AdmittancePU()

Admittance in per unit on the device base; converted to the system base (the
inverse of the [`ImpedancePU`](@ref) factor).
"""
struct AdmittancePU <: PUType end

"""
    RatingType

Supertype of device rating kinds declared with [`ModelRating`](@ref). A rating field
supplies the device base used by the `*PU` conversions ([`PUType`](@ref)); a device
without a rating uses the system base power and the bus base voltage.
"""
abstract type RatingType end

"""
    NominalPower()

Device power rating `Sn` (MVA), for example a generator's `mbase`.
"""
struct NominalPower <: RatingType end

"""
    NominalVoltage()

Device voltage rating `Vn` (kV).
"""
struct NominalVoltage <: RatingType end

"""
    NominalCurrent()

Device current rating.
"""
struct NominalCurrent <: RatingType end

"""
    ReferenceType

Supertype of the kinds of external parameters ([`ExtParam`](@ref)) a model reads from
a model it references, and of the base quantities a model can own or inherit (see
[`base_quantity_role`](@ref)).
"""
abstract type ReferenceType end

"""
    BaseVoltage()

Base voltage (kV) of the referenced bus, for example the Bus field `basekv`.
"""
struct BaseVoltage <: ReferenceType end

"""
    BasePower()

System base power (MVA).
"""
struct BasePower <: ReferenceType end

"""
    BaseCurrent()

Base current.
"""
struct BaseCurrent <: ReferenceType end

"""
    Dimensionless()

External parameter that is not a base quantity and is passed through without
conversion, such as a bus type code, area or zone.
"""
struct Dimensionless <: ReferenceType end

"""
    RetrievalPhase

Supertype of the phases in which an external parameter ([`ExtParam`](@ref)) is
retrieved: [`PrePerUnit`](@ref) or [`PostPerUnit`](@ref).
"""
abstract type RetrievalPhase end

"""
    PrePerUnit()

Retrieve the external parameter before per-unit conversion. Used for quantities
the conversion itself needs (base voltages) and for unitless data (types, indices).
"""
struct PrePerUnit <: RetrievalPhase end

"""
    PostPerUnit()

Retrieve the external parameter after per-unit conversion, so it is read in per
unit.
"""
struct PostPerUnit <: RetrievalPhase end

"""
    ModelParam

Metadata for a model parameter field (e.g., resistance, reactance).

- `name`: Field name in model struct
- `pu_type`: Per-unit type for conversion (e.g., `Impedance`, `Power`)
"""
struct ModelParam
    name::Symbol
    pu_type::PUType
    description::String
    units::String
end

ModelParam(name::Symbol, ::Type{T}; description::String="", units::String="") where T<:PUType =
    ModelParam(name, T(), description, units)

"""
    ModelObservable

Metadata for computed/observable output variables.
These are values computed from state variables and parameters, not stored in DAE.
Examples: power flows in branches, reactive power from PV generators.

# Fields
- `name::Symbol`: Variable name (e.g., :p_from, :q_gen)
- `pu_type::PUType`: Per-unit type for inverse conversion to physical units
- `description::String`: Human-readable description
- `units::String`: Units in per-unit system (e.g., "pu", "rad")
- `default_output::Bool`: Whether to include in reports by default (default: true)
"""
struct ModelObservable
    name::Symbol
    pu_type::PUType
    description::String
    units::String
    default_output::Bool
end

ModelObservable(name::Symbol, ::Type{T}; description::String="", units::String="", default_output::Bool=true) where T<:PUType =
    ModelObservable(name, T(), description, units, default_output)

"""
    RatingSource

Where a device rating ([`ModelRating`](@ref)) comes from: [`OwnField`](@ref), a field of the
model, or [`LinkedDevice`](@ref), the rating of the static device the model is linked to.
"""
abstract type RatingSource end

"""
    OwnField()

The rating is a field of the model (the default), e.g. the `mbase` of a PV generator.
"""
struct OwnField <: RatingSource end

"""
    LinkedDevice()

The rating is that of the device the model replaces, found through its links
(`get_replacement_map`, `ModelLinks`): the `field` of the rating names the field of the
linked device that holds it, and the linked device declares a rating of the same
[`RatingType`](@ref) in it. A converter model, whose parameters are per unit on the
converter's `MBASE` (a field of the generator record it replaces, not of its own record),
declares `ModelRating(:mbase, NominalPower; from = LinkedDevice())`.
"""
struct LinkedDevice <: RatingSource end

"""
    ModelRating

Declares which field provides a device rating: a field of the model itself, or (`from =
LinkedDevice()`) of the device it is linked to.
"""
struct ModelRating
    field::Symbol
    rating_type::RatingType
    description::String
    source::RatingSource
end

ModelRating(field::Symbol, rating_type::RatingType, description::String) =
    ModelRating(field, rating_type, description, OwnField())

ModelRating(field::Symbol, ::Type{T}; description::String="", from::RatingSource=OwnField()) where T<:RatingType =
    ModelRating(field, T(), description, from)

"""
    RatingRatio(name, NominalPower; description = "")

Declares a derived parameter of a model: `name` is the ratio of the model's resolved power
rating `Sn` (its [`ModelRating`](@ref), from its own field or from the device it replaces,
[`LinkedDevice`](@ref)) to the system base `Sb`, `Sn/Sb`. It is not a field of the model and
has no DYR entry; Powerful computes it with the per-unit coefficients (`calc_pu_coefficients`,
which keeps it in `sys.pu_coeffs[model][name]`) and the kernels read it as
`external.name[i]` (the external parameters of the model), and initialization as
`gather_init_data(ctx, M).name`.

A model that declares it multiplies by it the quantities it holds on its own (device) base
where they enter the network, which is on the system base: a power or a current into the
bus is `ratio` times the device-base value, and the system-base power read from the network
is divided by it. The model must declare a `NominalPower` rating (a derived ratio without one
is an error, not a ratio of 1).
"""
struct RatingRatio
    name::Symbol
    rating_type::RatingType
    description::String
    function RatingRatio(name::Symbol, rating_type::RatingType, description::String)
        rating_type isa NominalPower || throw(ArgumentError(
            "RatingRatio(:$name, $(nameof(typeof(rating_type)))): only a NominalPower ratio (Sn/Sb) is defined"))
        new(name, rating_type, description)
    end
end

RatingRatio(name::Symbol, ::Type{T}; description::String="") where T<:RatingType =
    RatingRatio(name, T(), description)

"""
    Terminal

Declares that a model field references a vertex (Bus, Node, etc.) or another model
in the network graph.

Used for both topology analysis and external parameter retrieval (ExtParam).
This is the single source of truth for model connectivity.

# Fields
- `vertex_model::Symbol`: Target model type (`:Bus`, `:Node`, `:GENCLS`, etc.) or group (`:bus_group`, `:generators`)
- `ref_key::Union{Symbol, Vector{Symbol}}`: Field name(s) in this model used to look up the target
  - Single field: `:i` (simple key like bus number)
  - Composite key: `[:BUS, :ID]` (for matching dynamic models)

# Examples
```julia
# Branch connects two buses (single-field keys)
Terminal(:Bus, :i)  # "from" bus at field :i
Terminal(:Bus, :j)  # "to" bus at field :j

# TGOV1 controls a generator (composite key)
Terminal(:generators, [:BUS, :ID])  # Match by bus number + machine ID
```
"""
struct Terminal
    vertex_model::Symbol
    ref_key::Union{Symbol, Vector{Symbol}}
end

"""
    ExtParam

Declares how to retrieve external parameters from other models.

ExtParam uses `Terminal` as base type to avoid duplicating
connectivity information. Terminal specifies WHAT vertex is referenced, ExtParam adds
WHICH FIELD to fetch from that vertex.

# Fields
- `name::Symbol`: Local parameter name in this model (e.g., `:i_basekv`, `:basekv`)
- `terminal::Terminal`: Vertex connection (contains vertex_model and field)
- `ref_field::Symbol`: Field to retrieve from the referenced vertex (e.g., `:basekv`, `:ide`)
- `ref_type::ReferenceType`: Type of reference (BaseVoltage, BasePower, etc.)
- `phase::RetrievalPhase`: When to retrieve (PrePerUnit or PostPerUnit)
- `required::Bool`: Whether this reference is required
- `description::String`: Human-readable description

# Examples
```julia
# Branch: Get base voltage from buses at both terminals
const BRANCH_TERMINALS = (Terminal(:Bus, :i), Terminal(:Bus, :j))

const BRANCH_EXT_PARAMS = [
    ExtParam(:i_basekv, BRANCH_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Base voltage at bus i"),
    ExtParam(:j_basekv, BRANCH_TERMINALS[2], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Base voltage at bus j"),
    ExtParam(:i_ide, BRANCH_TERMINALS[1], :ide, Dimensionless(), PrePerUnit(),
             required=true, description="Bus type at bus i")
]

# PQ Load: Get voltage from single attached bus
const PQ_TERMINALS = (Terminal(:Bus, :i),)

const PQ_EXT_PARAMS = [
    ExtParam(:basekv, PV_TERMINALS[1], :basekv, BaseVoltage(), PrePerUnit(),
             required=true, description="Connected bus base voltage")
]
```
"""
struct ExtParam
    name::Symbol           # Local parameter name
    terminal::Terminal     # Vertex reference (contains vertex_model + field)
    ref_field::Symbol      # Field to fetch from that vertex
    ref_type::ReferenceType
    phase::RetrievalPhase
    required::Bool
    description::String
end

ExtParam(
    name::Symbol,
    terminal::Terminal,
    ref_field::Symbol,
    ref_type::ReferenceType,
    phase::RetrievalPhase = PrePerUnit();
    required::Bool = true,
    description::String = ""
) = ExtParam(name, terminal, ref_field, ref_type, phase, required, description)

"""
    InitParam(name::Symbol, source_var::Symbol; description::String="")
    InitParam(name::Symbol; description::String="")

Metadata for initialized parameters (setpoints/constants computed during initialization).

These are constant values computed once during initialization and stored in `sys.init_params`.
They are used by models for holding equations (when no controller is present) or as references.

# Two-argument form: Variable → Setpoint copy
Use when copying a variable's initial value to a differently-named setpoint:
```julia
InitParam(:P_m0, :P_m; description="Mechanical power setpoint")
# Copies initial value of variable :P_m to setpoint :P_m0
```

# One-argument form: Computed setpoint
Use when the setpoint is computed directly in initialization (not copied from a variable):
```julia
InitParam(:pfaref0; description="Initial power factor angle")
# Setpoint :pfaref0 is provided directly by initialize_dynamic_states
```

# Storage
Values are stored in `sys.init_params[:ModelName].param_name` as vectors.

# Usage
Models access these via `ws.setpoints.P_m0[i]` in residual calculations.
"""
struct InitParam
    name::Symbol           # Parameter name (e.g., :P_m0)
    source_var::Symbol     # Variable to copy from (e.g., :P_m), or same as name for computed setpoints
    description::String
end

# Two-argument form: explicit variable → setpoint (must have different names)
function InitParam(name::Symbol, source_var::Symbol; description::String="")
    if name === source_var
        error("InitParam: use single-argument form InitParam(:$name) for computed setpoints. " *
              "Two-argument form is for copying a variable to a differently-named setpoint, " *
              "e.g., InitParam(:P_m0, :P_m) copies P_m → P_m0.")
    end
    InitParam(name, source_var, description)
end

# One-argument form: computed setpoint (same name in init result and setpoints)
InitParam(name::Symbol; description::String="") = InitParam(name, name, description)

"""
    InitInput

Declares an input required for dynamic model initialization.

Supports three input kinds:
1. **:from_observation**: Computed value from linked power flow model (e.g., P, Q from PV/Slack)
2. **:from_solution**: Direct variable from DAE solution (e.g., V, θ from Bus)
3. **:from_group_init**: State from any model in a group (polymorphic - e.g., P_m from :generators group)

# Fields
- `name::Symbol`: Local name in initialization context (e.g., :P_gen, :V_bus)
- `kind::Symbol`: Input kind - :from_observation, :from_solution, or :from_group_init
- `from_model::Symbol`: Source model type or GROUP name (e.g., :Generator, :Bus, or :generators)
- `field::Symbol`: Field/variable name in source (canonical group name for :from_group_init)
- `description::String`: Human-readable description

# Examples
```julia
# GENCLS needs P and Q from linked PV/Slack generators
InitInput(:P_gen, :from_observation, :Generator, :p_gen;
    description="Active power from linked generator")

# Bus voltage and angle from DAE solution
InitInput(:V_bus, :from_solution, :Bus, :v;
    description="Bus voltage magnitude")

# TGOV1 uses group reference for polymorphism (works with GENCLS, GENROU, etc.)
InitInput(:tm0, :from_group_init, :generators, :P_m;
    description="Mechanical power from any generator in group")
```
"""
struct InitInput
    name::Symbol                      # Local name for use in initialization
    kind::Symbol                      # :from_observation, :from_solution, or :from_group_init
    from_model::Symbol                # Source model type OR group name
    field::Symbol                     # Field name (canonical group name for :from_group_init)
    description::String
end

# Constructor with keyword arguments for convenience
function InitInput(
    name::Symbol,
    kind::Symbol,
    from_model::Symbol,
    field::Symbol;
    description::String = ""
)
    kind in (:from_observation, :from_solution, :from_group_init) ||
        error("kind must be :from_observation, :from_solution, or :from_group_init, got :$kind")

    InitInput(name, kind, from_model, field, description)
end

"""
    Disturbance

Abstract base type for all TDS disturbance events.

Disturbance events are discrete changes applied during time-domain simulation,
such as faults, trips, load changes, or setpoint modifications.

Concrete types:
- `FieldModification`: Universal primitive for modifying model fields
- Event types that also reset states (e.g. a fault clearing that zeroes fault currents)

Each concrete type must implement a `to_callback()` method to integrate
with DifferentialEquations.jl.
"""
abstract type Disturbance end

"""
    to_callback(event, sys, ...) -> Callback

Convert a Disturbance event to a DifferentialEquations.jl callback.

This is a generic function that specialized event types extend with their own methods.
The signature varies by event type:
- `FieldModification`: `to_callback(event, sys, DiffEq; initializealg=nothing)`
- Events that reset states: `to_callback(event, sys, am, adapter, reevaluate_fn)`; they choose their own initialization
"""
function to_callback end

# ============================================================================
# Variable Limits
# ============================================================================

"""
    VariableLimit(variable, lower, upper; when = nothing)

The bounds of a model variable, declared apart from how they are enforced: by an
automaton that holds the variable at them (anti-windup lowered to an automaton), see
[`variable_limits`](@ref). Each bound is `nothing` (no bound at
that end), the name of a model field, or a function of the model arrays computing it for
every instance (or one number for all), as a `MassCoefficient` is; a function of two
arguments, the model arrays and the model's initialized variables (a `NamedTuple` of the
variables' initial values), computes a bound that depends on the operating point. The initialization
checks the initial value of the variable against them, and moves a violated bound that
is a field to it (`enforce_initial_limits!`, `TDSConfig(adjust_limits)`); a computed bound
cannot be moved, and its violation is an error.

`when` (a function of the model arrays and the model's initialized variables, as a bound
of two arguments is, computing for every instance whether the limit applies: a `Bool`
vector or one `Bool`) restricts the check to the instances where it holds: a variable that
equals what the limits bound only under some flags or operating points (a reference that is
the reactive power only where a controller acts on it). `nothing` applies it everywhere.
Where a bound is moved to the initial value the start is computed again from the moved data
(`enforce_initial_limits!`).

# Example
```julia
VariableLimit(:LAG_y, :Vmin, :Vmax)
VariableLimit(:LAG_y, :Vmin, m -> 2 .* m.Vmax)
VariableLimit(:PI_xi, (m, initial) -> .-initial.Imax, (m, initial) -> initial.Imax)  # bounds read from the initialized Imax
VariableLimit(:Qref, :QMin, :QMax; when = (m, x) -> m.VFLAG .== 1)                   # only where VFLAG = 1
```
"""
struct VariableLimit{L<:Union{Nothing, Symbol, Function}, U<:Union{Nothing, Symbol, Function},
                     W<:Union{Nothing, Function}}
    variable::Symbol
    lower::L
    upper::U
    when::W
end

VariableLimit(variable::Symbol, lower, upper; when = nothing) = VariableLimit(variable, lower, upper, when)

"""
    InjectionForms(natural; power, current, voltage = Symbol[])

The forms of the injection a model delivers to the bus group of a terminal: the residuals
that carry its power (`p`, `q`), those that carry its current (`ip`, `iq`: in phase with and in
quadrature to the bus voltage, positive into the bus; `p = v ip`, `q = v iq`), and, for a model
that is an electromotive force E behind an impedance Z (a machine behind its subtransient
reactance), those that carry its short-circuit current (`jr`, `ji`: J = E/Z on the system base,
real and imaginary parts in the network frame), `voltage`. The model declares one form
(`natural`, `:power`, `:current` or `:voltage`) and computes the others from it; a connection is
delivered through exactly one of them (`injection_routes`): the declared form where the bus
implements it, otherwise, for `:voltage`, the current form and then the power form, and for the
other two the other one.

The voltage form is realized by the bus. A bus that implements it (a bridging device)
gives the model the voltage behind which its delivered current I flows, V_port = E - Z I (the
group variables `v_port`, `theta_port`), so the model's current (E - V_port)/Z is the current
the network receives; every other bus gives V_port = V, the bus voltage, which realizes the same
form algebraically, I = (E - V)/Z.
"""
struct InjectionForms
    natural::Symbol
    power::Vector{Symbol}
    current::Vector{Symbol}
    voltage::Vector{Symbol}

    function InjectionForms(natural::Symbol; power::Vector{Symbol}, current::Vector{Symbol},
                            voltage::Vector{Symbol} = Symbol[])
        natural in (:power, :current, :voltage) ||
            error("InjectionForms natural form must be :power, :current or :voltage, got :$natural")
        natural === :voltage && isempty(voltage) &&
            error("InjectionForms with the natural form :voltage needs the residuals of that form")
        new(natural, power, current, voltage)
    end
end

export InjectionForms

"""
    ModelMetadata

Complete metadata for a power system model.

# Formulation Support
Models can support multiple formulations (PowerInjection, CurrentInjection).
The `vars` and `residuals` fields define the default formulation (PowerInjection).
Alternative formulations can override these via `vars_by_formulation` and `residuals_by_formulation`.

# Automata
`automata` and `memory` declare the discrete behaviour of the model (see [`Automaton`](@ref)
and [`MemoryQuantity`](@ref)), as tuples; both are empty for a model without it. A model that
declares automata computes the values they name in `calc_automaton_values!`, and its
residual kernels read the discrete quantities and memory of its workspace
(`ws.automata.discrete`, `ws.automata.memory`, see `AutomatonState`).

# Example
```julia
function model_metadata(::Type{Branch})
    base = ModelMetadata(
        name = :Branch,
        vars = BRANCH_VARS,  # PowerInjection variables (θ, V)
        residuals = BRANCH_RESIDUALS,  # Power balance
        # ...
    )

    # Add CurrentInjection support
    base.vars_by_formulation[CurrentInjection] = BRANCH_VARS_CURRENT  # (Vr, Vi)
    base.residuals_by_formulation[CurrentInjection] = BRANCH_RESIDUALS_CURRENT  # Current balance

    return base
end
```
"""
Base.@kwdef mutable struct ModelMetadata
    name::Symbol
    vars::Vector{ModelVar}
    observable_vars::Vector{Symbol}
    residuals::Vector{ModelResidual}
    params::Vector{ModelParam} = ModelParam[]
    observations::Vector{ModelObservable} = ModelObservable[]
    ratings::Vector{ModelRating} = ModelRating[]
    derived_params::Vector{RatingRatio} = RatingRatio[]  # Ratios of the rating to the system base (Sn/Sb), not fields
    ext_params::Vector{ExtParam} = ExtParam[]
    initialized_params::Vector{InitParam} = InitParam[]
    init_inputs::Vector{InitInput} = InitInput[]  # Data needed from power flow for initialization
    variable_limits::Vector{VariableLimit} = VariableLimit[]  # Bounds enforced by other means (automata)
    # Discrete behaviour (automata.jl): declarations, as tuples (immutable, and the empty
    # default allocates nothing when metadata is built, which it is at every query)
    automata::Tuple{Vararg{Automaton}} = ()  # Modes, transitions, starts
    # Whether any mode of `automata` governs a state (`governs_states(automata)`), cached as a concrete
    # `Bool` at construction (round 19a, Phase D): `automata`'s declared field type, `Tuple{Vararg{Automaton}}`,
    # is not concrete (it admits every tuple length), so a hot-path call to `governs_states(metadata.automata)`
    # is a dynamic dispatch on an abstractly-typed value, not a branch the compiler can specialize away even
    # for a model that declares no automata at all. `automata` is set once and never reassigned after
    # construction, so caching this here is exact, not an approximation; read this field on every hot path
    # instead of recomputing `governs_states(metadata.automata)`.
    governs_any_state::Bool = governs_states(automata)
    memory::Tuple{Vararg{MemoryQuantity}} = ()  # Values kept between transitions, set by them
    status_field::Union{Symbol, Nothing} = nothing  # Field name for device status (:status, a Bool column)
    injection::Union{InjectionForms, Nothing} = nothing  # The forms of the bus injection (power, current), routed per connection

    # Evaluation pass: controls when model residuals are evaluated
    # Pass 0 (default): Standard accumulation (branches, loads)
    # Pass 1 (later): Post-processing that reads accumulated residuals (PV in CurrentInjection)
    evaluation_pass::Int = 0

    # Formulation-specific overrides (mutable for post-construction modification)
    vars_by_formulation::Dict{DataType, Vector{ModelVar}} = Dict{DataType, Vector{ModelVar}}()
    residuals_by_formulation::Dict{DataType, Vector{ModelResidual}} = Dict{DataType, Vector{ModelResidual}}()
end

# Load model types - now identified via LoadModel interface (see model_traits.jl)
# No abstract type needed - all loads inherit from AbstractModel directly

export ModelMetadata
export PUType, NoPU
export PowerNominal, ActivePower, ReactivePower, VoltageNominal, CurrentNominal, ImpedanceNominal, AdmittanceNominal, AngleDegrees
export PowerPU, VoltagePU, CurrentPU, ImpedancePU, AdmittancePU
export RatingType, NominalPower, NominalVoltage, NominalCurrent
export RatingSource, OwnField, LinkedDevice, RatingRatio
export ReferenceType, BaseVoltage, BasePower, BaseCurrent, Dimensionless
export RetrievalPhase, PrePerUnit, PostPerUnit
export ModelParam, ModelObservable, ModelRating, Terminal, ExtParam, InitParam, InitInput
export VariableLimit
export Disturbance, to_callback

### ============= End Model Metadata and Traits ============= ###

### ============= System Types ============= ###

"""
    SystemModel

Power system network model containing all components (buses, branches, loads, generators, etc.).

A pure data container storing per-unit models and index mappings. Formulation-specific infrastructure
(AddressManager, DAE vectors) is created at problem construction time, not stored here.

# Fields
$TYPEDFIELDS

# Constructors
```julia
SystemModel(case::ParsedSystemData)      # From parsed RAW file
SystemModel(components::NamedTuple)      # Direct from model StructArrays
```

# Example
```julia
using Powerful, PowerfulCases

# Load and parse a case file
case = parse_network(load("ieee14").raw)
sys = SystemModel(case)

# Access models
sys.models.Bus      # StructArray of all buses
sys.models.PQ       # StructArray of PQ loads
sys.models.Branch   # StructArray of branches

# Check system size
length(sys.models.Bus)  # Number of buses

# Create problem and solve
prob = ACPowerFlowProblem(sys)
sol = solve(prob, Newton())
```

# See Also
- [`ACPowerFlowProblem`](@ref PowerfulPF.ACPowerFlowProblem): Create power flow problem from system
- [`TDSProblem`](@ref PowerfulTDS.TDSProblem): Create time-domain simulation problem
- `parse_network`: Parse PSS/E RAW file
"""
mutable struct SystemModel
    "A named tuple of models (per-unit for computation)"
    models::NamedTuple

    "External parameters retrieved from other models (topology-derived references)"
    external::NamedTuple

    "Index mappings for external-to-internal conversion (MATPOWER ext2int pattern)"
    order::Union{IndexOrder, Nothing}

    "Array type for computational backend (Vector for CPU, CuVector for GPU)"
    array_type::Type{<:AbstractVector}

    "Optional: original parsed data in physical units (from ParsedSystemData)"
    parsed_data::Union{Nothing, Any}  # Any to avoid circular dependency with ParsedSystemData

    "Optional: per-unit conversion coefficients for inverse conversion (pu → physical)"
    pu_coeffs::Union{Nothing, NamedTuple}  # Coefficients from convert_system_to_per_unit

    "Optional: cached topology analysis results (connectivity, islands)"
    topology::Union{Nothing, Any}  # TopologyInfo, defined in topology.jl

    "System-wide settings (frequency and MVA base) used by every analysis of the system"
    settings::SystemConfig{Float64}

    # Inner constructor to validate component types
    # Note: external, order and settings are REQUIRED - use the outer constructor in
    # system_builder.jl for auto-population via ext2int; systems derived from another
    # system pass its settings on
    function SystemModel(components::NamedTuple,
                         external::NamedTuple,
                         order::Union{IndexOrder, Nothing};
                         settings::SystemConfig{Float64},
                         array_type::Type{<:AbstractVector}=Vector,
                         parsed_data::Union{Nothing, Any}=nothing,
                         pu_coeffs::Union{Nothing, NamedTuple}=nothing,
                         topology::Union{Nothing, Any}=nothing)
        # Ensure all elements are vectors/arrays of AbstractModel subtypes
        # Check eltype since StructArray wraps the base type
        all(V -> eltype(V) <: AbstractModel, values(components)) ||
            error("All components must be AbstractModel subtypes")

        return new(components, external, order, array_type, parsed_data, pu_coeffs, topology, settings)
    end
end

# Outer constructor - forwards to helper in system_builder.jl
# IMPORTANT: Must be defined here (next to struct) for Julia LSP to find it.
# LSP has a bug with re-exported functions where it can't trace methods
# defined in submodules. See: https://github.com/julia-vscode/julia-vscode/issues/1576
SystemModel(components::NamedTuple; settings::SystemConfig{Float64} = SystemConfig()) =
    _build_system_from_namedtuple(components; settings)

export SystemModel

### ========================================= ###


### ============= Data Format Types ============= ###

"""
    DataFormat

Supertype of input data formats. A model declares which formats it can be built from
with [`supports_format`](@ref) and implements [`parse_model`](@ref) for them.
"""
abstract type DataFormat end

"""
    PSSE <: DataFormat

PSS/E RAW power flow data format.
"""
struct PSSE <: DataFormat end

"""
    PSSEDYR <: DataFormat

PSS/E DYR dynamic data format.
"""
struct PSSEDYR <: DataFormat end
struct MATPOWER <: DataFormat end

export DataFormat, PSSE, PSSEDYR

### ========================================= ###

### ============= Dynamic Model Linking ============= ###

"""
    ModelLinks

Struct-of-arrays for linking dynamic models to power flow models.

Used in time-domain simulation to map DYR file models (GENCLS, GENROU, etc.)
to their corresponding power flow models (PV, Slack).

# Fields
- `source_type::Vector{Symbol}` - Dynamic model types (e.g., [:GENCLS, :GENCLS, ...])
- `source_idx::Vector{Int}` - Indices in dynamic model arrays
- `target_type::Vector{Symbol}` - Static model types (e.g., [:PV, :Slack, :PV, ...])
- `target_idx::Vector{Int}` - Indices in static model arrays

# Example
```julia
links = ModelLinks(
    [:GENCLS, :GENCLS, :GENROU],
    [1, 2, 1],
    [:PV, :Slack, :PV],
    [3, 1, 5]
)
# GENCLS[1] → PV[3], GENCLS[2] → Slack[1], GENROU[1] → PV[5]
```

# Filtering by Source Model
```julia
mask = links.source_type .== :GENCLS
target_types = links.target_type[mask]
target_indices = links.target_idx[mask]
```
"""
struct ModelLinks
    source_type::Vector{Symbol}
    source_idx::Vector{Int}
    target_type::Vector{Symbol}
    target_idx::Vector{Int}

    function ModelLinks(source_type, source_idx, target_type, target_idx)
        n = length(source_type)
        @assert length(source_idx) == n "source_idx length mismatch"
        @assert length(target_type) == n "target_type length mismatch"
        @assert length(target_idx) == n "target_idx length mismatch"
        new(source_type, source_idx, target_type, target_idx)
    end
end

# Constructor from empty
ModelLinks() = ModelLinks(Symbol[], Int[], Symbol[], Int[])

# Length for iteration support
Base.length(links::ModelLinks) = length(links.source_type)
Base.isempty(links::ModelLinks) = isempty(links.source_type)

export ModelLinks

### ========================================= ###
