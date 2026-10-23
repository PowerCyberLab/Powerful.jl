"""
Per-unit conversion utilities for power system models.

Provides generic, metadata-driven conversion from physical units to per-unit values.
"""

export extract_device_ratings, extract_external_bases
export linked_device_ratings, declares_linked_rating, convert_linked_models
export calc_pu_coefficients, convert_to_per_unit
export convert_system_to_per_unit
export ParsedSystemData

"""
    ParsedSystemData

Simple container for preserving original parsed data (physical units) before per-unit conversion.

This is an **optional** container - you only need it if you want to preserve the original
physical unit data alongside the per-unit computational model.

# Fields
- `models::NamedTuple`: Original model StructArrays in physical units (MW, MVAR, kV)
- `system_base::Float64`: System base power in MVA
- `metadata::Dict{String,Any}`: Optional metadata (case info, timestamps, etc.)

# Basic Workflow (no container needed)
```julia
models = load_models(case, all_models())
models_pu = convert_system_to_per_unit(models, 100.0)
sys = SystemModel(am, models_pu)
```

# Advanced Workflow (preserve original data)
```julia
# Parse and store original
models = load_models(case, all_models())
parsed = ParsedSystemData(models, 100.0)

# Convert and create system
models_pu = convert_system_to_per_unit(parsed.models, parsed.system_base)
sys = SystemModel(am, models_pu; parsed_data=parsed)

# Access both
parsed.models.PQ.pl[1]  # → 21.7 MW (original)
sys.models.PQ.pl[1]     # → 0.217 pu (computational)
```
"""
struct ParsedSystemData
    models::NamedTuple          # Original models in physical units
    system_base::Float64        # System base MVA
    metadata::Dict{String,Any}  # Optional metadata
end

# Simple constructor with default metadata
ParsedSystemData(models::NamedTuple, system_base::Float64; metadata=Dict{String,Any}()) =
    ParsedSystemData(models, system_base, metadata)

"""
    extract_device_ratings(model, metadata; linked = nothing)

Extract device ratings from model using formal rating declarations.

Returns a Dict mapping RatingType to vectors of values. A rating of the model's own field
(`OwnField`) is read from the model. One that comes from the linked device (`LinkedDevice`)
is taken from `linked`, the ratings of the linked devices (`linked_device_ratings`); an error
if there are none.
"""
function extract_device_ratings(
    model::StructArray{M},
    metadata::ModelMetadata;
    linked::Union{Nothing, AbstractDict} = nothing
) where M <: AbstractModel

    ratings = Dict{Type{<:RatingType}, Vector{Float64}}()

    for rating_decl in metadata.ratings
        rating_type = typeof(rating_decl.rating_type)
        if rating_decl.source isa LinkedDevice
            (linked === nothing || !haskey(linked, rating_type)) && error(
                "The $(nameof(rating_type)) of $(metadata.name) comes from the device it is linked to " *
                "(ModelRating(:$(rating_decl.field), $(nameof(rating_type)); from = LinkedDevice())), " *
                "but no linked device ratings were given for it")
            ratings[rating_type] = linked[rating_type]
        elseif hasfield(M, rating_decl.field)
            ratings[rating_type] = getproperty(model, rating_decl.field)
        else
            error("Model $(M) missing declared rating field: $(rating_decl.field)")
        end
    end

    return ratings
end

"""
    declares_linked_rating(metadata) -> Bool

Whether `metadata` declares a rating that comes from the linked device (`LinkedDevice`): the
model's parameters are per unit on that device's base.
"""
declares_linked_rating(metadata::ModelMetadata) = any(r -> r.source isa LinkedDevice, metadata.ratings)

# `Name(9, "1")`: the instance key as a user writes it (strings not wrapped in their type)
_key_text(name::Symbol, key::Tuple) =
    "$name($(join((x isa AbstractString ? repr(String(x)) : string(x) for x in key), ", ")))"

"""
    linked_device_ratings(sys, name, models, metadata, links) -> Dict{Type{<:RatingType}, Vector{Float64}}

The ratings of the instances of `models` (model `name`) that `metadata` declares to come from
the linked device (`LinkedDevice`), from the devices of `sys` that `links` (`ModelLinks`, as
`create_model_links` builds them) link the instances to. An `ArgumentError` names the
instance, the model and the reason where the rating cannot be resolved: the instance is
linked to no device, the model of the linked device declares no rating of that type, or it
declares it in another field than the one the model names.
"""
function linked_device_ratings(sys::SystemModel, name::Symbol, models::StructArray, metadata::ModelMetadata,
                               links::ModelLinks)
    n = length(models)
    key_fields = get_key_name(get_base_model_type(models))
    # (the instance named only where it is refused)
    key_of(i) = _key_text(name, instance_key(models, key_fields, i))
    link_of = Dict{Int, Int}()
    for l in 1:length(links)
        links.source_type[l] === name && (link_of[links.source_idx[l]] = l)
    end

    metadata_of = Dict{Symbol, ModelMetadata}()   # (building a model's metadata is not free)
    ratings = Dict{Type{<:RatingType}, Vector{Float64}}()
    for decl in metadata.ratings
        decl.source isa LinkedDevice || continue
        rating_type = typeof(decl.rating_type)
        values = Vector{Float64}(undef, n)
        for i in 1:n
            l = get(link_of, i, 0)
            l == 0 && throw(ArgumentError(
                "$(key_of(i)) takes its $(nameof(rating_type)) ($(decl.field)) from the device it is " *
                "linked to, but it is linked to none (its replacement map found no device)"))
            target = links.target_type[l]
            haskey(sys.models, target) || throw(ArgumentError(
                "$(key_of(i)) is linked to $target, which is not a model of the system"))
            target_models = getfield(sys.models, target)
            target_meta = get!(() -> model_metadata(get_base_model_type(target_models)), metadata_of, target)
            declared = findfirst(r -> r.source isa OwnField && typeof(r.rating_type) === rating_type, target_meta.ratings)
            declared === nothing && throw(ArgumentError(
                "$(key_of(i)) takes its $(nameof(rating_type)) ($(decl.field)) from the device it is " *
                "linked to, $target, which declares no $(nameof(rating_type))"))
            target_meta.ratings[declared].field === decl.field || throw(ArgumentError(
                "$(key_of(i)) takes its $(nameof(rating_type)) from the field $(decl.field) of the device " *
                "it is linked to, $target, which declares it in $(target_meta.ratings[declared].field)"))
            values[i] = getproperty(target_models, decl.field)[links.target_idx[l]]
        end
        ratings[rating_type] = values
    end
    return ratings
end

"""
    extract_external_bases(model, metadata, all_models)

Extract external base values using formal reference declarations, organized by terminal.

# Returns
NamedTuple keyed by terminal names, where each terminal has a Dict of reference types.

# Examples
- Single-terminal device (PQ): `(terminal = Dict(BaseVoltage => [13.8, ...]), )`
- Multi-terminal device (Branch): `(i = Dict(BaseVoltage => [138, ...]), j = Dict(BaseVoltage => [230, ...]))`
"""
function extract_external_bases(
    model_name::Symbol,
    metadata::ModelMetadata,
    external::NamedTuple
)
    # Group external references by terminal
    terminal_bases = Dict{Symbol, Dict{Type{<:ReferenceType}, Vector{Float64}}}()

    # Get this model's external params
    if !haskey(external, model_name)
        return NamedTuple()  # No external params for this model
    end

    model_external = getfield(external, model_name)

    for ext_ref in metadata.ext_params
        # Only process PrePerUnit phase params (basekv, etc.) for per-unit conversion
        if typeof(ext_ref.phase) != PrePerUnit
            continue
        end

        # Check if this external param was retrieved
        if !haskey(model_external, ext_ref.name)
            if ext_ref.required
                error("Required external parameter $(ext_ref.name) not found for model $model_name")
            else
                continue
            end
        end

        param_values = getfield(model_external, ext_ref.name)

        # Group by terminal ref_key (Symbol identifier from Terminal)
        terminal_field = ext_ref.terminal.ref_key

        # Initialize terminal dict if needed
        if !haskey(terminal_bases, terminal_field)
            terminal_bases[terminal_field] = Dict{Type{<:ReferenceType}, Vector{Float64}}()
        end

        # Store by terminal and reference type
        terminal_bases[terminal_field][typeof(ext_ref.ref_type)] = param_values
    end

    # Convert to NamedTuple for type stability
    return NamedTuple{Tuple(keys(terminal_bases))}(values(terminal_bases))
end

"""
    calc_pu_coefficients(metadata, model, model_name, system_base, external)

Calculate per-unit conversion coefficients using formal metadata interface.

# Arguments
- `metadata::ModelMetadata`: Model metadata with params, ratings, ext_params
- `model::StructArray`: Model instances
- `model_name::Symbol`: Name of this model in the system
- `system_base::Float64`: System base power (MVA)
- `external::NamedTuple`: Pre-retrieved external parameters from retrieve_external_params(PrePerUnit)
- `linked_ratings`: the ratings of the linked devices, for a model that declares
  ratings from them (`linked_device_ratings`)

# Returns
Dict{Symbol, Vector{Float64}} mapping parameter names to conversion coefficients, and the names of the
model's derived parameters (`RatingRatio`) to their values (Sn/Sb)
"""
function calc_pu_coefficients(
    metadata::ModelMetadata,
    model::StructArray{M},
    model_name::Symbol,
    system_base::Float64,
    external::NamedTuple;
    linked_ratings::Union{Nothing, AbstractDict} = nothing
) where M <: AbstractModel

    n = length(model)
    Sb = system_base

    # Extract device ratings (Sn, Vn) from this model, or from the device it is linked to
    device_ratings = extract_device_ratings(model, metadata; linked = linked_ratings)

    # Extract external bases (Vb from Bus) from pre-retrieved external params
    external_bases = extract_external_bases(model_name, metadata, external)

    # A parameter on the device's own base (Sn ./ Sb, or an impedance/admittance/current through
    # Zn = Vn^2/Sn) needs a device power rating; unlike the voltage base below, Sn never cancels
    if !haskey(device_ratings, NominalPower)
        for param in metadata.params
            _needs_device_rating(param.pu_type) &&
                error("Per-unit conversion of $model_name.$(param.name) ($(nameof(typeof(param.pu_type)))) " *
                      "needs a device power rating (NominalPower), but $model_name has none")
        end
    end

    # Get device power rating (every parameter that needs it has one, checked above)
    Sn = get(device_ratings, NominalPower, fill(Sb, n))

    # Get bus base voltage from primary terminal
    # Convention: single-terminal uses :terminal, multi-terminal uses :i
    primary = haskey(external_bases, :terminal) ? external_bases.terminal :
              haskey(external_bases, :i) ? external_bases.i : nothing
    Vb = isnothing(primary) ? nothing : get(primary, BaseVoltage, nothing)
    Vn_rating = get(device_ratings, NominalVoltage, nothing)

    if isnothing(Vb)
        # Without a bus base voltage, only conversions that do not depend on it are valid
        for param in metadata.params
            _needs_base_voltage(param.pu_type, Vn_rating) &&
                error("Per-unit conversion of $model_name.$(param.name) ($(nameof(typeof(param.pu_type)))) " *
                      "needs a bus base voltage, but $model_name has no BaseVoltage external parameter")
        end
        # No conversion below depends on Vb: with no voltage rating Vn = Vb cancels
        Vb = ones(n)
    end

    # Get device voltage rating (default to bus base if not provided)
    Vn = isnothing(Vn_rating) ? Vb : Vn_rating

    # Calculate base impedances
    Zn = Vn.^2 ./ Sn
    Zb = Vb.^2 ./ Sb

    # Coefficient lookup table
    # Nominal (physical units) → per-unit conversion
    coeff_types = Dict(
        PowerNominal() => fill(1.0 / Sb, n),  # MW → pu: divide by system base
        VoltageNominal() => 1.0 ./ Vb,  # kV → pu: divide by bus base voltage
        CurrentNominal() => 1.0 ./ (Sb ./ Vb),  # A → pu
        ImpedanceNominal() => Sb ./ (Vb.^2),  # Ω → pu
        AdmittanceNominal() => (Vb.^2) ./ Sb,  # S → pu
        AngleDegrees() => fill(deg2rad(1.0), n),  # degrees → radians: use deg2rad conversion

        # Already per-unit → base conversion (rare)
        PowerPU() => Sn ./ Sb,  # pu (device base) → pu (system base)
        VoltagePU() => Vn ./ Vb,  # pu (device base) → pu (system base)
        CurrentPU() => (Sn ./ Vn) ./ (Sb ./ Vb),
        ImpedancePU() => Zn ./ Zb,
        AdmittancePU() => Zb ./ Zn,

        NoPU() => ones(n),  # No conversion
    )

    # Build coefficient dictionary
    coeffs = Dict{Symbol, Vector{Float64}}()
    for param in metadata.params
        coeffs[param.name] = coeff_types[param.pu_type]
    end

    # The derived ratios (`RatingRatio`) are kept with the coefficients: Sn/Sb of the same rating
    for derived in metadata.derived_params
        haskey(device_ratings, typeof(derived.rating_type)) || error(
            "$model_name declares the derived parameter $(derived.name) (the ratio of its " *
            "$(nameof(typeof(derived.rating_type))) to the system base) but no such rating")
        haskey(coeffs, derived.name) && error("$model_name declares $(derived.name) as a parameter and as a derived parameter")
        coeffs[derived.name] = Sn ./ Sb
    end

    return coeffs
end

"""
    _needs_base_voltage(pu_type, Vn_rating) -> Bool

Whether converting a parameter of `pu_type` depends on the bus base voltage. Nominal
voltage, current, impedance and admittance always do; their per-unit (device base)
counterparts do only when the device declares a voltage rating (otherwise the device
base voltage is the bus base voltage and cancels).
"""
_needs_base_voltage(::Union{VoltageNominal, CurrentNominal, ImpedanceNominal, AdmittanceNominal}, _) = true
_needs_base_voltage(::Union{VoltagePU, CurrentPU, ImpedancePU, AdmittancePU}, Vn_rating) = !isnothing(Vn_rating)
_needs_base_voltage(::PUType, _) = false

"""
    _needs_device_rating(pu_type) -> Bool

Whether converting a parameter of `pu_type` depends on the device's own power rating `Sn`:
its per-unit-on-own-base forms (power, current, impedance, admittance), whose coefficients
are `Sn / Sb` or go through `Zn = Vn^2 / Sn`. Unlike the voltage base, `Sn` never cancels:
a device without the rating would be converted as if its own base were the system base.
"""
_needs_device_rating(::Union{PowerPU, CurrentPU, ImpedancePU, AdmittancePU}) = true
_needs_device_rating(::PUType) = false

"""
    convert_to_per_unit(model, coeffs)

Apply per-unit conversion coefficients to model parameters.

Returns a new StructArray with converted values.
"""
function convert_to_per_unit(
    model::StructArray{M},
    coeffs::Dict{Symbol, Vector{Float64}}
) where M <: AbstractModel

    # Extract all fields
    field_data = NamedTuple{fieldnames(M)}(
        getproperty(model, fname) for fname in fieldnames(M)
    )

    # Apply coefficients to parameters that have them
    converted_pairs = [fname => (haskey(coeffs, fname) ? field_values .* coeffs[fname] : field_values)
                       for (fname, field_values) in pairs(field_data)]

    # Return new StructArray
    return StructArray{M}(NamedTuple(converted_pairs))
end

"""
    convert_system_to_per_unit(models::NamedTuple, system_base::Float64)

Convert all models in a system to per-unit values using metadata-driven conversion.

This is the **high-level API** for per-unit conversion. It automatically:
- Extracts metadata for each model type
- Calculates conversion coefficients based on device ratings and external references
- Converts models that have parameters requiring conversion (PQ, PV, Transformer2W)
- Passes through models already in per-unit (Bus, Branch)

# Arguments
- `models::NamedTuple`: Named tuple of model StructArrays (e.g., from `load_models`)
- `system_base::Float64`: System base power in MVA
- `external::NamedTuple`: Pre-retrieved external parameters from retrieve_external_params(PrePerUnit)

# Returns
- Named tuple with fields:
  - `models::NamedTuple`: Models converted to per-unit
  - `coeffs::NamedTuple`: Per-unit conversion coefficients for inverse conversion (pu → physical)

# Example
```julia
using PowerfulCases
using Powerful

# 1. Parse case (physical units: MW, MVAR, kV)
case = parse_network(load("ieee14").raw)
models = load_models(case, all_models())

# 2. Convert to per-unit (on 100 MVA base); `external` comes from
#    retrieve_external_params(..., PrePerUnit), see `SystemModel(case)`
result = convert_system_to_per_unit(models, 100.0, external)
models_pu = result.models
coeffs = result.coeffs  # Stored for inverse conversion in reports
```

# Per-Unit Conversion Details

For PSS/E format:
- **PQ loads**: pl, ql, ip, iq (PowerPU), yp, yq (AdmittancePU) converted
- **PV generators**: pg, qg (PowerPU) converted using mbase rating
- **Transformers**: r1_2, x1_2 (ImpedancePU) converted using sbase1_2 rating
- **Branches**: Already in system per-unit, passed through unchanged
- **Buses**: No conversion needed (provides base voltage reference)
"""
function convert_system_to_per_unit(models::NamedTuple, system_base::Float64, external::NamedTuple)
    converted_pairs = []
    coeff_pairs = []

    for (name, model) in pairs(models)
        # Get base model type (unwraps aliases like PQNum → PQ)
        model_type = get_base_model_type(model)

        # Every model describes its parameters in model_metadata
        metadata = model_metadata(model_type)

        # Only convert if model has parameters requiring conversion
        if !isempty(metadata.params) || !isempty(metadata.derived_params)
            # Calculate conversion coefficients using pre-retrieved external params
            coeffs = calc_pu_coefficients(metadata, model, name, system_base, external)

            # Apply conversion
            model_pu = convert_to_per_unit(model, coeffs)
            push!(converted_pairs, name => model_pu)
            push!(coeff_pairs, name => coeffs)  # Store coefficients for inverse conversion
        else
            # No parameters to convert (e.g., Bus, Branch already in pu)
            push!(converted_pairs, name => model)
        end
    end

    models_pu = NamedTuple(converted_pairs)
    coeffs_nt = NamedTuple(coeff_pairs)

    return (models=models_pu, coeffs=coeffs_nt)
end

"""
    convert_linked_models(sys, new_models, links) -> (models = NamedTuple, coeffs = NamedTuple)

Per-unit conversion of the models among `new_models` that declare a rating from the device
they are linked to (`LinkedDevice`): their parameters are given per unit on that device's
base (for a converter, its `MBASE`), and are converted to the system base of `sys` through
the same coefficients as any model's (`calc_pu_coefficients`, with the linked device's
rating as `Sn`), which `coeffs` returns for `sys.pu_coeffs`. `links` (`ModelLinks`) link
their instances to the devices of `sys`. The other models are returned as they are (they
are per unit on the system base already), and have no entry in `coeffs`. A rating that
cannot be resolved is an `ArgumentError` (`linked_device_ratings`).
"""
function convert_linked_models(sys::SystemModel, new_models::NamedTuple, links::ModelLinks)
    converted = Pair{Symbol, Any}[]
    coeff_pairs = Pair{Symbol, Any}[]
    for (name, models) in pairs(new_models)
        metadata = model_metadata(get_base_model_type(models))
        if declares_linked_rating(metadata)
            ratings = linked_device_ratings(sys, name, models, metadata, links)
            coeffs = calc_pu_coefficients(metadata, models, name, Float64(sys.settings.mva_base), NamedTuple();
                                          linked_ratings = ratings)
            push!(converted, name => convert_to_per_unit(models, coeffs))
            push!(coeff_pairs, name => coeffs)
        else
            push!(converted, name => models)
        end
    end
    return (models = NamedTuple(converted), coeffs = NamedTuple(coeff_pairs))
end
