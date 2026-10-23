# Post-processing transformations for parsed model data.
#
# These functions transform raw parsed data into the final model structure.
# Each function has a single responsibility.

using StructArrays
using PowerfulCore: BUS_SLACK, parse_model, PSSE

export post_process_models

# ============================================================================
# Individual Post-Processors (Single Responsibility)
# ============================================================================

"""
    split_pv_and_slack(models::NamedTuple) -> NamedTuple

Split PV generators into PV and Slack based on connected bus type.

Generators at slack buses (ide == BUS_SLACK) are moved from PV to Slack model.
This enables clean separation of equation semantics:
- PV: Power injection + voltage control (SharedRes(:add), SharedRes(:set))
- Slack: Angle reference + voltage reference (SharedRes(:set), SharedRes(:set))

Returns models unchanged if PV or Bus not present.
"""
function split_pv_and_slack(models::NamedTuple)
    # Skip if prerequisites not present
    if !haskey(models, :PV) || !haskey(models, :Bus)
        return models
    end

    pv_generators = models.PV
    buses = models.Bus

    # Identify generators at slack buses
    bus_map = build_index_map(buses, :i)
    is_slack = map(gen_bus_num -> buses.ide[bus_map[(gen_bus_num,)]] == BUS_SLACK, pv_generators.i)

    # Create masks
    slack_mask = is_slack
    pv_mask = .!is_slack

    # Filter PV generators (non-slack)
    pv_only = pv_generators[pv_mask]

    # Construct Slack generators from filtered PV data
    slack_only = _construct_slack_from_pv(pv_generators, slack_mask)

    return merge(models, (PV = pv_only, Slack = slack_only))
end

"""
    _construct_slack_from_pv(pv_generators, mask) -> StructArray{Slack}

Construct Slack model instances from PV generator data.

Internal helper that handles the field-by-field construction.
"""
function _construct_slack_from_pv(pv_generators, mask)
    StructArray{Slack{Float64}}((
        i = pv_generators.i[mask],
        id = pv_generators.id[mask],
        pg = pv_generators.pg[mask],
        qg = pv_generators.qg[mask],
        qt = pv_generators.qt[mask],
        qb = pv_generators.qb[mask],
        vs = pv_generators.vs[mask],
        ireg = pv_generators.ireg[mask],
        mbase = pv_generators.mbase[mask],
        zr = pv_generators.zr[mask],
        zx = pv_generators.zx[mask],
        rt = pv_generators.rt[mask],
        xt = pv_generators.xt[mask],
        gtap = pv_generators.gtap[mask],
        status = pv_generators.status[mask],
        rmpct = pv_generators.rmpct[mask],
        pt = pv_generators.pt[mask],
        pb = pv_generators.pb[mask],
        o1 = pv_generators.o1[mask],
        f1 = pv_generators.f1[mask],
        o2 = pv_generators.o2[mask],
        f2 = pv_generators.f2[mask],
        o3 = pv_generators.o3[mask],
        f3 = pv_generators.f3[mask],
        o4 = pv_generators.o4[mask],
        f4 = pv_generators.f4[mask],
        wmod = pv_generators.wmod[mask],
        wpf = pv_generators.wpf[mask],
    ))
end

"""
    handle_3w_conversion(raw_data, models::NamedTuple, system_base::Float64) -> NamedTuple

Convert 3-winding transformers to 2-winding equivalents with virtual buses.

Parses Transformer3W from raw data, converts to star configuration,
and merges resulting 2W transformers and virtual buses into models.

Returns models unchanged if no 3W transformers present.
"""
function handle_3w_conversion(raw_data, models::NamedTuple, system_base::Float64)
    # Parse 3W transformer data (not a registered model, only used for conversion)
    t3w_data = parse_model(Transformer3W, raw_data)

    if isempty(t3w_data)
        return models
    end

    # Add 3W data temporarily for conversion
    models_with_3w = merge(models, (Transformer3W = t3w_data,))

    # Convert 3W to 2W + virtual buses
    converted_models, _ = convert_3w_to_star!(models_with_3w, system_base)

    # Remove Transformer3W from result (fully converted to 2W)
    return _remove_key(converted_models, :Transformer3W)
end

"""
    _remove_key(nt::NamedTuple, key::Symbol) -> NamedTuple

Remove a key from a NamedTuple, returning a new NamedTuple without that key.
"""
function _remove_key(nt::NamedTuple, key::Symbol)
    keys_to_keep = [k for k in keys(nt) if k != key]
    values_to_keep = [nt[k] for k in keys_to_keep]
    return NamedTuple{tuple(keys_to_keep...)}(tuple(values_to_keep...))
end

# ============================================================================
# Coordinator Function
# ============================================================================

"""
    post_process_models(raw_data, models::NamedTuple, model_types, system_base::Float64) -> NamedTuple

Apply all post-processing transformations to parsed models.

This is the single entry point for Formats module to call. It coordinates
individual post-processors based on which model types were requested.

# Arguments
- `raw_data`: Network (needed for 3W transformer parsing)
- `models`: NamedTuple of parsed models from load_models
- `model_types`: Requested model types (to determine which processors to run)
- `system_base`: System MVA base (for 3W conversion)

# Returns
NamedTuple with all transformations applied
"""
function post_process_models(raw_data, models::NamedTuple, model_types, system_base::Float64)
    result = models

    # 3W transformer conversion (if Transformer2W requested)
    if _type_requested(model_types, :Transformer2W)
        result = handle_3w_conversion(raw_data, result, system_base)
    end

    # PV/Slack split (if either requested)
    if _type_requested(model_types, :PV) || _type_requested(model_types, :Slack)
        result = split_pv_and_slack(result)
    end

    return result
end

"""
    _type_requested(model_types, name::Symbol) -> Bool

Check if a model type with given name is in the requested types.
"""
function _type_requested(model_types, name::Symbol)
    any(T -> nameof(T) == name, model_types)
end
