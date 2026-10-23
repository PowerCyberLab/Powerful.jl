# PSS/E RAW format loading
#
# This module handles loading models from Network (PSS/E RAW format).
# Individual models define their own parse_model methods; this file
# provides the top-level load_models function that orchestrates them.

# Import load_models from PowerfulCore to extend it (stub defined there)
import PowerfulCore: load_models

"""
    RAW_TABLE_POWER

The class of every record table of a PSS/E RAW case (the `PowerFlowData.Records` fields of a
`Network`), stated once. A table whose records carry power (they change the power flow) maps
to the function giving, for a parsed case, which of its records do, by their status and
setpoint; a table of data only maps to `nothing`.

- Carry power: buses (the network's nodes), loads, fixed shunts, generators, branches and
  transformers in service; switched shunts in service with a nonzero initial admittance `binit`;
  two-terminal DC lines not blocked (`mdc ≠ 0`) with a nonzero setpoint `setvl`; VSC DC lines in
  service (`mdc ≠ 0`); multi-terminal DC lines not blocked (`mdc ≠ 0`); FACTS devices in service
  (`mode ≠ 0`); impedance correction tables that an in-service transformer references (`tab1`,
  `tab2`, `tab3`).
- Data only: area interchanges (the schedules of area interchange control), multi-section line
  groups (groups of branches read one by one), zones, inter-area transfers and owners.

`check_unread_records` refuses a case with power-carrying records that no model reads, and
warns of the rest.
"""
const RAW_TABLE_POWER = (
    buses = net -> trues(length(net.buses)),
    loads = net -> net.loads.status,
    fixed_shunts = net -> net.fixed_shunts.status,
    generators = net -> net.generators.stat,
    branches = net -> net.branches.st,
    transformers = net -> net.transformers.stat .!= 0,
    area_interchanges = nothing,
    two_terminal_dc = net -> (net.two_terminal_dc.mdc .!= 0) .& (net.two_terminal_dc.setvl .!= 0),
    vsc_dc = net -> net.vsc_dc.mdc .!= 0,
    switched_shunts = net -> net.switched_shunts.stat .& (net.switched_shunts.binit .!= 0),
    impedance_corrections = net -> referenced_impedance_tables(net),
    multi_terminal_dc = net -> [line.line_id.mdc != 0 for line in net.multi_terminal_dc.lines],
    multi_section_lines = nothing,
    zones = nothing,
    area_transfers = nothing,
    owners = nothing,
    facts = net -> net.facts.mode .!= 0,
)

# The impedance correction tables an in-service transformer references (by its windings' tab1, tab2, tab3)
function referenced_impedance_tables(net::Network)
    t = net.transformers
    used = Set{Int}()
    for k in eachindex(t.i)
        t.stat[k] == 0 && continue
        for tab in (t.tab1[k], t.tab2[k], t.tab3[k])
            ismissing(tab) || tab == 0 || push!(used, tab)
        end
    end
    return [Int(n) in used for n in net.impedance_corrections.i]
end

# How a refusal names record k of an unread table
function record_label(net::Network, table::Symbol, k::Int)
    records = getfield(net, table)
    table === :multi_terminal_dc && return strip(String(records.lines[k].line_id.name))
    table in (:two_terminal_dc, :vsc_dc, :facts) && return strip(String(records.name[k]))
    table === :switched_shunts && return "bus $(records.i[k])"
    table === :impedance_corrections && return "table $(records.i[k])"
    table in (:loads, :fixed_shunts, :generators) && return "bus $(records.i[k]) id $(strip(String(records.id[k])))"
    table in (:branches, :transformers) && return "$(records.i[k])-$(records.j[k]) ckt $(strip(String(records.ckt[k])))"
    table === :buses && return "bus $(records.i[k])"
    return "record $k"
end

"""
    check_unread_records(raw_data::Network, models)

Check the record tables of `raw_data` with records that none of `models` reads
(`format_records`), by their class (`RAW_TABLE_POWER`):
- records that carry power refuse the case, naming each table with its counts and records:
  its power flow would be solved without them (e.g. the DC lines of a case, which no model
  reads);
- the records of a power-carrying table that carry none (out of service, or at a zero
  setpoint) and the records of a data-only table are warned of, with their counts and why,
  as the DYR parser warns of the models it ignores;
- a table whose class is not stated is refused.
Refuses a model of `models` that does not declare the records it reads.
"""
function check_unread_records(raw_data::Network, models)
    undeclared = [nameof(M) for M in models if !declares_format_records(M, PSSE)]
    isempty(undeclared) || throw(ArgumentError(
        "These models read PSS/E RAW data (supports_format) but do not declare the records they read: " *
        join(undeclared, ", ") * " (a method of PowerfulCore.format_records(::Type{M}, ::Type{PSSE}) " *
        "returning the record tables of a PowerFlowData.Network the model reads, e.g. (:buses,))"))
    read = Set{Symbol}(r for M in models for r in format_records(M, PSSE))
    unread = sort([name => n for (name, n) in record_counts(raw_data) if n > 0 && name ∉ read]; by = p -> -p.second)

    refused, no_power, data_only = String[], String[], String[]
    for (name, n) in unread
        if !haskey(RAW_TABLE_POWER, name)
            push!(refused, "$name ($n records; the class of this table, power-carrying or data only, is not stated in RAW_TABLE_POWER)")
            continue
        end
        carries = RAW_TABLE_POWER[name]
        if carries === nothing
            push!(data_only, "$name ($n)")
            continue
        end
        power = findall(carries(raw_data))
        if !isempty(power)
            labels = [record_label(raw_data, name, k) for k in power[1:min(end, 10)]]
            push!(refused, "$name ($(length(power)) of $n records carry power: " * join(labels, ", ") *
                           (length(power) > 10 ? ", ..." : "") * ")")
        end
        length(power) < n && push!(no_power, "$name ($(n - length(power)))")
    end
    isempty(refused) || throw(ArgumentError(
        "This RAW case has records that carry power and that no model reads, so its power flow would be " *
        "solved without them: " * join(refused, "; ")))
    isempty(no_power) || @warn "Unsupported RAW records ignored, as they carry no power (out of service, " *
        "or at a zero setpoint): " * join(no_power, ", ")
    isempty(data_only) || @warn "Unsupported RAW records ignored (data only): " * join(data_only, ", ")
    return nothing
end

"""
$(SIGNATURES)

Load models from a Network for specified model types. Refuses a case with power-carrying
records that none of the models reads, and warns of the other unread records
(`check_unread_records`).

Automatically handles:
- 3-winding transformer conversion to 2-winding transformers with virtual buses
- PV/Slack generator splitting based on bus type

Post-processing logic is delegated to `post_process_models` for single responsibility.

# Arguments
- `raw_data::Network`: Parsed network data (from `parse_network`)
- `model_types`: Vector or Tuple of model types to load (e.g., `all_models()`)

# Returns
- `NamedTuple` of model StructArrays, keyed by model name

# Example
```julia
using PowerfulPF

case = parse_network("cases/ieee14.raw")
models = load_models(case, all_models())
```
"""
function load_models(raw_data::Network, model_types)
    # The models read the PSS/E revision 33 record layouts. PowerFlowData parses
    # revision 33 as such and every other revision with the revision 30 layout, which
    # does not match revisions 31 and 32 (e.g. their bus records have no GL, BL).
    raw_data.version == 33 || throw(ArgumentError(
        "PSS/E RAW revision $(raw_data.version) is not supported: Powerful reads " *
        "revision 33 files (convert the case to revision 33)"))

    # Extract system base MVA from case data (single source of truth)
    system_base = raw_data.caseid.sbase

    # Filter to only models that support PSSE format, and check the records none of them reads
    supported = filter(T -> is_format_supported(T, PSSE), model_types)
    check_unread_records(raw_data, supported)

    # Parse all supported models
    models = NamedTuple(nameof(T) => parse_model(T, raw_data) for T in supported)

    # Delegate all post-processing to Models module
    models = post_process_models(raw_data, models, model_types, system_base)

    # Filter out empty models (no instances) for cleaner downstream handling
    non_empty_models = filter(pair -> !isempty(pair.second), pairs(models))

    return NamedTuple(non_empty_models)
end
