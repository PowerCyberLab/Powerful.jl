"""
Text formatting for power flow reports.

This module provides text-based rendering of power flow analysis results.
The core design separates data collection (PowerFlowReport) from presentation (format_text),
enabling support for multiple output formats (text, markdown, HTML, LaTeX) via PrettyTables.jl.
"""

using OrderedCollections
using Printf

# ============================================================================
# PowerFlowReport Structure
# ============================================================================

"""
    PowerFlowReport

Structured power flow report data, format-agnostic.

Contains all information needed to generate a report in any format (text, markdown, JSON, etc.).
The report captures system metadata, convergence information, and model result tables.
It covers every model of the system that has outputs (see [`get_default_outputs`](@ref)).

# Fields
- `counts::OrderedDict{Symbol, Int}`: Instances of each reported model, in presentation order
- `converged::Bool`: Whether the power flow converged
- `iterations::Int`: Number of solver iterations
- `final_residual::Float64`: Final residual norm (infinity norm)
- `sections::OrderedDict{Symbol, NamedTuple}`: Model data tables in presentation order

# Example
```julia
# Build report from solution
report = build_powerflow_report(sol, config)

# Render as text
text = format_text(report)
println(text)

# Future: other formats
# markdown = format_markdown(report)
# json = format_json(report)
```
"""
struct PowerFlowReport
    # Instances of each reported model (in presentation order)
    counts::OrderedDict{Symbol, Int}

    # Convergence info
    converged::Bool
    iterations::Int
    final_residual::Float64

    # Model data tables (in presentation order)
    sections::OrderedDict{Symbol, NamedTuple}
end

# ============================================================================
# Report Construction
# ============================================================================

"""
    build_powerflow_report(sol, config::ReportConfig) -> PowerFlowReport

Build a structured power flow report from a solution.

Extracts system metadata, convergence statistics, and model result tables
according to the provided configuration. The resulting report is format-agnostic
and can be rendered in multiple ways.

# Arguments
- `sol`: Power flow solution from solve(): a `PowerFlowSolution` in any formulation, or
  a `MultiIslandSolution` (its tables merge the islands, and it converged when every
  island did)
- `config::ReportConfig`: Report configuration (units, formatting, columns)

# Returns
`PowerFlowReport` containing all report data

# Example
```julia
prob = ACPowerFlowProblem(sys)
sol = solve(prob, RobustNR())

# Build report with custom config
config = ReportConfig(units=UnitPreferences(voltage=:kv))
report = build_powerflow_report(sol, config)

# Render as text
text = format_text(report)
```
"""
function build_powerflow_report(sol::Union{PowerFlowSolution, MultiIslandSolution}, config::ReportConfig)
    systems = _report_systems(sol)

    # One section per reported model, with its instance count over every island
    counts = OrderedDict{Symbol, Int}()
    sections = OrderedDict{Symbol, NamedTuple}()
    for name in _report_models(systems)
        counts[name] = sum(sys -> haskey(sys.models, name) ? length(sys.models[name]) : 0, systems)
        table = collect_model_table(sol, name; config=config)
        isnothing(table) || (sections[name] = table)
    end

    return PowerFlowReport(counts, sol.converged, sol.iterations, sol.residual_norm, sections)
end

# Models a report covers, in presentation order: every registered model with instances
# in the system and default outputs, the buses first, then the devices at one bus
# (generators, loads, shunts), then those between two buses (branches, transformers)
function _report_models(systems)
    present(M) = any(sys -> haskey(sys.models, nameof(M)) && !isempty(sys.models[nameof(M)]), systems)
    models = [M for M in all_models() if present(M) && !isempty(get_default_outputs(M))]
    n_bus_terminals(M) = count(t -> t.vertex_model === :Bus, get_terminals(M))
    sort!(models; by = n_bus_terminals, alg = MergeSort)  # stable: registration order within a group
    return Symbol[nameof(M) for M in models]
end

# Systems a report covers: the solved system, or every island
_report_systems(sol::PowerFlowSolution) = [sol.prob.sys]
_report_systems(sol::MultiIslandSolution) = sol.islands

# ============================================================================
# Table Formatting (inlined - no separate function needed)
# ============================================================================

# ============================================================================
# Report Rendering
# ============================================================================

"""
    format_text(report::PowerFlowReport; config=Default, backend=:text) -> String

Render a PowerFlowReport in various formats using PrettyTables.jl.

Generates a complete power flow analysis report with:
- System summary (instances of each reported model)
- Convergence information (status, iterations, residual)
- Model result tables (buses, generators, loads, shunts, branches, transformers)

# Arguments
- `report::PowerFlowReport`: Report data structure
- `config::ReportConfig`: Report configuration (for formatting precision)
- `backend::Symbol`: Output format (:text, :markdown, :html, :latex)

# Returns
Formatted string suitable for console or file output

# Example
```julia
config = ReportConfig(formatting=FormattingPreferences(float_precision=2))
report = build_powerflow_report(sol, config)

# Text output (default)
text = format_text(report; config=config)
println(text)

# Markdown output (for documentation)
md = format_text(report; config=config, backend=:markdown)
write("report.md", md)

# HTML output (for web reports)
html = format_text(report; config=config, backend=:html)
write("report.html", html)
```
"""
function format_text(report::PowerFlowReport; config::ReportConfig=Default, backend::Symbol=:text)
    io = IOBuffer()

    println(io, "=" ^ 80)
    println(io, "POWER FLOW ANALYSIS REPORT")
    println(io, "=" ^ 80)
    println(io)

    # System summary: instances of each reported model
    println(io, "System Summary:")
    for (name, n) in report.counts
        println(io, "  ", name, ": ", n)
    end
    println(io)

    # Convergence info
    println(io, "Convergence:")
    println(io, "  Converged: ", report.converged ? "Yes" : "No")
    println(io, "  Iterations: ", report.iterations)
    println(io, "  Final residual: ", @sprintf("%.3e", report.final_residual))
    println(io)

    # Model sections
    for (section_name, table) in report.sections
        println(io, "-" ^ 80)
        println(io, uppercase(string(section_name)), " RESULTS")
        println(io, "-" ^ 80)
        println(io)

        # Delegate table formatting to shared helper (from reports.jl)
        if !isnothing(table) && !isempty(keys(table)) && length(table[first(keys(table))]) > 0
            format_table_pretty(io, table; config=config, backend=backend)
        end

        println(io)
    end

    println(io, "=" ^ 80)

    return String(take!(io))
end

# ============================================================================
# Exports
# ============================================================================

export PowerFlowReport
export build_powerflow_report
export format_text
