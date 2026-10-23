using Powerful
using Powerful.Models
using PowerfulCore
using PowerfulPF
using PowerfulTDS
using PowerfulCases
using Documenter
using TOML

# Standalone packages describe their own documentation pages in docs/pages/*.toml:
# which package to load, which guide/developer/API pages they add, where (an `order`
# relative to the others), and the text for the generated lists below. make.jl
# discovers whichever files are present (readdir, not a list of names) and never
# names a package or feature itself; a package some checkouts don't have is simply a
# file this checkout doesn't have.
const PAGES_DIR = joinpath(@__DIR__, "pages")
const FEATURES = sort(
    [TOML.parsefile(joinpath(PAGES_DIR, f)) for f in readdir(PAGES_DIR) if endswith(f, ".toml")];
    by = f -> f["order"],
)

for feature in FEATURES
    @eval using $(Symbol(feature["package"]))
end

DocMeta.setdocmeta!(Powerful, :DocTestSetup, :(using Powerful); recursive=true)

# --- Build from a temporary copy of docs/src ---------------------------------------
#
# Several pages are filled in at build time (see below) and the notebook pages are
# regenerated on every build. None of that may touch the tracked docs/src/, so the
# whole tree is copied once into a gitignored directory and every read/write below
# targets the copy. The copy is nested under docs/src/ so its path relative to
# BUILD_ROOT_DIR is "docs/src/...", same as the tracked files' path relative to the
# repository root; combined with the `remotes` entry below, source/edit links still
# point at the real docs/src/ files in the real repository, not at the gitignored copy.
const BUILD_ROOT_DIR = joinpath(@__DIR__, ".build_src")
const BUILD_SRC_DIR = joinpath(BUILD_ROOT_DIR, "docs", "src")
rm(BUILD_ROOT_DIR; recursive = true, force = true)
mkpath(dirname(BUILD_SRC_DIR))
cp(joinpath(@__DIR__, "src"), BUILD_SRC_DIR)

# Notebook pages: regenerated on every build from the Literate.jl sources
# notebooks/NN_*.jl into <build_src>/notebooks/. Documenter runs their code while
# building, so outputs always come from the current code. A notebook that needs
# hardware this builder doesn't have declares so itself (see notebooks/build.jl) and
# is shown as non-executed code.
ENV["GKSwstype"] = "100"  # headless GR backend for the Plots used by the notebooks
include(joinpath(@__DIR__, "..", "notebooks", "build.jl"))
notebook_pages = generate_notebook_pages(BUILD_SRC_DIR)

# --- Generated passages in shared pages -------------------------------------------
#
# guides/index.md, index.md and api.md list guides, packages and API pages by name.
# Each list ends in one `<!--GENERATED:key-->` placeholder; this build appends
# whatever the discovered features contribute (sorted by `order`), so the committed
# source never names a package or feature itself, and a feature that declares nothing
# for a given list (e.g. one not meant to appear in a table) simply contributes
# nothing to it.

join_and(items::Vector{String}) =
    isempty(items) ? "" : length(items) == 1 ? items[1] : join(items[1:end-1], ", ") * " and " * items[end]

function render_generated!(relpath::AbstractString, substitutions::Dict{String, String})
    path = joinpath(BUILD_SRC_DIR, relpath)
    content = read(path, String)
    for (key, replacement) in substitutions
        placeholder = "<!--GENERATED:$(key)-->"
        content = isempty(replacement) ? replace(content, placeholder * "\n" => "") :
                                          replace(content, placeholder => replacement)
    end
    write(path, content)
    return nothing
end

guide_relpath(file::AbstractString) = replace(file, r"^guides/" => "")

function render_guide_list(list_name::AbstractString)
    lines = String[]
    for feature in FEATURES
        guide = get(feature, "guide", nothing)
        guide === nothing && continue
        get(guide, "list", nothing) == list_name || continue
        push!(lines, "- **[$(guide["title"])]($(guide_relpath(guide["file"])))** - $(guide["list_summary"])")
    end
    return join(lines, "\n")
end

render_generated!("guides/index.md", Dict(
    "guides-studies" => render_guide_list("studies"),
    "guides-topics" => render_guide_list("topics"),
))

let
    study_labels = [feature["guide"]["study_label"] for feature in FEATURES
                    if get(feature, "guide", nothing) !== nothing && haskey(feature["guide"], "study_label")]
    more_studies = isempty(study_labels) ? "" : ", " * join_and(study_labels)

    example_bodies = [feature["example"]["body"] for feature in FEATURES if haskey(feature, "example")]

    package_rows = String[]
    for feature in FEATURES
        guide = get(feature, "guide", nothing)
        guide === nothing && continue
        haskey(guide, "package_summary") || continue
        push!(package_rows, "| `$(feature["package"])` | $(guide["package_summary"]), loaded separately |")
    end

    render_generated!("index.md", Dict(
        "more-studies" => more_studies,
        "more-examples" => join(example_bodies, "\n"),
        "packages-table" => join(package_rows, "\n"),
    ))
end

let
    # The same "does this feature have a package_summary" test as the index.md
    # packages table selects which features are "standalone packages" for this
    # sentence, so the two generated lists always agree with each other.
    standalone = [feature["package"] for feature in FEATURES
                  if get(feature, "guide", nothing) !== nothing && haskey(feature["guide"], "package_summary")]
    push!(standalone, "PowerfulCases")

    api_rows = String[]
    for feature in FEATURES
        api = get(feature, "api", nothing)
        api === nothing && continue
        for page in api["pages"]
            haskey(page, "summary") || continue
            title = get(page, "title", api["title"])
            push!(api_rows, "| [$title]($(page["file"])) | $(page["summary"]) ($(feature["package"])) |")
        end
    end

    render_generated!("api.md", Dict(
        "standalone-intro" => join_and(standalone) * " are standalone packages that are loaded explicitly.",
        "api-table" => join(api_rows, "\n"),
    ))
end

# --- Pages array ---------------------------------------------------------------
#
# Each entry below carries the same `order` a feature's TOML does, so the core pages
# and the discovered features merge with one sort instead of a branch per feature.
# A feature's own pages are inserted as a block (developer pages after "TDS Event
# System", API pages after "Time-Domain Simulation") rather than finely interleaved
# with the core pages at that feature's exact original position, which keeps the
# insertion itself generic instead of position-per-feature.

core_guides = [
    (order = 100, title = "Overview", file = "guides/index.md"),
    (order = 200, title = "Getting Started", file = "guides/getting_started.md"),
    (order = 300, title = "Power Flow Analysis", file = "guides/powerflow.md"),
    (order = 400, title = "Multi-Island Systems", file = "guides/multi_island.md"),
    (order = 500, title = "DC Analysis", file = "guides/dc_analysis.md"),
    (order = 700, title = "Time-Domain Simulation", file = "guides/tds.md"),
    (order = 900, title = "Dynamic Models", file = "guides/models.md"),
    (order = 1200, title = "Eigenvalue Analysis", file = "guides/eigenvalue.md"),
    (order = 1300, title = "Frequency Response", file = "guides/bode.md"),
    (order = 1400, title = "Automatic Differentiation", file = "guides/autodiff.md"),
    (order = 1600, title = "Configuration", file = "guides/configuration.md"),
    (order = 1700, title = "Accessing Results", file = "guides/accessing_variables.md"),
    (order = 1800, title = "Report Customization", file = "guides/reports.md"),
]
feature_guides = [(order = feature["order"], title = feature["guide"]["title"], file = feature["guide"]["file"])
                   for feature in FEATURES if haskey(feature, "guide")]
guide_pages = [g.title => g.file for g in sort(vcat(core_guides, feature_guides); by = g -> g.order)]

pages_array = [
    "Home" => "index.md",
    "Guides" => guide_pages,
]

push!(pages_array, "Notebooks" => notebook_pages)

developer_pages = [
    "Overview" => "developer/overview.md",
    "Three-Layer Indexing" => "developer/indexing.md",
    "Multi-Physics Formulations" => "developer/formulations.md",
    "Computation Architecture" => "developer/workspaces.md",
    "Implementing Models" => "developer/models.md",
    "TDS Architecture" => "developer/tds.md",
    "TDS Event System" => "developer/events.md",
]
for feature in FEATURES
    for page in get(feature, "developer", [])
        push!(developer_pages, page["title"] => page["file"])
    end
end
append!(developer_pages, [
    "Per-Unit & Multi-Island" => "developer/miscellaneous.md",
    "Docstring Guidelines" => "developer/docstrings.md",
])

api_pages = [
    "Overview" => "api.md",
    "Powerful" => "api/powerful.md",
    "Core" => [
        "Types and Variables" => "api/foundation.md",
        "Metadata and Traits" => "api/metadata.md",
        "Formulations, Units and Registry" => "api/formulations.md",
        "DAE and Addressing" => "api/dae.md",
        "Jacobians and Newton" => "api/ad.md",
        "Time Series and Injection Recording" => "api/recording.md",
        "System" => "api/system.md",
        "Islands and Outages" => "api/outages.md",
    ],
    "Power Flow" => [
        "Network Models" => "api/models.md",
        "Power Flow" => "api/powerflow.md",
        "Results and Reports" => "api/results.md",
    ],
    "Time-Domain Simulation" => [
        "Dynamic Models" => "api/dynamic_models.md",
        "Simulation" => "api/tds.md",
        "Validation and Small-Signal Analysis" => "api/small_signal.md",
    ],
]
for feature in FEATURES
    api = get(feature, "api", nothing)
    api === nothing && continue
    pages = api["pages"]
    entry = length(pages) == 1 ? api["title"] => pages[1]["file"] :
                                  api["title"] => [get(p, "title", api["title"]) => p["file"] for p in pages]
    push!(api_pages, entry)
end
push!(api_pages, "Test Cases" => "api/cases.md")

append!(pages_array, [
    "Developer Documentation" => developer_pages,
    "API Reference" => api_pages,
    "Index" => "index_page.md",
])

# index_page.md (the alphabetical `@index`) lists every API page that has a table
# summary (the same ones api.md's table does); regenerated here so it never names a
# page a checkout doesn't have.
let
    api_index_pages = [
        "api/powerful.md", "api/foundation.md", "api/metadata.md", "api/formulations.md",
        "api/dae.md", "api/ad.md", "api/recording.md", "api/system.md", "api/outages.md",
        "api/models.md", "api/powerflow.md", "api/results.md", "api/dynamic_models.md",
        "api/tds.md", "api/small_signal.md",
    ]
    for feature in FEATURES
        api = get(feature, "api", nothing)
        api === nothing && continue
        for page in api["pages"]
            haskey(page, "summary") && push!(api_index_pages, page["file"])
        end
    end
    push!(api_index_pages, "api/cases.md")

    pages_list = join(("    \"$p\"," for p in api_index_pages), "\n")
    write(joinpath(BUILD_SRC_DIR, "index_page.md"), """
        # Index

        Alphabetical index of all documented types, functions, macros and constants.

        ```@index
        Pages = [
        $pages_list
        ]
        ```
        """)
end

# Modules passed to `checkdocs=:exports`: the core packages, plus any loaded feature.
modules = Module[Powerful, PowerfulCore, PowerfulPF, PowerfulTDS, PowerfulCases]
for feature in FEATURES
    push!(modules, getfield(Main, Symbol(feature["package"])))
end

# Build documentation
makedocs(;
    modules=modules,
    authors="Hantao Cui <cuihantao@gmail.com>",
    repo="https://github.com/PowerCyberLab/Powerful.jl/blob/{commit}{path}#{line}",
    remotes=Dict(
        # Without this, source/edit links would point at the gitignored build copy
        # (docs/.build_src/docs/src/...), which doesn't exist in the real repository.
        BUILD_ROOT_DIR => (Documenter.Remotes.URL("https://github.com/PowerCyberLab/Powerful.jl/blob/{commit}{path}#{line}"), "main"),
        # PowerfulCases is a separate repository; where it's installed from a package
        # URL rather than a git checkout, Documenter cannot auto-detect its remote.
        pkgdir(PowerfulCases) => (Documenter.Remotes.URL("https://github.com/cuihantao/PowerfulCases/blob/{commit}{path}#{line}"), "main"),
    ),
    sitename="Powerful.jl",
    source=joinpath(".build_src", "docs", "src"),
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://powercyberlab.github.io/Powerful.jl",
        edit_link = "main",
        repolink = "https://github.com/PowerCyberLab/Powerful.jl",
        assets=String[],
        collapselevel = 1,  # Collapse submenu items by default
        inventory_version = string(pkgversion(Powerful)),
        size_threshold_ignore = ["index_page.md"],  # alphabetical index of every documented name
        # The API reference documents every exported name of the loaded packages, which
        # gives a search index of about 1 MiB in the full checkout
        search_size_threshold_warn = 2^21,
    ),
    pages=pages_array,
    checkdocs=:exports,
)

deploydocs(;
    repo = "github.com/PowerCyberLab/Powerful.jl.git",
    devbranch = "main",
)
