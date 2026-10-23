# Build helpers for the tutorial notebooks.
#
# The Literate.jl sources `notebooks/NN_*.jl` are the single source of truth. From them
# - `generate_ipynb()` writes the Jupyter notebooks `notebooks/NN_*.ipynb` (without
#   outputs), run by `convert_notebooks.jl`;
# - `generate_notebook_pages(src_dir)` writes the documentation pages into
#   `<src_dir>/notebooks/` on every docs build (called by `docs/make.jl`). Documenter
#   executes their code blocks during the build, so the published outputs always come
#   from the current code.
#
# A notebook that needs hardware the documentation builder does not have declares so
# itself, in a `# requires: gpu` comment on the line right after its title (see
# `notebook_requirement`). Its page then shows the code as a plain, non-executed Julia
# block instead of running it.

using Literate

const NOTEBOOK_DIR = @__DIR__

"""
    notebook_sources() -> Vector{String}

Absolute paths of the notebook sources `NN_*.jl`, sorted by number.
"""
function notebook_sources()
    names = sort(filter(f -> occursin(r"^\d{2}_.+\.jl$", f), readdir(NOTEBOOK_DIR)))
    isempty(names) && error("No notebook sources NN_*.jl found in $NOTEBOOK_DIR")
    return [joinpath(NOTEBOOK_DIR, name) for name in names]
end

"""
    notebook_requirement(src) -> Union{Nothing, String}

The notebook's declared requirement, read from a `# requires: ...` comment on the
line right after its title (`# # Title`), or `nothing` if it declares none.
"""
function notebook_requirement(src::AbstractString)
    lines = readlines(src)
    length(lines) < 2 && return nothing
    m = match(r"^# requires: (.+)$", lines[2])
    return m === nothing ? nothing : strip(m[1])
end

requires_gpu(src::AbstractString) = notebook_requirement(src) == "gpu"

"""
    strip_requires_directive(content) -> String

Remove the `# requires: ...` directive line, if present right after the title, before
handing the source to Literate: it is metadata for `notebook_requirement`, not prose,
and must not appear in a generated notebook page or `.ipynb` file.
"""
function strip_requires_directive(content::AbstractString)
    return replace(content, r"^(# #[^\n]*\n)# requires: [^\n]*\n" => s"\1"; count = 1)
end

"""
    notebook_title(src) -> String

Title of a notebook: its first line, which must be the level-1 heading `# # Title`.
"""
function notebook_title(src::AbstractString)
    first_line = readline(src)
    m = match(r"^# # (.+)$", first_line)
    m === nothing && error("First line of $src must be the notebook title `# # Title`, got: $first_line")
    return strip(m[1])
end

# Documentation sources on GitHub, for links from the .ipynb files
const DOCS_SOURCE_URL = "https://github.com/PowerCyberLab/Powerful.jl/blob/main/docs/src"

"""
    notebook_links(content) -> String

Rewrite the documentation-relative links of a notebook source for the `.ipynb` file:
`../guides/x.md` points to the documentation source on GitHub, `NN_name.md` (another
notebook page) to `NN_name.ipynb`.
"""
function notebook_links(content::AbstractString)
    content = replace(content, r"\]\(\.\./([^)]+\.md)\)" => SubstitutionString("]($DOCS_SOURCE_URL/\\1)"))
    return replace(content, r"\]\((\d{2}_[^)/]+)\.md\)" => s"](\1.ipynb)")
end

"""
    generate_ipynb()

Write `NN_*.ipynb` next to each source, without executing the code (no outputs).
"""
function generate_ipynb()
    for src in notebook_sources()
        Literate.notebook(src, NOTEBOOK_DIR; execute = false, credit = false,
                          preprocess = content -> notebook_links(strip_requires_directive(content)))
    end
    return nothing
end

"""
    generate_notebook_pages(src_dir) -> Vector{Pair{String, String}}

Regenerate the documentation pages of all notebooks in `<src_dir>/notebooks/`
(removing whatever was there before) and return the Documenter `pages` entries.
"""
function generate_notebook_pages(src_dir::AbstractString)
    out_dir = joinpath(src_dir, "notebooks")
    rm(out_dir; recursive = true, force = true)
    mkpath(out_dir)

    pages = ["Overview" => "notebooks/index.md"]
    index_items = String[]
    for src in notebook_sources()
        name = splitext(basename(src))[1]
        title = notebook_title(src)
        if requires_gpu(src)
            ## Plain Julia code blocks: shown, not executed by Documenter
            Literate.markdown(src, out_dir; codefence = "```julia" => "```", credit = false,
                              preprocess = strip_requires_directive)
            push!(index_items, "- [$title]($name.md) (requires an NVIDIA GPU; shown without outputs)")
        else
            ## `@example` blocks: executed by Documenter during the build
            Literate.markdown(src, out_dir; credit = false, preprocess = strip_requires_directive)
            push!(index_items, "- [$title]($name.md)")
        end
        push!(pages, title => "notebooks/$name.md")
    end

    write(joinpath(out_dir, "index.md"), """
        # Notebooks

        Tutorials that walk through Powerful.jl step by step. The outputs on these pages
        are produced by running the code when the documentation is built. Pages marked
        as requiring a GPU are shown without outputs.

        Each tutorial is also available as a Julia script (`.jl`) and a Jupyter notebook
        (`.ipynb`) in the `notebooks/` directory of the repository.

        $(join(index_items, "\n"))
        """)
    return pages
end
