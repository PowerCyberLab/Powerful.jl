using Powerful
using Documenter
using JSON

# Include the NotebookConverter module
include("utils/NotebookConverter.jl")
using .NotebookConverter

DocMeta.setdocmeta!(Powerful, :DocTestSetup, :(using Powerful); recursive=true)

# Convert notebooks before building docs
output_dir = NotebookConverter.convert_notebooks()

# Define the pages
makedocs(;
    modules=[Powerful],
    authors="Hantao Cui <cuihantao@gmail.com>",
    repo="https://github.com/sensl/Powerful.jl/blob/{commit}{path}#{line}",
    sitename="Powerful.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sensl.gitlab.io/Powerful.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Index" => "index_page.md",
        "Notebooks" => "$output_dir/index.md",
        "API Reference" => "api.md",
    ],
    checkdocs=:exports,
)
