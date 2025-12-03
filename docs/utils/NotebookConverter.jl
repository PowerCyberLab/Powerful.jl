"""
Notebook to Markdown conversion utilities for Documenter integration.
"""
module NotebookConverter

using JSON

export convert_notebooks

# Configuration constants
const DEFAULT_NOTEBOOK_DIR = "dev_notebooks"
const DEFAULT_OUTPUT_DIR = "notebooks"
const CODE_BLOCK_START = "```julia"
const CODE_BLOCK_END = "```"
const OUTPUT_BLOCK_START = "```"
const OUTPUT_BLOCK_END = "```"

# Regex for ANSI escape codes used for colors and formatting in terminal output
const ANSI_ESCAPE_REGEX = r"\x1b\[[0-9;]*[a-zA-Z]|\[[0-9]+m|\[[0-9]+"

"""
    strip_ansi_codes(text)

Remove ANSI escape sequences (color codes, formatting) from text.
"""
function strip_ansi_codes(text)
    return replace(text, ANSI_ESCAPE_REGEX => "")
end

"""
    extract_title_from_notebook(notebook_data)

Extract a title from notebook, either from the first markdown cell or generate from filename.
"""
function extract_title_from_notebook(notebook_data, fallback_title)
    # Get all cells and find first markdown cell with a title
    cells = get(notebook_data, "cells", [])

    for cell in cells
        if get(cell, "cell_type", "") == "markdown"
            source = get(cell, "source", [])
            if !isempty(source) && startswith(source[1], "# ")
                return strip(source[1][3:end])
            end
        end
    end

    # No title found, return fallback
    return fallback_title
end

"""
    process_markdown_cell(io, cell)

Process a markdown cell and write it to the output stream.
"""
function process_markdown_cell(io, cell)
    for line in get(cell, "source", [])
        write(io, line)
    end
    write(io, "\n\n")
end

"""
    process_code_cell(io, cell)

Process a code cell and write it with its outputs to the output stream.
"""
function process_code_cell(io, cell)
    # Single admonition block for the code cell and its output
    write(io, """
    !!! note "Code Cell"

    """)

    # Write label for input section
    write(io, "    **Input:**\n\n")

    # Write code content with indentation for the admonition block
    write(io, "    $CODE_BLOCK_START\n")
    for line in get(cell, "source", [])
        # Add 4 spaces for proper admonition block indentation
        write(io, "    $(line)")
    end
    write(io, "\n    $CODE_BLOCK_END\n\n")

    # Process outputs
    outputs = get(cell, "outputs", [])
    if !isempty(outputs)
        # Add output label within the same block
        write(io, "    **Output:**\n\n")

        # Process outputs with proper indentation
        process_outputs_indented(io, outputs)
    end

    # Add extra newline after the cell
    write(io, "\n")
end

"""
    process_outputs(io, outputs)

Process cell outputs and write them to the output stream (non-indented).
"""
function process_outputs(io, outputs)
    # First try text/plain format (usually most readable)
    for output in outputs
        if haskey(output, "data") && haskey(output["data"], "text/plain")
            write(io, "$OUTPUT_BLOCK_START\n")
            for line in output["data"]["text/plain"]
                write(io, strip_ansi_codes(line))
            end
            write(io, "\n$OUTPUT_BLOCK_END\n\n")
            return  # Exit after finding first suitable output
        end
    end

    # Try stream output format if no text/plain
    for output in outputs
        if haskey(output, "text")
            write(io, "$OUTPUT_BLOCK_START\n")
            for line in output["text"]
                write(io, strip_ansi_codes(line))
            end
            write(io, "\n$OUTPUT_BLOCK_END\n\n")
            return
        end
    end
end

"""
    process_outputs_indented(io, outputs)

Process cell outputs with proper indentation for admonition blocks.
"""
function process_outputs_indented(io, outputs)
    # First try text/plain format (usually most readable)
    for output in outputs
        if haskey(output, "data") && haskey(output["data"], "text/plain")
            write(io, "    $OUTPUT_BLOCK_START\n")
            for line in output["data"]["text/plain"]
                # Clean ANSI escape codes and add indentation
                write(io, "    $(strip_ansi_codes(line))")
            end
            write(io, "\n    $OUTPUT_BLOCK_END\n")
            return  # Exit after finding first suitable output
        end
    end

    # Try stream output format if no text/plain
    for output in outputs
        if haskey(output, "text")
            write(io, "    $OUTPUT_BLOCK_START\n")
            for line in output["text"]
                # Clean ANSI escape codes and add indentation
                write(io, "    $(strip_ansi_codes(line))")
            end
            write(io, "\n    $OUTPUT_BLOCK_END\n")
            return
        end
    end
end

"""
    create_notebook_index(output_dir, notebooks_info)

Create an index page for the notebooks.
"""
function create_notebook_index(output_dir, notebooks_info)
    index_path = joinpath(output_dir, "index.md")

    open(index_path, "w") do io
        write(io, """
        # Development Notebooks

        This section contains Jupyter notebooks that demonstrate the usage and features of Powerful.jl.

        ## Available Notebooks

        """)

        for (filename, title, description) in notebooks_info
            md_name = replace(filename, ".ipynb" => ".md")
            write(io, "- [**$title**]($(md_name)) - $description\n")
        end
    end
end

"""
    convert_notebook(notebook_path, output_path)

Convert a single Jupyter notebook to markdown format.
"""
function convert_notebook(notebook_path, output_path)
    # Parse notebook JSON
    notebook_data = JSON.parsefile(notebook_path)
    notebook_filename = basename(notebook_path)

    # Generate title from filename as fallback
    fallback_title = replace(notebook_filename, ".ipynb" => "")
    title = extract_title_from_notebook(notebook_data, fallback_title)

    # Track if we found a title in the notebook cells
    found_title_in_cell = false
    title_cell_index = -1

    # Find if first cell contains the title
    cells = get(notebook_data, "cells", [])
    for (i, cell) in enumerate(cells)
        if get(cell, "cell_type", "") == "markdown"
            source = get(cell, "source", [])
            if !isempty(source)
                first_line = source[1]
                if startswith(first_line, "# ")
                    found_title_in_cell = true
                    title_cell_index = i
                    break
                end
            end
        end
    end

    # Extract first paragraph for description
    description = ""
    for cell in cells
        if get(cell, "cell_type", "") == "markdown"
            source = get(cell, "source", [])
            if length(source) > 1
                # Skip header and get next paragraph
                for i in 2:length(source)
                    line = source[i]
                    if !isempty(strip(line))
                        description = strip(line)
                        if length(description) > 100
                            description = description[1:97] * "..."
                        end
                        break
                    end
                end
            end
            if !isempty(description)
                break
            end
        end
    end

    # If no description found, use generic one
    if isempty(description)
        description = "Jupyter notebook example"
    end

    # Create markdown file
    open(output_path, "w") do io
        # Only write title if not found in cells
        if !found_title_in_cell
            write(io, "# $title\n\n")
        end

        # Process each cell
        for (i, cell) in enumerate(cells)
            # Skip the title cell if it was processed above
            if i == title_cell_index
                # For title cells, we write the title but skip any other content
                write(io, "# $title\n\n")
                continue
            end

            cell_type = get(cell, "cell_type", "")
            if cell_type == "markdown"
                process_markdown_cell(io, cell)
            elseif cell_type == "code"
                process_code_cell(io, cell)
            end
        end
    end

    return (fallback_title, title, description)
end

"""
    convert_notebooks(;notebooks_dir=DEFAULT_NOTEBOOK_DIR, output_dir=nothing)

Convert all Jupyter notebooks in a directory to markdown files for Documenter.

# Arguments
- `notebooks_dir`: Directory containing Jupyter notebooks
- `output_dir`: Output directory (relative to docs/src), defaults to "notebooks"
"""
function convert_notebooks(;notebooks_dir=DEFAULT_NOTEBOOK_DIR, output_dir=nothing)
    # Get the docs directory (parent of utils)
    docs_dir = dirname(dirname(@__FILE__))

    # Resolve paths
    base_dir = dirname(docs_dir)
    notebooks_path = isabspath(notebooks_dir) ? notebooks_dir : joinpath(base_dir, notebooks_dir)

    # Create output directory
    output_dirname = isnothing(output_dir) ? DEFAULT_OUTPUT_DIR : output_dir
    output_path = joinpath(docs_dir, "src", output_dirname)
    isdir(output_path) || mkdir(output_path)

    # Find all notebooks
    notebooks = filter(f -> endswith(f, ".ipynb"), readdir(notebooks_path))

    # Track notebook info for index creation
    notebooks_info = []

    # Process each notebook
    for notebook in notebooks
        nb_path = joinpath(notebooks_path, notebook)
        out_path = joinpath(output_path, replace(notebook, ".ipynb" => ".md"))

        try
            info = convert_notebook(nb_path, out_path)
            push!(notebooks_info, (notebook, info...))
            @info "Converted $notebook to markdown"
        catch e
            @warn "Failed to convert $notebook: $e"
            # Create fallback markdown file
            open(out_path, "w") do io
                write(io, "# $(replace(notebook, ".ipynb" => ""))\n\n")
                write(io, "This notebook couldn't be automatically converted.")
            end

            # Add to notebook info with placeholder
            push!(notebooks_info, (notebook,
                                  replace(notebook, ".ipynb" => ""),
                                  "Notebook conversion failed"))
        end
    end

    # Create index
    create_notebook_index(output_path, notebooks_info)

    return output_dirname
end

end # module