# Notebooks

Step-by-step tutorials for Powerful.jl.

| Source | Topic | Executed in docs build |
|--------|-------|------------------------|
| `01_getting_started.jl` | Load a case, solve power flow, inspect and modify results | yes |
| `02_formulations_and_solvers.jl` | AC/DC formulations, solver algorithms, KLU orderings, Jacobian structure, convergence | yes |
| `03_multi_island_power_flow.jl` | Power flow on a network with disconnected islands | yes |
| `04_tds_initialization_validation.jl` | Checking the initial conditions of a time-domain simulation | yes |

## Files

- `NN_*.jl`: notebook sources in [Literate.jl](https://fredrikekre.github.io/Literate.jl/v2/fileformat/)
  format. **These are the single source of truth; edit only these.**
- `NN_*.ipynb`: Jupyter notebooks generated from the sources, committed **without outputs**
  (run them in Jupyter to see results). Regenerate them after editing a source. In the
  `.ipynb` files, links to documentation pages (`../guides/x.md`) point to the
  documentation sources on GitHub, and links to other notebooks point to their `.ipynb`.
- `build.jl`: shared helpers: the notebook list, titles, `.ipynb` and docs page generation.
- `convert_notebooks.jl`: regenerates the `.ipynb` files.

The page title of each notebook is its first line, `# # Title`, and the numeric prefix
sets the order. A new `NN_*.jl` file is picked up automatically by both the `.ipynb`
generation and the docs build.

## Running a Notebook

The notebooks use Plots.jl and DataFrames.jl, which are dependencies of the `docs`
environment. From the repository root:

```bash
julia --project=docs notebooks/01_getting_started.jl
```

or open the `.ipynb` file in Jupyter with a Julia kernel whose environment provides
Powerful.jl, PowerfulCases, Plots and DataFrames.

A notebook that needs extra hardware or packages names them in its own header (see
`notebook_requirement` in `build.jl`) and may need a different project environment to
run standalone, e.g. the repository root environment instead of `docs`.

## Regenerating the `.ipynb` Files

After editing a `.jl` source, regenerate the notebooks and commit both files:

```bash
julia --project=docs notebooks/convert_notebooks.jl
```

## Notebooks in the Documentation

`docs/make.jl` regenerates the notebook pages into `docs/src/notebooks/` (gitignored,
cleared first) on every documentation build, so local and CI builds use the same
sources and never pick up stale pages. The code blocks are Documenter `@example`
blocks, executed during the build, so the published outputs always come from the
current code.

A notebook that needs hardware the documentation builder does not have declares so in
a `# requires: ...` comment on the line right after its title (`notebook_requirement`
in `build.jl`). Its page then shows the code as a plain Julia block without outputs,
with a note saying so.

## Authoring Notes

### Comments Inside Code Blocks

Literate.jl treats every line starting with `# ` as markdown, even inside a function,
which splits the code into separate cells and breaks the function. Use `##` for
comments inside code; Literate.jl turns them into regular `#` comments:

```julia
# BAD - breaks the function into multiple cells
function my_function()
    x = 1
    # This comment creates a new cell!
    y = 2
end

# GOOD - keeps function intact
function my_function()
    x = 1
    ## This is a regular Julia comment
    y = 2
end
```

Trailing comments after code on the same line (`x = 1  # note`) are fine.

### Headings

Use exactly one level-1 heading (`# # Title`, the first line), level 2 (`# ## Part`)
for major sections and level 3 or deeper for subsections. Documenter shows level 1
and 2 in the sidebar; multiple level-1 headings flatten the table of contents.

### Documentation-Only Content

Lines prefixed with `#md` appear only on the documentation page, not in the `.ipynb`
file — useful for a `!!! note` admonition about a requirement the notebook declares
in its header.

### Printed Output Next to a Returned Value

On a documentation page, Documenter shows what a code block prints only when the
block's last expression returns `nothing`. If a block prints and then ends with a
value (a plot, a table), put a `#-` line after the prints to split it into two blocks.

### Side Effects

Notebooks run during the documentation build. Do not write files into the repository;
use `mktempdir()` for anything written to disk (see `04_tds_initialization_validation.jl`).
