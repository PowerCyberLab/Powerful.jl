#!/usr/bin/env julia
#
# Regenerate the Jupyter notebooks `notebooks/NN_*.ipynb` from the Literate.jl sources
# `notebooks/NN_*.jl`. The notebooks are written without outputs; run them in Jupyter
# to see the results. Run this after editing a source, and commit both files.
#
# Usage (from the repository root):
#   julia --project=docs notebooks/convert_notebooks.jl
#
# The documentation pages are generated separately, on every docs build, by
# `docs/make.jl` (see `build.jl`).

include(joinpath(@__DIR__, "build.jl"))

generate_ipynb()
