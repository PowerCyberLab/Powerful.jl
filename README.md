# Powerful





<!-- [![Build Status](https://github.com/sensl/Powerful.jl/badges/main/pipeline.svg)](https://github.com/sensl/Powerful.jl/pipelines) -->


Yet another power systems simulation package meant to be powerful.


## Introduction

### Installation

Install `Powerful.jl` with

```julia
using Pkg; Pkg.add("Powerful")
```

## Documentation

### Building Documentation with Notebooks

The documentation includes converted Jupyter notebooks from the `dev_notebooks` directory.
To build the documentation:

```julia
cd docs
julia --project -e 'using Pkg; Pkg.instantiate(); include("make.jl")'
```

The notebooks are automatically converted to markdown and included in the documentation under the "Notebooks" section.
