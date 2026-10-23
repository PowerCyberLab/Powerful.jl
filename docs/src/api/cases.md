# Test Cases (PowerfulCases)

Access to the bundled and downloadable test cases. Case files are always loaded
through PowerfulCases: `load("ieee14")` returns a `CaseBundle` whose `raw` and `dyr`
fields are file paths.

## Loading Cases

```@autodocs
Modules = [PowerfulCases]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/PowerfulCases.jl"]
Private = false
```

## Case Manifests

```@autodocs
Modules = [PowerfulCases]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/manifest.jl"]
Private = false
```

## Remote Cases and Cache

```@autodocs
Modules = [PowerfulCases]
Order = [:module, :constant, :type, :function, :macro]
Pages = ["src/registry.jl", "src/cache.jl"]
Private = false
```
