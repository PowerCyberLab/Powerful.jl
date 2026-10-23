#!/usr/bin/env julia
# Live documentation server with automatic rebuilding
# Usage: julia --project=docs docs/serve.jl [--port PORT] [--host HOST]

using Pkg
Pkg.activate(@__DIR__)

# Always regenerate manifest to stay in sync with Powerful's dependencies
manifest_path = joinpath(@__DIR__, "Manifest.toml")
if isfile(manifest_path)
    rm(manifest_path)
end
Pkg.develop(path=joinpath(dirname(@__DIR__), "Powerful"))
Pkg.instantiate()

using LiveServer
using Logging

# Parse command-line arguments
host = "127.0.0.1"
port = 9000

for i in eachindex(ARGS)
    if ARGS[i] == "--port" && i < length(ARGS)
        global port = parse(Int, ARGS[i+1])
    elseif ARGS[i] == "--host" && i < length(ARGS)
        global host = ARGS[i+1]
    end
end

println("=" ^ 70)
println("Powerful.jl Documentation Server")
println("=" ^ 70)
println()

# Build documentation initially
println("Building documentation...")
println("-" ^ 70)
flush(stdout)

try
    include("make.jl")
    println()
    println("-" ^ 70)
    println("✓ Initial documentation build complete")
    println()
catch e
    @error "Failed to build documentation" exception=(e, catch_backtrace())
    println()
    println("Continuing with server startup...")
    println()
end

# Start server
println("Starting live server...")
println()
println("  Documentation: http://$host:$port")
println("  Build directory: $(joinpath(@__DIR__, "build"))")
println()
println("Watching for changes in:")
println("  • docs/src/ (documentation source)")
println("  • docs/make.jl (build script)")
println("  • dev_notebooks/ (Jupyter notebooks)")
println("  • ../src/ (source code docstrings)")
println()
println("Press Ctrl+C to stop")
println("=" ^ 70)
println()
flush(stdout)

# Use servedocs which automatically rebuilds on file changes
try
    servedocs(
        foldername=".",  # We're already in docs/, look in current dir
        host=host,
        port=port,
        doc_env=true,  # Use docs environment
        skip_dirs=["build", "src/notebooks"],  # Don't watch generated files
        literate="",  # We handle notebooks via make.jl
        verbose=false,  # Reduce noise
        launch_browser=false  # Don't auto-open browser
    )
catch e
    if isa(e, InterruptException)
        println()
        println("=" ^ 70)
        println("Server stopped")
        println("=" ^ 70)
    else
        rethrow(e)
    end
end
