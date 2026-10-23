# Multi-Island Systems

A network can split into several electrically disconnected parts (islands), for
example after branch outages. Powerful.jl detects islands when a power flow problem
is created and solves each island separately.

This guide uses the two-island variant of the IEEE 14-bus case, in which buses 1-6
and buses 7-14 are not connected.

## Detecting Islands

```@example multi_island
using Powerful
using PowerfulCases

raw = file(load("ieee14"), :raw; variant = "island")
sys = SystemModel(parse_network(raw))

topo = analyze_topology(sys)
topo.n_islands
```

`analyze_topology` returns a `TopologyInfo` with the connectivity of the network,
including the number of islands, the island of each bus (`topo.island_map`, by
internal bus index), and the buses of each island (`topo.components`).

`split_islands` builds one `SystemModel` per island. For a connected system it
returns the original `SystemModel`:

```@example multi_island
islands = split_islands(sys)

for (k, island) in enumerate(islands)
    println("Island $k: buses ", Int.(sort(island.models.Bus.i)),
            ", slack bus ", Int.(island.models.Slack.i))
end
```

Each island needs its own slack bus (bus type 3) to fix its voltage angle and supply
its losses. An island without one is reported with a warning and does not converge.

## Solving

Create the problem and solve as for a connected system:

```@example multi_island
prob = ACPowerFlowProblem(sys)
```

```@example multi_island
sol = solve(prob, Newton())
```

`ACPowerFlowProblem` detects the islands and `solve` solves each one. For a system
with more than one island, the result is a `MultiIslandSolution`:

- `sol.converged` is `true` only if every island converged
- `sol.iterations` and `sol.residual_norm` are the largest over the islands
- `sol.islands` holds the island systems and `sol.solutions` their
  `PowerFlowSolution`s

```@example multi_island
[(converged = s.converged, iterations = s.iterations) for s in sol.solutions]
```

## Accessing Results

`get_var`, `get_all_vars`, and `get_vars_dict` work on a `MultiIslandSolution` as on a
single solution. `get_var` finds the island that holds the key:

```@example multi_island
v1 = get_var(sol, :Bus, :v, (1,))          # Bus 1, in island 1
v14 = get_var(sol, :Bus, :v, (14,))        # Bus 14, in island 2
(v1, v14)
```

```@example multi_island
θ = get_vars_dict(sol, :Bus, :theta)       # All buses of all islands
rad2deg(θ[(1,)]), rad2deg(θ[(8,)])
```

Each island has its own angle reference: the angle of its slack bus (buses 1 and 8
here) is fixed at the value in the case file.

`get_all_vars` concatenates the values of the islands in island order:

```@example multi_island
length(get_all_vars(sol, :Bus, :v))
```

Report tables merge all islands:

```@example multi_island
using DataFrames
DataFrame(collect_model_table(sol, :Bus))
```

## Solving Islands Individually

For control over each island, for example to pass an initial guess, build one
problem per island:

```@example multi_island
solutions = [solve(ACPowerFlowProblem(island), Newton()) for island in islands]

get_var(solutions[2], :Bus, :v, (14,))
```

An initial guess (`initial_guess` keyword) is supported only for a single island,
since each island has its own variable layout.

## Disabling Island Detection

With `auto_split = false`, the whole system is solved as one problem without
splitting:

```@example multi_island
prob_whole = ACPowerFlowProblem(sys; auto_split = false)
sol_whole = solve(prob_whole, Newton())
get_var(sol_whole, :Bus, :v, (14,))
```

This requires every island to have a slack bus, as the solution above does.

## Further Reading

- [Power Flow Analysis](powerflow.md) - Formulations and algorithms
- [Accessing Results](accessing_variables.md) - Reading solutions
- [API Reference: Core Infrastructure](../api/system.md) - Topology and system functions
