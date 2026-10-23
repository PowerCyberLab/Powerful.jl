# User Guides

Task-oriented guides for using Powerful.jl.

## Getting Started

New to Powerful.jl? Start here:

- **[Getting Started](getting_started.md)** - The basic workflow and key concepts
- **[Power Flow Analysis](powerflow.md)** - Formulations, algorithms, and options
- **[Accessing Results](accessing_variables.md)** - Reading solution variables and result tables
- **[Report Customization](reports.md)** - Units, columns, and presets for result tables
- **[Multi-Island Systems](multi_island.md)** - Detecting and solving disconnected networks

## Studies

- **[Time-Domain Simulation](tds.md)** - Dynamic simulation and transient analysis
- **[Eigenvalue Analysis](eigenvalue.md)** - Small-signal stability and mode analysis
- **[Frequency Response](bode.md)** - Linearization and Bode plots
- **[DC Analysis](dc_analysis.md)** - DC power flow, PTDF/LODF, and DC contingency screening
<!--GENERATED:guides-studies-->

## Additional Topics

- **[Models](models.md)** - Available static and dynamic models
- **[Configuration](configuration.md)** - System, model, and solver settings via TOML or code
- **[Automatic Differentiation](autodiff.md)** - Derivatives of power flow and simulation results
<!--GENERATED:guides-topics-->

## Learning Path

1. **Beginners**: Read [Getting Started](getting_started.md) to understand `SystemModel`, problems, and solutions
2. **Intermediate**: Study [Power Flow Analysis](powerflow.md) and [Accessing Results](accessing_variables.md)
3. **Dynamics**: Learn [Time-Domain Simulation](tds.md), then [Eigenvalue Analysis](eigenvalue.md) and [Frequency Response](bode.md)
4. **Security**: Use [Contingency Analysis](accc.md) and [DC Analysis](dc_analysis.md) for N-1/N-k assessment
5. **Reference**: Use the [API Reference](../api.md) for complete function documentation

---

**Note**: These guides focus on *how to accomplish tasks*. For complete API details, see the [API Reference](../api.md). For implementation details, see the [Developer Documentation](../developer/overview.md).
