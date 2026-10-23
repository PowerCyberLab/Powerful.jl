"""
PowerfulYbus - Y-Bus Network Formulation for Transient Stability

This package provides a sparse admittance-matrix (Y-bus) network formulation for
transient stability simulation: the network is the sparse admittance matrix of each
pre-computed topology, with the generators represented as Norton-equivalent current
sources behind their internal voltage.

# Key Features
- Full sparse Y matrix per network topology, pre-computed so a fault or line trip is a
  cheap switch to an already-factored matrix rather than a new factorization
- Generator admittances `1/(ra + j*xd')` and constant-impedance loads in Y;
  constant-power and constant-current loads and current sources as injections
- Network solve at every bus by sparse Newton iteration; one linear solve when all
  loads are constant impedance
- Solvers: a fixed-step partitioned trapezoidal integrator, and OrdinaryDiffEq.jl with
  the network equations in the DAE (`SimultaneousNetwork`, sparse Jacobian) or solved
  in every evaluation (`EliminatedNetwork`, dense Jacobian)

# Usage
```julia
using PowerfulYbus

# Create Y-bus TDS problem from standard TDS problem
prob_ybus = YBusTDSProblem(tds_prob; events=events)
sol = solve(prob_ybus)                    # DiffEq.jl, SimultaneousNetwork()
sol = solve_trapezoidal(prob_ybus; dt=0.001)
```

# See Also
- PowerfulTDS for standard residual-based TDS
"""
module PowerfulYbus

using DocStringExtensions
using StructArrays
using SparseArrays
using LinearAlgebra
using KLU: KLU, klu, klu!

# PowerfulTDS, PowerfulPF, PowerfulCore
using PowerfulTDS
using PowerfulTDS: fault_admittance, pf_snapshot_positions, working_copy, init_tds_buffers!
# The automaton runtime of the phasor simulation, run on the Y-bus DAE (`dae_residuals!`)
using PowerfulTDS: AutomatonRuntime, automaton_runtime, start!, build_automaton_callback, settle_instant!,
    rows_without_their_variable,
    erased, HostReinit, FixedStepAutomata, step_automata!, settle_fixed!, retrace_automata!, _event_models,
    FixedInstant, StartResolve
import PowerfulTDS: dae_residuals!, dual_workspace, algebraic_equations_changed!, rate_dual_type, _start_state_rates,
    _start_algebraic_rate
using PowerfulPF
using PowerfulCore

# Import core types needed for trait implementations
import PowerfulCore: AbstractModel, BUS_PV, BUS_SLACK
using PowerfulCore: has_automata, AutomatonState
import PowerfulCore: automaton_states, sync_automaton_memory!, automaton_memory_synced  # extended for YBusSystemBuffers

# ============================================================================
# Formulation Traits
# ============================================================================

# ============================================================================
# Solution Strategies (Two-Node vs Terminal-Node)
# ============================================================================

# ============================================================================
# Traits - Network Partition Classification
# ============================================================================
include("traits/network_partition_traits.jl")
# exports: NetworkPartitionRole, FixedNetworkElement, MachineWithInternalVoltage,
#          NonConformingLoad, StaticGenerator, NotInYBusFormulation, UnsupportedInYBus, network_partition_role

include("traits/nc_load_traits.jl")
# exports: NCLoadType, VoltageDependentNCLoad, CurrentSourceNCLoad, nc_load_type, declares_nc_load_type

include("traits/ybus_interface.jl")
# exports: Interface functions for Y matrix contributions

# ============================================================================
# Y-Bus Construction (needed by network solver)
# ============================================================================
include("ybus/ybus_builder.jl")
# exports: FullYMatrix, build_full_ybus


include("ybus/topology_cache.jl")
# exports: TopologyKey, TopologyCache, switch_topology!, current_Y, add_topology!

include("network/bus_voltage_inputs.jl")
# exports: BusVoltageInput, bus_voltage_inputs, set_bus_voltages!

# ============================================================================
# Network Solver (needs Y-bus, NC load workspace)
# ============================================================================
include("network/internal_voltage.jl")
# exports: InternalVoltageBuffer, gather_internal_voltages!, get_internal_voltage

include("network/current_source.jl")
# exports: CurrentSourceBuffer, setup_current_source_buffer!, gather_source_currents!

include("network/network_injections.jl")
# exports: injects_nonlinear_current, nc_load_current, nc_load_current_jacobian,
#          source_current, source_current_jacobian

include("network/static_generators.jl")
# exports: StaticGeneratorBuses, static_generator_buses


include("network/network_solver.jl")
# exports: NetworkCache, NetworkSolution, RetainedNetwork, solve_network!,
#          similar_network_cache, retained_network!, network_failure_message


# ============================================================================
# Y-Bus Residuals Trait (needs NetworkSolution from network_solver.jl)
# ============================================================================
include("traits/ybus_residuals.jl")
# exports: calc_ybus_residuals!

# ============================================================================
# TDS Integration
# ============================================================================
include("integration/ybus_system_buffers.jl")
# exports: YBusSystemBuffers, evaluate_ybus_residuals!

# Event System (pre-computed topologies)
include("integration/algebraic_newton.jl")
include("integration/ybus_event_schedule.jl")
# exports: YBusEventSchedule, TimeSegment, get_topology_at, has_topology_change, build_time_segments

include("integration/ybus_event_processor.jl")
# exports: process_ybus_events

include("integration/ybus_tds_problem.jl")
# exports: YBusTDSProblem, get_network_cache, get_topology_cache


include("integration/ybus_residual.jl")
# exports: NetworkTreatment, EliminatedNetwork, SimultaneousNetwork, YBusResidual,
#          consistent_algebraics!

include("integration/ybus_solve.jl")
# exports: solve (for YBusTDSProblem)

include("integration/trapezoidal_solver.jl")
# exports: solve_trapezoidal (returns TDSSolution from PowerfulCore)

end # module
