# TDS Initialization Validation

Learn how to validate time-domain simulation (TDS) initialization and diagnose steady-state issues.

## Setup


```julia
using Powerful
using Powerful.PowerFlowData
using Powerful.Models: GENCLS
using StructArrays
using DataFrames
```

## Why Validation Matters

Time-domain simulations start from a steady-state operating point. If the initial
conditions don't satisfy the DAE (differential-algebraic equation) residuals f=0 and g=0,
the simulation will have undesirable initial transients before applying disturbances.

**Common issues:**
- Incorrect generator initialization
- Missing dynamic models for critical buses
- Inconsistent power flow and dynamic model data

**Solution:** Use validation tools to catch these errors before time-stepping begins.

## 1. Load System and Solve Power Flow


```julia
# Load IEEE 14-bus test case
case = PowerFlowData.parse_network(joinpath(pkgdir(Powerful), "cases", "ieee14.raw"))
sys_pf = SystemModel(case)

# Solve power flow (establishes steady-state)
prob_pf = ACPowerFlowProblem(sys_pf)
sol_pf = solve(prob_pf, NewtonRaphson())

println("Power flow converged: ", sol_pf.converged)
println("Iterations: ", sol_pf.iterations)
```

    Power flow converged: true
    Iterations: 4


## 2. Identify Generator Buses

Dynamic models replace static generator models (Slack, PV) at selected buses.


```julia
# Get generator bus numbers
slack_models = sys_pf.models[:Slack]
pv_models = sys_pf.models[:PV]

slack_bus = slack_models.i[1]  ## Slack bus number
pv_buses = pv_models.i         ## PV bus numbers

println("Slack bus: $slack_bus")
println("PV buses: $pv_buses")
println("Total generators: ", 1 + length(pv_buses))
```

    Slack bus: 1
    PV buses: Int32[2, 3, 6, 8]
    Total generators: 5


## 3. Create Dynamic Models

Replace all generators with classical machine models (GENCLS).
This is the recommended approach for TDS.


```julia
# Combine all generator buses
gen_buses = vcat(slack_bus, pv_buses)

# Create GENCLS models for each generator
gencls_array = [
    GENCLS(
        BUS = gen_buses[i],
        ID = "1",
        H = 5.0,      ## Inertia constant (seconds)
        D = 0.0,      ## Damping coefficient
        xd_p = 0.3,   ## Transient reactance (pu)
        status = true
    )
    for i in eachindex(gen_buses)
]

# Convert to StructArray (required format)
gencls_models = StructArray(gencls_array)
dynamic_models = (GENCLS = gencls_models,)

println("Created $(length(gencls_models)) GENCLS models")
```

    Created 5 GENCLS models


## 4. Create TDS Problem

The TDSProblem initializes dynamic states based on the power flow solution.


```julia
prob = TDSProblem(sys_pf, sol_pf, dynamic_models; tspan=(0.0, 10.0))

println("TDS problem created for $(length(gencls_models)) dynamic generators")
```

    TDS problem created for 5 dynamic generators


## 5. Validate Initialization

Check that all residuals are zero at t=0 (steady-state).


```julia
result = validate_tds_initialization(prob; atol=1e-6)

println("\nValidation result:")
println("  Valid: ", result.is_valid)
println("  Max differential residual: ", result.max_diff_residual)
println("  Max algebraic residual: ", result.max_alg_residual)
```

    
    Validation result:
      Valid: true
      Max differential residual: 5.551115123125783e-17
      Max algebraic residual: 3.2540081740251026e-12


## 6. Inspect Validation Report

Print detailed diagnostics if validation fails.


```julia
print_validation_report(result)
```

    
    ======================================================================
    TDS Initialization Validation Report
    ======================================================================
    ✓ All residuals within tolerance
      Max differential residual: 5.55e-17
      Max algebraic residual:    3.25e-12


## 7. Understanding Residuals

**Differential residuals (f):** Time derivatives of state variables (e.g., δ̇, ω̇)
- Should be zero at steady-state (no acceleration)
- Non-zero indicates initial transient

**Algebraic residuals (g):** Power balance equations
- Should be zero at all times
- Non-zero indicates power flow mismatch


```julia
if result.is_valid
    println("\n✓ Initialization passed - ready for time-domain simulation")
else
    println("\n✗ Initialization failed - check report above")
end
```

    
    ✓ Initialization passed - ready for time-domain simulation


## 8. What If Validation Fails?

**Common causes:**
1. **Missing dynamic models**: Some generators still use static (Slack/PV) models
2. **Incorrect parameters**: Generator reactance or inertia values mismatch
3. **Power flow mismatch**: Dynamic model injections ≠ static model injections

**Example: Intentional failure** (replace only one generator)


```julia
# Create problem with only 1 GENCLS (others remain static)
gen_buses_partial = [slack_bus]  ## Only replace slack
gencls_partial = StructArray([
    GENCLS(BUS = slack_bus, ID = "1", H = 5.0, D = 0.0, xd_p = 0.3, status = true)
])
dynamic_models_partial = (GENCLS = gencls_partial,)

prob_partial = TDSProblem(sys_pf, sol_pf, dynamic_models_partial; tspan=(0.0, 10.0))
result_partial = validate_tds_initialization(prob_partial; atol=1e-6)

println("\nPartial replacement validation:")
println("  Valid: ", result_partial.is_valid)
println("  (Expected: true - network models handle remaining generators)")
```

    
    Partial replacement validation:
      Valid: true
      (Expected: true - network models handle remaining generators)


## 9. Save Validation Report

Export diagnostics to file for documentation or debugging.


```julia
write_validation_report(result, "validation_report.txt")
println("\n✓ Validation report saved to: validation_report.txt")
```

    [ Info: Validation report written to validation_report.txt
    
    ✓ Validation report saved to: validation_report.txt


## 10. Inspect Initial States

View the computed initial conditions for dynamic models in a formatted table.


```julia
# Print all instances and all state variables
print_initialization_summary(prob)

# Filter by instance keys (show only specific generators)
println("\n--- Filtered: First two generators ---")
print_initialization_summary(prob; keys = Dict(:GENCLS => [(1, "1"), (2, "1")]))

# Filter by state variables (show only delta and omega)
println("\n--- Filtered: Only delta and omega ---")
print_initialization_summary(prob; vars = [:delta, :omega])

# Include output variables (P_e, Q_e) - requires residual evaluation
println("\n--- With output variables (P_e, Q_e) ---")
print_initialization_summary(prob; vars = [:delta, :omega, :P_e, :Q_e], include_output_vars = true)
```

    
    ======================================================================
    TDS Initialization Summary
    ======================================================================
    
    GENCLS: 5 instances
    ----------------------------------------------------------------------
    ┌─────┬────┬───────────┬──────────┬───────────┬──────────┬──────────┬───────────┬──────────┬──────────┐
    │ BUS │ ID │     delta │    omega │       I_d │      I_q │    psi_d │     psi_q │      P_m │      E_q │
    ├─────┼────┼───────────┼──────────┼───────────┼──────────┼──────────┼───────────┼──────────┼──────────┤
    │   1 │  1 │  0.240265 │ 1.000000 │ -0.001447 │ 0.816995 │ 1.000413 │ -0.245099 │ 0.816978 │ 0.999979 │
    │   2 │  1 │  0.071603 │ 1.000000 │  0.380370 │ 0.351276 │ 1.024595 │ -0.105383 │ 0.400000 │ 1.138706 │
    │   3 │  1 │  0.050109 │ 1.000000 │  0.196981 │ 0.376378 │ 1.003669 │ -0.112913 │ 0.400000 │ 1.062763 │
    │   6 │  1 │ -0.041486 │ 1.000000 │  0.400352 │ 0.261526 │ 1.027007 │ -0.078458 │ 0.300000 │ 1.147113 │
    │   8 │  1 │  0.072279 │ 1.000000 │  0.195314 │ 0.322869 │ 1.025436 │ -0.096861 │ 0.350000 │ 1.084030 │
    └─────┴────┴───────────┴──────────┴───────────┴──────────┴──────────┴───────────┴──────────┴──────────┘
    
    ======================================================================
    
    --- Filtered: First two generators ---
    
    ======================================================================
    TDS Initialization Summary
    ======================================================================
    
    GENCLS: 2 instances
    ----------------------------------------------------------------------
    ┌─────┬────┬──────────┬──────────┬───────────┬──────────┬──────────┬───────────┬──────────┬──────────┐
    │ BUS │ ID │    delta │    omega │       I_d │      I_q │    psi_d │     psi_q │      P_m │      E_q │
    ├─────┼────┼──────────┼──────────┼───────────┼──────────┼──────────┼───────────┼──────────┼──────────┤
    │   1 │  1 │ 0.240265 │ 1.000000 │ -0.001447 │ 0.816995 │ 1.000413 │ -0.245099 │ 0.816978 │ 0.999979 │
    │   2 │  1 │ 0.071603 │ 1.000000 │  0.380370 │ 0.351276 │ 1.024595 │ -0.105383 │ 0.400000 │ 1.138706 │
    └─────┴────┴──────────┴──────────┴───────────┴──────────┴──────────┴───────────┴──────────┴──────────┘
    
    ======================================================================
    
    --- Filtered: Only delta and omega ---
    
    ======================================================================
    TDS Initialization Summary
    ======================================================================
    
    GENCLS: 5 instances
    ----------------------------------------------------------------------
    ┌─────┬────┬───────────┬──────────┐
    │ BUS │ ID │     delta │    omega │
    ├─────┼────┼───────────┼──────────┤
    │   1 │  1 │  0.240265 │ 1.000000 │
    │   2 │  1 │  0.071603 │ 1.000000 │
    │   3 │  1 │  0.050109 │ 1.000000 │
    │   6 │  1 │ -0.041486 │ 1.000000 │
    │   8 │  1 │  0.072279 │ 1.000000 │
    └─────┴────┴───────────┴──────────┘
    
    ======================================================================
    
    --- With output variables (P_e, Q_e) ---
    
    ======================================================================
    TDS Initialization Summary
    ======================================================================
    
    GENCLS: 5 instances
    ----------------------------------------------------------------------
    ┌─────┬────┬───────────┬──────────┬──────────┬───────────┐
    │ BUS │ ID │     delta │    omega │      P_e │       Q_e │
    ├─────┼────┼───────────┼──────────┼──────────┼───────────┤
    │   1 │  1 │  0.240265 │ 1.000000 │ 0.816978 │ -0.201692 │
    │   2 │  1 │  0.071603 │ 1.000000 │ 0.400000 │  0.352707 │
    │   3 │  1 │  0.050109 │ 1.000000 │ 0.400000 │  0.155205 │
    │   6 │  1 │ -0.041486 │ 1.000000 │ 0.300000 │  0.390645 │
    │   8 │  1 │  0.072279 │ 1.000000 │ 0.350000 │  0.169009 │
    └─────┴────┴───────────┴──────────┴──────────┴───────────┘
    
    ======================================================================


## 11. Verify System Composition

After dynamic model replacement, verify which models remain in the system.


```julia
# Create ADCache to access model workspaces
adcache = ADCache(prob.address_manager, prob.sys, prob.formulation)

println("\nTDS model composition:")
for (model_name, workspace) in pairs(adcache.float_ws.model_workspaces)
    n = length(prob.sys.models[model_name])
    if n > 0
        println("  $model_name: $n instances")
    end
end
```

    
    TDS model composition:
      Branch: 16 instances
      Transformer2W: 4 instances
      Bus: 14 instances
      PQ: 11 instances
      GENCLS: 5 instances


## 12. Compare Generator Counts

Show how dynamic models replace static generators.


```julia
println("\nTDS system generator composition:")
println("  Slack: $(length(prob.sys.models[:Slack])) (replaced)")
println("  PV: $(length(prob.sys.models[:PV])) (replaced)")
println("  GENCLS: $(length(prob.sys.models[:GENCLS])) (dynamic)")
```

    
    TDS system generator composition:
      Slack: 0 (replaced)
      PV: 0 (replaced)
      GENCLS: 5 (dynamic)


## Summary

**Key takeaways:**
1. Always validate TDS initialization before running simulations
2. Use `validate_tds_initialization(prob)` to check residuals
3. Aim for `max_residual < 1e-6` for accurate simulations
4. Replace all generators with dynamic models for consistent results
5. Export reports for documentation and debugging

**Next steps:**
- Run time-domain simulations (see TDS guide)
- Apply disturbances (line trips, load changes)
- Analyze transient stability

---

*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
