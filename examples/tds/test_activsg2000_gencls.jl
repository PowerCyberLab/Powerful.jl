#!/usr/bin/env julia
"""Test TDS initialization with all GENCLS for ACTIVSg2000."""

using Powerful
using PowerfulCases
using PowerfulCore: create_adcache, ForwardDiffBackend
using PowerfulTDS: TDSAdapter, ODEInterface, evaluate!
using PowerfulCore: AlgebRes

println("=== Testing ACTIVSg2000 with all GENCLS ===")

# Load case using collection/case format to avoid ambiguity
case = load("synthetic/ACTIVSg2000")

# Parse and solve power flow
println("1. Parsing RAW and solving power flow...")
sys_pf = SystemModel(parse_network(case.raw))
pf_sol = solve(ACPowerFlowProblem(sys_pf), Newton())
println("   Power flow converged: $(pf_sol.converged)")

# Parse DYR
println("2. Parsing DYR file...")
dyr = parse_dynamics(case.dyr, sys_pf)
println("   GENCLS: $(haskey(sys_pf.models, :GENCLS) ? length(sys_pf.models.GENCLS) : 0)")
println("   SEXS: $(haskey(sys_pf.models, :SEXS) ? length(sys_pf.models.SEXS) : 0)")
println("   TGOV1: $(haskey(sys_pf.models, :TGOV1) ? length(sys_pf.models.TGOV1) : 0)")

# Create TDS problem
println("3. Creating TDS problem...")
init = setup_tds(pf_sol, dyr)
prob = TDSProblem(init; tspan=(0.0, 0.1))

# Evaluate initial residuals
println("4. Evaluating initial residuals...")
adcache = create_adcache(prob.address_manager, prob.sys, prob.formulation, ForwardDiffBackend())
adapter = TDSAdapter(adcache, prob)
u0 = vcat(adapter.x_buf, adapter.y_buf)
du = similar(u0)
evaluate!(du, adapter, u0, 0.0, ODEInterface())

n_diff = adapter.n_diff
dx = du[1:n_diff]
g = du[n_diff+1:end]

println("\n=== Results ===")
println("Max |dx|: $(maximum(abs, dx))")
println("Max |g|: $(maximum(abs, g))")

# Check problematic residuals
am = prob.address_manager
buses = prob.sys.models.Bus
p_addrs = am.addresses[(AlgebRes(), :Bus, :p)]
q_addrs = am.addresses[(AlgebRes(), :Bus, :q)]

tol = 1e-3
n_p_issues = sum(abs(g[addr]) > tol for addr in p_addrs)
n_q_issues = sum(abs(g[addr]) > tol for addr in q_addrs)

println("Bus:p issues (>$tol): $n_p_issues")
println("Bus:q issues (>$tol): $n_q_issues")

if n_p_issues + n_q_issues == 0
    println("\n✓ SUCCESS: TDS initialization at steady state!")
else
    println("\n✗ Some issues remain - showing first 15:")
    let shown = 0
        for (i, (p_addr, q_addr)) in enumerate(zip(p_addrs, q_addrs))
            if abs(g[p_addr]) > tol || abs(g[q_addr]) > tol
                shown += 1
                shown > 15 && break
                println("  Bus $(buses.i[i]) (type=$(buses.ide[i])): p=$(round(g[p_addr], digits=4)), q=$(round(g[q_addr], digits=4))")
            end
        end
    end
end
