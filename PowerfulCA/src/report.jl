# ACCC Report Functions

"""
    print_accc_report(result; vmin=0.95, vmax=1.05, rating=:rate_a)

Print ACCC report showing voltage and branch flow violations.

Instances a contingency switched off (`out_of_service` of its result) are not
checked: an outaged branch carries no flow, and the buses of an island the
contingency de-energized are listed as isolated instead of as voltage violations.

# Arguments
- `result`: Output from `solve_accc`
- `vmin`: Low voltage limit in pu (default 0.95)
- `vmax`: High voltage limit in pu (default 1.05)
- `rating`: Rating set to use for branch limits (:rate_a, :rate_b, :rate_c)

# Example
```julia
result = solve_accc(sys, contingencies)
print_accc_report(result; vmin=0.9, vmax=1.1)
```
"""
function print_accc_report(result; vmin::Float64=0.95, vmax::Float64=1.05, rating::Symbol=:rate_a)
    base_sol = result.base
    sys = base_sol.prob.sys
    am = base_sol.prob.address_manager
    sbase = sys.parsed_data.system_base

    # Get addresses for voltage extraction
    v_addrs = get_addresses(am, :Bus, :v)
    buses = sys.models.Bus
    branches = sys.models.Branch
    n_bus = length(buses)
    n_branches = length(branches)

    # Base case voltages
    v_base = base_sol.dae.y[v_addrs]

    # Get branch ratings (MVA)
    branch_ratings = getproperty(branches, rating)

    # Compute base case branch flows
    p_from_base = compute_observations(base_sol, Branch, :p_from) .* sbase
    q_from_base = compute_observations(base_sol, Branch, :q_from) .* sbase
    s_from_base = sqrt.(p_from_base.^2 .+ q_from_base.^2)

    # Header
    println("=" ^ 80)
    println("AC CONTINGENCY ANALYSIS REPORT")
    println("=" ^ 80)
    println("Rating set: $(rating)    Voltage limits: $(vmin) - $(vmax) pu")
    println("System base: $(sbase) MVA")
    println()

    # Base case summary
    println("-" ^ 80)
    println("BASE CASE: $(base_sol.converged ? "CONVERGED" : "FAILED") in $(base_sol.iterations) iterations")
    println("-" ^ 80)
    println()

    # Collect all violations
    voltage_violations = NamedTuple[]
    branch_violations = NamedTuple[]
    isolated_buses = NamedTuple[]
    nonconverged = NamedTuple[]

    for r in result.contingencies
        if !r.converged
            push!(nonconverged, (name=r.name, iterations=r.iterations, residual=r.residual_norm,
                                 singular=r.singular))
            continue
        end

        dead_buses = get(r.out_of_service, :Bus, Int[])
        out_branches = get(r.out_of_service, :Branch, Int[])
        for i in dead_buses
            push!(isolated_buses, (contingency = r.name, bus = Int(buses.i[i]), basekv = buses.basekv[i]))
        end

        # Extract voltages from stored y vector
        v_cont = r.y[v_addrs]

        # Check voltage violations of the energized buses
        for i in 1:n_bus
            i in dead_buses && continue
            v = v_cont[i]
            if v < vmin || v > vmax
                push!(voltage_violations, (
                    contingency = r.name,
                    bus = Int(buses.i[i]),
                    basekv = buses.basekv[i],
                    v_cont = v,
                    v_base = v_base[i],
                    v_min = vmin,
                    v_max = vmax,
                    violation_type = v < vmin ? "LOW" : "HIGH",
                ))
            end
        end

        # Compute branch flows for this contingency using reusable physics kernel
        theta_i_addrs = get_addresses(am, :Branch, :theta_i)
        theta_j_addrs = get_addresses(am, :Branch, :theta_j)
        vi_addrs = get_addresses(am, :Branch, :v_i)
        vj_addrs = get_addresses(am, :Branch, :v_j)

        theta_i = r.y[theta_i_addrs]
        theta_j = r.y[theta_j_addrs]
        v_i = r.y[vi_addrs]
        v_j = r.y[vj_addrs]

        # Check branch overloads of the in-service branches
        for k in 1:n_branches
            k in out_branches && continue
            rate = branch_ratings[k]
            if rate <= 0
                continue  # Skip branches with no rating
            end

            # Use proper physics kernel for branch flow calculation
            flows = calc_branch_power_flows(
                v_i[k], theta_i[k], v_j[k], theta_j[k],
                branches.r[k], branches.x[k], branches.b[k],
                branches.gi[k], branches.bi[k], branches.gj[k], branches.bj[k]
            )
            s_i = sqrt(flows.p_i^2 + flows.q_i^2) * sbase  # MVA

            pct_loading = 100.0 * s_i / rate

            if pct_loading > 100.0
                push!(branch_violations, (
                    contingency = r.name,
                    from_bus = Int(branches.i[k]),
                    to_bus = Int(branches.j[k]),
                    circuit = String(branches.ckt[k]),
                    flow_mva = s_i,
                    rating_mva = rate,
                    pct_loading = pct_loading,
                    base_flow_mva = s_from_base[k],
                ))
            end
        end
    end

    # Print voltage violations
    if !isempty(voltage_violations)
        println("-" ^ 80)
        println("VOLTAGE VIOLATIONS")
        println("-" ^ 80)
        println(rpad("CONTINGENCY", 16), rpad("BUS", 8), rpad("BASE_KV", 10),
                rpad("V-CONT", 10), rpad("V-BASE", 10), rpad("V-MIN", 8), rpad("V-MAX", 8), "TYPE")
        println("-" ^ 80)
        for v in voltage_violations
            println(rpad(v.contingency, 16), rpad(v.bus, 8), rpad(round(v.basekv, digits=2), 10),
                    rpad(round(v.v_cont, digits=5), 10), rpad(round(v.v_base, digits=5), 10),
                    rpad(v.v_min, 8), rpad(v.v_max, 8), v.violation_type)
        end
        println()
    end

    # Print branch violations
    if !isempty(branch_violations)
        println("-" ^ 80)
        println("BRANCH OVERLOADS")
        println("-" ^ 80)
        println(rpad("CONTINGENCY", 16), rpad("FROM", 8), rpad("TO", 8), rpad("CKT", 5),
                rpad("FLOW", 10), rpad("RATING", 10), rpad("%LOAD", 8), "BASE_FLOW")
        println("-" ^ 80)
        for b in branch_violations
            println(rpad(b.contingency, 16), rpad(b.from_bus, 8), rpad(b.to_bus, 8), rpad(b.circuit, 5),
                    rpad(round(b.flow_mva, digits=1), 10), rpad(round(b.rating_mva, digits=1), 10),
                    rpad(round(b.pct_loading, digits=1), 8), round(b.base_flow_mva, digits=1))
        end
        println()
    end

    # Print buses isolated by contingencies
    if !isempty(isolated_buses)
        println("-" ^ 80)
        println("ISOLATED BUSES (de-energized: no path to a slack bus)")
        println("-" ^ 80)
        println(rpad("CONTINGENCY", 16), rpad("BUS", 8), "BASE_KV")
        println("-" ^ 80)
        for b in isolated_buses
            println(rpad(b.contingency, 16), rpad(b.bus, 8), round(b.basekv, digits=2))
        end
        println()
    end

    # Print non-converged contingencies
    if !isempty(nonconverged)
        println("-" ^ 80)
        println("NON-CONVERGED CONTINGENCIES")
        println("-" ^ 80)
        println(rpad("CONTINGENCY", 20), rpad("ITERATIONS", 12), "RESIDUAL")
        println("-" ^ 80)
        for nc in nonconverged
            println(rpad(nc.name, 20), rpad(nc.iterations, 12),
                    nc.singular ? "singular Jacobian" : @sprintf("%.2e", nc.residual))
        end
        println()
    end

    # Summary
    println("=" ^ 80)
    println("SUMMARY")
    println("=" ^ 80)
    n_cont = length(result.contingencies)
    n_conv = count(r -> r.converged, result.contingencies)
    println("Total contingencies:     $(n_cont)")
    println("Converged:               $(n_conv)")
    println("Non-converged:           $(n_cont - n_conv)")
    println("Voltage violations:      $(length(voltage_violations))")
    println("Branch overloads:        $(length(branch_violations))")
    println("Isolated buses:          $(length(isolated_buses))")
    println("=" ^ 80)

    return (
        voltage_violations = voltage_violations,
        branch_violations = branch_violations,
        isolated_buses = isolated_buses,
        nonconverged = nonconverged,
    )
end
