# Precompilation workload for PowerfulTDS
#
# Compiles TDS-specific components (DYR parsing, model kernels, functors)
# and TDS infrastructure (ADCache, TDSAdapter, solve setup).
#
# Uses a small IEEE 14-bus test case with GENROU, SEXS, TGOV1 models.
# Skips DiffEq.solve() since Julia can't cache cross-package specializations.
#
# To disable: Comment out the `include("precompile.jl")` line in PowerfulTDS.jl

using PrecompileTools

# ============================================================================
# Test Data: IEEE 14-bus system with GENROU, SEXS, TGOV1 models
# ============================================================================

const _PRECOMPILE_RAW = """
0,   100.00, 33, 0, 1, 60.00     / PSS(R)E 33 RAW created by rawd33
 IEEE 14 BUS TEST CASE
 03/06/14 CONTO           100.0  1962 W
   1,'BUS1        ',  69.0000,3,   1,   1,   1,1.03000,   0.0000,1.10000,0.90000,1.10000,0.90000
   2,'BUS2        ',  69.0000,2,   1,   1,   1,1.01970,  -1.6032,1.10000,0.90000,1.10000,0.90000
   3,'BUS3        ',  69.0000,1,   1,   1,   1,1.00042,  -3.4433,1.10000,0.90000,1.10000,0.90000
   4,'BUS4        ',  69.0000,1,   1,   1,   1,0.99858,  -4.2812,1.10000,0.90000,1.10000,0.90000
   5,'BUS5        ',  69.0000,1,   1,   1,   1,1.00443,  -3.6850,1.10000,0.90000,1.10000,0.90000
0 / END OF BUS DATA, BEGIN LOAD DATA
   3,'1 ',1,   1,   1,    50.000,    25.000,     0.000,     0.000,     0.000,     0.000,   1,1,0
   4,'1 ',1,   1,   1,    47.800,    10.000,     0.000,     0.000,     0.000,     0.000,   1,1,0
   5,'1 ',1,   1,   1,     7.600,     1.600,     0.000,     0.000,     0.000,     0.000,   1,1,0
0 / END OF LOAD DATA, BEGIN FIXED SHUNT DATA
0 / END OF FIXED SHUNT DATA, BEGIN GENERATOR DATA
   1,'1 ',    81.442,     1.962,   100.000,   -50.000,1.03000,     0,   100.000, 0.00000E+0, 2.30000E-1, 0.00000E+0, 0.00000E+0,1.00000,1,  100.0,   200.000,    50.000,   1,1.0000
   2,'1 ',    40.000,    15.000,    15.000,   -40.000,1.03000,     0,   100.000, 0.00000E+0, 1.30000E-1, 0.00000E+0, 0.00000E+0,1.00000,1,  100.0,    50.000,    10.000,   1,1.0000
0 / END OF GENERATOR DATA, BEGIN BRANCH DATA
   1,     2,'1 ', 1.93800E-2, 5.91700E-2,   0.05280,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
   1,     5,'1 ', 5.40300E-2, 2.23040E-1,   0.04920,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
   2,     3,'1 ', 4.69900E-2, 1.97970E-1,   0.04380,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
   2,     4,'1 ', 5.81100E-2, 1.76320E-1,   0.03400,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
   2,     5,'1 ', 5.69500E-2, 1.73880E-1,   0.03460,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
   3,     4,'1 ', 6.70100E-2, 1.71030E-1,   0.01280,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
   4,     5,'1 ', 1.33500E-2, 4.21100E-2,   0.00000,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
0 / END OF BRANCH DATA, BEGIN TRANSFORMER DATA
0 / END OF TRANSFORMER DATA, BEGIN AREA DATA
 1,     1,     0.000,   999.990,'AREA1       '
0 / END OF AREA DATA, BEGIN TWO-TERMINAL DC DATA
0 / END OF TWO-TERMINAL DC DATA, BEGIN VSC DC LINE DATA
0 / END OF VSC DC LINE DATA, BEGIN IMPEDANCE CORRECTION DATA
0 / END OF IMPEDANCE CORRECTION DATA, BEGIN MULTI-TERMINAL DC DATA
0 / END OF MULTI-TERMINAL DC DATA, BEGIN MULTI-SECTION LINE DATA
0 / END OF MULTI-SECTION LINE DATA, BEGIN ZONE DATA
 1,'ZONE1       '
0 / END OF ZONE DATA, BEGIN INTER-AREA TRANSFER DATA
0 / END OF INTER-AREA TRANSFER DATA, BEGIN OWNER DATA
 1,'OWNER1      '
0 / END OF OWNER DATA, BEGIN FACTS DEVICE DATA
0 / END OF FACTS DEVICE DATA, BEGIN SWITCHED SHUNT DATA
0 / END OF SWITCHED SHUNT DATA, BEGIN GNE DATA
0 / END OF GNE DATA, BEGIN INDUCTION MACHINE DATA
0 / END OF INDUCTION MACHINE DATA
Q
"""

const _PRECOMPILE_DYR = """
// IEEE 5-bus with GENROU generators, SEXS exciters, and TGOV1 governors

// GENROU: BUS 'GENROU' ID Td10 Td20 Tq10 Tq20
//         H D Xd Xq Xd1
//         Xq1 Xd2 Xl S10 S12 /
1 'GENROU' 1  6.5  0.06  0.2  0.05
            4.0  0.0  1.8  1.75  0.6
            0.8  0.23  0.15  0.09  0.38 /
2 'GENROU' 1  6.5  0.06  0.2  0.05
            6.5  0.0  1.8  1.75  0.6
            0.8  0.28  0.15  0.09  0.38 /

// SEXS: BUS 'SEXS' ID TATB TB K TE EMIN EMAX /
1 'SEXS' 1  0.1  10.0  100.0  0.1  0.0  3.0 /
2 'SEXS' 1  0.1  10.0  100.0  0.1  0.0  3.0 /

// TGOV1: BUS 'TGOV1' ID R T1 Vmax Vmin T2 T3 Dt /
1 'TGOV1' 1  0.05  0.1  1.0  0.0  0.5  2.1  0.1 /
2 'TGOV1' 1  0.05  0.1  1.0  0.0  0.5  2.1  0.1 /
"""

# ============================================================================
# Workload Helpers
# ============================================================================

"""
Create system with power flow solution for TDS precompilation.
"""
function _precompile_system(raw::String)
    # Parse RAW data and solve power flow
    case = parse_network(IOBuffer(raw))
    sys = SystemModel(case)
    prob = ACPowerFlowProblem(sys)
    sol = solve(prob, Newton())
    return sys, sol
end

"""Precompile DYR parsing of `dyr` text against `sys`."""
function _precompile_parse_dynamics(sys, dyr::String)
    dyr_path = tempname() * ".dyr"
    try
        write(dyr_path, dyr)
        return parse_dynamics(dyr_path, sys)
    finally
        rm(dyr_path; force=true)
    end
end

"""Precompile TDS model metadata access."""
function _precompile_model_metadata()
    for M in (GENCLS, GENROU, GENTPJ1, GENSAL, TGOV1, SEXS)
        meta = model_metadata(M)
        _ = meta.vars
        _ = meta.residuals
        _ = meta.initialized_params
    end
end

"""Precompile events and config construction."""
function _precompile_events_and_callbacks()
    events = [BusFault(0.1, 0.2, 1, "F1"; xf = 0.01)]
    config = SimConfig(solver = (TDS = TDSConfig(
        abstol = 1e-8,
        reltol = 1e-6,
        convert_pq_to_zip = false,
    ),))
    return events, config
end

"""Create TDSProblem for precompilation."""
function _precompile_tds_problem(pf_sol, dyr_data, events, config)
    init = setup_tds(pf_sol, dyr_data; config = config)
    return TDSProblem(init;
        tspan = (0.0, 0.3),
        events = events
    )
end

# ============================================================================
# Main Workload
# ============================================================================

@setup_workload begin
    # Register this package's models during setup
    register_pending_models!(_PENDING_MODELS)

    @compile_workload begin
        # 1. System with power flow solution
        sys, pf_sol = _precompile_system(_PRECOMPILE_RAW)

        # 2. DYR parsing
        dyr_data = _precompile_parse_dynamics(sys, _PRECOMPILE_DYR)

        # 3. Model metadata access
        _precompile_model_metadata()

        # 4. Events and config
        events, config = _precompile_events_and_callbacks()

        # 5. Create TDSProblem (compiles TDSAdapter, ADCache, etc.)
        prob = _precompile_tds_problem(pf_sol, dyr_data, events, config)

        # 6. Precompile TDS infrastructure (ForwardDiff backend - default)
        adcache = create_adcache(prob)
        adapter = TDSAdapter(adcache, prob)
        create_diffeq_problem(adapter, prob.tspan)
        _simulation_callbacks(prob, adapter, _default_initializealg(prob); abstol = 1e-8)

        # 7. Compile ALL calc_residuals! kernels for Float64
        # This compiles: gather_variables!, gather_residuals!, calc_residuals!, scatter_residuals!
        # for EVERY model type in the system (Bus, Branch, PQ, PV, GENROU, SEXS, TGOV1, etc.)
        evaluate_system_residuals!(adapter.sys_ws, adapter.f_buf, adapter.g_buf,
                                   adapter.x_buf, adapter.y_buf, adapter.status, prob.tspan[1])

        # 8. Precompile Enzyme backend infrastructure
        # This significantly reduces first-run time when using EnzymeBackend()
        enzyme_adcache = create_adcache(prob; backend = EnzymeBackend())
        enzyme_adapter = TDSAdapter(enzyme_adcache, prob)
        create_diffeq_problem(enzyme_adapter, prob.tspan)

        # Note: We skip DiffEq.solve() because Julia can't cache cross-package
        # method specializations. Skipping saves significant precompile time.
    end
end
