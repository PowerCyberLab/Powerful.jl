# Precompilation workload for PowerfulPF
#
# Precompiles the most common code paths:
# - parse_network (parsing RAW files)
# - SystemModel construction (per-unit conversion, indexing)
# - ACPowerFlowProblem creation (address management)
# - Newton solver (sparse Jacobian, KLU factorization)
# - Result extraction (get_var, collect_model_table)

using PrecompileTools

@setup_workload begin
    # Register this package's models (dependencies already registered via their __init__)
    register_pending_models!(_PENDING_MODELS)

    # Minimal 2-bus system in PSS/E RAW format (v33)
    # Bus 1: Slack (type 3), Bus 2: PV (type 2) with load
    # One branch connecting them
    MINIMAL_RAW = """
0,   100.00, 33, 0, 1, 60.00     / PSS(R)E-33 Minimal Test
 Precompile workload case
 2-bus system
    1,'Slack       ', 100.0000,3,   1,   1,   1,1.00000,   0.0000,1.10000,0.90000,1.10000,0.90000
    2,'GenBus      ', 100.0000,2,   1,   1,   1,1.00000,   0.0000,1.10000,0.90000,1.10000,0.90000
0 / END OF BUS DATA, BEGIN LOAD DATA
    2,'1 ',1,   1,   1,    50.000,    20.000,     0.000,     0.000,     0.000,     0.000,   1,1,0
0 / END OF LOAD DATA, BEGIN FIXED SHUNT DATA
0 / END OF FIXED SHUNT DATA, BEGIN GENERATOR DATA
    1,'1 ',     0.000,     0.000,   100.000,  -100.000,1.00000,     0,   100.000, 0.00000E+0, 1.00000E-1, 0.00000E+0, 0.00000E+0,1.00000,1,  100.0,   100.000,     0.000,   1,1.0000
    2,'1 ',    50.000,     0.000,   100.000,   -50.000,1.00000,     0,   100.000, 0.00000E+0, 2.00000E-1, 0.00000E+0, 0.00000E+0,1.00000,1,  100.0,   100.000,     0.000,   1,1.0000
0 / END OF GENERATOR DATA, BEGIN BRANCH DATA
     1,     2,'1 ', 1.00000E-2, 1.00000E-1,   0.00000,  100.00,  100.00,    0.00,  0.00000,  0.00000,  0.00000,  0.00000,1,1,   0.00,   1,1.0000
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

    @compile_workload begin
        # Parse RAW data from string (IOBuffer)
        case = parse_network(IOBuffer(MINIMAL_RAW))

        # Build SystemModel (per-unit conversion, ext2int indexing)
        sys = SystemModel(case)

        # Create power flow problem (address management)
        prob = ACPowerFlowProblem(sys)

        # Solve with Newton (sparse Jacobian, KLU)
        sol = solve(prob, Newton())

        # Result extraction
        _ = get_var(sol, :Bus, :v, (1,))
        _ = get_var(sol, :Bus, :theta, (2,))
        _ = collect_model_table(sol, :Bus)
    end
end
