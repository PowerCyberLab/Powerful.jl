# ============================================================================
# Standard Group Definitions
# ============================================================================
#
# This file is AUTO-GENERATED from groups.toml.
# DO NOT EDIT MANUALLY.
#
# ============================================================================

"""
    register_standard_groups!()

Register the standard dynamic model groups.

Call this function during module initialization to register the built-in
group interfaces (generators, governors, exciters, stabilizers).
"""
function register_standard_groups!()

    # -------------------------------------------------------------------------
    # Bus_Group Group
    # -------------------------------------------------------------------------
    # Bus-like models that provide voltage interface (Bus, or other bus models)
    register_group!(GroupInterface(
        name = :bus_group,
        required_vars = [
            GroupVar(:v, :algebraic; description="Voltage magnitude (pu)"),
            GroupVar(:theta, :algebraic; description="Voltage angle (rad)"),
        ],
        required_residuals = [
            GroupRes(:p, :algebraic, :shared; description="Active power balance equation"),
            GroupRes(:q, :algebraic, :shared; description="Reactive power balance equation"),
        ],
        optional_residuals = [
            GroupRes(:ip, :algebraic, :shared; description="Active current injection, in phase with the bus voltage (pu), positive into the bus"),
            GroupRes(:iq, :algebraic, :shared; description="Reactive current injection, in quadrature to the bus voltage (pu), positive into the bus; q = v iq"),
        ],
        key_fields = [:i]
    ))

    # -------------------------------------------------------------------------
    # Generators Group
    # -------------------------------------------------------------------------
    # Synchronous machines (GENCLS, GENROU, GENSAL, etc.)
    register_group!(GroupInterface(
        name = :generators,
        required_vars = [
            GroupVar(:omega, :state; description="Rotor speed (pu)"),
            GroupVar(:delta, :state; description="Rotor angle (rad)"),
            GroupVar(:P_m, :algebraic; description="Mechanical power input (pu)"),
            GroupVar(:E_fd, :algebraic; description="Field voltage input (pu)"),
        ],
        required_residuals = [
            GroupRes(:g_P_m, :algebraic, :shared; description="Mechanical power input residual"),
            GroupRes(:g_E_fd, :algebraic, :shared; description="Field voltage input residual"),
        ],
        key_fields = [:BUS, :ID]
    ))

    # -------------------------------------------------------------------------
    # Governors Group
    # -------------------------------------------------------------------------
    # Turbine-governors (TGOV1, HYGOV, IEEEG1, etc.)
    register_group!(GroupInterface(
        name = :governors,
        required_vars = [
            GroupVar(:pout, :algebraic; description="Turbine mechanical power output (pu)"),
        ],
        required_residuals = GroupRes[],
        key_fields = [:BUS, :ID]
    ))

    # -------------------------------------------------------------------------
    # Exciters Group
    # -------------------------------------------------------------------------
    # Excitation systems (SEXS, ESST1A, EXST1, etc.)
    register_group!(GroupInterface(
        name = :exciters,
        required_vars = [
            GroupVar(:E_fd, :algebraic; description="Field voltage output (pu)"),
        ],
        required_residuals = GroupRes[],
        key_fields = [:BUS, :ID]
    ))

    # -------------------------------------------------------------------------
    # Stabilizers Group
    # -------------------------------------------------------------------------
    # Power system stabilizers (PSS2A, STAB1, etc.)
    register_group!(GroupInterface(
        name = :stabilizers,
        required_vars = [
            GroupVar(:vs, :algebraic; description="Stabilizer output signal (pu)"),
        ],
        required_residuals = GroupRes[],
        key_fields = [:BUS, :ID]
    ))

    # -------------------------------------------------------------------------
    # Ren_Gen Group
    # -------------------------------------------------------------------------
    # Renewable generator/converter models (REGCA1, etc.)
    register_group!(GroupInterface(
        name = :ren_gen,
        required_vars = [
            GroupVar(:Ipcmd, :algebraic; description="Active current command (pu, Mbase)"),
            GroupVar(:Iqcmd, :algebraic; description="Reactive current command (pu, Mbase); positive absorbs reactive power"),
            GroupVar(:Pe, :algebraic; description="Active power injected into the bus (pu, Mbase), positive into the grid"),
            GroupVar(:Qe, :algebraic; description="Reactive power injected into the bus (pu, Mbase), positive into the grid"),
        ],
        required_residuals = [
            GroupRes(:g_Ipcmd, :algebraic, :shared; description="Active current command residual"),
            GroupRes(:g_Iqcmd, :algebraic, :shared; description="Reactive current command residual"),
        ],
        key_fields = [:BUS, :ID]
    ))

    # -------------------------------------------------------------------------
    # Ren_Exciter Group
    # -------------------------------------------------------------------------
    # Renewable electrical controls (REECA1, etc.)
    register_group!(GroupInterface(
        name = :ren_exciter,
        required_vars = [
            GroupVar(:Ipcmd_out, :algebraic; description="Active current command output (pu)"),
            GroupVar(:Iqcmd_out, :algebraic; description="Reactive current command output (pu); positive injects reactive power"),
        ],
        required_residuals = GroupRes[],
        key_fields = [:BUS, :ID]
    ))

    # -------------------------------------------------------------------------
    # Ren_Gfm Group
    # -------------------------------------------------------------------------
    # Grid-forming converter models (REGFM_A1, etc.)
    register_group!(GroupInterface(
        name = :ren_gfm,
        required_vars = [
            GroupVar(:Pref, :algebraic; description="Active power reference, converter base (pu)"),
            GroupVar(:Qref, :algebraic; description="Reactive power reference, converter base (pu)"),
            GroupVar(:Vref, :algebraic; description="Voltage set point (pu)"),
            GroupVar(:Pe, :algebraic; description="Active power injected into the bus (pu), positive into the grid"),
            GroupVar(:Qe, :algebraic; description="Reactive power injected into the bus (pu), positive into the grid"),
        ],
        required_residuals = [
            GroupRes(:g_Pref, :algebraic, :shared; description="Active power reference residual"),
            GroupRes(:g_Qref, :algebraic, :shared; description="Reactive power reference residual"),
            GroupRes(:g_Vref, :algebraic, :shared; description="Voltage set point residual"),
        ],
        key_fields = [:BUS, :ID]
    ))

    return nothing
end
