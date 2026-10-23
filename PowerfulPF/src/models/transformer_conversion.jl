"""
Three-winding transformer conversion to star configuration with virtual buses.

Handles PSS/E data format conversion modes (CZ, CM) properly before computing
equivalent star impedances.
"""

using StructArrays
using InlineStrings

export convert_3w_to_star!

# ============================================================================
# 3-Winding Transformer Status Helpers
# ============================================================================

"""
    _is_winding_in_service(stat::Int8, winding_num::Int) -> Bool

Determine if a specific winding is in service based on 3WT status code.

# PSS/E 3WT Status Codes:
- 0: All windings out of service
- 1: All windings in service
- 2: Winding 2 out of service
- 3: Winding 3 out of service
- 4: Only winding 1 in service (windings 2 and 3 out)
"""
function _is_winding_in_service(stat::Int8, winding_num::Int)::Bool
    stat == 0 && return false  # All off
    stat == 1 && return true   # All on
    stat == 2 && return winding_num != 2  # W2 out
    stat == 3 && return winding_num != 3  # W3 out
    stat == 4 && return winding_num == 1  # Only W1 on
    return true  # Default: in service (unknown status treated as on)
end

# ============================================================================
# PSS/E Data Format Conversion Helpers
# ============================================================================

"""
    cz3_to_winding_pu(load_loss_w, z_mag, sbase_winding) -> (r, x)

Impedance given in PSS/E CZ=3 form (R = 3-phase load loss in W at rated flow,
X = |Z| in pu on the winding base) as (r, x) in pu on the winding base:
R_pu = P_loss_MW / S_winding_MVA and X = sqrt(|Z|² - R²). If R > |Z| (invalid
data), |Z| is used as X.
"""
function cz3_to_winding_pu(load_loss_w, z_mag, sbase_winding)
    r = load_loss_w / 1e6 / sbase_winding
    x = z_mag^2 >= r^2 ? sqrt(z_mag^2 - r^2) : z_mag
    return (r, x)
end

"""
    tap_ratio_pu(windv, nomv, cw, bus_basekv) -> Float64

Off-nominal tap ratio in per unit of the winding's bus base voltage `bus_basekv` (kV),
from the winding voltage `windv` according to the RAW winding data code CW:
- CW=1: `windv` is already in pu of the bus base voltage
- CW=2: `windv` is the winding voltage in kV: `windv / bus_basekv`
- CW=3: `windv` is in pu of the winding nominal voltage `nomv` (kV, where 0 means the
  bus base voltage): `windv * nomv / bus_basekv`

`windv = 0` is taken as nominal (1.0): RAW files converted from MATPOWER carry its
"ratio 0 = nominal" convention (e.g. the phase shifters of case89pegase, which then
reproduce the solved voltages stored in the case file). Errors for negative `windv`,
other CW codes, or CW=2/3 with a non-positive bus base voltage.
"""
function tap_ratio_pu(windv, nomv, cw, bus_basekv)
    cw in (1, 2, 3) || error("Transformer winding data code CW=$cw is not supported (only CW=1, 2 and 3)")
    windv >= 0.0 || error("Transformer winding voltage WINDV must not be negative, got $windv")
    windv == 0.0 && return 1.0
    cw == 1 && return windv
    bus_basekv > 0.0 || error("CW=$cw tap conversion needs a positive bus base voltage, got $bus_basekv kV")
    cw == 2 && return windv / bus_basekv
    nomv_kv = nomv > 0.0 ? nomv : bus_basekv  # NOMV = 0: nominal voltage is the bus base voltage
    return windv * nomv_kv / bus_basekv
end

"""
    convert_impedance_to_pu(r, x, cz::Int8, sbase_winding::Float64, sbase_system::Float64) -> (Float64, Float64)

Convert transformer impedance from PSS/E format to per-unit on system base.

# PSS/E CZ Codes:
- 1: Impedance in pu on system MVA base (SBASE)
- 2: Impedance in pu on specified winding MVA base (SBASE_WINDING)
- 3: Load loss in watts and impedance magnitude in pu on winding MVA base

# Arguments
- `r`: Resistance value (interpretation depends on `cz`)
- `x`: Reactance value (interpretation depends on `cz`)
- `cz`: PSS/E impedance data I/O code
- `sbase_winding`: Winding MVA base (from SBASE1_2, SBASE2_3, or SBASE3_1)
- `sbase_system`: System MVA base

# Returns
- `(r_pu, x_pu)`: Resistance and reactance in pu on system MVA base
"""
function convert_impedance_to_pu(
    r::T,
    x::T,
    cz::Int8,
    sbase_winding::T,
    sbase_system::T
) where {T <: AbstractFloat}
    if cz == 1
        # Already in pu on system base
        return (r, x)

    elseif cz == 2
        # In pu on winding base - convert to system base
        # Z_system = Z_winding * (S_system / S_winding)
        ratio = sbase_system / sbase_winding
        return (r * ratio, x * ratio)

    elseif cz == 3
        # r = load loss in watts, x = impedance magnitude in pu on winding base
        #
        # PSS/E CZ=3 format:
        # - R1-2: Total 3-phase load loss in watts at rated power flow
        # - X1-2: Impedance magnitude |Z| in pu on winding 1 base
        #
        # Convert load loss from watts to per-unit resistance:
        # At rated current I_rated (corresponding to S_winding), the I²R losses = P_loss
        # In per-unit, I_rated = 1, so R_pu = P_loss / S_winding² in per-unit power
        # Equivalently: R_pu = P_loss_MW / S_winding_MVA
        r_pu_winding, x_pu_winding = cz3_to_winding_pu(r, x, sbase_winding)

        # Convert from winding base to system base
        # Z_system = Z_winding * (S_system / S_winding)
        ratio = sbase_system / sbase_winding
        r_pu = r_pu_winding * ratio
        x_pu = x_pu_winding * ratio

        return (r_pu, x_pu)

    else
        error("Invalid CZ code: $cz (must be 1, 2, or 3)")
    end
end

"""
    convert_magnetizing_to_pu(mag1, mag2, cm, sbase_winding, sbase_system, bus_basekv, nomv1) -> (g, b)

Magnetizing admittance from RAW transformer data, in pu on the system MVA base and the
winding 1 bus base voltage `bus_basekv` (kV). It is a shunt at the winding 1 bus.

# RAW CM codes
- CM=1: `mag1`, `mag2` are the conductance and susceptance already in pu on the system
  MVA base and winding 1 bus voltage base
- CM=2: `mag1` is the no-load loss in watts and `mag2` the exciting current in pu on the
  winding 1-2 MVA base `sbase_winding` and nominal voltage `nomv1` (kV; 0 means the bus
  base voltage). On that base g = P_NL / S and |y| = I_exc, and the magnetizing branch is
  inductive: b = -sqrt(|y|² - g²). Both are then moved to the system base by
  (sbase_winding / sbase_system) * (bus_basekv / nomv1)².

Errors for other CM codes, for CM=2 with a non-positive winding MVA base or bus base
voltage, or for CM=2 data whose exciting current is smaller than the no-load loss
conductance (no real susceptance).
"""
function convert_magnetizing_to_pu(mag1, mag2, cm, sbase_winding, sbase_system, bus_basekv, nomv1)
    cm == 1 && return (mag1, mag2)
    cm == 2 || error("Invalid CM code: $cm (must be 1 or 2)")
    sbase_winding > 0.0 || error("CM=2 magnetizing data needs a positive winding MVA base, got $sbase_winding")
    bus_basekv > 0.0 || error("CM=2 magnetizing data needs a positive bus base voltage, got $bus_basekv kV")

    # On the transformer base (sbase_winding, nominal winding 1 voltage)
    g_t = (mag1 / 1e6) / sbase_winding
    y_t = mag2
    y_t >= g_t || error("CM=2 exciting current $y_t pu is smaller than the no-load loss " *
                        "conductance $g_t pu on the winding base; no susceptance exists")
    b_t = -sqrt(y_t * y_t - g_t * g_t)

    # To the system MVA base and the winding 1 bus base voltage
    nomv_kv = nomv1 > 0.0 ? nomv1 : bus_basekv
    scale = (sbase_winding / sbase_system) * (bus_basekv / nomv_kv)^2
    return (g_t * scale, b_t * scale)
end

# ============================================================================
# Star Conversion Functions
# ============================================================================

"""
    compute_star_impedances(t3w, sbase_system::Float64) -> (ComplexF64, ComplexF64, ComplexF64)

Compute star (wye) equivalent impedances from 3-winding transformer delta impedances.

# Delta-to-Star Transformation:
Given delta impedances Z₁₂, Z₂₃, Z₃₁ (between windings), compute star impedances:
- Z₁ = (Z₁₂ + Z₃₁ - Z₂₃) / 2  (from bus i to star center)
- Z₂ = (Z₁₂ + Z₂₃ - Z₃₁) / 2  (from bus j to star center)
- Z₃ = (Z₂₃ + Z₃₁ - Z₁₂) / 2  (from bus k to star center)

# Arguments
- `t3w`: Transformer3W instance
- `sbase_system`: System MVA base

# Returns
- `(Z1, Z2, Z3)`: Star impedances as complex numbers (in pu on system base)
"""
function compute_star_impedances(t3w, sbase_system::Float64)
    # Convert each delta impedance to pu on system base
    r12_pu, x12_pu = convert_impedance_to_pu(
        t3w.r1_2, t3w.x1_2, t3w.cz, t3w.sbase1_2, sbase_system
    )

    r23_pu, x23_pu = convert_impedance_to_pu(
        t3w.r2_3, t3w.x2_3, t3w.cz, t3w.sbase2_3, sbase_system
    )

    r31_pu, x31_pu = convert_impedance_to_pu(
        t3w.r3_1, t3w.x3_1, t3w.cz, t3w.sbase3_1, sbase_system
    )

    # Form complex impedances
    Z12 = complex(r12_pu, x12_pu)
    Z23 = complex(r23_pu, x23_pu)
    Z31 = complex(r31_pu, x31_pu)

    # Delta to star conversion
    Z1 = (Z12 + Z31 - Z23) / 2
    Z2 = (Z12 + Z23 - Z31) / 2
    Z3 = (Z23 + Z31 - Z12) / 2

    # Validate results
    if real(Z1) < -1e-6 || real(Z2) < -1e-6 || real(Z3) < -1e-6
        @debug "Negative resistance in star impedance conversion (unusual but possible)" Z1 Z2 Z3 maxlog=3
    end

    if abs(Z1) < 1e-10 || abs(Z2) < 1e-10 || abs(Z3) < 1e-10
        @debug "Very small star impedance may cause numerical issues" Z1 Z2 Z3 maxlog=3
    end

    return (Z1, Z2, Z3)
end

"""
    create_virtual_bus(bus_num::Int32, all_buses, bus_idx, t3w) -> Bus

Create a virtual bus for 3-winding transformer star center.

The virtual bus inherits properties from the reference bus (winding 1 / bus i)
and is configured as a PQ bus with zero load. Its voltage is the star point voltage of
the record (`vmstar`, `anstar`), as every other bus holds the voltage of its bus record:
the power flow starts each bus from the voltages stored in the case.

# Arguments
- `bus_num`: Virtual bus number (from VIRTUAL_BUS_BASE range)
- `all_buses`: StructArray of all existing buses
- `bus_idx`: `build_index_map(all_buses, :i)`
- `t3w`: Transformer3W instance

# Returns
- Virtual `Bus` instance
"""
function create_virtual_bus(bus_num::Int32, all_buses, bus_idx, t3w)
    # Find all three connected buses
    terminal_keys = ((t3w.i,), (t3w.j,), (t3w.k,))
    all(k -> haskey(bus_idx, k), terminal_keys) ||
        error("Connected buses not found for 3W transformer: $(t3w.i), $(t3w.j), $(t3w.k)")
    # Use bus i as reference for other properties
    ref_bus = all_buses[bus_idx[terminal_keys[1]]]

    return Bus(;
        i = bus_num,
        name = String15("VBUS_$(t3w.ckt)"),  # Virtual Bus marker (InlineString15)
        basekv = ref_bus.basekv,  # Use bus i's base voltage
        ide = BUS_PQ,  # PQ type
        area = ref_bus.area,
        zone = ref_bus.zone,
        owner = ref_bus.owner,
        vm = t3w.vmstar,  # Star point voltage magnitude of the record (pu)
        va = t3w.anstar,  # Star point angle of the record (degrees, as the bus records before per unit)
        nvhi = 1.1,     # Normal voltage high limit
        nvlo = 0.9,     # Normal voltage low limit
        evhi = 1.2,     # Emergency voltage high limit
        evlo = 0.8,     # Emergency voltage low limit
    )
end

"""
    create_2w_from_star(t3w, from_bus::Int32, to_bus::Int32, Z_star::ComplexF64,
                        winding_num::Int, sbase_system::Float64, all_buses) -> Transformer2W

Create a 2-winding transformer from a 3-winding transformer star leg.

# Arguments
- `t3w`: Original Transformer3W instance
- `from_bus`: From bus number (original winding bus)
- `to_bus`: To bus number (virtual bus at star center)
- `Z_star`: Star impedance (complex, in pu on system base)
- `winding_num`: Which winding (1, 2, or 3)
- `sbase_system`: System MVA base
- `all_buses`: All buses (for looking up base voltages)
- `bus_idx`: `build_index_map(all_buses, :i)`

# Returns
- `Transformer2W` instance
"""
function create_2w_from_star(
    t3w,
    from_bus::Int32,
    to_bus::Int32,
    Z_star::ComplexF64,
    winding_num::Int,
    sbase_system::Float64,
    all_buses,
    bus_idx
)
    # Select winding parameters based on winding number
    if winding_num == 1
        windv_raw = t3w.windv1
        nomv = t3w.nomv1
        ang = t3w.ang1
    elseif winding_num == 2
        windv_raw = t3w.windv2
        nomv = t3w.nomv2
        ang = t3w.ang2
    elseif winding_num == 3
        windv_raw = t3w.windv3
        nomv = t3w.nomv3
        ang = t3w.ang3
    else
        error("Invalid winding number: $winding_num")
    end

    # Tap ratio in pu of the winding bus base voltage (see tap_ratio_pu)
    windv = tap_ratio_pu(windv_raw, nomv, t3w.cw, all_buses.basekv[bus_idx[(from_bus,)]])
    cw_out = Int8(1)  # windv is now in pu of the bus base voltage

    # Extract star impedance components (already in pu on system base)
    r_star = real(Z_star)
    x_star = imag(Z_star)

    # Magnetizing admittance placement:
    # The RAW model places it at the winding 1 bus, so only the winding 1 two-winding
    # transformer (from bus i) carries it; it is a shunt at bus i, outside the tap.
    local mag1_val::Float64, mag2_val::Float64
    if winding_num == 1
        mag1_val, mag2_val = convert_magnetizing_to_pu(
            t3w.mag1, t3w.mag2, t3w.cm, t3w.sbase1_2, sbase_system,
            all_buses.basekv[bus_idx[(from_bus,)]], t3w.nomv1
        )
    else
        # Secondary/tertiary windings: zero magnetizing admittance
        mag1_val = 0.0
        mag2_val = 0.0
    end

    # Create unique circuit ID within String3 (3 char) limit
    # Format: original_ckt + winding_num, intelligently truncated
    original_ckt_str = String(t3w.ckt)
    if length(original_ckt_str) <= 1
        # Original is 1 char or empty: append "_W{num}" → e.g., "1" → "1_1"
        ckt_id = String3(original_ckt_str * "_" * string(winding_num))
    elseif length(original_ckt_str) == 2
        # Original is 2 chars: append winding num → e.g., "AB" → "AB1"
        ckt_id = String3(original_ckt_str * string(winding_num))
    else
        # Original is 3+ chars: take first char + winding num → e.g., "ABC" → "A1"
        ckt_id = String3(original_ckt_str[1:1] * string(winding_num))
    end

    # Create 2W transformer
    return Transformer2W(;
        i = from_bus,
        j = to_bus,
        ckt = ckt_id,  # Unique circuit ID within 3-char limit
        cw = cw_out,    # CW=2/3 converted to CW=1
        cz = Int8(1),   # Impedance now in pu on system base (after conversion)
        cm = Int8(1),   # Magnetizing admittance now in pu on system base
        mag1 = mag1_val,  # Only winding 1 has magnetizing admittance
        mag2 = mag2_val,  # Zero for windings 2 and 3
        nmetr = t3w.nmetr,
        name = String15("$(t3w.name)_W$winding_num"),  # InlineString15 required
        # 3WT stat: 0=all off, 1=all on, 2=W2 out, 3=W3 out, 4=only W1 on
        status = _is_winding_in_service(t3w.stat, winding_num),
        o1 = t3w.o1,
        f1 = t3w.f1,
        o2 = t3w.o2,
        f2 = t3w.f2,
        o3 = t3w.o3,
        f3 = t3w.f3,
        o4 = t3w.o4,
        f4 = t3w.f4,
        vecgrp = t3w.vecgrp,
        r1_2 = r_star,
        x1_2 = x_star,
        sbase1_2 = sbase_system,  # Now on system base
        windv1 = windv,
        nomv1 = nomv,
        ang1 = ang,     # Preserve original winding phase shift
        windv2 = 1.0,   # Virtual bus side: unity tap
        nomv2 = nomv,   # Note: Transformer2W has no ang2 field
        # Inherit ratings from original transformer
        rata1 = winding_num == 1 ? t3w.rata1 : winding_num == 2 ? t3w.rata2 : t3w.rata3,
        ratb1 = winding_num == 1 ? t3w.ratb1 : winding_num == 2 ? t3w.ratb2 : t3w.ratb3,
        ratc1 = winding_num == 1 ? t3w.ratc1 : winding_num == 2 ? t3w.ratc2 : t3w.ratc3,
        # Control parameters (typically not used for virtual transformer legs)
        cod1 = Int8(0),
        cont1 = Int32(0),
        rma1 = 1.5,
        rmi1 = 0.5,
        vma1 = 1.5,
        vmi1 = 0.5,
        ntp1 = Int16(33),
        tab1 = Int64(0),  # Note: tab1 is Int64 in Transformer2W struct
        cr1 = 0.0,
        cx1 = 0.0,
        cnxa1 = 0.0,
    )
end

# ============================================================================
# Main Conversion Function
# ============================================================================

"""
    convert_3w_to_star!(models::NamedTuple, sbase_system::Float64) -> (NamedTuple, Union{Nothing, Dict})

Convert 3-winding transformers to star configuration with virtual buses.

This function:
1. Creates virtual buses for each 3W transformer
2. Converts delta impedances to star impedances (properly handling CZ/CM codes)
3. Creates three 2W transformers for each 3W transformer
4. Removes Transformer3W from models tuple
5. Returns metadata for tracking virtual buses

# Arguments
- `models`: NamedTuple of model vectors (must include Transformer3W if present)
- `sbase_system`: System MVA base (for impedance conversions)

# Returns
- `(models_new, metadata)`: New models tuple without Transformer3W, and conversion metadata dictionary

# Metadata Dictionary Keys
- `"original_buses"`: Vector{Tuple{Int32, Int32, Int32}} - Original 3W transformer buses (i, j, k)
- `"virtual_buses"`: Vector{Int32} - Virtual bus numbers created
- `"circuit_ids"`: Vector{String} - Circuit IDs of original 3W transformers
"""
function convert_3w_to_star!(models::NamedTuple, sbase_system::Float64)
    # Check if 3W transformers exist
    if !haskey(models, :Transformer3W) || isempty(models.Transformer3W)
        return models, nothing
    end

    t3w = models.Transformer3W
    n_3w = length(t3w)

    @debug "Converting $n_3w three-winding transformer(s) to star configuration..."

    # Pre-allocate tracking arrays
    original_buses = Vector{Tuple{Int32, Int32, Int32}}(undef, n_3w)
    virtual_buses = Vector{Int32}(undef, n_3w)
    circuit_ids = Vector{String}(undef, n_3w)

    # Bus number → index, built once for all transformers
    bus_idx = build_index_map(models.Bus, :i)

    # Pre-allocate arrays for new components
    new_buses = Vector{Bus{Float64}}(undef, n_3w)
    new_t2w = Vector{Transformer2W{Float64}}(undef, 3 * n_3w)

    # Convert each 3W transformer
    for idx in 1:n_3w
        # Generate virtual bus number from safe range
        virtual_bus_num = Int32(VIRTUAL_BUS_BASE + idx)

        # Track conversion metadata
        original_buses[idx] = (t3w.i[idx], t3w.j[idx], t3w.k[idx])
        virtual_buses[idx] = virtual_bus_num
        circuit_ids[idx] = String(t3w.ckt[idx])

        # Create virtual bus
        new_buses[idx] = create_virtual_bus(virtual_bus_num, models.Bus, bus_idx, t3w[idx])

        # Compute star impedances (with proper CZ conversion)
        Z1, Z2, Z3 = compute_star_impedances(t3w[idx], sbase_system)

        # Create three 2W transformers (with proper CM conversion)
        new_t2w[3*(idx-1) + 1] = create_2w_from_star(t3w[idx], t3w.i[idx], virtual_bus_num, Z1, 1, sbase_system, models.Bus, bus_idx)
        new_t2w[3*(idx-1) + 2] = create_2w_from_star(t3w[idx], t3w.j[idx], virtual_bus_num, Z2, 2, sbase_system, models.Bus, bus_idx)
        new_t2w[3*(idx-1) + 3] = create_2w_from_star(t3w[idx], t3w.k[idx], virtual_bus_num, Z3, 3, sbase_system, models.Bus, bus_idx)

        @debug("  3W Transformer '$(t3w.ckt[idx])' ($(t3w.i[idx])-$(t3w.j[idx])-$(t3w.k[idx])) " *
                "→ Virtual bus $virtual_bus_num + 3 x 2W transformers")
    end

    # Combine buses and transformers
    bus_combined = vcat(models.Bus, StructArray(new_buses))
    t2w_combined = vcat(models.Transformer2W, StructArray(new_t2w))

    # Build new models tuple (without Transformer3W)
    new_models_vec = Pair{Symbol, Any}[]
    for name in keys(models)
        if name == :Bus
            push!(new_models_vec, name => bus_combined)
        elseif name == :Transformer2W
            push!(new_models_vec, name => t2w_combined)
        elseif name == :Transformer3W
            # Skip - fully converted
            continue
        else
            push!(new_models_vec, name => models[name])
        end
    end

    # Create conversion metadata as dictionary
    t3w_metadata = Dict{String,Any}(
        "original_buses" => original_buses,
        "virtual_buses" => virtual_buses,
        "circuit_ids" => circuit_ids
    )

    return NamedTuple(new_models_vec), t3w_metadata
end

