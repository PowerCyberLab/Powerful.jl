export SwitchedShunt, SwitchedShuntVec

"""
    $TYPEDEF

Switched shunt: blocks of capacitors and reactors at a bus, switched in steps by a
control (`modsw`), with its initial admittance `binit`, the solved admittance in a solved
case. The power flow and the simulations hold every switched shunt at `binit`, as a fixed
susceptance with the physics of [`Shunt`](@ref) (`AbstractShunt`); its control is not
modeled, and the RAW reader warns of the in-service switched shunts with a control mode
(`modsw ≠ 0`). The record's control data and blocks are kept for a control model.

The fields are those of the PSS/E RAW revision 33 record (PowerFlowData `SwitchedShunts33`),
which has one switched shunt per bus: the bus is the key.

# Fields

$TYPEDFIELDS
"""
@kwdef mutable struct SwitchedShunt{Tv} <: AbstractShunt
    "Bus number (the key: one switched shunt per bus)"
    i::Int32

    "Status: true = in service (STAT)"
    status::Bool

    """
    Control mode (MODSW): 0 fixed; 1 discrete and 2 continuous adjustment, controlling the
    voltage at bus `i` or `swrem`; 3 discrete, controlling the reactive power of the plant at
    `swrem`; 4 discrete, controlling the reactive power of the VSC DC line converter at
    `swrem`; 5 discrete, controlling the admittance of the switched shunt at `swrem`.
    Not modeled: the shunt is held at `binit`.
    """
    modsw::Int8

    "Adjustment method (ADJM): false = blocks in input order, true = nearest total admittance"
    adjm::Bool

    "Upper limit of the controlled voltage (modsw 1, 2) or reactive power range (modsw 3–5), pu"
    vswhi::Float64

    "Lower limit of the controlled voltage (modsw 1, 2) or reactive power range (modsw 3–5), pu"
    vswlo::Float64

    "Bus whose voltage or equipment reactive power is controlled (SWREM); 0 = bus `i`"
    swrem::Int32

    "Percent of the total Mvar for the voltage at a remotely controlled bus (RMPCT)"
    rmpct::Float64

    "Name of the VSC DC line whose converter it controls (RMIDNT, modsw 4)"
    rmidnt::String15

    """
    Initial admittance (BINIT), Mvar at 1.0 pu voltage: the susceptance it is held at.
    Positive = capacitive (supplies reactive power), negative = inductive.
    """
    binit::Tv

    "Number of steps of block 1 (N1); the first zero N or B ends the blocks"
    n1::Int32
    "Admittance of each step of block 1 (B1), Mvar at 1.0 pu voltage"
    b1::Tv
    "Number of steps of block 2"
    n2::Int32
    "Admittance of each step of block 2, Mvar at 1.0 pu voltage"
    b2::Tv
    "Number of steps of block 3"
    n3::Int32
    "Admittance of each step of block 3, Mvar at 1.0 pu voltage"
    b3::Tv
    "Number of steps of block 4"
    n4::Int32
    "Admittance of each step of block 4, Mvar at 1.0 pu voltage"
    b4::Tv
    "Number of steps of block 5"
    n5::Int32
    "Admittance of each step of block 5, Mvar at 1.0 pu voltage"
    b5::Tv
    "Number of steps of block 6"
    n6::Int32
    "Admittance of each step of block 6, Mvar at 1.0 pu voltage"
    b6::Tv
    "Number of steps of block 7"
    n7::Int32
    "Admittance of each step of block 7, Mvar at 1.0 pu voltage"
    b7::Tv
    "Number of steps of block 8"
    n8::Int32
    "Admittance of each step of block 8, Mvar at 1.0 pu voltage"
    b8::Tv
end

@register_model SwitchedShunt

"""
    SwitchedShuntVec{T}

Alias for `StructArray{SwitchedShunt{T}}`, the struct-of-arrays container holding all
[`SwitchedShunt`](@ref) instances of a system.
"""
const SwitchedShuntVec{T} = StructArray{SwitchedShunt{T}}

get_key_name(::Type{T}) where T<:SwitchedShunt = [:i]

# Admittance held: the susceptance binit (a switched shunt has no conductance)
shunt_admittance_params(shunt::StructArray{<:SwitchedShunt}) = (binit = shunt.binit,)
Base.@propagate_inbounds shunt_admittance(params::NamedTuple{(:binit,)}, i) = (zero(eltype(params.binit)), params.binit[i])

# binit and the block admittances are in Mvar at 1.0 pu voltage
const SWITCHED_SHUNT_PARAMS = [
    ModelParam(:binit, PowerNominal, description="Initial admittance (held)", units="MVAR");
    [ModelParam(Symbol(:b, k), PowerNominal, description="Admittance of each step of block $k", units="MVAR")
     for k in 1:8]
]

model_metadata(::Type{SwitchedShunt}) = shunt_metadata(:SwitchedShunt, SWITCHED_SHUNT_PARAMS)

### === Model Format Adapters ===

supports_format(::Type{SwitchedShunt}, ::Type{PSSE}) = FormatSupport{PSSE}()
format_records(::Type{SwitchedShunt}, ::Type{PSSE}) = (:switched_shunts,)

"""
    parse_model(::Type{SwitchedShunt}, raw, ::FormatSupport{PSSE})

The switched shunts of a RAW case (revision 33 records). Refuses a bus with more than
one switched shunt record (the bus is the key), and warns of the in-service switched
shunts with a control mode (`modsw ≠ 0`), which are held at `binit`.
"""
function parse_model(::Type{SwitchedShunt}, raw, ::FormatSupport{PSSE})
    ss = raw.switched_shunts
    records = Dict{Int32, Int}()
    for b in ss.i
        records[b] = get(records, b, 0) + 1
    end
    repeated = sort([b for (b, n) in records if n > 1])
    isempty(repeated) || throw(ArgumentError(
        "More than one switched shunt record at bus " * join(repeated, ", ") *
        ": a RAW revision 33 case has one switched shunt per bus, its key"))

    controlled = count(k -> ss.stat[k] && ss.modsw[k] != 0, eachindex(ss.i))
    controlled > 0 && @warn "Switched shunts in service with a control mode (modsw ≠ 0): $controlled. " *
        "Each is held at its initial admittance binit: switched shunt control is not modeled, " *
        "which matters once the operating point moves from the case's (contingencies, simulations)"

    return StructArray{SwitchedShunt{Float64}}((
        i = ss.i, status = ss.stat, modsw = ss.modsw, adjm = ss.adjm,
        vswhi = ss.vswhi, vswlo = ss.vswlo, swrem = ss.swrem, rmpct = ss.rmpct,
        rmidnt = ss.rmidnt, binit = ss.binit,
        n1 = ss.n1, b1 = ss.b1, n2 = ss.n2, b2 = ss.b2, n3 = ss.n3, b3 = ss.b3, n4 = ss.n4, b4 = ss.b4,
        n5 = ss.n5, b5 = ss.b5, n6 = ss.n6, b6 = ss.b6, n7 = ss.n7, b7 = ss.b7, n8 = ss.n8, b8 = ss.b8,
    ))
end
