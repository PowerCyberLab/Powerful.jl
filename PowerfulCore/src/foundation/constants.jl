"""
    BUS_PQ

Bus type code 1 (PSS/E `IDE`): load bus, with specified active and reactive power
injection.

See also [`BUS_PV`](@ref), [`BUS_SLACK`](@ref), [`BUS_ISOLATED`](@ref).
"""
const BUS_PQ = Int8(1)

"""
    BUS_PV

Bus type code 2 (PSS/E `IDE`): generator bus with specified active power and
controlled voltage magnitude.
"""
const BUS_PV = Int8(2)

"""
    BUS_SLACK

Bus type code 3 (PSS/E `IDE`): reference (swing) bus with specified voltage magnitude
and angle.
"""
const BUS_SLACK = Int8(3)

"""
    BUS_ISOLATED

Bus type code 4 (PSS/E `IDE`): isolated (disconnected) bus.
"""
const BUS_ISOLATED = Int8(4)

export BUS_PQ, BUS_PV, BUS_SLACK, BUS_ISOLATED
