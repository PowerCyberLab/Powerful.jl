# ============================================================================
# Nonlinear Network Current Injections
# ============================================================================
#
# Current injections that are not constant admittances in Y: the constant-power
# and constant-current parts of loads, and current sources (REGCA1, ...). They
# make the network equations nonlinear and are solved by the network Newton
# iteration (network_solver.jl), in rectangular (V_r, V_i) coordinates because
# they depend on |V| and so are not holomorphic in V.
#
# ============================================================================

using PowerfulPF: PQ, calc_load_current, calc_load_current_jacobian

"""
    injects_nonlinear_current(model::StructArray, idx::Int) -> Bool

Whether instance `idx` of a `NonConformingLoad` model has a current injection
outside the admittance matrix. A PQ load does when any of its constant-power or
constant-current parts (`pl`, `ql`, `ip`, `iq`) is nonzero; its constant-impedance
part (`yp`, `yq`) is always in Y. Current sources always do.
"""
injects_nonlinear_current(model::StructArray, idx::Int) = true
injects_nonlinear_current(model::StructArray{<:PQ}, idx::Int) =
    !(iszero(model.pl[idx]) && iszero(model.ql[idx]) &&
      iszero(model.ip[idx]) && iszero(model.iq[idx]))

"""
    nc_load_current(model::StructArray{<:PQ}, idx, v_r, v_i, vmin) -> (i_r=..., i_i=...)

Current drawn by the constant-power and constant-current parts of PQ load `idx` at
the rectangular voltage (`calc_load_current` without the constant-impedance part,
which is in Y). Consumption convention: the injection into the bus is the negative.
"""
@inline function nc_load_current(model::StructArray{<:PQ}, idx::Int, v_r, v_i, vmin)
    return calc_load_current(v_r, v_i, model.pl[idx], model.ql[idx], model.ip[idx],
                             model.iq[idx], zero(v_r), zero(v_r), vmin)
end

"""
    nc_load_current_jacobian(model::StructArray{<:PQ}, idx, v_r, v_i, vmin) -> (rr, ri, ir, ii)

Rectangular derivative of [`nc_load_current`](@ref) (`calc_load_current_jacobian`
without the constant-impedance part).
"""
@inline function nc_load_current_jacobian(model::StructArray{<:PQ}, idx::Int, v_r, v_i, vmin)
    return calc_load_current_jacobian(v_r, v_i, model.pl[idx], model.ql[idx], model.ip[idx],
                                      model.iq[idx], zero(v_r), zero(v_r), vmin)
end

"""
    source_current(v_r, v_i, Ip, Iq) -> (i_r=..., i_i=...)

Current injected by a current source with active and reactive components `Ip`, `Iq`
(`Iq > 0` generates reactive power) aligned with the terminal voltage:
I = (Ip - j Iq) V / |V|, so S = V I* = (Ip + j Iq)|V|.
"""
@inline function source_current(v_r, v_i, Ip, Iq)
    inv_v = one(v_r) / sqrt(v_r * v_r + v_i * v_i)
    return (i_r=(Ip * v_r + Iq * v_i) * inv_v, i_i=(Ip * v_i - Iq * v_r) * inv_v)
end

"""
    source_current_jacobian(v_r, v_i, Ip, Iq) -> (rr, ri, ir, ii)

Exact rectangular derivative of [`source_current`](@ref) (derived with SymPy). With
m = |V|, c_r + j c_i = (Ip - j Iq) V:

    ∂i_r/∂v_r =  Ip/m - v_r c_r/m³      ∂i_r/∂v_i =  Iq/m - v_i c_r/m³
    ∂i_i/∂v_r = -Iq/m - v_r c_i/m³      ∂i_i/∂v_i =  Ip/m - v_i c_i/m³
"""
@inline function source_current_jacobian(v_r, v_i, Ip, Iq)
    inv_v = one(v_r) / sqrt(v_r * v_r + v_i * v_i)
    inv_v3 = inv_v * inv_v * inv_v
    c_r = Ip * v_r + Iq * v_i
    c_i = Ip * v_i - Iq * v_r
    return (rr=Ip * inv_v - v_r * c_r * inv_v3, ri=Iq * inv_v - v_i * c_r * inv_v3,
            ir=-Iq * inv_v - v_r * c_i * inv_v3, ii=Ip * inv_v - v_i * c_i * inv_v3)
end

"""
    nonlinear_instances(sys, model_name::Symbol) -> Vector{Int}

In-service instances of `model_name` with a current injection outside Y
(`injects_nonlinear_current`).
"""
function nonlinear_instances(sys, model_name::Symbol)
    model = sys.models[model_name]
    active = in_service(model)
    return [idx for idx in eachindex(model) if active[idx] && injects_nonlinear_current(model, idx)]
end

"""
    network_load_models(sys) -> Tuple

The voltage-dependent load models whose current is iterated by the network Newton
solve: one `(model, instances, buses)` entry per `NonConformingLoad` model that is
not a current source (current sources go through `CurrentSourceBuffer`), where
`instances` are the `nonlinear_instances` and `buses` their internal bus indices. A
tuple, so iteration over it is type-stable.
"""
function network_load_models(sys)
    entries = Any[]
    for model_name in keys(sys.models)
        model = sys.models[model_name]
        network_partition_role(eltype(model)) isa NonConformingLoad || continue
        nc_load_type(eltype(model)) isa CurrentSourceNCLoad && continue
        instances = nonlinear_instances(sys, model_name)
        isempty(instances) && continue
        push!(entries, (model, instances, device_bus_indices(sys, model_name)[instances]))
    end
    return Tuple(entries)
end

# --- Load current and Jacobian accumulation over the load-model tuple --------

@inline _add_load_currents!(I_inj, ::Tuple{}, V, vmin) = nothing

@inline function _add_load_currents!(I_inj, loads::Tuple, V, vmin)
    model, instances, buses = first(loads)
    @inbounds for k in eachindex(instances)
        pos = buses[k]
        v = V[pos]
        i = nc_load_current(model, instances[k], real(v), imag(v), vmin)
        I_inj[pos] -= complex(i.i_r, i.i_i)  # consumption → injection
    end
    _add_load_currents!(I_inj, Base.tail(loads), V, vmin)
end

@inline _add_load_jacobians!(blocks, ::Tuple{}, V, vmin) = nothing

@inline function _add_load_jacobians!(blocks, loads::Tuple, V, vmin)
    model, instances, buses = first(loads)
    @inbounds for k in eachindex(instances)
        pos = buses[k]
        v = V[pos]
        d = nc_load_current_jacobian(model, instances[k], real(v), imag(v), vmin)
        # Mismatch = Y V - I_source - I_inj with I_inj = -I_load: ∂/∂V gains +∂I_load/∂V
        b = blocks[pos]
        blocks[pos] = (b[1] + d.rr, b[2] + d.ri, b[3] + d.ir, b[4] + d.ii)
    end
    _add_load_jacobians!(blocks, Base.tail(loads), V, vmin)
end

export injects_nonlinear_current, nc_load_current, nc_load_current_jacobian
export source_current, source_current_jacobian
