# ============================================================================
# Y-Bus Residual Evaluation Trait
# ============================================================================
#
# Residuals of the machines with an internal voltage at the Y-bus network
# solution (their terminal voltages).
#
# ============================================================================

using PowerfulCore: ResidualBuffers, FormulationTrait, calc_residuals!

"""
    calc_ybus_residuals!(ws::ResidualBuffers, model::StructArray, external::NamedTuple,
                         network_sol::NetworkSolution, formulation::FormulationTrait,
                         model_to_gen_idx::Vector{Int}; active_mask)

Residuals of a machine with an internal voltage (`network_partition_role(M) isa
MachineWithInternalVoltage`) at the Y-bus network solution: the terminal voltages `v`,
`theta` of its in-service instances from `network_sol` (`gather_terminal_voltages!`), then
the model's own residuals at them (`calc_residuals!`, those of the phasor TDS: flux and
swing dynamics, stator current constraints, control inputs and power outputs).

`model_to_gen_idx` maps each model array index to its index among the network generators
(0 when out of service). A machine model provides its Norton source to the network with
[`get_machine_impedance`](@ref) and [`get_internal_voltage`](@ref).
"""
function calc_ybus_residuals!(
    ws::ResidualBuffers,
    model::StructArray,
    external::NamedTuple,
    network_sol::NetworkSolution,
    formulation::FormulationTrait,
    model_to_gen_idx::Vector{Int};
    active_mask::AbstractVector{Bool}
)
    gather_terminal_voltages!(ws, network_sol, model_to_gen_idx, active_mask)
    calc_residuals!(ws, model, external, formulation; active_mask=active_mask)
    return nothing
end

export calc_ybus_residuals!

"""
    gather_terminal_voltages!(ws, network_sol, model_to_gen_idx, active_mask)

Set `ws.var_vals.v` and `ws.var_vals.theta` of each active generator from its terminal
voltage `network_sol.V_terminal[model_to_gen_idx[i]]`. Inactive generators are not in
the network generator set (`model_to_gen_idx[i] == 0`) and their residuals are masked,
so they are skipped. Errors if an active generator is not in the set (e.g. it was put
in service after the network was built).
"""
function gather_terminal_voltages!(ws, network_sol, model_to_gen_idx::Vector{Int},
                                   active_mask::AbstractVector{Bool})
    @inbounds for i in eachindex(model_to_gen_idx)
        active_mask[i] || continue
        k = model_to_gen_idx[i]
        k > 0 || error("Generator instance $i is in service but was not in service when the " *
                       "Y-bus network was built")
        V_t = network_sol.V_terminal[k]
        ws.var_vals.v[i] = abs(V_t)
        ws.var_vals.theta[i] = angle(V_t)
    end
    return nothing
end
