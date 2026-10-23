"""
Utility functions for working with model types and instances.
"""

using StructArrays: StructArray


#==============================================================================
# Algebraic initial guess helpers
==============================================================================#

"""
    set_algebraic_guess!(y, am, sys, formulation)

Populate an algebraic guess vector in-place using per-model hooks.
Dispatches to `set_algebraic_guess!(::AbstractVector, ::AddressManager, ::SystemModel,
::Type{<:AbstractModel}, ::FormulationTrait)` for each model in the system.
"""
function set_algebraic_guess!(
    y::AbstractVector,
    am::AddressManager,
    sys::SystemModel,
    formulation::FormulationTrait
)
    for model_vec in values(sys.models)
        model_type = get_base_model_type(model_vec)
        set_algebraic_guess!(y, am, sys, model_type, formulation)
    end
    return nothing
end

# Default no-op for models without a custom guess
function set_algebraic_guess!(
    y::AbstractVector,
    am::AddressManager,
    sys::SystemModel,
    ::Type{<:AbstractModel},
    formulation::FormulationTrait
)
    (; y, am, sys, formulation)  # suppress unused warnings
    return nothing
end
