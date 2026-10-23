"""
Foundation module - Core type system for Powerful.jl

This module contains the foundational types, traits, and metadata system
that rarely change. It has minimal external dependencies (only StructArrays).

# Contents
- Core types: AbstractModel, ModelVar, ModelResidual, AddressManager, etc.
- Traits: FormulationTrait, BaseQuantityRole, FormatSupport
- Metadata: ModelMetadata, model_metadata interface
- Constants: Bus type codes
- Utilities: Coordinate conversions, metadata queries
"""
module Foundation

using StructArrays
using DocStringExtensions

# ============================================================================
# Core Types (must come first - defines AbstractModel, ModelMetadata, etc.)
# ============================================================================
include("system_config.jl")  # SystemConfig, owned by SystemModel (core_types.jl)
export SystemConfig
include("automata.jl")  # Automaton metadata, read by ModelMetadata (core_types.jl)
include("core_types.jl")
export VIRTUAL_BUS_BASE
export AbstractModel
export PrivateRes
export ModelRetrievedAddresses
export GroupLinkRef, is_heterogeneous

# ============================================================================
# Constants
# ============================================================================
include("constants.jl")

# ============================================================================
# Format Traits
# ============================================================================
include("format_traits.jl")
export NoFormatSupport

# ============================================================================
# Formulation Traits
# ============================================================================
include("formulation_traits.jl")

# ============================================================================
# Model Variable/Residual Constructors
# ============================================================================
include("model_var.jl")
export get_property, get_description, get_units, get_mass_coeff
export external_vars

# ============================================================================
# Model Metadata Interface
# ============================================================================
include("model_metadata.jl")

# ============================================================================
# Metadata Query Predicates
# ============================================================================
include("metadata_queries.jl")

# ============================================================================
# Base Quantity Traits
# ============================================================================
include("base_quantity_traits.jl")
export BaseQuantityRole, OwnsBase, InheritsBase, NoBase
export base_quantity_role, owns_base, inherits_base, needs_base
export get_base_quantity_vector

# ============================================================================
# Group System (Dynamic Model Interfaces)
# ============================================================================
include("groups.jl")
export get_controlled_group

# ============================================================================
# Model Traits (TDS/Dynamic Model Support)
# ============================================================================
include("model_traits.jl")
export PowerSystemComponent, PowerFlowModel, DynamicModel
export GeneratorModel, LoadModel, NetworkElement
export requires_initialization, get_dyr_model_name
export check_parameter_domains
export get_replacement_map, get_data_inheritance
export find_control_terminal, parse_dynamic_model

# ============================================================================
# Solution Interface (Abstract interface for steady-state solutions)
# ============================================================================
include("solution_interface.jl")

end # module
