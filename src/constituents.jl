"""Physical basis carried by a constituent field."""
abstract type ConcentrationBasis end

"""Moles of constituent per mole of dry air [mol mol⁻¹ dry]."""
struct DryMoleFraction <: ConcentrationBasis end

"""Moles of constituent per mole of moist air [mol mol⁻¹ wet]."""
struct WetMoleFraction <: ConcentrationBasis end

"""Constituent mass per mass of carrier air [kg kg⁻¹]."""
struct MassMixingRatio <: ConcentrationBasis end

"""Number per volume [m⁻³]."""
struct NumberDensity <: ConcentrationBasis end

"""Mass per volume [kg m⁻³]."""
struct MassDensity <: ConcentrationBasis end

"""Moles contained in each layer [mol]."""
struct LayerMoles <: ConcentrationBasis end

"""Mass contained in each layer [kg]."""
struct LayerMass <: ConcentrationBasis end

"""Number of particles contained in each layer [1]."""
struct LayerParticleNumber <: ConcentrationBasis end

"""
    ConstituentField(values, basis; metadata=(;))

Values for a gas or aerosol component with an explicit physical basis. Metadata is for
provenance and molecular/material properties; algorithms dispatch on `basis`, never on a
unit string in metadata.
"""
struct ConstituentField{B <: ConcentrationBasis, V, M}
    values::V
    basis::B
    metadata::M
end

ConstituentField(values, basis::ConcentrationBasis; metadata=(;)) =
    ConstituentField{typeof(basis), typeof(values), typeof(metadata)}(
        values, basis, metadata)

concentration_basis(field::ConstituentField) = field.basis
constituent_values(field::ConstituentField) = field.values
