"""Coordinate used for a sectional aerosol size grid."""
abstract type ParticleSizeCoordinate end

"""Particle diameter [m]."""
struct ParticleDiameter <: ParticleSizeCoordinate end

"""Particle radius [m]."""
struct ParticleRadius <: ParticleSizeCoordinate end

"""Dry particle mass [kg particle⁻¹]."""
struct ParticleMass <: ParticleSizeCoordinate end

"""
    SectionalSizeGrid(edges, coordinate; metadata=(;))

Section edges in SI units. A grid with `B+1` edges describes `B` bins.
"""
struct SectionalSizeGrid{E, C <: ParticleSizeCoordinate, M}
    edges::E
    coordinate::C
    metadata::M
end

function SectionalSizeGrid(edges, coordinate::ParticleSizeCoordinate; metadata=(;))
    length(edges) >= 2 || throw(ArgumentError("sectional size grid needs at least one bin"))
    _allfinite(edges) || throw(ArgumentError("sectional size edges must be finite"))
    all(>(zero(eltype(edges))), edges) || throw(ArgumentError(
        "sectional size edges must be positive SI quantities"))
    all(>(zero(eltype(edges))), diff(edges)) || throw(ArgumentError(
        "sectional size edges must increase strictly"))
    return SectionalSizeGrid{typeof(edges), typeof(coordinate), typeof(metadata)}(
        edges, coordinate, metadata)
end

Base.length(grid::SectionalSizeGrid) = length(grid.edges) - 1

"""Sentinel for a column without an aerosol payload."""
struct NoAerosols end

"""
    SectionalAerosolState(size_grid, number, components; metadata=(;))

Layer-resolved sectional aerosol state. Arrays use `(bin, layer)` order. `number` and
every named component are `ConstituentField`s, so GCHP-native dry mole fractions,
number densities, or already-converted extensive amounts cannot be confused.

Optical depth, single-scattering albedo, and phase moments do not belong here: they are
derived optical properties and may be computed with different microphysical assumptions.
"""
struct SectionalAerosolState{G, N <: ConstituentField, C <: NamedTuple, M}
    size_grid::G
    number::N
    components::C
    metadata::M
end

function SectionalAerosolState(size_grid::SectionalSizeGrid,
                               number::ConstituentField,
                               components::NamedTuple;
                               metadata=(;))
    state = SectionalAerosolState{typeof(size_grid), typeof(number),
                                  typeof(components), typeof(metadata)}(
        size_grid, number, components, metadata)
    validate_sectional_aerosols(state)
    return state
end

aerosol_size_grid(state::SectionalAerosolState) = state.size_grid
aerosol_number(state::SectionalAerosolState) = state.number
aerosol_components(state::SectionalAerosolState) = state.components

function validate_sectional_aerosols(state::SectionalAerosolState;
                                     expected_layers=nothing)
    number = constituent_values(state.number)
    ndims(number) == 2 || throw(DimensionMismatch(
        "sectional aerosol number must use (bin, layer) storage"))
    size(number, 1) == length(state.size_grid) || throw(DimensionMismatch(
        "aerosol number bin count does not match the size grid"))
    expected_layers === nothing || size(number, 2) == expected_layers ||
        throw(DimensionMismatch("aerosol layer count does not match the atmospheric column"))
    _allfinite(number) || throw(ArgumentError("aerosol number values must be finite"))

    for name in keys(state.components)
        component = getproperty(state.components, name)
        component isa ConstituentField || throw(ArgumentError(
            "aerosol component $name must be a ConstituentField"))
        size(constituent_values(component)) == size(number) || throw(DimensionMismatch(
            "aerosol component $name must match (bin, layer) number storage"))
        _allfinite(constituent_values(component)) || throw(ArgumentError(
            "aerosol component $name contains non-finite values"))
    end
    return state
end
