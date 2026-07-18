"""
    EarthStateInterface

Small, solver-independent contracts for exchanging physically described Earth-system
state. The initial API focuses on one atmospheric column, including thermodynamics,
trace gases, and optional sectional aerosols.

The package deliberately does not perform radiative transfer, transport, meteorological
I/O, or aerosol microphysics. Those packages own their computations and implement this
interface at their boundaries.
"""
module EarthStateInterface

include("common.jl")
include("vertical_coordinates.jl")
include("constituents.jl")
include("aerosols.jl")
include("atmospheric_columns.jl")
include("validation.jl")
include("layer_reduction.jl")
include("optical_reduction.jl")
include("gchp_schema.jl")

export Unavailable, isavailable

export VerticalOrientation, TopToBottom, BottomToTop
export PressureCoordinate, pressure_interfaces, pressure_centers
export vertical_orientation, nlayers
export RepresentativePressure, ArithmeticMeanPressure, LogMeanPressure
export representative_pressure

export ConcentrationBasis, DryMoleFraction, WetMoleFraction
export MassMixingRatio, NumberDensity, MassDensity
export LayerMoles, LayerMass, LayerParticleNumber
export ConstituentField, concentration_basis, constituent_values

export ParticleSizeCoordinate, ParticleDiameter, ParticleRadius, ParticleMass
export SectionalSizeGrid, SectionalAerosolState, NoAerosols
export aerosol_number, aerosol_components, aerosol_size_grid

export AtmosphericColumn, atmospheric_column
export pressure_coordinate, temperature, specific_humidity
export dry_air_mass, dry_air_column_moles, layer_volume
export trace_gases, aerosols, metadata
export species_names, constituent, has_constituent

export validate_column, validate_sectional_aerosols

export LayerPartition, partition_by_count, partition_at_interfaces
export StateReduction, DryAirMassWeightedState
export merge_column

export OpticalReductionRoute, StateThenOptics, OpticsThenMerge
export merge_optical_depth, effective_cross_section

export GCHPColumnSchema, GCHPLocation, GCHPSectionalAerosolSchema
export GCHPColumnSource, open_gchp, read_column
export gchp_gas_variable, gchp_tomas_variable

end
