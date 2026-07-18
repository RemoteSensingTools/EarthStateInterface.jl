"""
    GCHPColumnSchema(; kwargs...)

Names and physical conventions needed by an optional GCHP/GEOS-Chem reader. This type is
kept free of NetCDF dependencies so the core interface remains lightweight.

Defaults match common GCHP History output: model storage is surface-to-TOA, pressure
thickness and surface pressure are in hPa, humidity may be g kg⁻¹, and
`SpeciesConcVV_*` gases/aerosols are dry-air mole fractions.
"""
Base.@kwdef struct GCHPColumnSchema
    pressure_thickness::String = "Met_DELP"
    dry_pressure_thickness::String = "Met_DELPDRY"
    surface_pressure::String = "Met_PS2WET"
    temperature::String = "Met_T"
    specific_humidity::String = "Met_SPHU"
    dry_air_mass::String = "Met_AD"
    layer_volume::String = "Met_AIRVOL"
    gas_prefix::String = "SpeciesConcVV_"
    tomas_number_prefix::String = "SpeciesConcVV_NK"
    stored_orientation::VerticalOrientation = BottomToTop()
end

"""Cubed-sphere cell and time index used by the optional GCHP reader."""
Base.@kwdef struct GCHPLocation
    x::Int
    y::Int
    face::Int
    time::Int = 1
end

"""
    GCHPSectionalAerosolSchema(size_grid; components, air_molar_mass=...)

Mapping for a sectional aerosol payload in GCHP History output. The size grid is required
explicitly: TOMAS mass bins do not imply unique physical diameter edges without a density
or composition assumption, so this interface does not invent a placeholder diameter.
"""
struct GCHPSectionalAerosolSchema{G, C, FT}
    size_grid::G
    components::C
    air_molar_mass::FT
end

function GCHPSectionalAerosolSchema(size_grid::SectionalSizeGrid;
                                    components=(:SF, :SS, :DUST, :OCOB, :OCIL,
                                                :ECOB, :ECIL, :AW),
                                    air_molar_mass=28.9644e-3)
    return GCHPSectionalAerosolSchema{typeof(size_grid), typeof(components),
                                      typeof(air_molar_mass)}(
        size_grid, components, air_molar_mass)
end

"""Open collection of GCHP datasets. Construct with [`open_gchp`](@ref)."""
struct GCHPColumnSource{D, S <: GCHPColumnSchema, A, P}
    datasets::D
    schema::S
    aerosol_schema::A
    paths::P
end

"""Open one or more GCHP History files. Implemented by the NCDatasets extension."""
function open_gchp end

"""Read one canonical atmospheric column from a column source."""
function read_column end

gchp_gas_variable(schema::GCHPColumnSchema, species::Symbol) =
    string(schema.gas_prefix, species)

function gchp_tomas_variable(schema::GCHPColumnSchema,
                             component::Symbol,
                             bin::Integer)
    bin > 0 || throw(ArgumentError("TOMAS bin index must be positive"))
    prefix = component === :NK ? schema.tomas_number_prefix :
        string(schema.gas_prefix, component)
    return string(prefix, lpad(bin, 2, '0'))
end
