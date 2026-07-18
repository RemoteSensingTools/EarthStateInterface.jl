# Interface functions are declared here so external packages can implement the contract
# for their own storage without inheriting from a package-owned abstract type.
function pressure_coordinate end
function temperature end
function specific_humidity end
function dry_air_mass end
function dry_air_column_moles end
function layer_volume end
function trace_gases end
function aerosols end
function metadata end

"""
    AtmosphericColumn(; pressure, temperature, specific_humidity,
                        dry_air_mass=Unavailable(),
                        dry_air_column_moles=Unavailable(),
                        layer_volume=Unavailable(),
                        trace_gases=(;), aerosols=NoAerosols(), metadata=(;))

Canonical single-column interchange container.

Required fields use SI units and share one vertical orientation:

- pressure interfaces: Pa
- temperature: K
- specific humidity: kg kg⁻¹ moist air

Trace gases and aerosols carry their own explicit concentration bases. Dry-air cell mass
and dry-air column moles are alternative carrier amounts: transport/GCHP commonly has the
former, while spectroscopy commonly has the latter. Layer volume is optional generally,
but required for conservative number- or mass-density reduction.
"""
struct AtmosphericColumn{P <: PressureCoordinate, T, Q, D, C, V,
                         G <: NamedTuple, A, M}
    pressure::P
    temperature::T
    specific_humidity::Q
    dry_air_mass::D
    dry_air_column_moles::C
    layer_volume::V
    trace_gases::G
    aerosols::A
    metadata::M
end

function AtmosphericColumn(; pressure::PressureCoordinate,
                           temperature,
                           specific_humidity,
                           dry_air_mass=Unavailable(),
                           dry_air_column_moles=Unavailable(),
                           layer_volume=Unavailable(),
                           trace_gases::NamedTuple=(;),
                           aerosols=NoAerosols(),
                           metadata=(;))
    column = AtmosphericColumn{typeof(pressure), typeof(temperature),
                               typeof(specific_humidity), typeof(dry_air_mass),
                               typeof(dry_air_column_moles), typeof(layer_volume),
                               typeof(trace_gases),
                               typeof(aerosols), typeof(metadata)}(
        pressure, temperature, specific_humidity, dry_air_mass,
        dry_air_column_moles, layer_volume, trace_gases, aerosols, metadata)
    validate_column(column)
    return column
end

atmospheric_column(; kwargs...) = AtmosphericColumn(; kwargs...)

pressure_coordinate(column::AtmosphericColumn) = column.pressure
temperature(column::AtmosphericColumn) = column.temperature
specific_humidity(column::AtmosphericColumn) = column.specific_humidity
dry_air_mass(column::AtmosphericColumn) = column.dry_air_mass
dry_air_column_moles(column::AtmosphericColumn) = column.dry_air_column_moles
layer_volume(column::AtmosphericColumn) = column.layer_volume
trace_gases(column::AtmosphericColumn) = column.trace_gases
aerosols(column::AtmosphericColumn) = column.aerosols
metadata(column::AtmosphericColumn) = column.metadata

pressure_interfaces(column::AtmosphericColumn) =
    pressure_interfaces(pressure_coordinate(column))
pressure_centers(column::AtmosphericColumn) =
    pressure_centers(pressure_coordinate(column))
vertical_orientation(column::AtmosphericColumn) =
    vertical_orientation(pressure_coordinate(column))
nlayers(column::AtmosphericColumn) = nlayers(pressure_coordinate(column))

species_names(column) = keys(trace_gases(column))
has_constituent(column, name::Symbol) = name in species_names(column)

function constituent(column, name::Symbol)
    has_constituent(column, name) || throw(KeyError(name))
    return getproperty(trace_gases(column), name)
end
