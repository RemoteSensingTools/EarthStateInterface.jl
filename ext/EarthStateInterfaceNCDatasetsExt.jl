module EarthStateInterfaceNCDatasetsExt

using EarthStateInterface
using NCDatasets

import Base: close
import EarthStateInterface: open_gchp, read_column

function open_gchp(paths::Union{AbstractString, AbstractVector{<:AbstractString}};
                   schema::GCHPColumnSchema=GCHPColumnSchema(),
                   aerosol_schema=nothing)
    normalized_paths = paths isa AbstractString ? [String(paths)] : String.(paths)
    isempty(normalized_paths) && throw(ArgumentError("at least one GCHP path is required"))
    opened = NCDataset[]
    try
        append!(opened, (NCDataset(path, "r") for path in normalized_paths))
        datasets = Tuple(opened)
        return GCHPColumnSource(datasets, schema, aerosol_schema, normalized_paths)
    catch
        foreach(close, opened)
        rethrow()
    end
end

function open_gchp(f::Function,
                   paths::Union{AbstractString, AbstractVector{<:AbstractString}};
                   kwargs...)
    source = open_gchp(paths; kwargs...)
    try
        return f(source)
    finally
        close(source)
    end
end

function close(source::GCHPColumnSource)
    foreach(close, source.datasets)
    return nothing
end

function _dataset_with(source::GCHPColumnSource, variable_name::AbstractString)
    for dataset in source.datasets
        haskey(dataset, variable_name) && return dataset
    end
    return nothing
end

function _require_variable(source::GCHPColumnSource, variable_name::AbstractString)
    dataset = _dataset_with(source, variable_name)
    dataset === nothing && throw(ArgumentError(
        "GCHP variable $variable_name was not found in $(source.paths)"))
    return dataset[variable_name]
end

function _dimension_index(dimension_name::AbstractString,
                          location::GCHPLocation,
                          variable,
                          dimension::Integer)
    lowercase_name = lowercase(dimension_name)
    if lowercase_name in ("xdim", "x", "lon", "longitude")
        return location.x
    elseif lowercase_name in ("ydim", "y", "lat", "latitude")
        return location.y
    elseif lowercase_name in ("nf", "face", "tile")
        return location.face
    elseif lowercase_name in ("time", "datetime")
        return location.time
    elseif lowercase_name in ("lev", "layer", "layers")
        return Colon()
    elseif size(variable, dimension) == 1
        return 1
    end
    throw(ArgumentError(
        "cannot map dimension '$dimension_name' of a GCHP variable; " *
        "configure or pre-slice the dataset"))
end

function _read_variable(source::GCHPColumnSource,
                        variable_name::AbstractString,
                        location::GCHPLocation)
    variable = _require_variable(source, variable_name)
    indices = Tuple(_dimension_index(String(dimension_name), location, variable, dimension)
                    for (dimension, dimension_name) in enumerate(dimnames(variable)))
    values = variable[indices...]
    return values isa AbstractArray ? Array(values) : fill(values), variable
end

function _read_profile(source, variable_name, location)
    values, variable = _read_variable(source, variable_name, location)
    return vec(values), variable
end

function _read_scalar(source, variable_name, location)
    values, variable = _read_variable(source, variable_name, location)
    length(values) == 1 || throw(DimensionMismatch(
        "$variable_name must be scalar at one GCHP location"))
    return only(values), variable
end

_units(variable) = haskey(variable.attrib, "units") ? String(variable.attrib["units"]) : ""

function _pressure_to_pa(values, variable_name, units)
    normalized = lowercase(strip(units))
    if normalized in ("pa", "pascal", "pascals")
        return values
    elseif normalized in ("hpa", "mb", "mbar", "millibar")
        return values .* 100
    end
    throw(ArgumentError(
        "$variable_name units '$units' are unsupported; expected Pa, hPa, or mb"))
end

function _humidity_to_kgkg(values, units)
    normalized = lowercase(replace(strip(units), " " => ""))
    if normalized in ("kgkg-1", "kg/kg", "1", "")
        return values
    elseif normalized in ("gkg-1", "g/kg")
        return values ./ 1000
    end
    throw(ArgumentError(
        "specific-humidity units '$units' are unsupported; expected kg kg-1 or g kg-1"))
end

function _to_top_down(values, orientation::TopToBottom)
    return collect(values)
end

function _to_top_down(values, orientation::BottomToTop)
    return reverse(collect(values))
end

function _trace_gases(source, location, gas_species, orientation)
    names = Symbol[]
    fields = ConstituentField[]
    for species in gas_species
        name = Symbol(species)
        variable_name = gchp_gas_variable(source.schema, name)
        _dataset_with(source, variable_name) === nothing && continue
        values, variable = _read_profile(source, variable_name, location)
        units = lowercase(replace(_units(variable), " " => ""))
        units in ("molmol-1dry", "mol/mol-dry", "1", "") || throw(ArgumentError(
            "$variable_name units '$(_units(variable))' are not dry mole fraction"))
        push!(names, name)
        push!(fields, ConstituentField(_to_top_down(values, orientation),
                                      DryMoleFraction();
                                      metadata=(source_variable=variable_name,)))
    end
    return NamedTuple{Tuple(names)}(Tuple(fields))
end

function _sectional_aerosols(source,
                             location,
                             ::Nothing,
                             orientation,
                             dry_mass,
                             volume)
    return NoAerosols()
end

function _sectional_aerosols(source,
                             location,
                             mapping::GCHPSectionalAerosolSchema,
                             orientation,
                             dry_mass,
                             volume)
    nbin = length(mapping.size_grid)
    nlay = length(dry_mass)
    number = Matrix{promote_type(eltype(dry_mass), eltype(volume))}(undef, nbin, nlay)
    dry_air_moles = dry_mass ./ mapping.air_molar_mass

    for bin in 1:nbin
        variable_name = gchp_tomas_variable(source.schema, :NK, bin)
        native, _ = _read_profile(source, variable_name, location)
        values = _to_top_down(native, orientation)
        # GCHP TOMAS NK: 1000 * particles / mol dry air.
        number[bin, :] .= (values ./ 1000) .* dry_air_moles ./ volume
    end

    component_names = Tuple(Symbol.(mapping.components))
    component_fields = map(component_names) do component
        values = similar(number)
        for bin in 1:nbin
            variable_name = gchp_tomas_variable(source.schema, component, bin)
            native, variable = _read_profile(source, variable_name, location)
            units = lowercase(replace(_units(variable), " " => ""))
            units in ("molmol-1dry", "mol/mol-dry", "1", "") || throw(ArgumentError(
                "$variable_name units '$(_units(variable))' are not dry mole fraction"))
            values[bin, :] .= _to_top_down(native, orientation)
        end
        ConstituentField(values, DryMoleFraction();
                         metadata=(source=:gchp_tomas, component=component))
    end
    components = NamedTuple{component_names}(Tuple(component_fields))
    return SectionalAerosolState(mapping.size_grid,
        ConstituentField(number, NumberDensity();
                         metadata=(source=:gchp_tomas, native_units="1000 particles/mol dry air")),
        components;
        metadata=(source=:gchp_tomas,),
    )
end

function read_column(source::GCHPColumnSource,
                     location::GCHPLocation;
                     gases=(:CO2, :CH4, :N2O, :CO, :O3, :C2H6),
                     FT::Type{<:AbstractFloat}=Float64)
    schema = source.schema
    native_dp, dp_variable = _read_profile(
        source, schema.pressure_thickness, location)
    native_surface_pressure, ps_variable = _read_scalar(
        source, schema.surface_pressure, location)
    dp = FT.(_pressure_to_pa(native_dp, schema.pressure_thickness,
                            _units(dp_variable)))
    surface_pressure = FT(_pressure_to_pa(native_surface_pressure,
                                          schema.surface_pressure,
                                          _units(ps_variable)))

    # Build interfaces from surface toward TOA regardless of source layer order,
    # then return the package's canonical top-to-bottom pressure coordinate.
    bottom_up_dp = schema.stored_orientation isa BottomToTop ? dp : reverse(dp)
    native_interfaces = [surface_pressure;
                         surface_pressure .- cumsum(bottom_up_dp)]
    tolerance = 100 * eps(FT) * max(surface_pressure, one(FT))
    abs(native_interfaces[end]) <= tolerance && (native_interfaces[end] = zero(FT))
    interfaces = reverse(native_interfaces)

    native_temperature, temperature_variable = _read_profile(
        source, schema.temperature, location)
    lowercase(strip(_units(temperature_variable))) in ("k", "kelvin") ||
        throw(ArgumentError("$(schema.temperature) must use K"))
    temperature_values = FT.(_to_top_down(native_temperature,
                                          schema.stored_orientation))

    native_humidity, humidity_variable = _read_profile(
        source, schema.specific_humidity, location)
    humidity_values = FT.(_to_top_down(
        _humidity_to_kgkg(native_humidity, _units(humidity_variable)),
        schema.stored_orientation))

    native_dry_mass, dry_mass_variable = _read_profile(
        source, schema.dry_air_mass, location)
    lowercase(replace(_units(dry_mass_variable), " " => "")) in ("kg", "") ||
        throw(ArgumentError("$(schema.dry_air_mass) must use kg"))
    dry_mass_values = FT.(_to_top_down(native_dry_mass,
                                      schema.stored_orientation))

    native_volume, volume_variable = _read_profile(
        source, schema.layer_volume, location)
    lowercase(replace(_units(volume_variable), " " => "")) in ("m3", "m^3", "") ||
        throw(ArgumentError("$(schema.layer_volume) must use m3"))
    volume_values = FT.(_to_top_down(native_volume,
                                    schema.stored_orientation))

    dry_column_moles = if _dataset_with(source, schema.dry_pressure_thickness) === nothing
        Unavailable()
    else
        native_dry_dp, dry_dp_variable = _read_profile(
            source, schema.dry_pressure_thickness, location)
        dry_dp = FT.(_pressure_to_pa(native_dry_dp,
                                     schema.dry_pressure_thickness,
                                     _units(dry_dp_variable)))
        dry_dp_top_down = _to_top_down(dry_dp, schema.stored_orientation)
        gravity = FT(9.80665)
        air_molar_mass = FT(28.9644e-3)
        dry_dp_top_down ./ (gravity * air_molar_mass)
    end

    gas_fields = _trace_gases(source, location, gases, schema.stored_orientation)
    aerosol_state = _sectional_aerosols(source, location, source.aerosol_schema,
                                        schema.stored_orientation,
                                        dry_mass_values, volume_values)

    return AtmosphericColumn(
        pressure=PressureCoordinate(interfaces; orientation=TopToBottom()),
        temperature=temperature_values,
        specific_humidity=humidity_values,
        dry_air_mass=dry_mass_values,
        dry_air_column_moles=dry_column_moles,
        layer_volume=volume_values,
        trace_gases=gas_fields,
        aerosols=aerosol_state,
        metadata=(source=:gchp, paths=source.paths, location=location),
    )
end

end
