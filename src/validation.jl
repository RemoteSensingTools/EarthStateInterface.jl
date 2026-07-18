"""
    validate_column(column)

Validate the atmospheric-column interface and return `column`. This function also acts as
the contract test entry point for third-party column types.
"""
function validate_column(column)
    coordinate = pressure_coordinate(column)
    coordinate isa PressureCoordinate || throw(ArgumentError(
        "pressure_coordinate(column) must return a PressureCoordinate"))
    _validate_pressure_coordinate(coordinate)

    n = nlayers(coordinate)
    length(temperature(column)) == n || throw(DimensionMismatch(
        "temperature must contain one value per layer"))
    length(specific_humidity(column)) == n || throw(DimensionMismatch(
        "specific humidity must contain one value per layer"))
    _allfinite(temperature(column)) || throw(ArgumentError(
        "temperature must be finite"))
    all(T -> isfinite(T) && T > zero(T), temperature(column)) ||
        throw(ArgumentError("temperature must be positive [K]"))
    all(q -> isfinite(q) && zero(q) <= q < one(q), specific_humidity(column)) ||
        throw(ArgumentError("specific humidity must satisfy 0 ≤ q < 1 [kg kg⁻¹]"))

    if isavailable(dry_air_mass(column))
        length(dry_air_mass(column)) == n || throw(DimensionMismatch(
            "dry-air mass must contain one value per layer"))
        all(m -> isfinite(m) && m > zero(m), dry_air_mass(column)) ||
            throw(ArgumentError("dry-air mass must be positive [kg]"))
    end
    if isavailable(dry_air_column_moles(column))
        length(dry_air_column_moles(column)) == n || throw(DimensionMismatch(
            "dry-air column moles must contain one value per layer"))
        all(amount -> isfinite(amount) && amount > zero(amount),
            dry_air_column_moles(column)) || throw(ArgumentError(
            "dry-air column moles must be positive [mol m⁻²]"))
    end
    if isavailable(layer_volume(column))
        length(layer_volume(column)) == n || throw(DimensionMismatch(
            "layer volume must contain one value per layer"))
        all(v -> isfinite(v) && v > zero(v), layer_volume(column)) ||
            throw(ArgumentError("layer volume must be positive [m³]"))
    end

    gases = trace_gases(column)
    gases isa NamedTuple || throw(ArgumentError("trace_gases(column) must be a NamedTuple"))
    for name in keys(gases)
        field = getproperty(gases, name)
        field isa ConstituentField || throw(ArgumentError(
            "trace gas $name must be a ConstituentField"))
        length(constituent_values(field)) == n || throw(DimensionMismatch(
            "trace gas $name must contain one value per layer"))
        _allfinite(constituent_values(field)) || throw(ArgumentError(
            "trace gas $name contains non-finite values"))
    end

    aerosol_state = aerosols(column)
    aerosol_state isa NoAerosols ||
        validate_sectional_aerosols(aerosol_state; expected_layers=n)
    return column
end
