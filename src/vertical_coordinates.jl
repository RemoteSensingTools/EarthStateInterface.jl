"""Orientation of layer arrays along their vertical dimension."""
abstract type VerticalOrientation end

"""Layer index 1 is nearest the top of atmosphere."""
struct TopToBottom <: VerticalOrientation end

"""Layer index 1 is nearest the surface."""
struct BottomToTop <: VerticalOrientation end

"""
    PressureCoordinate(interfaces; centers=Unavailable(), orientation=TopToBottom())

Pressure coordinate in SI units. `interfaces` has `N+1` entries for `N` layers.
Optional `centers` contains source-provided representative layer pressures; it is not
silently synthesized because different applications require different definitions.
"""
struct PressureCoordinate{I, C, O <: VerticalOrientation}
    interfaces::I
    centers::C
    orientation::O
end

function PressureCoordinate(interfaces;
                            centers=Unavailable(),
                            orientation::VerticalOrientation=TopToBottom())
    coordinate = PressureCoordinate{typeof(interfaces), typeof(centers), typeof(orientation)}(
        interfaces, centers, orientation)
    _validate_pressure_coordinate(coordinate)
    return coordinate
end

pressure_interfaces(coordinate::PressureCoordinate) = coordinate.interfaces

function pressure_centers(coordinate::PressureCoordinate)
    return _require_available(coordinate.centers, "source-provided pressure centers")
end

vertical_orientation(coordinate::PressureCoordinate) = coordinate.orientation
nlayers(coordinate::PressureCoordinate) = length(coordinate.interfaces) - 1

function _validate_pressure_coordinate(coordinate::PressureCoordinate)
    interfaces = coordinate.interfaces
    length(interfaces) >= 2 || throw(ArgumentError(
        "pressure coordinate needs at least two interfaces"))
    _allfinite(interfaces) || throw(ArgumentError("pressure interfaces must be finite"))
    all(p -> p >= zero(p), interfaces) || throw(ArgumentError(
        "pressure interfaces must be non-negative SI pressures [Pa]"))

    differences = diff(interfaces)
    if coordinate.orientation isa TopToBottom
        all(>(zero(eltype(differences))), differences) || throw(ArgumentError(
            "TopToBottom pressure interfaces must increase strictly toward the surface"))
    else
        all(<(zero(eltype(differences))), differences) || throw(ArgumentError(
            "BottomToTop pressure interfaces must decrease strictly away from the surface"))
    end

    if isavailable(coordinate.centers)
        length(coordinate.centers) == nlayers(coordinate) || throw(DimensionMismatch(
            "pressure centers must contain one value per layer"))
        _allfinite(coordinate.centers) || throw(ArgumentError(
            "pressure centers must be finite"))
        for layer in 1:nlayers(coordinate)
            lo, hi = extrema((interfaces[layer], interfaces[layer + 1]))
            lo <= coordinate.centers[layer] <= hi || throw(ArgumentError(
                "pressure center $layer lies outside its layer interfaces"))
        end
    end
    return coordinate
end

"""Strategy for diagnosing one representative pressure per layer."""
abstract type RepresentativePressure end

"""Arithmetic mean of bounding interface pressures."""
struct ArithmeticMeanPressure <: RepresentativePressure end

"""Logarithmic pressure mean, with an arithmetic fallback at a zero-pressure TOA."""
struct LogMeanPressure <: RepresentativePressure end

function representative_pressure(coordinate::PressureCoordinate,
                                 ::ArithmeticMeanPressure)
    interfaces = pressure_interfaces(coordinate)
    return (interfaces[1:end-1] .+ interfaces[2:end]) ./ 2
end

function representative_pressure(coordinate::PressureCoordinate,
                                 ::LogMeanPressure)
    interfaces = pressure_interfaces(coordinate)
    return map(interfaces[1:end-1], interfaces[2:end]) do first, second
        lo, hi = extrema((first, second))
        lo == zero(lo) ? (lo + hi) / 2 : (hi - lo) / log(hi / lo)
    end
end
