"""Model used to evaluate gravitational acceleration for hydrostatic air columns."""
abstract type GravityModel end

"""
    ConstantGravity([g=9.80665])

Use one gravitational acceleration `g` [m s⁻²] in every layer. This matches the
usual reanalysis convention and deliberately ignores latitude and altitude.
"""
struct ConstantGravity{T <: Real} <: GravityModel
    g::T

    function ConstantGravity(g::T) where {T <: Real}
        isfinite(g) && g > zero(g) || throw(ArgumentError(
            "gravitational acceleration must be finite and positive"))
        return new{T}(g)
    end
end

ConstantGravity() = ConstantGravity(9.80665)

"""
    HelmertLatitudeGravity([reference=9.80665], [amplitude=0.02586])

Latitude-only approximation used by RRTMGP:
`g(φ) = reference - amplitude*cos(2φ)`, where `φ` is geodetic latitude.
Altitude is deliberately ignored.
"""
struct HelmertLatitudeGravity{T <: Real} <: GravityModel
    reference::T
    amplitude::T

    function HelmertLatitudeGravity(reference::T, amplitude::T) where {T <: Real}
        isfinite(reference) && reference > zero(reference) || throw(ArgumentError(
            "reference gravity must be finite and positive"))
        isfinite(amplitude) && zero(amplitude) <= amplitude < reference ||
            throw(ArgumentError("gravity amplitude must satisfy 0 ≤ amplitude < reference"))
        return new{T}(reference, amplitude)
    end
end

HelmertLatitudeGravity() = HelmertLatitudeGravity(9.80665, 0.02586)
function HelmertLatitudeGravity(reference::Real, amplitude::Real)
    values = promote(reference, amplitude)
    return HelmertLatitudeGravity(values...)
end

"""
    SphericalAltitudeGravity([surface_gravity=9.80665], [earth_radius=6_371_000])

Altitude-only inverse-square model,
`g(h) = surface_gravity * (earth_radius / (earth_radius + h))²`.
Latitude is deliberately ignored. Heights are geometric layer-center altitudes [m].
"""
struct SphericalAltitudeGravity{T <: Real} <: GravityModel
    surface_gravity::T
    earth_radius::T

    function SphericalAltitudeGravity(surface_gravity::T,
                                      earth_radius::T) where {T <: Real}
        isfinite(surface_gravity) && surface_gravity > zero(surface_gravity) ||
            throw(ArgumentError("surface gravity must be finite and positive"))
        isfinite(earth_radius) && earth_radius > zero(earth_radius) ||
            throw(ArgumentError("Earth radius must be finite and positive"))
        return new{T}(surface_gravity, earth_radius)
    end
end

SphericalAltitudeGravity() = SphericalAltitudeGravity(9.80665, 6_371_000.0)
function SphericalAltitudeGravity(surface_gravity::Real, earth_radius::Real)
    values = promote(surface_gravity, earth_radius)
    return SphericalAltitudeGravity(values...)
end

"""
    SomiglianaAltitudeGravity()

Latitude-and-altitude gravity model. Surface gravity follows the WGS 84
Somigliana formula. Its altitude dependence is the transparent inverse-square
continuation `g(φ, h) = g(φ, 0) * (a / (a + h))²`.

This is suitable for sensitivity studies and atmospheric retrieval columns, but it is
not a high-order geopotential model. Heights are geometric layer-center altitudes [m]
above the reference ellipsoid, not geopotential height.
"""
struct SomiglianaAltitudeGravity{T <: Real} <: GravityModel
    equatorial_gravity::T
    somigliana_constant::T
    eccentricity_squared::T
    semi_major_axis::T

    function SomiglianaAltitudeGravity(equatorial_gravity::T,
                                       somigliana_constant::T,
                                       eccentricity_squared::T,
                                       semi_major_axis::T) where {T <: Real}
        all(isfinite, (equatorial_gravity, somigliana_constant,
                       eccentricity_squared, semi_major_axis)) ||
            throw(ArgumentError("gravity-model constants must be finite"))
        equatorial_gravity > zero(equatorial_gravity) || throw(ArgumentError(
            "equatorial gravity must be positive"))
        zero(eccentricity_squared) <= eccentricity_squared < one(eccentricity_squared) ||
            throw(ArgumentError("eccentricity squared must satisfy 0 ≤ e² < 1"))
        semi_major_axis > zero(semi_major_axis) || throw(ArgumentError(
            "semi-major axis must be positive"))
        return new{T}(equatorial_gravity, somigliana_constant,
                      eccentricity_squared, semi_major_axis)
    end
end

SomiglianaAltitudeGravity() = SomiglianaAltitudeGravity(
    9.7803253359, 0.00193185265241, 0.00669437999013, 6_378_137.0)
function SomiglianaAltitudeGravity(equatorial_gravity::Real,
                                   somigliana_constant::Real,
                                   eccentricity_squared::Real,
                                   semi_major_axis::Real)
    values = promote(equatorial_gravity, somigliana_constant,
                     eccentricity_squared, semi_major_axis)
    return SomiglianaAltitudeGravity(values...)
end

function _check_latitude(latitude)
    isfinite(latitude) && -90 <= latitude <= 90 || throw(ArgumentError(
        "geodetic latitude must be finite and lie in [-90, 90] degrees"))
    return latitude
end

function _check_altitude(altitude, radius)
    isfinite(altitude) && altitude > -radius || throw(ArgumentError(
        "geometric altitude must be finite and greater than -Earth radius"))
    return altitude
end

"""
    gravity(model; latitude=nothing, altitude=nothing)

Evaluate gravitational acceleration [m s⁻²]. Each concrete model requires only the
coordinates represented by its name; unused coordinates may be supplied so a common call
can compare several models.
"""
gravity(model::ConstantGravity; latitude=nothing, altitude=nothing) = model.g

function gravity(model::HelmertLatitudeGravity; latitude=nothing, altitude=nothing)
    latitude === nothing && throw(ArgumentError(
        "HelmertLatitudeGravity requires geodetic latitude [degrees]"))
    φ = deg2rad(_check_latitude(latitude))
    return model.reference - model.amplitude * cos(2φ)
end

function gravity(model::SphericalAltitudeGravity; latitude=nothing, altitude=nothing)
    altitude === nothing && throw(ArgumentError(
        "SphericalAltitudeGravity requires geometric altitude [m]"))
    h = _check_altitude(altitude, model.earth_radius)
    return model.surface_gravity * (model.earth_radius / (model.earth_radius + h))^2
end


function gravity(model::SomiglianaAltitudeGravity; latitude=nothing, altitude=nothing)
    latitude === nothing && throw(ArgumentError(
        "SomiglianaAltitudeGravity requires geodetic latitude [degrees]"))
    altitude === nothing && throw(ArgumentError(
        "SomiglianaAltitudeGravity requires geometric altitude [m]"))
    φ = deg2rad(_check_latitude(latitude))
    h = _check_altitude(altitude, model.semi_major_axis)
    sin²φ = sin(φ)^2
    surface_gravity = model.equatorial_gravity *
        (1 + model.somigliana_constant * sin²φ) /
        sqrt(1 - model.eccentricity_squared * sin²φ)
    return surface_gravity * (model.semi_major_axis /
                              (model.semi_major_axis + h))^2
end

"""
    AirColumnAmounts

Per-layer hydrostatic amounts. `dry_air_moles`, `water_vapor_moles`, and
`wet_air_moles` are in mol m⁻². Here *wet air* means the complete moist mixture, so
`wet_air_moles == dry_air_moles + water_vapor_moles`. `layer_gravity` is in m s⁻².
"""
struct AirColumnAmounts{D, W, M, G}
    dry_air_moles::D
    water_vapor_moles::W
    wet_air_moles::M
    layer_gravity::G
end

dry_air_moles(amounts::AirColumnAmounts) = amounts.dry_air_moles
water_vapor_moles(amounts::AirColumnAmounts) = amounts.water_vapor_moles
wet_air_moles(amounts::AirColumnAmounts) = amounts.wet_air_moles
layer_gravity(amounts::AirColumnAmounts) = amounts.layer_gravity

function _layer_altitudes(altitude::Real, n)
    return fill(altitude, n)
end

function _layer_altitudes(altitude, n)
    altitude === nothing && return fill(nothing, n)
    length(altitude) == n || throw(DimensionMismatch(
        "altitude must be one geometric layer-center height per atmospheric layer"))
    return altitude
end

function _layer_gravity(column, model::GravityModel, latitude, altitude)
    altitudes = _layer_altitudes(altitude, nlayers(column))
    return [gravity(model; latitude, altitude=h) for h in altitudes]
end

"""
    air_column_amounts(column, gravity_model=ConstantGravity();
                       latitude=nothing, altitude=nothing,
                       dry_air_molar_mass=0.0289644,
                       water_molar_mass=0.01801528)

Compute dry-air, water-vapor, and total wet-air vertical columns in every layer
[mol m⁻²]. Pressure differences are converted to hydrostatic moist-air mass using
`dm = |dp| / g`; specific humidity `q` then gives water mass `q*dm` and dry mass
`(1-q)*dm` without an intermediate mixing-ratio convention.

`latitude` is one geodetic latitude [degrees] for the column. `altitude` may be one
geometric altitude [m] used for all layers or one layer-center altitude per layer. Models
that do not represent latitude or altitude ignore that coordinate. The pressure coordinate
may use either vertical orientation.
"""
function air_column_amounts(column,
                            gravity_model::GravityModel=ConstantGravity();
                            latitude=nothing,
                            altitude=nothing,
                            dry_air_molar_mass=0.0289644,
                            water_molar_mass=0.01801528)
    validate_column(column)
    isfinite(dry_air_molar_mass) && dry_air_molar_mass > 0 || throw(ArgumentError(
        "dry-air molar mass must be positive [kg mol⁻¹]"))
    isfinite(water_molar_mass) && water_molar_mass > 0 || throw(ArgumentError(
        "water molar mass must be positive [kg mol⁻¹]"))

    Δp = abs.(diff(pressure_interfaces(column)))
    g = _layer_gravity(column, gravity_model, latitude, altitude)
    moist_mass = Δp ./ g
    water_mass = specific_humidity(column) .* moist_mass
    dry_mass = moist_mass .- water_mass
    dry_moles = dry_mass ./ dry_air_molar_mass
    water_moles = water_mass ./ water_molar_mass
    return AirColumnAmounts(dry_moles, water_moles,
                            dry_moles .+ water_moles, g)
end
