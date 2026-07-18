using EarthStateInterface

"""
    from_vsmartmom_profile(profile; pressure_scale=100, include_h2o=false)

Translate the current vSmartMOM `AtmosphericProfile` into the shared contract without
making EarthStateInterface depend on vSmartMOM. The defaults reflect vSmartMOM's current
hPa pressure and molecules cm⁻² vertical-column conventions.

This is a migration adapter, not the desired long-term vSmartMOM API.
"""
function from_vsmartmom_profile(profile;
                                pressure_scale=100,
                                include_h2o::Bool=false)
    gas_keys = sort!(collect(keys(profile.vmr)))
    include_h2o || filter!(key -> uppercase(String(key)) != "H2O", gas_keys)
    gas_names = Tuple(Symbol.(gas_keys))
    gas_fields = map(gas_keys) do key
        values = profile.vmr[key]
        layer_values = values isa Real ? fill(values, length(profile.T)) : values
        ConstituentField(layer_values, DryMoleFraction();
                         metadata=(source=:vsmartmom, source_name=String(key)))
    end
    gases = NamedTuple{gas_names}(Tuple(gas_fields))

    avogadro = 6.02214076e23
    dry_column_moles = profile.vcd_dry .* (1e4 / avogadro)
    pressure = PressureCoordinate(
        profile.p_half .* pressure_scale;
        centers=profile.p_full .* pressure_scale,
        orientation=TopToBottom(),
    )

    return AtmosphericColumn(
        pressure=pressure,
        temperature=profile.T,
        specific_humidity=profile.q,
        dry_air_column_moles=dry_column_moles,
        trace_gases=gases,
        metadata=(source=:vsmartmom_atmospheric_profile,),
    )
end
