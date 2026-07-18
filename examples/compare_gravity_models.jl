using EarthStateInterface

# A compact idealized sea-level profile. Altitudes are geometric layer centers and
# correspond approximately to the pressure layers in top-to-bottom storage order.
column = AtmosphericColumn(
    pressure=PressureCoordinate(
        [100.0, 1_000.0, 5_000.0, 15_000.0, 30_000.0,
         50_000.0, 70_000.0, 85_000.0, 95_000.0, 101_325.0],
    ),
    temperature=[220.0, 225.0, 235.0, 250.0, 265.0, 275.0, 283.0, 288.0, 292.0],
    specific_humidity=[0.0, 0.0, 1e-5, 1e-4, 5e-4, 0.002, 0.005, 0.008, 0.010],
)
altitude = [44_000.0, 32_000.0, 22_000.0, 14_000.0, 8_000.0,
            4_000.0, 2_000.0, 900.0, 250.0]

constant = air_column_amounts(column, ConstantGravity())
reference = sum(dry_air_moles(constant))

function report(label, amounts)
    difference_percent = 100 * (sum(dry_air_moles(amounts)) / reference - 1)
    println(rpad(label, 34), round(difference_percent; digits=4), " %")
end

println("Dry-column change relative to constant g = 9.80665 m s⁻²")
report("constant", constant)
for latitude in (0.0, 45.0, 90.0)
    report("Helmert latitude $(latitude)°",
           air_column_amounts(column, HelmertLatitudeGravity(); latitude))
end
report("altitude only",
       air_column_amounts(column, SphericalAltitudeGravity(); altitude))
for latitude in (0.0, 45.0, 90.0)
    report("latitude + altitude $(latitude)°",
           air_column_amounts(column, SomiglianaAltitudeGravity(); latitude, altitude))
end
