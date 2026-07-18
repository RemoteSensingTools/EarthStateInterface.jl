@testset "Gravity-aware air columns" begin
    pressure = PressureCoordinate([0.0, 20_000.0, 60_000.0, 100_000.0])
    column = AtmosphericColumn(
        pressure=pressure,
        temperature=[220.0, 260.0, 290.0],
        specific_humidity=[0.0, 0.005, 0.01],
    )

    constant = air_column_amounts(column)
    expected_mass = diff(pressure_interfaces(column)) ./ 9.80665
    @test dry_air_moles(constant) ≈
          (1 .- specific_humidity(column)) .* expected_mass ./ 0.0289644
    @test water_vapor_moles(constant) ≈
          specific_humidity(column) .* expected_mass ./ 0.01801528
    @test wet_air_moles(constant) ≈
          dry_air_moles(constant) .+ water_vapor_moles(constant)
    @test layer_gravity(constant) == fill(9.80665, 3)

    helmert_equator = gravity(HelmertLatitudeGravity(); latitude=0.0)
    helmert_pole = gravity(HelmertLatitudeGravity(); latitude=90.0)
    @test helmert_equator ≈ 9.80665 - 0.02586
    @test helmert_pole ≈ 9.80665 + 0.02586
    @test HelmertLatitudeGravity(9.8f0, 0.02).reference isa Float64
    @test_throws ArgumentError gravity(HelmertLatitudeGravity())
    @test_throws ArgumentError gravity(HelmertLatitudeGravity(); latitude=91)

    spherical = SphericalAltitudeGravity()
    @test gravity(spherical; altitude=0.0) ≈ 9.80665
    @test gravity(spherical; altitude=5_000.0) < gravity(spherical; altitude=0.0)
    @test_throws ArgumentError gravity(spherical)

    combined = SomiglianaAltitudeGravity()
    @test gravity(combined; latitude=0.0, altitude=0.0) ≈ 9.7803253359
    @test gravity(combined; latitude=90.0, altitude=0.0) ≈ 9.8321849379 rtol=1e-9
    @test gravity(combined; latitude=45.0, altitude=3_000.0) <
          gravity(combined; latitude=45.0, altitude=0.0)
    @test_throws ArgumentError gravity(combined; latitude=45.0)

    varying = air_column_amounts(
        column, combined;
        latitude=35.0,
        altitude=[20_000.0, 8_000.0, 1_000.0],
    )
    @test length(layer_gravity(varying)) == nlayers(column)
    @test layer_gravity(varying)[1] < layer_gravity(varying)[end]
    @test sum(dry_air_moles(varying)) != sum(dry_air_moles(constant))
    @test_throws DimensionMismatch air_column_amounts(
        column, combined; latitude=35.0, altitude=[1_000.0, 0.0])

    bottom_up = AtmosphericColumn(
        pressure=PressureCoordinate([100_000.0, 60_000.0, 20_000.0, 0.0];
                                    orientation=BottomToTop()),
        temperature=reverse(temperature(column)),
        specific_humidity=reverse(specific_humidity(column)),
    )
    @test reverse(dry_air_moles(air_column_amounts(bottom_up))) ≈
          dry_air_moles(constant)
end
