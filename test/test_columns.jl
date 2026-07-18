@testset "Atmospheric column contract" begin
    pressure = PressureCoordinate(
        Float32[0, 10_000, 30_000, 60_000, 100_000];
        orientation=TopToBottom(),
    )
    gases = (
        CO2=ConstituentField(Float32[400, 401, 402, 403] .* 1f-6,
                             DryMoleFraction()),
        O3=ConstituentField(Float32[3, 2, 1, 0.1] .* 1f-6,
                            DryMoleFraction()),
    )
    column = AtmosphericColumn(
        pressure=pressure,
        temperature=Float32[220, 240, 270, 290],
        specific_humidity=Float32[0, 0.001, 0.004, 0.01],
        dry_air_mass=Float32[1, 2, 3, 4],
        dry_air_column_moles=Float32[10, 20, 30, 40],
        layer_volume=Float32[10, 9, 8, 7],
        trace_gases=gases,
        metadata=(source=:idealized,),
    )

    @test validate_column(column) === column
    @test nlayers(column) == 4
    @test species_names(column) == (:CO2, :O3)
    @test constituent(column, :CO2) === gases.CO2
    @test dry_air_column_moles(column) == Float32[10, 20, 30, 40]
    @test representative_pressure(pressure, ArithmeticMeanPressure()) ==
          Float32[5_000, 20_000, 45_000, 80_000]
    @test representative_pressure(pressure, LogMeanPressure())[1] == 5_000
    @test_throws ArgumentError pressure_centers(column)

    @test_throws ArgumentError PressureCoordinate(Float32[100_000, 10_000, 0];
                                                   orientation=TopToBottom())
    bottom_up = PressureCoordinate(Float32[100_000, 10_000, 0];
                                   orientation=BottomToTop())
    @test nlayers(bottom_up) == 2
end
