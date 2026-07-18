@testset "Conservative state reduction" begin
    pressure = PressureCoordinate([0.0, 10_000, 30_000, 60_000, 100_000])
    size_grid = SectionalSizeGrid([1e-8, 1e-7, 1e-6], ParticleDiameter())
    aerosol = SectionalAerosolState(
        size_grid,
        ConstituentField([1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0], NumberDensity()),
        (
            sulfate=ConstituentField(
                [0.1 0.2 0.3 0.4; 0.2 0.3 0.4 0.5], DryMoleFraction()),
        ),
    )
    column = AtmosphericColumn(
        pressure=pressure,
        temperature=[220.0, 240, 270, 290],
        specific_humidity=[0.0, 0.001, 0.004, 0.01],
        dry_air_mass=[1.0, 2, 3, 4],
        dry_air_column_moles=[10.0, 20, 30, 40],
        layer_volume=[10.0, 9, 8, 7],
        trace_gases=(
            CO2=ConstituentField([400.0, 402, 404, 406] .* 1e-6,
                                  DryMoleFraction()),
            tagged_mass=ConstituentField([1.0, 2, 3, 4], LayerMass()),
        ),
        aerosols=aerosol,
    )

    partition = LayerPartition([1:2, 3:4], 4)
    merged = merge_column(column, partition)

    @test pressure_interfaces(merged) == [0.0, 30_000, 100_000]
    @test dry_air_mass(merged) == [3.0, 7.0]
    @test dry_air_column_moles(merged) == [30.0, 70.0]
    @test layer_volume(merged) == [19.0, 15.0]
    @test temperature(merged) ≈ [(220 + 2 * 240) / 3, (3 * 270 + 4 * 290) / 7]
    @test constituent_values(constituent(merged, :CO2)) ≈
          [(400 + 2 * 402) / 3, (3 * 404 + 4 * 406) / 7] .* 1e-6
    @test constituent_values(constituent(merged, :tagged_mass)) == [3.0, 7.0]

    merged_number = constituent_values(aerosol_number(aerosols(merged)))
    @test merged_number[:, 1] ≈ [28 / 19, 104 / 19]
    @test merged_number[:, 2] ≈ [52 / 15, 112 / 15]

    by_count = partition_by_count(7, 3)
    @test by_count.groups == [1:2, 3:4, 5:7]
    by_pressure = partition_at_interfaces(pressure, [0.0, 30_000, 100_000])
    @test by_pressure.groups == [1:2, 3:4]

    column_amount_only = AtmosphericColumn(
        pressure=pressure,
        temperature=[220.0, 240, 270, 290],
        specific_humidity=[0.0, 0.001, 0.004, 0.01],
        dry_air_column_moles=[10.0, 20, 30, 40],
        trace_gases=(CO2=ConstituentField(fill(400e-6, 4), DryMoleFraction()),),
    )
    merged_amount_only = merge_column(column_amount_only, partition)
    @test dry_air_mass(merged_amount_only) isa Unavailable
    @test dry_air_column_moles(merged_amount_only) == [30.0, 70.0]
end
