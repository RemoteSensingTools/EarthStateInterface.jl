using NCDatasets

function write_gchp_fixture(path)
    NCDataset(path, "c") do dataset
        defDim(dataset, "Xdim", 1)
        defDim(dataset, "Ydim", 1)
        defDim(dataset, "nf", 1)
        defDim(dataset, "lev", 3)
        defDim(dataset, "time", 1)

        layer_dims = ("Xdim", "Ydim", "nf", "lev", "time")
        surface_dims = ("Xdim", "Ydim", "nf", "time")

        function layer_variable(name, values, units)
            variable = defVar(dataset, name, Float64, layer_dims)
            variable.attrib["units"] = units
            variable[:, :, :, :, :] = reshape(values, 1, 1, 1, 3, 1)
        end

        surface_pressure = defVar(dataset, "Met_PS2WET", Float64, surface_dims)
        surface_pressure.attrib["units"] = "hPa"
        surface_pressure[:, :, :, :] .= 1000

        layer_variable("Met_DELP", [400, 350, 250], "hPa")
        layer_variable("Met_DELPDRY", [395, 348, 250], "hPa")
        layer_variable("Met_T", [290, 260, 220], "K")
        layer_variable("Met_SPHU", [10, 2, 0], "g kg-1")
        layer_variable("Met_AD", [4, 3, 2], "kg")
        layer_variable("Met_AIRVOL", [8, 9, 10], "m3")
        layer_variable("SpeciesConcVV_CO2", [420e-6, 410e-6, 400e-6],
                       "mol mol-1 dry")

        for bin in 1:2
            suffix = lpad(bin, 2, '0')
            layer_variable("SpeciesConcVV_NK$suffix", fill(1000.0 * bin, 3), "1")
            layer_variable("SpeciesConcVV_SF$suffix", fill(1e-9 * bin, 3),
                           "mol mol-1 dry")
        end
    end
    return path
end

@testset "Optional GCHP NetCDF reader" begin
    mktempdir() do directory
        path = write_gchp_fixture(joinpath(directory, "gchp.nc4"))
        size_grid = SectionalSizeGrid([1e-8, 1e-7, 1e-6], ParticleDiameter())
        aerosol_schema = GCHPSectionalAerosolSchema(
            size_grid; components=(:SF,))

        column = open_gchp(path; aerosol_schema=aerosol_schema) do source
            read_column(source, GCHPLocation(x=1, y=1, face=1); gases=(:CO2,))
        end

        @test pressure_interfaces(column) ≈ [0, 25_000, 60_000, 100_000]
        @test temperature(column) == [220, 260, 290]
        @test specific_humidity(column) == [0, 0.002, 0.01]
        @test dry_air_mass(column) == [2, 3, 4]
        @test dry_air_column_moles(column) ≈
              reverse([395, 348, 250]) .* 100 ./ (9.80665 * 28.9644e-3)
        @test constituent_values(constituent(column, :CO2)) ==
              [400e-6, 410e-6, 420e-6]
        @test aerosol_number(aerosols(column)).basis isa NumberDensity
        @test size(constituent_values(aerosol_number(aerosols(column)))) == (2, 3)
        @test validate_column(column) === column
    end
end
