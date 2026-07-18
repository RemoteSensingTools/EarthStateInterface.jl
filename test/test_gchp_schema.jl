@testset "GCHP variable schema" begin
    schema = GCHPColumnSchema()
    @test schema.stored_orientation isa BottomToTop
    @test gchp_gas_variable(schema, :CO2) == "SpeciesConcVV_CO2"
    @test gchp_tomas_variable(schema, :NK, 1) == "SpeciesConcVV_NK01"
    @test gchp_tomas_variable(schema, :DUST, 15) == "SpeciesConcVV_DUST15"
end
