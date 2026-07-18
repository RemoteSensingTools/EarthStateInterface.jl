@testset "Optical reduction routes" begin
    partition = LayerPartition([1:2, 3:4], 4)
    optical_depth = [1.0 10.0; 2.0 20.0; 3.0 30.0; 4.0 40.0]
    @test merge_optical_depth(optical_depth, partition) == [3.0 30.0; 7.0 70.0]
    @test merge_optical_depth([1.0, 2, 3, 4], partition) == [3.0, 7.0]

    amount = [1.0, 2, 3, 4]
    cross_section = [2.0, 3, 5, 7]
    effective = effective_cross_section(cross_section, amount, partition)
    @test effective ≈ [8 / 3, 43 / 7]
    @test effective .* [sum(amount[1:2]), sum(amount[3:4])] ≈
          merge_optical_depth(cross_section .* amount, partition)

    spectral_cross_section = [2.0 4.0; 3.0 6.0; 5.0 10.0; 7.0 14.0]
    spectral_effective = effective_cross_section(
        spectral_cross_section, amount, partition)
    @test spectral_effective[:, 2] ≈ 2 .* spectral_effective[:, 1]
    @test_throws ArgumentError effective_cross_section(
        cross_section, zeros(4), partition)
end
