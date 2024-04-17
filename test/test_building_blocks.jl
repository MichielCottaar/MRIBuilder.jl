@testset "test_building_blocks.jl" begin
    @testset "edge_times" begin
        t_axis = 0:0.1:6
        pulse = GenericPulse(t_axis, t_axis .^ 2)
        bb = BuildingBlock([(0, [1., 0., 0.]), (7, [1., 0., 0.])], [(0., pulse)])
        @test all(edge_times(bb) .≈ [0., 6., 7.])
    end
end