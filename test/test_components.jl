@testset "test_components.jl" begin
    @testset "probing GenericPulse" begin
        gp = GenericPulse([0., 1., 3., 4.], [1., 2., 3., 0.], [0., 90., 180., 180.])
        @test amplitude(gp, 0.) ≈ 1.
        @test amplitude(gp, 1.) ≈ 2.
        @test amplitude(gp, 0.25) ≈ 1.25
        @test amplitude(gp, 2.5) ≈ 2.75
        @test amplitude(gp, 3.5) ≈ 1.5
        @test amplitude(gp, 4.) ≈ 0.

        @test phase(gp, 0.) ≈ 0.
        @test phase(gp, 1.) ≈ 90.
        @test phase(gp, 0.25) ≈ 90/4
        @test phase(gp, 2.5) ≈ 90 * 7/4
        @test phase(gp, 3.5) ≈ 180.
        @test phase(gp, 4.) ≈ 180.

        @test frequency(gp, 0.) ≈ 1/4
        @test frequency(gp, 1.) ≈ 1/6
        @test frequency(gp, 0.25) ≈ 1/4
        @test frequency(gp, 2.5) ≈ 1/8
        @test frequency(gp, 3.) ≈ 1/12
        @test frequency(gp, 3.5) ≈ 0.
        @test frequency(gp, 4.) ≈ 0.
    end

end
