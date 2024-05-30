@testset "test_components.jl" begin
    @testset "probing GenericPulse" begin
        gp = GenericPulse([0., 1., 3., 4.], [1., 2., 3., 0.], [0., 90., 180., 180.])
        @test variables.amplitude(gp, 0.) ≈ 1.
        @test variables.amplitude(gp, 1.) ≈ 2.
        @test variables.amplitude(gp, 0.25) ≈ 1.25
        @test variables.amplitude(gp, 2.5) ≈ 2.75
        @test variables.amplitude(gp, 3.5) ≈ 1.5
        @test variables.amplitude(gp, 4.) ≈ 0.

        @test variables.phase(gp, 0.) ≈ 0.
        @test variables.phase(gp, 1.) ≈ 90.
        @test variables.phase(gp, 0.25) ≈ 90/4
        @test variables.phase(gp, 2.5) ≈ 90 * 7/4
        @test variables.phase(gp, 3.5) ≈ 180.
        @test variables.phase(gp, 4.) ≈ 180.

        @test variables.frequency(gp, 0.) ≈ 1/4
        @test variables.frequency(gp, 1.) ≈ 1/6
        @test variables.frequency(gp, 0.25) ≈ 1/4
        @test variables.frequency(gp, 2.5) ≈ 1/8
        @test variables.frequency(gp, 3.) ≈ 1/12
        @test variables.frequency(gp, 3.5) ≈ 0.
        @test variables.frequency(gp, 4.) ≈ 0.
    end
    @testset "SincPulse" begin
        sp = SincPulse(amplitude=1., frequency=2., phase=0., lobe_duration=10., Nzeros=(2, 1))
        
        @test variables.duration(sp) ≈ 30.

        @test variables.amplitude(sp, 20) ≈ 1.
        @test variables.amplitude(sp, 0) ≈ 0. atol=1e-8
        @test variables.amplitude(sp, 10) ≈ 0. atol=1e-8

        @test variables.phase(sp, 20) ≈ 0. atol=1e-8
        @test variables.phase(sp, 0) ≈ -360 * 40
        @test variables.phase(sp, 10) ≈ -360 * 20
        @test variables.phase(sp, 30) ≈ 360 * 20

        @test variables.frequency(sp, π) ≈ 2.
    end

end
