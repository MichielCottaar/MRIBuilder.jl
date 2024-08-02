@testset "test_post_hoc.jl" begin
    import Rotations: RotationVec

    @testset "adjust different components" begin
        @testset "finite gradients" begin
            dwi = DWI(bval=1., TE=:min)
            @test variables.bval(dwi) ≈ 1.
            qval_orig = variables.qval(dwi[:gradient])
            @test all(variables.qvec(dwi[:gradient]) .≈ [qval_orig, 0., 0.])

            @testset "scale and change orientation" begin
                new_dwi = adjust(dwi, diffusion=(scale=0.5, orientation=[0., 1., 0.]))
                @test variables.bval(new_dwi) ≈ 0.25
                @test all(variables.qvec(new_dwi[:gradient]) .≈ [0., qval_orig/2, 0.])
            end
            @testset "Rotate gradient" begin
                new_dwi = adjust(dwi, gradient=(rotation=RotationVec(0., 0., π/4), ))
                @test variables.bval(new_dwi) ≈ 1.
                @test all(variables.qvec(new_dwi[:gradient]) .≈ [qval_orig/√2, qval_orig/√2, 0.])
            end
        end
        @testset "instant gradients" begin
            dwi = DWI(bval=1., TE=80, Δ=40, gradient=(type=:instant, ))
            @test variables.bval(dwi) ≈ 1.
            qval_orig = variables.qval(dwi[:gradient])
            @test all(variables.qvec(dwi[:gradient]) .≈ [qval_orig, 0., 0.])

            @testset "scale and change orientation" begin
                new_dwi = adjust(dwi, diffusion=(scale=0.5, orientation=[0., 1., 0.]))
                @test variables.bval(new_dwi) ≈ 0.25
                @test all(variables.qvec(new_dwi[:gradient]) .≈ [0., qval_orig/2, 0.])

                @testset "multiple adjustments" begin
                    new_dwi = adjust(dwi, diffusion=(scale=[0.5, 1.], orientation=[0., 1., 0.]), merge=false)
                    @test length(new_dwi) == 2
                    @test variables.bval(new_dwi[1]) ≈ 0.25
                    @test variables.bval(new_dwi[2]) ≈ 1.
                    @test all(variables.qvec(new_dwi[1][:gradient]) .≈ [0., qval_orig/2, 0.])
                    @test all(variables.qvec(new_dwi[2][:gradient]) .≈ [0., qval_orig, 0.])
                    
                    new_dwi = adjust(dwi, diffusion=(scale=[0.5, 1.], orientation=[0., 1., 0.]))
                    @test variables.duration(new_dwi) ≈ 160
                    @test length(new_dwi) == 2

                    new_dwi = adjust(dwi, diffusion=(scale=[0.5, 1.], orientation=[0., 1., 0.]), merge=(wait_time=10, ))
                    @test variables.duration(new_dwi) ≈ 170
                    @test length(new_dwi) == 3
                end
            end
            @testset "Rotate gradient" begin
                new_dwi = adjust(dwi, gradient=(rotation=RotationVec(0., 0., π/4), ))
                @test variables.bval(new_dwi) ≈ 1.
                @test all(variables.qvec(new_dwi[:gradient]) .≈ [qval_orig/√2, qval_orig/√2, 0.])
            end
        end
        @testset "Slice selection" begin
            dwi = DWI(bval=1., slice_thickness=2., refocus=(spoiler=0.5, ))
            @test variables.frequency(dwi) == (excitation=0., refocus=0.)
            new_dwi = adjust(dwi, slice=(shift=2., ))
            @test variables.frequency(new_dwi).excitation ≈ variables.bandwidth(dwi).excitation
            @test variables.frequency(new_dwi).refocus ≈ variables.bandwidth(dwi).refocus
            another_dwi = adjust(new_dwi, slice=(shift=-3., ))
            @test variables.frequency(another_dwi).excitation ≈ -variables.bandwidth(dwi).excitation * 0.5
            @test variables.frequency(another_dwi).refocus ≈ -variables.bandwidth(dwi).refocus * 0.5
        end
    end
end
