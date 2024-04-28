@testset "test_post_hoc.jl" begin
    import Rotations: RotationVec

    @testset "adjust different components" begin
        @testset "finite gradients" begin
            dwi = DWI(bval=1., TE=:min)
            @test bval(dwi) ≈ 1.
            qval_orig = qval(dwi[:gradient])
            @test all(qval3(dwi[:gradient]) .≈ [qval_orig, 0., 0.])

            @testset "scale and change orientation" begin
                new_dwi = adjust(dwi, DWI=(scale=0.5, orientation=[0., 1., 0.]))
                @test bval(new_dwi) ≈ 0.25
                @test all(qval3(new_dwi[:gradient]) .≈ [0., qval_orig/2, 0.])
            end
            @testset "Rotate gradient" begin
                new_dwi = adjust(dwi, gradient=(rotation=RotationVec(0., 0., π/4), ))
                @test bval(new_dwi) ≈ 1.
                @test all(qval3(new_dwi[:gradient]) .≈ [qval_orig/√2, qval_orig/√2, 0.])
            end
        end
        @testset "instant gradients" begin
            dwi = DWI(bval=1., TE=80, Δ=40, gradient=(type=:instant, ))
            @test bval(dwi) ≈ 1.
            qval_orig = qval(dwi[:gradient])
            @test all(qval3(dwi[:gradient]) .≈ [qval_orig, 0., 0.])

            @testset "scale and change orientation" begin
                new_dwi = adjust(dwi, DWI=(scale=0.5, orientation=[0., 1., 0.]))
                @test bval(new_dwi) ≈ 0.25
                @test all(qval3(new_dwi[:gradient]) .≈ [0., qval_orig/2, 0.])

                @testset "multiple adjustments" begin
                    new_dwi = adjust(dwi, DWI=(scale=[0.5, 1.], orientation=[0., 1., 0.]), merge=false)
                    @test length(new_dwi) == 2
                    @test bval(new_dwi[1]) ≈ 0.25
                    @test bval(new_dwi[2]) ≈ 1.
                    @test all(qval3(new_dwi[1][:gradient]) .≈ [0., qval_orig/2, 0.])
                    @test all(qval3(new_dwi[2][:gradient]) .≈ [0., qval_orig, 0.])
                    
                    new_dwi = adjust(dwi, DWI=(scale=[0.5, 1.], orientation=[0., 1., 0.]))
                    @test duration(new_dwi) ≈ 160

                    new_dwi = adjust(dwi, DWI=(scale=[0.5, 1.], orientation=[0., 1., 0.]), merge=(wait_time=10, ))
                    @test duration(new_dwi) ≈ 170
                end
            end
            @testset "Rotate gradient" begin
                new_dwi = adjust(dwi, gradient=(rotation=RotationVec(0., 0., π/4), ))
                @test bval(new_dwi) ≈ 1.
                @test all(qval3(new_dwi[:gradient]) .≈ [qval_orig/√2, qval_orig/√2, 0.])

            end
        end
    end
end
