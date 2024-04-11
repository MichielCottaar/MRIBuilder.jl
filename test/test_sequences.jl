@testset "test_sequences.jl" begin
    
    @testset "DW-SE" begin
        @testset "Instant pulse & readout" begin
            @testset "Minimise TE" begin
                seq = DiffusionSpinEcho(TE=:min, bval=1.)
                @test length(seq) == 9
                grad_duration = TE(seq) / 2
                @test all(isapprox.(duration.(seq), [0., 0., grad_duration, 0., 0., 0., grad_duration, 0., 0.], atol=1e-6))
                @test length([iter_instant_pulses(seq)...]) == 2
                @test length([iter_instant_gradients(seq)...]) == 0.
                @test bval(seq) ≈ 1.
                @test 40. < TE(seq) < 50.
                @test TR(seq) ≈ TE(seq)
            end
            @testset "Maximise b-value" begin
                seq = DiffusionSpinEcho(TE=80., bval=:max)
                @test length(seq) == 9
                @test all(isapprox.(duration.(seq), [0., 0., 40., 0., 0., 0., 40., 0., 0.], atol=1e-4, rtol=1e-4))
                @test length([iter_instant_pulses(seq)...]) == 2
                @test length([iter_instant_gradients(seq)...]) == 0.
                @test TE(seq) ≈ 80.
                @test TR(seq) ≈ 80.
                @test 4.8 < bval(seq) < 4.9

                # can also maximise q-value
                seq2 = DiffusionSpinEcho(TE=80., qval=:max)
                @test all(isapprox.(duration.(seq), duration.(seq2), atol=1e-4, rtol=1e-4))
                @test TE(seq) ≈ TE(seq2) atol=1e-4 rtol=1e-4
                @test bval(seq) ≈ bval(seq2) atol=1e-4 rtol=1e-4
            end
            @testset "Set diffusion time Δ" begin
                seq = DiffusionSpinEcho(TE=80., Δ=70., qval=:max)
                @test all(isapprox.(duration.(seq), [0., 0., 10., 30., 0., 30., 10., 0., 0.], atol=1e-4, rtol=1e-4))
                @test Δ(seq) ≈ 70.
                @test TE(seq) ≈ 80.
                @test TR(seq) ≈ 80.
                @test 0.72 < bval(seq) < 0.73
            end
            @testset "Set gradient duration" begin
                seq = DiffusionSpinEcho(TE=80., gradient=(duration=10., ), bval=:max)
                @test all(isapprox.(duration.(seq), [0., 0., 10., 30., 0., 30., 10., 0., 0.], atol=1e-4, rtol=1e-4))
                @test Δ(seq) ≈ 70.
                @test TE(seq) ≈ 80.
                @test TR(seq) ≈ 80.
                @test 0.72 < bval(seq) < 0.73
            end
        end

    end
end