@testset "test_sequences.jl" begin
    @testset "GradientEcho" begin
        seq = GradientEcho(TE=40)
        @test length(seq) == 3
        @test duration(seq) == 40
        @test all(isapprox.(duration.(seq), [0., 40., 0.], atol=1e-6))
        @test length(collect(iter_instant_pulses(seq))) == 1
        @test length(collect(iter_instant_gradients(seq))) == 0.
        @test duration_dephase(seq) ≈ 40.
        @test duration_transverse(seq) ≈ 40.
    end
    @testset "SpinEcho" begin
        seq = SpinEcho(TE=40)
        @test length(seq) == 5
        @test duration(seq) == 40
        @test all(isapprox.(duration.(seq), [0., 20., 0., 20., 0.], atol=1e-6))
        @test length(collect(iter_instant_pulses(seq))) == 2
        @test length(collect(iter_instant_gradients(seq))) == 0.
        @test duration_dephase(seq) ≈ 0. atol=1e-4
        @test duration_transverse(seq) ≈ 40.
    end
    
    @testset "DW-SE" begin
        min_rise_time = Default_Scanner.gradient / Default_Scanner.slew_rate
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
                @test rise_time(seq[:gradient]) ≈ min_rise_time rtol=1e-4
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
                @test rise_time(seq[:gradient]) ≈ min_rise_time rtol=1e-4
                @test all(isapprox.(edge_times(seq, tol=1e-3), [0., min_rise_time, 40. - min_rise_time, 40, 40 + min_rise_time, 80 - min_rise_time, 80.], atol=1e-4))

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
                @test readout_times(seq)[1] ≈ TE(seq)
                @test rise_time(seq[:gradient]) ≈ min_rise_time rtol=1e-4
                @test all(isapprox.(edge_times(seq), [0., min_rise_time, 10. - min_rise_time, 10., 40, 70, 70 + min_rise_time, 80 - min_rise_time, 80.], atol=1e-4))
            end
            @testset "Set gradient duration" begin
                seq = DiffusionSpinEcho(TE=80., gradient=(duration=10., ), bval=:max)
                @test all(isapprox.(duration.(seq), [0., 0., 10., 30., 0., 30., 10., 0., 0.], atol=1e-4, rtol=1e-4))
                @test Δ(seq) ≈ 70. rtol=1e-4
                @test TE(seq) ≈ 80.
                @test TR(seq) ≈ 80.
                @test 0.72 < bval(seq) < 0.73
                @test readout_times(seq)[1] ≈ TE(seq)
                @test rise_time(seq[:gradient]) ≈ min_rise_time rtol=1e-4
            end
        end
        @testset "DW-SE with finite RF pulses" begin
            @testset "slice-select DW-SE" begin
                seq = DiffusionSpinEcho(TR=:min, bval=2., slice_thickness=2.)
                @test length(seq) == 9
                @test duration(seq[1]) > 1.
                for index in 1:9
                    if index in [2, 4, 8, 9]
                        @test abs(duration(seq[index])) < 1e-6
                    else
                        @test abs(duration(seq[index])) > 1
                    end
                end
                @test duration(seq[:gradient]) ≈ duration(seq[:gradient2])
                @test bval(seq) ≈ 2.
                @test length(readout_times(seq)) == 1
                @test readout_times(seq)[1] > TE(seq)
            end
            @testset "voxel-wise DW-SE" begin
                seq = DiffusionSpinEcho(TR=:min, bval=2., voxel_size=2., fov=(20, 20, 20))
                @test length(seq) == 9
                @test duration(seq[1]) > 1.
                for index in 1:9
                    if index in [2, 6, 8]
                        @test abs(duration(seq[index])) < 1e-6
                    else
                        @test abs(duration(seq[index])) > 1
                    end
                end
                @test duration(seq[:gradient]) ≈ duration(seq[:gradient2])
                @test bval(seq) ≈ 2.
                @test length(readout_times(seq)) > 50
                @test 67 < TE(seq) < 68
                @test 72 < TR(seq) < 73

                (pulse, t_pulse) = get_pulse(seq, 1.)
                @test 1. - t_pulse ≈ min_rise_time rtol=1e-4
                @test flip_angle(pulse) ≈ 90.
                @test iszero(phase(pulse))
                @test isnothing(get_pulse(seq, 10.))

                (pulse, t_pulse) = get_pulse(seq, 35.)
                @test 1.1 < t_pulse < 1.2
                @test flip_angle(pulse) ≈ 180.
                @test iszero(phase(pulse))
            end
        end

    end
end