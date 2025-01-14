@testset "test_IO.jl" begin
    directory = joinpath(pwd(), "example_pulseseq")
    @testset "read v1.3.1 fiddisp.seq file" begin
        seq = read_sequence(joinpath(directory, "fiddisp_v131.seq"))

        @test length(seq) == 3

        # starting with RF pulse
        @test length(events(seq[1])) == 1
        (index, pulse) = events(seq[1])[1]

        @test start_time(seq, 1, index) == 0.1
        @test end_time(seq, 1, index) == 0.22

        for (time, ampl, phase_check) in [
            (0.09, 0., isnan),
            (0.11, 2.5, iszero),
            (0.19, 2.5, iszero),
            (0.21, 0., iszero),
            (0.23, 0., isnan),
        ]
            @test phase_check(variables.phase(seq, time))
            @test variables.amplitude(seq, time) ≈ ampl
        end

        # Single ADC event
        @test length(variables.readout_times(seq)) == 1024
        @test variables.readout_times(seq)[1] ≈ 0.22 + 5 + 0.02 + 0.5 * 0.3125
        @test variables.readout_times(seq)[end] ≈ 0.22 + 5 + 0.02 + 1023.5 * 0.3125

        @test variables.duration(seq) ≈ 0.22 + 5 + 0.02 + 1024 * 0.3125
    end
    @testset "read all v1.4.0 files in 01_from_FID_to_PRESS" begin
        path = joinpath(directory, "01_from_FID_to_PRESS_v140")
        for fn in readdir(path)
            if !endswith(fn, ".seq")
                continue
            end
            full_filename = joinpath(path, fn)
            seq = read_sequence(full_filename)
        end
    end
    @testset "check 01_from_FID_to_PRESS/01_FID.seq" begin
        fn = joinpath(directory, "01_from_FID_to_PRESS_v140", "01_FID.seq")
        seq = read_sequence(fn)
        @test length(seq) == 2
        (index, pulse) = events(seq[1])[1]
        @test variables.flip_angle(pulse) ≈ 90
        @test start_time(seq, 1, index) ≈ 0.1  # determined by RF dead time
        @test end_time(seq, 1, index) ≈ 0.6  # RF dead time + RF duration
        @test variables.phase(seq, 0.3) == 0
        @test variables.amplitude(seq, 0.3) ≈ 0.5

        start_adc = (
            0.6 +  # end of RF pulse
            0.02 +  # RF ringdown time (added by block duration)
            29.730  # Delay until ADC start (to get start at ADC at TE=30)
        )
        @test length(variables.readout_times(seq)) == 8192
        @test variables.readout_times(seq)[1] ≈ start_adc + 0.5 * 0.03125
        @test variables.readout_times(seq)[end] ≈ start_adc + 8191.5 * 0.03125
    end
    @testset "check 01_from_FID_to_PRESS/06_PRESS_center.seq" begin
        fn = joinpath(directory, "01_from_FID_to_PRESS_v140", "06_PRESS_center.seq")
        seq = read_sequence(fn)
        @test length(seq) == 7

        (index, excitation) = events(seq[1])[1]
        #@test flip_angle(excitation) ≈ 90
        @test start_time(seq, 1, index) ≈ 0.1  # determined by RF dead time
        @test end_time(seq, 1, index) ≈ 3.1  # RF dead time + RF duration
        @test variables.phase(seq, 0.3) ≈ 0 + rad2deg(0.5)
        @test variables.phase(seq, 1.6) ≈ 0
        @test length(variables.readout_times(seq)) == 4096

        refocus_pulses = (
            events(seq[3])[1][2],
            events(seq[5])[1][2],
        )
        for p in refocus_pulses
            @test variables.duration(p) ≈ 3  # should have been 4 for refocus pulses, but there is an error in the matlab generation
            @test variables.phase(p, 0.) ≈ 90 rtol=1e-5
            @test variables.phase(p, 0.38) ≈ 90 + rad2deg(0.5) rtol=1e-5
            @test variables.phase(p, 1.5) ≈ 90 rtol=1e-5
        end
    end
    for format in (:pulseq, :serialize)
        @testset "check that format $(format) works for all v1.4.0 files in 01_from_FID_to_PRESS" begin
            path = joinpath(directory, "01_from_FID_to_PRESS_v140")
            for fn in readdir(path)
                @testset "checking $fn" begin
                    if !endswith(fn, ".seq")
                        continue
                    end
                    full_filename = joinpath(path, fn)
                    seq_orig = read_sequence(full_filename)
                    
                    io = IOBuffer()
                    write_sequence(io, seq_orig, format=format)
                    seek(io, 0)
                    seq_json = read_sequence(io, format=format)
                    @test variables.duration(seq_orig) ≈ variables.duration(seq_json)
                    @test length(seq_orig) == length(seq_json)
                    @test all(variables.duration.(seq_orig) .≈ variables.duration.(seq_json))
                    @test iszero(length(iter_instant_gradients(seq_json)))
                    @test iszero(length(iter_instant_pulses(seq_json)))
                    @test all(variables.readout_times(seq_json) .≈ variables.readout_times(seq_orig))
                end
            end
        end
        @testset "check that format $(format) works for some sequences with instant pulses/gradients/readouts" begin
            for seq_orig in [
                DWI(TE=80., bval=2., scanner=Siemens_Connectom),
                DWI(TE=80., bval=2., scanner=Siemens_Connectom, fov=(10, 10, 10), voxel_size=2.),
                DWI(TE=80., bval=2., scanner=Siemens_Terra, gradient=(type=:instant, )),
                SpinEcho(TE=30.),
            ]
                io = IOBuffer()
                write_sequence(io, seq_orig; format=format)
                seek(io, 0)
                seq_json = read_sequence(io; format=format)
                @test isapprox(variables.duration(seq_orig), variables.duration(seq_json), atol=1e-3)
                @test length(iter_blocks(seq_orig)) == length(iter_blocks(seq_json))
                @test all(isapprox.(variables.duration.(seq_orig), variables.duration.(seq_json), atol=1e-3))
                @test length(iter_instant_gradients(seq_json)) == length(iter_instant_gradients(seq_json))
                @test length(iter_instant_pulses(seq_json)) == length(iter_instant_pulses(seq_json))
                @test length(variables.readout_times(seq_orig)) == length(variables.readout_times(seq_json))
                @test all(isapprox.(variables.readout_times(seq_orig), variables.readout_times(seq_json), atol=1e-3))
            end
        end
    end
    @testset "Test Pulseq shape compression" begin
        import MRIBuilder.SequenceIO.PulseqIO.Parsers: compress, uncompress
        import MRIBuilder.SequenceIO.PulseqIO.Types: PulseqShape

        for (amplitudes, result) in [
            ([0., 0.1, 0.25, 0.5, ones(7)..., 0.75, 0.5, 0.25, 0.], [0., 0.1, 0.15, 0.25, 0.5, 0., 0., 4, -0.25, -0.25, 2]),
            (zeros(100), [0., 0., 98]),
            (ones(100), [1., 0., 0., 97]),
        ]
            shape = PulseqShape(amplitudes)
            res = compress(shape)
            @test res.num == length(amplitudes)
            @test all(res.samples == result)
            @test all(uncompress(res).samples == amplitudes)
        end
    end
end