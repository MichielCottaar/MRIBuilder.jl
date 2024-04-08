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
            @test phase_check(phase(seq, time))
            @test amplitude(seq, time) ≈ ampl
        end

        # Single ADC event
        @test length(readout_times(seq)) == 1024
        @test readout_times(seq)[1] ≈ 0.22 + 5 + 0.02 + 0.5 * 0.3125
        @test readout_times(seq)[end] ≈ 0.22 + 5 + 0.02 + 1023.5 * 0.3125

        @test TR(seq) ≈ 0.22 + 5 + 0.02 + 1024 * 0.3125
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
        @test flip_angle(pulse) ≈ 90
        @test start_time(seq, 1, index) ≈ 0.1  # determined by RF dead time
        @test end_time(seq, 1, index) ≈ 0.6  # RF dead time + RF duration
        @test phase(seq, 0.3) == 0
        @test amplitude(seq, 0.3) ≈ 0.5

        start_adc = (
            0.6 +  # end of RF pulse
            0.02 +  # RF ringdown time (added by block duration)
            29.730  # Delay until ADC start (to get start at ADC at TE=30)
        )
        @test length(readout_times(seq)) == 8192
        @test readout_times(seq)[1] ≈ start_adc + 0.5 * 0.03125
        @test readout_times(seq)[end] ≈ start_adc + 8191.5 * 0.03125
    end
    @testset "check 01_from_FID_to_PRESS/06_PRESS_center.seq" begin
        fn = joinpath(directory, "01_from_FID_to_PRESS_v140", "06_PRESS_center.seq")
        seq = read_sequence(fn)
        @test length(seq) == 7

        (index, excitation) = events(seq[1])[1]
        #@test flip_angle(excitation) ≈ 90
        @test start_time(seq, 1, index) ≈ 0.1  # determined by RF dead time
        @test end_time(seq, 1, index) ≈ 3.1  # RF dead time + RF duration
        @test phase(seq, 0.3) ≈ 0 + rad2deg(0.5)
        @test phase(seq, 1.6) ≈ 0
        @test length(readout_times(seq)) == 4096

        refocus_pulses = (
            events(seq[3])[1][2],
            events(seq[5])[1][2],
        )
        for p in refocus_pulses
            @test duration(p) ≈ 3  # should have been 4 for refocus pulses, but there is an error in the matlab generation
            @test phase(p, 0.) ≈ 90 rtol=1e-5
            @test phase(p, 0.38) ≈ 90 + rad2deg(0.5) rtol=1e-5
            @test phase(p, 1.5) ≈ 90 rtol=1e-5
        end
    end
    if false
        # JSON encoding has not been implemented yet
        @testset "check that JSON encoding works for all v1.4.0 files in 01_from_FID_to_PRESS" begin
            path = joinpath(directory, "01_from_FID_to_PRESS_v140")
            for fn in readdir(path)
                @testset "checking $fn" begin
                    if !endswith(fn, ".seq")
                        continue
                    end
                    full_filename = joinpath(path, fn)
                    seq_orig = read_sequence(full_filename)

                    
                    io = IOBuffer()
                    write_sequence(io, seq_orig)
                    s = String(io.data)
                    seq_json = read_sequence(s)
                    @test seq_orig.TR == seq_json.TR
                    @test length(seq_orig.gradients) == length(seq_json.gradients)
                    for (g1, g2) in zip(seq_orig.gradients, seq_json.gradients)
                        @test all(g1.origin .== g2.origin)
                        @test all(g1.shape.times .== g2.shape.times)
                        @test all(g1.shape.amplitudes .== g2.shape.amplitudes)
                    end
                    @test length(seq_orig.pulses) == length(seq_json.pulses)
                    for (p1, p2) in zip(seq_orig.pulses, seq_json.pulses)
                        @test p1.max_amplitude == p2.max_amplitude
                        @test all(p1.amplitude.times .== p2.amplitude.times)
                        @test all(p1.amplitude.amplitudes .== p2.amplitude.amplitudes)
                        @test all(p1.phase.times .== p2.phase.times)
                        @test all(p1.phase.amplitudes .== p2.phase.amplitudes)
                    end
                    @test iszero(length(seq_json.instants))
                    @test length(seq_json.readout_times) > 0
                    @test all(seq_json.readout_times .== seq_orig.readout_times)
                end
            end
        end
        @testset "check that JSON encoding works for some sequences with instant pulses/gradients" begin
            for seq_orig in [
                DWI(TE=80., bval=2., TR=100., scanner=Siemens_Connectom),
                DWI(TE=80., bval=2., gradient_duration=0., TR=100, scanner=Siemens_Terra),
                SpinEcho(TE=30., TR=100., scanner=Scanner(B0=1.5)),
            ]
                io = IOBuffer()
                write_sequence(io, seq_orig)
                s = String(io.data)
                seq_json = read_sequence(s)
                @test seq_json.TR == 100.
                @test seq_json.scanner.B0 == seq_orig.scanner.B0
                @test seq_json.scanner.gradient == seq_orig.scanner.gradient
                @test seq_json.scanner.slew_rate == seq_orig.scanner.slew_rate
                @test length(seq_orig.gradients) == length(seq_json.gradients)
                for (g1, g2) in zip(seq_orig.gradients, seq_json.gradients)
                    @test all(g1.origin .== g2.origin)
                    @test all(g1.shape.times .== g2.shape.times)
                    @test all(g1.shape.amplitudes .== g2.shape.amplitudes)
                end
                @test iszero(length(seq_json.pulses))
                @test length(seq_json.instants) == length(seq_orig.instants)
                for (i1, i2) in zip(seq_orig.instants, seq_json.instants)
                    @test i1.time == i2.time
                    if i1 isa mr.InstantRFPulse
                        @test i2 isa mr.InstantRFPulse
                        @test i1.flip_angle == i2.flip_angle
                        @test i1.phase == i2.phase
                    else
                        @test i1 isa mr.InstantGradient
                        @test i2 isa mr.InstantGradient
                        @test all(i1.origin .== i2.origin)
                        @test all(i1.qvec .== i2.qvec)
                    end
                end
                @test length(seq_json.readout_times) > 0
                @test all(seq_json.readout_times .== seq_orig.readout_times)
            end
        end
    end
end