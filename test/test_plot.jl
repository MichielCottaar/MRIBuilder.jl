using CairoMakie
using VisualRegressionTests
using Gtk
@testset "test_plot.jl" begin
    dir = @__DIR__
    isCI = get(ENV, "CI", "false") == "true"
    @testset "pulseseq example PRESS sequence" begin
        function plot_press(fname)
            fn = joinpath(dir, "example_pulseseq", "01_from_FID_to_PRESS_v140", "06_PRESS_center.seq")
            sequence = read_sequence(fn)
            f = Figure()
            Axis(f[1, 1])
            plot!(sequence)
            CairoMakie.save(fname, f)
        end

        @visualtest plot_press "$dir/plots/pulseq_press.png" !isCI
    end

    @testset "Instantaneous gradients & pulses" begin
        function plot_perfect_dwi(fname)
            sequence = DWI(TE=80., bval=2., Δ=40., gradient=(type=:instant, ))
            f = plot_sequence(sequence)
            CairoMakie.save(fname, f)
        end
    
        @visualtest plot_perfect_dwi "$dir/plots/perfect_dwi.png" !isCI
    end
    
    @testset "Finite RF pulses" begin
        function plot_finite_rf(fname)
            times = -2:0.1:2
            raw_amp = sinc.(times)
            sequence = build_sequence() do 
                Sequence(
                    GenericPulse(times .+ 2, raw_amp),
                    InstantPulse(flip_angle=60., phase=0.),
                    GenericPulse(times .+ 2, raw_amp .* 2),
                    InstantPulse(flip_angle=120., phase=0.),
                ) 
            end
            f = Figure()
            Axis(f[1, 1])
            plot_sequence(sequence)
            CairoMakie.save(fname, f)
        end
    
        @visualtest plot_finite_rf "$dir/plots/finite_rf.png" !isCI
    end
    
    @testset "Finite gradients" begin
        function plot_finite_dwi(fname, single_gradient=false)
            sequence = DWI(bval=2., TE=80, TR=100, orientation=[0, -1, 1], scanner=Siemens_Prisma)
            f = Figure()
            Axis(f[1, 1])
            plot!(sequence; single_gradient=single_gradient)
            CairoMakie.save(fname, f)
        end
    
        @visualtest fn->plot_finite_dwi(fn, true) "$dir/single_grad_dwi.png" !isCI
        @visualtest fn->plot_finite_dwi(fn, false) "$dir/multi_grad_dwi.png" !isCI
    end
    
        

end