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

end