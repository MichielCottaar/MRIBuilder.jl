module MakieMRIBuilder
using Makie
import MRIBuilder.Plot: SequenceDiagram, range_full, normalise, plot_sequence
import MRIBuilder: BaseSequence, BaseBuildingBlock, variables

@recipe Plot_Sequence (sequence,) begin
    color = @inherit textcolor
    linecolor = Makie.automatic
    linewidth = 1.5
    instant_width = 3.
    textcolor = Makie.automatic
    font = @inherit font
    fonts = @inherit fonts
    fontsize = @inherit fontsize
    Makie.mixin_generic_plot_attributes()...
end

argument_names(::Type{<: Plot_Sequence}, N) = (:sequence, )


function Makie.plot!(scene:: Plot_Sequence)

    map!(scene.attributes, [:linewidth, :instant_width], :linewidth_instant) do l, i
        l * i
    end

    map!(scene.attributes, [:textcolor, :color], :textcolor_final) do t, c
        t === Makie.automatic ? c : t
    end

    map!(scene.attributes, [:linecolor, :color], :linecolor_final) do l, c
        l === Makie.automatic ? c : l
    end

    Makie.register_computation!(scene.attributes, [:sequence], [:text, :text_x, :text_y, :times, :amplitudes, :event_times, :event_amplitudes, :fake_x, :fake_y]) do inputs, changed, cached
        sequence = inputs.sequence
        sequence_diagram = normalise(SequenceDiagram(sequence))

        text = String[]
        text_x = Float64[]
        text_y = Float64[]
        times = Float64[]
        amplitudes = Float64[]

        event_times = Float64[]
        event_amplitudes = Float64[]

        current_y = 0.
        for label in (:ADC, :Gz, :Gy, :Gx, :G, :RFy, :RFx)
            line = getproperty(sequence_diagram, label)
            (lower, upper) = range_full(line)
            if lower ≈ upper
                continue
            end
            shift = current_y - lower
            push!(text, string(label))
            push!(text_x, -variables.duration(sequence) / 40)
            push!(text_y, shift)
            append!(times, line.times)
            append!(amplitudes, line.amplitudes .+ shift)
            push!(times, NaN)
            push!(amplitudes, NaN)
            for (time, amplitude) in zip(line.event_times, line.event_amplitudes)
                append!(event_times, [time, time, NaN])
                append!(event_amplitudes, [shift, amplitude + shift, NaN])
            end
            current_y += (upper - lower) + 0.1
        end

        fake_x = [-variables.duration(sequence) / 10., 0.]
        fake_y = [current_y/2, current_y/2.]

        return (text, text_x, text_y, times, amplitudes, event_times, event_amplitudes, fake_x, fake_y)
    end

    Makie.text!(scene, scene.attributes, scene[:text_x], scene[:text_y]; align=(:right, :center), color=scene[:textcolor_final])
    Makie.lines!(scene, scene.attributes, scene[:times], scene[:amplitudes]; color=scene[:linecolor_final])
    Makie.lines!(scene, scene.attributes, scene[:event_times], scene[:event_amplitudes]; linewidth=scene[:linewidth_instant], color=scene[:linecolor_final])

    Makie.lines!(scene, scene[:fake_x], scene[:fake_y], color=(:black, 0.))
end

Makie.plottype(::Union{BaseBuildingBlock, BaseSequence}) = Plot_Sequence

function plot_sequence(sequence; figure=(), axis=(xgridvisible=false, ygridvisible=false), kwargs...)
    f = Figure(; figure...)
    ax = Axis(f[1, 1]; axis...)
    p = plot!(ax, sequence; kwargs...)
    ax.xlabel[] = "Time (ms)"
    hideydecorations!(ax)
    hidespines!(ax, :l, :r, :t)
    return Makie.FigureAxisPlot(f, ax, p)
end


end