module MakieMRIBuilder
using Makie
import MakieCore
import MRIBuilder.Plot: Plot_Sequence, SequenceDiagram, range_line
import MRIBuilder: BaseSequence, BaseBuildingBlock

function Makie.plot!(scene:: Plot_Sequence)
    kwargs = Dict([
        key => scene[key] for key in [
            :visible, :overdraw, :transparency, :fxaa, :inspectable, :depth_shift, :model, :space
        ]
    ])
    text_kwargs = Dict([
        key => scene[key] for key in [
            :textcolor, :font, :fonts, :fontsize
        ]
    ])
    text_kwargs[:color] = map((a, c) -> a === MakieCore.automatic ? c : a, scene[:textcolor], scene[:color])
    line_color = map((a, c) -> a === MakieCore.automatic ? c : a, scene[:linecolor], scene[:color])
    instant_width = map((a, c) -> a * c, scene[:linewidth], scene[:instant_width])

    lift(scene[:sequence]) do sequence
        sequence_diagram = SequenceDiagram(sequence)
        current_y = 0.
        for label in (:ADC, :Gz, :Gy, :Gx, :G, :RFy, :RFx)
            line = getproperty(sequence_diagram, label)
            full_range = range_line(line)
            norm = max(abs.(full_range)...)
            (lower, upper) = full_range ./ norm
            if norm > 0
                shift = current_y - lower
                Makie.text!(scene, string(label) * " "; position=(0., shift), align=(:right, :center), text_kwargs..., kwargs...)
                Makie.lines!(scene, line.times, line.amplitudes ./ norm .+ shift; color=line_color, linewidth=scene[:linewidth], kwargs...)
                for (time, amplitude) in zip(line.event_times, line.event_amplitudes)
                    Makie.lines!([time, time], [0., amplitude / norm] .+ shift; color=line_color, linewidth=instant_width, kwargs...)
                end
                current_y += (upper - lower) + 0.1
            end
        end
    end
    
end

Makie.plottype(::Union{BaseBuildingBlock, BaseSequence}) = Plot_Sequence

end