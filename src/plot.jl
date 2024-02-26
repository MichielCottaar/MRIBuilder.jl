module Plot
import MakieCore: generic_plot_attributes!
import ..Containers: BaseBuildingBlock, BaseSequence, waveform, events, start_time, ndim_grad, waveform_sequence
import ..Variables: duration, flip_angle, phase, make_generic, gradient_orientation
import ..Components: RFPulseComponent, ADC, InstantPulse, NoGradient

"""
    SinglePlotLine(times, amplitudes, event_times, event_amplitudes)

A single line in a sequence diagram (e.g., RFx, Gy, ADC).
"""
struct SinglePlotLine
    times :: Vector{Float64}
    amplitudes :: Vector{Float64}
    event_times :: Vector{Float64}
    event_amplitudes :: Vector{Float64}
end

duration(pl::SinglePlotLine) = iszero(length(pl.times)) ? 0. : pl.times[end]
SinglePlotLine(times::Vector{Float64}, amplitudes::Vector{Float64}) = SinglePlotLine(times, amplitudes, Float64[], Float64[])

function SinglePlotLine(control_points::AbstractVector{<:Tuple}, duration::Number)
    times = [cp[1] for cp in control_points]
    amplitudes = [cp[2] for cp in control_points]
    if times[1] > 0
        pushfirst!(times, 0.)
        pushfirst!(amplitudes, 0.)
    end
    if times[end] < duration
        push!(times, duration)
        push!(amplitudes, 0.)
    end
    return SinglePlotLine(times, amplitudes, [], [])
end

SinglePlotLine(duration::Number) = SinglePlotLine([0., duration], [0., 0.])
function SinglePlotLine(plot_lines::SinglePlotLine...)
    times = Float64[]
    amplitudes = Float64[]
    event_times = Float64[]
    event_amplitudes = Float64[]
    current_time = 0.
    for pl in plot_lines
        append!(times, pl.times .+ current_time)
        append!(amplitudes, pl.amplitudes)
        append!(event_times, pl.event_times .+ current_time)
        append!(event_amplitudes, pl.event_amplitudes)
        current_time += duration(pl)
    end
    return SinglePlotLine(times, amplitudes, event_times, event_amplitudes)
end

"""
    range_line(single_plot_line)

Returns the minimum and maximum amplitude for a [`SinglePlotLine`](@ref)
"""
range_line(spl::SinglePlotLine) = (
    min(minimum(spl.amplitudes; init=0.), minimum(spl.event_amplitudes; init=0.)),
    max(maximum(spl.amplitudes; init=0.), maximum(spl.event_amplitudes; init=0.)),
)

normalise(spl::SinglePlotLine, value) = iszero(value) ? spl : SinglePlotLine(spl.times, spl.amplitudes ./ value, spl.event_times, spl.event_amplitudes ./ value)

"""
    SequenceDiagram(; RFx, RFy, Gx, Gy, Gz, ADC)

All the lines forming a sequence diagram.

Each parameter should be a [`SinglePlotLine`](@ref) if provided.
Any parameters not provided will be set to a [`SinglePlotLine`](@ref) with zero amplitude.
"""
struct SequenceDiagram
    RFx :: SinglePlotLine
    RFy :: SinglePlotLine
    G :: SinglePlotLine
    Gx :: SinglePlotLine
    Gy :: SinglePlotLine
    Gz :: SinglePlotLine
    ADC :: SinglePlotLine
end

function SequenceDiagram(actual_duration; kwargs...)
    durations = duration.([values(kwargs)...])
    @assert all(isapprox.(actual_duration, durations, rtol=1e-3))
    res = SinglePlotLine[]
    for symbol in (:RFx, :RFy, :G, :Gx, :Gy, :Gz, :ADC)
        push!(res, get(kwargs, symbol, SinglePlotLine(actual_duration)))
    end
    return SequenceDiagram(res...)
end


function SequenceDiagram(bbb::BaseBuildingBlock)
    kwargs = Dict{Symbol, SinglePlotLine}()

    if !all([wvs isa NoGradient for (_, wvs) in waveform_sequence(bbb)])
        orientation = ndim_grad(bbb) == 3 ? 1. : gradient_orientation(bbb)
        for (index, symbol) in enumerate((:Gx, :Gy, :Gz))
            kwargs[symbol] = SinglePlotLine([(time, (orientation .* amplitude)[index]) for (time, amplitude) in waveform(bbb)], duration(bbb))
        end
    end

    for (name, raw_event) in events(bbb)
        delay = start_time(bbb, name)
        event = make_generic(raw_event)
        if event isa RFPulseComponent
            if :RFx in keys(kwargs)
                error("Cannot plot a building block with more than 1 RF pulse.")
            end
            for (symbol, fn) in [(:RFx, cosd), (:RFy, sind)]
                if event isa InstantPulse
                    kwargs[symbol] = SinglePlotLine([0., duration(bbb)], [0., 0.], [delay], [flip_angle(event) * fn(phase(event))])
                else
                    points = [(t + delay, a * fn(p)) for (t, a, p) in zip(event.time, event.amplitude, event.phase)]
                    kwargs[symbol] = SinglePlotLine(points, duration(bbb))
                end
            end
        elseif event isa ADC
            if :ADC in keys(kwargs)
                error("Cannot plot a building block with more than 1 ADC event.")
            end
            if iszero(duration(event))
                kwargs[:ADC] = SinglePlotLine(
                    [0., duration(bbb)],
                    [0., 0.],
                    [delay], [1.]
                )
            else
                kwargs[:ADC] = SinglePlotLine(
                    [0., delay, delay + duration(event), duration(bbb)],
                    [0., 1., 1., 0.],
                )
            end
        end
    end
    return SequenceDiagram(duration(bbb); kwargs...)
end

function SequenceDiagram(seq::BaseSequence)
    as_lines = SequenceDiagram.(seq)
    return SequenceDiagram([
        SinglePlotLine([getproperty(line, symbol) for line in as_lines]...) 
        for symbol in [:RFx, :RFy, :G, :Gx, :Gy, :Gz, :ADC]]...)
end

function normalise(sd::SequenceDiagram)
    rf_norm = max(abs.(range_line(sd.RFx))..., abs.(range_line(sd.RFy))...)
    grad_norm = max(
        abs.(range_line(sd.G))...,
        abs.(range_line(sd.Gx))...,
        abs.(range_line(sd.Gy))...,
        abs.(range_line(sd.Gz))...,
    )
    SequenceDiagram(
        normalise(sd.RFx, rf_norm),
        normalise(sd.RFy, rf_norm),
        normalise(sd.G, grad_norm),
        normalise(sd.Gx, grad_norm),
        normalise(sd.Gy, grad_norm),
        normalise(sd.Gz, grad_norm),
        sd.ADC
    )
end

"""
    plot_sequence(sequence; figure=(), axis=(), attributes...)
    plot(sequence; attributes...)
    plot!([scene,] sequence; attributes...)

Plot the sequence diagram.

Calling `plot_sequence` will result in a much cleaner sequence diagram (recommended).
However, if you want to combine this diagram with other plots you will have to use `plot` or `plot!` instead.

If called as `plot_sequence` the user can also supply `Makie.Figure` (`figure=(...)`) and `Makie.Axis` (`axis=(...)`) keywords.
If called using the `plot` or `plot!` interface, only the attributes listed below can be supplied

This function will only work if [`Makie`](https://makie.org) is installed and imported.

## Attributes

### Line properties
- `linecolor` sets the color of the lines. If you want to set the text color to the same value, you can also use `color=...`.
- `linewidth=1.5` sets the width of the lines.
- `instant_width=3.` sets the width of any instant gradients or pulses with respect to the `linewidth`.

### Text properties
- `textcolor` sets the color of the text. If you want to set the line color to the same value, you can also use `color=...`.
- `font` sets whether the rendered text is :regular, :bold, or :italic.
- `fontsize`: set the size of each character.

$(Base.Docs.doc(generic_plot_attributes!))
"""
function plot_sequence end

end