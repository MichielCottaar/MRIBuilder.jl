module Plot
import MakieCore: @recipe, theme, generic_plot_attributes!, Attributes, automatic, shading_attributes!, colormap_attributes!
import ..Containers: BaseBuildingBlock, BaseSequence, waveform, events, start_time, ndim_grad
import ..Variables: duration, flip_angle, phase, make_generic
import ..Components: RFPulseComponent, ADC, InstantPulse

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
        push!(amplitudes, duration)
    end
    return SinglePlotLine(times, amplitudes, [], [])
end

SinglePlotLine(duration::Number) = SinglePlotLine([0., duration], [0., 0.])
function SinglePlotLine(plot_lines::SinglePlotLine...)
    times = Float64[]
    amplitudes = Float64[]
    current_time = 0.
    for pl in plot_lines
        append!(times, pl.times .+ current_time)
        append!(amplitudes, pl.amplitudes)
        current_time += duration(pl)
    end
    return SinglePlotLine(times, amplitudes)
end

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

function SequenceDiagram(; kwargs...)
    @assert length(kwargs) > 0
    durations = duration.([values(kwargs)...])
    @assert all(isapprox.(durations[1], durations, rtol=1e-3))
    res = SinglePlotLine[]
    for symbol in (:RFx, :RFy, :G, :Gx, :Gy, :Gz, :ADC)
        push!(res, get(kwargs, symbol, SinglePlotLine(durations[1])))
    end
    return SequenceDiagram(res...)
end


function SequenceDiagram(bbb::BaseBuildingBlock)
    kwargs = Dict{Symbol, SinglePlotLine}()

    if ndim_grad(bbb) == 3
        for (index, symbol) in enumerate((:Gx, :Gy, :Gz))
            kwargs[symbol] = SinglePlotLine([(time, amplitude[index]) for (time, amplitude) in waveform(bbb)], duration(bbb))
        end
    else
        kwargs[:G] = SinglePlotLine([(time, amplitude) for (time, amplitude) in waveform(bbb)], duration(bbb))
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
                    kwargs[symbol] = SinglePlotLine([0., duration(bbb)], [0., 0.], [delay], [flip_angle(event) * fn(phase(evenet))])
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
    return SequenceDiagram(; kwargs...)
end

function SequenceDiagram(seq::BaseSequence)
    as_lines = SequenceDiagram.(seq)
    return SequenceDiagram([
        SinglePlotLine([getproperty(line, symbol) for line in as_lines]...) 
        for symbol in [:RFx, :RFy, :G, :Gx, :Gy, :Gz, :ADC]]...)
end

end