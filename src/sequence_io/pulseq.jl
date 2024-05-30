"""
Module converting MRIBuilder sequences to and from sequences recognised by [`PulseqIO`](@ref).
"""
module Pulseq
import Interpolations: linear_interpolation
import ..PulseqIO.Types: PulseqSequence, PulseqBlock, PulseqTrapezoid, PulseqGradient, PulseqRFPulse, PulseqShape, PulseqADC, PulseqExtension
import ...Containers: Sequence, BuildingBlock, BaseBuildingBlock, events, waveform, iter_blocks
import ...Components: GenericPulse, ADC, RFPulseComponent
import ...Scanners: Scanner
import ...Variables: variables, make_generic

function Sequence(pulseq::PulseqSequence; scanner=nothing, B0=nothing)
    if isnothing(scanner)
        use_B0 = isnothing(B0) ? get(pulseq.definitions, :B0, 3.) : B0
        scanner = Scanner(B0=use_B0)
    end
    blocks = BuildingBlock.(pulseq.blocks; pulseq.definitions..., version=pulseq.version)
    return Sequence(blocks; name=Symbol(get(pulseq.definitions, :Name, "from_pulseq")), scanner=scanner)
end


function BuildingBlock(pulseq::PulseqBlock; version, BlockDurationRaster, RadiofrequencyRasterTime, GradientRasterTime, kwargs...)
    stated_duration = pulseq.duration * BlockDurationRaster * 1e3

    events = []
    if !isnothing(pulseq.rf)
        f(samples) = isnothing(pulseq.rf.time) ? [samples[1], samples..., samples[end]] : samples
        if isnothing(pulseq.rf.time)
            time = [0., ((1:length(pulseq.rf.magnitude.samples)) .- 0.5)..., length(pulseq.rf.magnitude.samples)] .* RadiofrequencyRasterTime .* 1e3
        else
            time = pulseq.rf.time.samples .* 1e3 * RadiofrequencyRasterTime
        end
        push!(events, (
            pulseq.rf.delay * 1e-3,
            GenericPulse(
                time, 
                f(pulseq.rf.magnitude.samples) * pulseq.rf.amplitude * 1e-3, 
                rad2deg.(f(pulseq.rf.phase.samples) .+ pulseq.rf.phase_offset .+ pulseq.rf.frequency .* time .* 2π)
            )
        ))
    end
    if !isnothing(pulseq.adc)
        dwell_time = pulseq.adc.dwell * 1e-6
        push!(events, (
            pulseq.adc.delay * 1e-3,
            ADC(
                pulseq.adc.num,
                dwell_time,
                dwell_time * pulseq.adc.num / 2,
                1.
            )
        ))
    end

    grads = [pulseq.gx, pulseq.gy, pulseq.gz]
    min_duration = max(
        maximum(e[1] + variables.duration(e[2]) for e in events; init=0.),
        maximum(vcat(_control_times.(grads, GradientRasterTime)...); init=0.)
    )

    if min_duration > stated_duration
        if version == v"1.3.1"
            stated_duration = min_duration
        else
            error("Minimum duration to play all RF/gradient/ADC events exceeds stated duration.")
        end
    end

    times = sort(unique(vcat([0., stated_duration], _control_times.(grads, GradientRasterTime)...)))
    waveform = [(t, _get_amplitude.(grads, t, GradientRasterTime)) for t in times]

    return BuildingBlock(waveform, events)
end

_control_times(::Nothing, ::Number) = Float64[]
_control_times(trap::PulseqTrapezoid, ::Number) = cumsum([trap.delay, trap.rise, trap.flat, trap.fall]) * 1e-3
function _control_times(grad::PulseqGradient, raster::Number)
    if isnothing(grad.time)
        return ((1:length(grad.shape.samples)) .- 0.5) .* (1e3 * raster)
    else
        return grad.time.samples .* (1e3 * raster)
    end
end

_get_amplitude(::Nothing, ::Number, ::Number) = 0.
function _get_amplitude(trap::PulseqTrapezoid, time::Number, raster::Number)
    amp = trap.amplitude * 1e-3
    edges = _control_times(trap, raster)
    if !(edges[1] < time < edges[end])
        return 0.
    elseif time < edges[2]
        return amp * (time - edges[1]) / (1e-3 * trap.rise)
    elseif time < edges[3]
        return amp
    else
        return amp * (edges[end] - time) / (1e-3 * trap.fall)
    end
end

function _get_amplitude(grad::PulseqGradient, time::Number, raster::Number)
    amp = grad.amplitude * 1e-3
    edges = _control_times(grad, raster)
    return amp * linear_interpolation(edges, grad.shape.samples, extrapolation_bc=0.)(time)
end


function PulseqSequence(seq::Sequence{S}) where {S}
    definitions = (
        Name=S,
        AdcRasterTime=1e-9,
        BlockDurationRaster=1e-9,
        RadiofrequencyRasterTime=1e-9,
        GradientRasterTime=1e-9,
        TotalDuration=duration(seq) * 1e-3,
        B0=seq.scanner.B0,
    )
    blocks = [PulseqBlock(block; definitions...) for (_, block) in iter_blocks(seq)]
    return PulseqSequence(
        v"1.4.0",
        definitions,
        blocks
    )
end

function PulseqBlock(block::BaseBuildingBlock; BlockDurationRaster, AdcRasterTime, kwargs...)
    rf = nothing
    adc = nothing
    for (key, event) in events(block)
        if event isa RFPulseComponent
            if !isnothing(rf)
                error("Pulseq does not support a single building block containing multiple RF pulses.")
            end
            gen = make_generic(event)
            rf = PulseqRFPulse(
                maximum(gen.amplitude) * 1e3,
                PulseqShape(gen.amplitude ./ maximum(gen.amplitude)),
                PulseqShape(deg2rad.(gen.phase)),
                PulseqShape(gen.time .* 1e-3),
                Int(div(start_time(block, key), 1e-3, RoundNearest)),
                0., 
                0.
            )
        elseif event isa ADC
            if !isnothing(rf)
                error("Pulseq does not support a single building block containing multiple ADC events.")
            end
            adc = PulseqADC(
                variables.nsamples(event),
                div(variables.dwell_time(event), AdcRasterTime, RoundNearest),
                Int(div(delay, 1e-3, RoundNearest)),
                0., 0.
            )
        else
            error("Cannot write $(typeof(event)) to Pulseq.")
        end
    end

    grads = []
    times = [t for (t, _) in waveform(block)]
    for dim in 1:3
        amplitudes = [g[dim] for (_, g) in waveform(block)]
        if iszero(maximum(abs.(amplitudes); init=0.))
            push!(grads, nothing)
        else
            push!(grads, PulseqGradient(
                maximum(amplitudes) * 1e3,
                PulseqShape(amplitudes ./ maximum(amplitudes)),
                PulseqShape(times .* 1e-3),
                0.,
            ))
        end
    end
    
    return PulseqBlock(
        Int(div(1e-3 * variables.duration(block), BlockDurationRaster, RoundNearest)),
        rf,
        grads...,
        adc,
        Tuple{PulseqExtension, Int}[]
    )
end

end