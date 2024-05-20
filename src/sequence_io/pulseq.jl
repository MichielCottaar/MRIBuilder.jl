"""
Module converting MRIBuilder sequences to and from sequences recognised by [`PulseqIO`](@ref).
"""
module Pulseq
import Interpolations: linear_interpolation
import ..PulseqIO.Types: PulseqSequence, PulseqBlock, PulseqTrapezoid, PulseqGradient
import ...Containers: Sequence, BuildingBlock
import ...Components: GenericPulse, ADC
import ...Scanners: Scanner

function Sequence(pulseq::PulseqSequence; scanner=nothing, B0=3)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    blocks = BuildingBlock.(pulseq.blocks; pulseq.definitions...)
    return Sequence(blocks; name=Symbol(get(pulseq.definitions, :name, "from_pulseq")), scanner=scanner)
end


function BuildingBlock(pulseq::PulseqBlock; BlockDurationRaster, RadiofrequencyRasterTime, GradientRasterTime, kwargs...)
    duration = pulseq.duration * BlockDurationRaster * 1e3

    events = []
    if !isnothing(pulseq.rf)
        if isnothing(pulseq.rf.time)
            time = ((1:length(pulseq.rf.magnitude.samples)) .- 1) .* RadiofrequencyRasterTime .* 1e3
        else
            time = pulseq.rf.time.samples .* 1e3
        end
        push!(events, (
            pulseq.rf.delay * 1e-3,
            GenericPulse(
                time, 
                pulseq.rf.magnitude.samples * pulseq.rf.amplitude * 1e-3, 
                rad2deg.(pulseq.rf.phase.samples .+ pulseq.rf.phase_offset .+ pulseq.rf.frequency .* time .* 2π)
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
    times = sort(unique(vcat([0., duration], _control_times.(grads, GradientRasterTime)...)))
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


end