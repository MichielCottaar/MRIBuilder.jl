"""
Module converting MRIBuilder sequences to and from sequences recognised by [`PulseqIO`](@ref).
"""
module Pulseq
import Interpolations: linear_interpolation
import ..PulseqIO.Types: PulseqSequence, PulseqBlock, PulseqTrapezoid, PulseqGradient, PulseqRFPulse, PulseqShape, PulseqADC
import ..PulseqIO.Extensions: parse_extension, get_extension_name, add_extension_definition!, PulseqExtension, PulseqExtensionDefinition
import ...Containers: Sequence, BuildingBlock, BaseBuildingBlock, events, waveform, iter_blocks, start_time
import ...Components: GenericPulse, ADC, RFPulseComponent, InstantPulse, InstantGradient, InstantGradient3D, SingleReadout
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
            (
                iszero(dwell_time) ?
                SingleReadout() :
                ADC(
                    pulseq.adc.num,
                    dwell_time,
                    dwell_time * pulseq.adc.num / 2,
                    1.
                )
            )
        ))
    end

    append!(events, pulseq.ext)

    grads = [pulseq.gx, pulseq.gy, pulseq.gz]
    min_duration = max(
        maximum(e[1] + variables.duration(e[2]) for e in events; init=0.),
        maximum(vcat(_control_times.(grads, GradientRasterTime)...); init=0.)
    )

    if min_duration > stated_duration
        if version == v"1.3.1" || isapprox(min_duration, stated_duration, rtol=1e-3)
            stated_duration = min_duration
        else
            error("Minimum duration to play all RF/gradient/ADC events exceeds stated duration.")
        end
    end

    times = sort(unique(vcat([0., stated_duration], _control_times.(grads, GradientRasterTime)...)))
    if length(times) == 1
        push!(times, 0.)
    end
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
    if time ≈ edges[1]
        return grad.shape.samples[1]
    elseif time ≈ edges[end]
        return grad.shape.samples[end]
    end
    return amp * linear_interpolation(edges, grad.shape.samples, extrapolation_bc=0.)(time)
end


function PulseqSequence(seq::Sequence{S}) where {S}
    definitions = (
        Name=S,
        AdcRasterTime=1e-9,
        BlockDurationRaster=1e-9,
        RadiofrequencyRasterTime=1e-9,
        GradientRasterTime=1e-9,
        TotalDuration=variables.duration(seq) * 1e-3,
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
    ext = []
    for (key, event) in events(block)
        gen = make_generic(event)
        if event isa Union{InstantPulse, InstantGradient}
            push!(ext, (start_time(block, key), event))
        elseif event isa RFPulseComponent
            if !isnothing(rf)
                error("Pulseq does not support a single building block containing multiple RF pulses.")
            end
            rf = PulseqRFPulse(
                maximum(gen.amplitude) * 1e3,
                PulseqShape(gen.amplitude ./ maximum(gen.amplitude)),
                PulseqShape(deg2rad.(gen.phase)),
                PulseqShape(gen.time .* 1e-3),
                Int(div(start_time(block, key), 1e-3, RoundNearest)),
                0., 
                0.
            )
        elseif gen isa ADC
            if !isnothing(rf)
                error("Pulseq does not support a single building block containing multiple ADC events.")
            end
            adc = PulseqADC(
                variables.nsamples(gen),
                div(variables.dwell_time(gen) * 1e-3, AdcRasterTime, RoundNearest),
                Int(div(start_time(block, key), 1e-3, RoundNearest)),
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
        max_ampl = maximum(abs.(amplitudes); init=0.)
        if iszero(max_ampl)
            push!(grads, nothing)
        else
            push!(grads, PulseqGradient(
                max_ampl * 1e3,
                PulseqShape(amplitudes ./ max_ampl),
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
        ext
    )
end


# I/O of InstantPulse
function parse_extension(ext::PulseqExtensionDefinition{:InstantPulse})
    mapping = Dict{Int, Tuple{Float64, InstantPulse}}()
    for line in ext.content
        (id, delay, flip_angle, phase) = parse.((Int, Float64, Float64, Float64), split(line))
        mapping[id] = (delay, InstantPulse(flip_angle, phase, nothing))
    end
    return mapping
end

get_extension_name(::Tuple{<:Number, InstantPulse}) = :InstantPulse

function add_extension_definition!(content::Vector{String}, obj::Tuple{Number, InstantPulse})
    to_store = (obj[1], obj[2].flip_angle, obj[2].phase)
    for line in content
        (id, this_line...) = parse.((Int, Float64, Float64, Float64), split(line))
        if all(to_store .≈ this_line)
            return id
        end
    end
    push!(content, "$(length(content) + 1) " * join(string.(to_store), " "))
    return length(content)
end


# I/O of InstantGradient
function parse_extension(ext::PulseqExtensionDefinition{:InstantGradient})
    mapping = Dict{Int, Tuple{Float64, InstantGradient3D}}()
    for line in ext.content
        (id, delay, qvec...) = parse.((Int, Float64, Float64, Float64, Float64), split(line))
        mapping[id] = (delay, InstantGradient3D([qvec...], nothing))
    end
    return mapping
end

get_extension_name(::Tuple{<:Number, <:InstantGradient}) = :InstantGradient

function add_extension_definition!(content::Vector{String}, obj::Tuple{<:Number, <:InstantGradient})
    to_store = (obj[1], variables.qvec(obj[2])...)
    for line in content
        (id, this_line...) = parse.((Int, Float64, Float64, Float64, Float64), split(line))
        if all(to_store .≈ this_line)
            return id
        end
    end
    push!(content, "$(length(content) + 1) " * join(string.(to_store), " "))
    return length(content)
end


end