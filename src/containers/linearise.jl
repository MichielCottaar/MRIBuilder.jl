module Linearise
import StaticArrays: SVector
import ...Components: GradientWaveform, split_timestep, InstantPulse, InstantGradient3D, InstantGradient1D, edge_times
import ...Variables: amplitude, phase, gradient_strength3, duration, qval
import ..Abstract: start_time, end_time, ContainerBlock, iter_instant_gradients, iter_instant_pulses
import ..BaseSequences: BaseSequence, Sequence
import ..BuildingBlocks: BaseBuildingBlock


struct LinearPart{T}
    start_value :: T
    end_value :: T
end

Base.iszero(lp::LinearPart) = iszero(lp.start_value) && iszero(lp.end_value)
Base.iszero(lp::LinearPart{<:AbstractVector}) = all(iszero.(lp.start_value)) && all(iszero.(lp.end_value))
(lp::LinearPart)(time::Number) = @. (1 - time) * lp.start_value + time * lp.end_value

"""
    SequencePart(sequence, time1, time2)

Represents the time between `time1` and `time2` of a larger [`LinearSequence`](@ref)

The gradient, RF amplitude, and RF phase are all be modeled as changing linearly during this time.
"""
struct SequencePart
    gradient :: LinearPart{SVector{3, Float64}}
    gradient_origin :: SVector{3, Float64}
    rf_amplitude :: LinearPart{Float64}
    rf_phase :: LinearPart{Float64}
    duration :: Float64
end


function SequencePart(sequence::BaseSequence{N}, time1::Number, time2::Number) where {N}
    tmean = (time1 + time2) / 2
    nTR = div(tmean, duration(sequence), RoundDown)
    time1 -= nTR * duration(sequence)
    time2 -= nTR * duration(sequence)
    tmean -= nTR * duration(sequence)
    if -1e-9 < time1 < 0.
        time1 = 0.
    end
    if time2 > duration(sequence) && time2 ≈ duration(sequence)
        time2 = duration(sequence)
    end
    if !(0 <= time1 <= time2 <= duration(sequence))
        error("Initial time ($time1) and final time ($time2) are not part of the same TR.")
    end
    for key in 1:N
        if (end_time(sequence, key) > tmean)
            return SequencePart(sequence[key], time1 - start_time(sequence, key), time2 - start_time(sequence, key))
        end
    end
end

function SequencePart(bb::BaseBuildingBlock, time1::Number, time2::Number)
    if -1e-9 < time1 < 0.
        time1 = 0.
    end
    if time2 > duration(bb) && time2 ≈ duration(bb)
        time2 = duration(bb)
    end
    if !(0 <= time1 <= time2 <= duration(bb))
        error("Sequence timings are out of bound")
    end

    tmean = (time1 + time2) / 2

    function fit_linear_part(fn)
        start = fn(bb, time1)
        final = fn(bb, time2)
        middle = fn(bb, tmean)
        mid_deviation = @. middle - (start + final) / 2
        return LinearPart(start .+ mid_deviation ./ 3, final .+ mid_deviation ./ 3)
    end

    return SequencePart(
        fit_linear_part(gradient_strength3),
        zero(SVector{3, Float64}),
        fit_linear_part(amplitude),
        fit_linear_part(phase),
        time2 - time1
    )
end

duration(sp::SequencePart) = sp.duration

"""
    split_times(sequence(s); precision=0.01, max_timestep=Inf)

Suggests at what times to split one or more sequence into linear parts (see [`linearise`](@ref)).

The split times will include any time when (for any of the provided sequences):
- the starting/end points of building blocks, gradients, RF pulses, or ADC readouts.
- a gradient or RF pulse is discontinuous in the first derivative
- the time of any instantaneous gradients, RF pulses, or readouts.

Continuous gradient waveforms or RF pulses might be split up further to ensure the linear approximations meet the required `precision` (see [`split_timestep`](@ref)).
"""
split_times(sequence::BaseSequence, args...; kwargs...) = split_times([sequence], args...; kwargs...)
split_times(sequences::AbstractVector{<:BaseSequence}; kwargs...) = split_times(sequences, 0., maximum(duration.(sequences)); kwargs...)

function split_times(sequences::AbstractVector{<:BaseSequence}, tstart::Number, tfinal::Number; precision=0.01, max_timestep=Inf)
    edges = Float64.([tstart, tfinal])
    for sequence in sequences
        raw_edges = edge_times(sequence)
        nTR_start = Int(div(tstart, duration(sequence), RoundDown))
        nTR_final = Int(div(tfinal, duration(sequence), RoundUp))
        for nTR in nTR_start:nTR_final
            for time in raw_edges .+ (nTR * duration(sequence))
                if tstart < time < tfinal
                    push!(edges, time)
                end
            end
        end
    end
    edges = sort(unique(edges))
    splits = Float64[]

    for (t1, t2) in zip(edges[1:end-1], edges[2:end])
        tmean = (t1 + t2) / 2

        timestep = Float64(max_timestep)
        for sequence in sequences
            (block_time, block) = sequence(mod(tmean, duration(sequence)))
            for key in keys(block)
                if !(start_time(block, key) < block_time < end_time(block, key))
                    continue
                end
                new_timestep = split_timestep(block[key], precision)
                if new_timestep < timestep
                    timestep = new_timestep
                end
            end
        end
        nsteps = isinf(timestep) ? 1 : Int(div(t2 - t1, timestep, RoundUp))
        append!(splits, range(t1, t2, length=nsteps+1))
    end
    return unique(splits)
end


"""
    LinearSequence(sequence, times)
    LinearSequence(sequence; max_timestep=Inf, precision=0.01)
    LinearSequence(sequence, time1, time2; max_timestep=Inf, precision=0.01)

A piece-wise linear approximation of a sequence.

By default it represents the sequence between `time1=0` and `time2=TR`

The gradient, RF amplitude, and RF phase are all be modeled as changing linearly during this time.

If multiple sequences are provided a vector of `LinearSequence` objects are returned, all of which are split at the same time.
"""
struct LinearSequence
    finite_parts :: Vector{SequencePart}
    instant_pulses :: Vector{Union{Nothing, InstantPulse}}
    instant_gradient :: Vector{Union{Nothing, InstantGradient3D}}
end

LinearSequence(sequence::BaseSequence; kwargs...) = LinearSequence(sequence, 0., duration(sequence); kwargs...)
LinearSequence(container::Union{BaseSequence, AbstractVector{<:BaseSequence}}, tstart::Number, tfinal::Number; kwargs...) = LinearSequence(container, split_times(container, tstart, tfinal; kwargs...))
LinearSequence(containers::AbstractVector{<:BaseSequence}, times::AbstractVector{<:Number}) = map(c -> LinearSequence(c, times), containers)
function LinearSequence(container::BaseSequence, times::AbstractVector{<:Number}) 
    parts = [SequencePart(container, t1, t2) for (t1, t2) in zip(times[1:end-1], times[2:end])]

    tstart = times[1]
    tfinal = times[end]

    pulses = Union{Nothing, InstantPulse}[nothing for _ in 1:length(times)]
    gradients = Union{Nothing, InstantGradient3D}[nothing for _ in 1:length(times)]
    for nTR in div(tstart, duration(container), RoundDown):div(tfinal, duration(container), RoundUp)
        for (to_store, func) in [
            (pulses, iter_instant_pulses),
            (gradients, iter_instant_gradients),
        ]
            for (time, pulse) in func(container)
                real_time = time + nTR * duration(container)
                if !(tstart <= real_time < tfinal)
                    continue
                end
                if pulse isa InstantGradient1D
                    pulse = InstantGradient3D(qval(pulse) .* pulse.orientation, pulse.group)
                end
                index = findmin(t -> abs(t - real_time), times)[2] 
                to_store[index] = pulse
            end
        end
    end
    return LinearSequence(parts, pulses, gradients)
end


end