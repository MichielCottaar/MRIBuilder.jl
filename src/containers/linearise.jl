module Linearise
import StaticArrays: SVector
import ...Components: GradientWaveform, split_timestep
import ...Variables: amplitude, phase, gradient_strength3, duration
import ..Abstract: edge_times, start_time, end_time, ContainerBlock
import ..BaseSequences: BaseSequence, Sequence
import ..BuildingBlocks: BaseBuildingBlock


struct LinearPart{T}
    start_value :: T
    end_value :: T
end


"""
    SequencePart(sequence, time1, time2)

Represents the time between `time1` and `time2` of a larger [`Sequence`](@ref)

The gradient, RF amplitude, and RF phase are all be modeled as changing linearly during this time.

See [`linearise`](@ref) to split a sequence into such linear parts.
"""
struct SequencePart
    gradient :: LinearPart{SVector{3, Float64}}
    gradient_origin :: SVector{3, Float64}
    rf_amplitude :: LinearPart{Float64}
    rf_phase :: LinearPart{Float64}
    duration :: Float64
end


function SequencePart(sequence::BaseSequence{N}, time1::Number, time2::Number) where {N}
    if -1e-9 < time1 < 0.
        time1 = 0.
    end
    if time2 > duration(sequence) && time2 ≈ duration(sequence)
        time2 = duration(sequence)
    end
    if !(0 <= time1 <= time2 <= duration(sequence))
        error("Sequence timings are out of bound")
    end
    tmean = (time1 + time2) / 2
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
        @show bb time1 time2 duration(bb)
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
split_times(sequence::BaseSequence; kwargs...) = split_times([sequence]; kwargs...)
split_times(sequences::AbstractVector{<:BaseSequence}; kwargs...) = split_times(sequences, 0., maximum(duration.(sequences)); kwargs...)

function split_times(sequences::AbstractVector{<:BaseSequence}, tstart::Number, tfinal::Number; precision=0.01, max_timestep=Inf)
    edges = [tstart, tfinal]
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
            (block_time, block) = sequence[mod(tmean, duration(sequence))]
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
    linearise(sequence(s), times)
    linearise(sequence(s), tstart, tfinal; max_timestep=Inf, precision=0.01)

Splits any [`Sequence`](@ref) into a series of [`SequencePart`](@ref) objects where the gradients/pulses are approximated to be linear.

If the `times` are not explicitly set they will be obtained from [`split_times`](@ref) (using the values of `max_timestep` and `precision`).
"""
linearise(container::Union{BaseSequence, AbstractVector{<:BaseSequence}}, tstart::Number, tfinal::Number; kwargs...) = linearise(container, split_times(container, tstart, tfinal; kwargs...))
linearise(containers::AbstractVector{<:BaseSequence}, times::AbstractVector{<:Number}) = [linearise(c, times) for c in containers]
linearise(container::BaseSequence, times::AbstractVector{<:Number}) = [SequencePart(container, t1, t2) for (t1, t2) in zip(times[1:end-1], times[2:end])]


end