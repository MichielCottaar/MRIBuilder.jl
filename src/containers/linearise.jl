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
    if !(0 <= time1 <= time2 <= duration(sequence))
        error("Sequence timings are out of bound")
    end
    for key in 1:N
        if (end_time(sequence, key) > time1)
            return SequencePart(sequence[key], time1 - start_time(sequence, key), time2 - start_time(sequence, key))
        end
    end
end

function SequencePart(bb::BaseBuildingBlock, time1::Number, time2::Number)
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
    split_times(sequence; precision=0.01)

Suggests at what times to split a sequence into linear parts (see [`linearise`](@ref)).

The split times will include any time when:
- the starting/end points of building blocks, gradients, RF pulses, or ADC readouts.
- a gradient or RF pulse is discontinuous in the first derivative
- the time of any instantaneous gradients, RF pulses, or readouts.

Continuous gradient waveforms or RF pulses might be split up further to ensure the linear approximations meet the required `precision` (see [`split_timestep`](@ref)).
"""
function split_times(sequence::BaseSequence; precision=0.01)
    splits = Float64[]
    start_time = 0.
    for block in sequence
        append!(splits, start_time .+ split_times(block; precision=precision))
        start_time += duration(block)
    end
    return sort(unique(splits))
end


function split_times(bb::BaseBuildingBlock; precision=0.01)
    splits = Float64[]
    edges = edge_times(bb)
    for (t1, t2) in zip(edges[1:end-1], edges[2:end])
        tmean = (t1 + t2) / 2

        # determine where to split this domain
        timestep = Inf
        for key in keys(bb)
            if !(start_time(bb, key) < tmean < end_time(bb, key))
                continue
            end
            new_timestep = split_timestep(bb[key], precision)
            if new_timestep < timestep
                timestep = new_timestep
            end
        end
        nsteps = (isinf(timestep) ? 1 : Int(div(t2 - t1, timestep, RoundUp))) + 1
        append!(splits, range(t1, t2, length=nsteps))
    end
    return sort(unique(splits))
end



"""
    linearise(sequence[, times])

Splits any [`BaseSequence`](@ref) or [`BaseBuildingBlock`](@ref) into a series of [`SequencePart`](@ref) objects where the gradients/pulses are approximated to be linear.

If the `times` are not explicitly set they will be obtained from [`split_times`](@ref).
"""
linearise(container::ContainerBlock) = linearise(container, split_times(container))
linearise(container::ContainerBlock, times::AbstractVector{<:Number}) = [SequencePart(container, t1, t2) for (t1, t2) in zip(times[1:end-1], times[2:end])]


end