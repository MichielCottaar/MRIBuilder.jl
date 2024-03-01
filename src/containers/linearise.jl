module Linearise
import LinA
import StaticArrays: SVector
import ...Components: GradientWaveform, split_times
import ...Variables: amplitude, phase, gradient_strength3, duration
import ..Abstract: edge_times, start_time, end_time
import ..BaseSequences: BaseSequence, Sequence
import ..BuildingBlocks: BaseBuildingBlock


struct LinearPart{T}
    start_value :: T
    end_value :: T
end


"""
Represents a small part of a larger [`Sequence`](@ref), where the gradient, RF amplitude, and RF phase can all be modeled as changing linearly.
"""
struct SequencePart
    gradient :: LinearPart{SVector{3, Float64}}
    gradient_origin :: SVector{3, Float64}
    rf_amplitude :: LinearPart{Float64}
    rf_phase :: LinearPart{Float64}
    duration :: Float64
end

duration(sp::SequencePart) = sp.duration


function linearise(seq::BaseSequence; kwargs...)
    parts = SequencePart[]
    for block in seq
        append!(parts, linearise(block; kwargs...))
    end
    return parts
end

function linearise(bb::BaseBuildingBlock; precision=0.01)
    parts = SequencePart[]
    if duration(bb) <= 0.
        return parts
    end
    edges = edge_times(bb)
    for (t1, t2) in zip(edges[1:end-1], edges[2:end])
        tmean = (t1 + t2) / 2

        # determine where to split this domain
        possible_splits = Any[[t1, t2], ]
        for key in keys(bb)
            if !(start_time(bb, key) < tmean < end_time(bb, key))
                continue
            end
            block = bb[key]
            append!(possible_splits, start_time(bb, key) .+ split_times(block, t1 - start_time(bb, key), t2 - start_time(bb, key), precision))
        end
        tsplits = argmax(length, possible_splits)

        for (t1b, t2b) in zip(tsplits[1:end-1], tsplits[2:end])
            tmeanb = (t1b + t2b) / 2

            function fit_linear_part(fn)
                start = fn(bb, t1b)
                final = fn(bb, t2b)
                middle = fn(bb, tmeanb)
                mid_deviation = @. middle - (start + final) / 2
                return LinearPart(start .+ mid_deviation ./ 3, final .+ mid_deviation ./ 3)
            end

            push!(parts, SequencePart(
                fit_linear_part(gradient_strength3),
                zero(SVector{3, Float64}),
                fit_linear_part(amplitude),
                fit_linear_part(phase),
                t2b - t1b
            ))
        end
    end
    return parts
end


end