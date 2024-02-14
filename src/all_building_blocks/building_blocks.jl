module BuildingBlocks
import LinearAlgebra: norm
import ..BaseBuildingBlocks: BaseBuildingBlock, events, waveform_sequence
import ...Variables: VariableType, duration, make_generic, get_pulse, get_readout
import ...Components: BaseComponent, DelayedEvent, RFPulseComponent, ReadoutComponent
import ...Readouts: InstantReadout, ADC
import ...Pulses: RFPulseBlock
import ...Gradients: GradientBlock
import ...BuildingBlocks: scanner_constraints!, make_generic

"""
    BuildingBlock(waveform, events; duration=nothing, orientation=nothing)

Generic [`BaseBuildingBlock`](@ref) that can capture any overlapping gradients, RF pulses, and/or readouts.
The gradients cannot contain any free variables.

## Arguments
- `waveform`: Sequence of 2-element tuples with (time, (Gx, Gy, Gz)). If `orientation` is set then the tuple is expected to look like (time, G). This cannot contain any free variables.
- `events`: Sequence of 2-element tuples with (index, pulse/readout). The start time of the pulse/readout at the start of the gradient waveform element with index `index` (use [`DelayedEvent`](@ref) to make this earlier or later).
- `duration`: duration of this `BuildingBlock`. If not set then it will be assumed to be the time of the last element in `waveform`.
- `orientation`: orientation of the gradients in the waveform. If not set, then the full gradient vector should be given explicitly.
"""
struct BuildingBlock <: BaseBuildingBlock
    parts :: Vector{<:BaseComponent}
    function BuildingBlock(parts::AbstractVector{<:BaseComponent})
        res = new(duration, parts)
        for (_, part) in waveform_sequence(parts)
            scanner_constraints!(part)
        end
        return res
    end
end

function BuildingBlock(waveform::AbstractVector, events::AbstractVector; duration=nothing, orientation=nothing)
    events = Any[events...]
    waveform = Any[waveform...]
    ndim = isnothing(orientat) ? 1 : 3
    zero_grad = isnothing(orientation) ? zeros(3) : 0.
    if length(waveform) == 0 || waveform[1][1] > 0.
        pushfirst!(waveform, (0., zero_grad))
        events = [(i+1, e) for (i, e) in events]
    end

    if isnothing(duration)
        duration = waveform[end][1]
    end
    if !(duration ≈ waveform[end][1])
        @assert duration > waveform[end][1]
        push!(waveform, (duration, zero_grad))
    end
    components = BaseComponent[]
    for (index_grad, ((prev_time, prev_grad), (time, grad))) in enumerate(zip(waveform[1:end-1], waveform[2:end]))
        duration = time - prev_time
        if norm(prev_grad) <= 1e-12 && norm(grad) <= 1e-12
            push!(components, NoGradientBlock{ndim}(duration))
        elseif norm(prev_grad) ≈ norm(grad)
            push!(components, ConstantGradientBlock(prev_grad, duration))
        else
            push!(components, ChangingGradientBlock(prev_grad, (grad .- prev_grad) ./ duration, duration))
        end
        while length(events) > 0 && index_grad == events[1][1]
            (_, event) = popfirst!(events)
            push!(components, event)
        end
    end
    return components
end

make_generic(other_block::BaseBuildingBlock) = BuildingBlock(duration(other_block), [other_block...])
Base.keys(bb::BuildingBlock) = 1:length(bb.parts)
Base.getindex(bb::BuildingBlock, i::Integer) = bb.parts[i]

duration(bb::BuildingBlock) = sum(duration, waveform_sequence(bb); init=0.)

function get_pulse(bb::BuildingBlock)
    pulses = [p for p in events(bb) if p isa RFPulseComponent]
    if length(pulses) == 0
        error("BuildingBlock does not contain any pulses.")
    end
    if length(pulses) == 1
        return pulses[1]
    end
    error("BuildingBlock contains more than one pulse. Not sure which one to return.")
end

function get_readout(bb::BuildingBlock)
    readouts = [r for r in events(bb) if r isa ReadoutComponent]
    if length(readouts) == 0
        error("BuildingBlock does not contain any readouts.")
    end
    if length(readouts) == 1
        return readouts[1]
    end
    error("BuildingBlock contains more than one readout. Not sure which one to return.")
end

end