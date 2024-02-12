module Generic
import ..Abstract: AbstractOverlapping, interruptions, duration, waveform_sequence
import ...Variables: VariableType
import ...Wait: WaitBlock
import ...Readouts: InstantReadout, ADC
import ...Pulses: RFPulseBlock
import ...Gradients: GradientBlock
import ...BuildingBlocks: scanner_constraints!, make_generic

"""
    GenericWaveform(duration, waveform, interruptions)

Generic description that can capture any overlapping gradients, RF pulses, and/or readouts.

Interruptions are stores as tuples with 3 fields:
- `index`: which part of the `waveform` is being interrupted (cannot be a free variable).
- `time`: time of the interruption relative to the start of the `waveform` part (between 0 and the length of this part). This can be a variable.
- `block`: [`RFPulseBlock`](@ref) or [`InstantReadout`](@ref)
"""
struct GenericOverlapping <: AbstractOverlapping
    duration :: VariableType
    waveform :: Vector{Union{WaitBlock, GradientBlock}}
    interruptions :: Vector{NamedTuple{(:index, :time, :object), Tuple{Int64, <:VariableType, <:Union{RFPulseBlock, InstantReadout, ADC}}}}
    function GenericOverlapping(duration::VariableType, waveform::AbstractVector{<:Union{WaitBlock, GradientBlock}}, interruptions::AbstractVector=[])
        res = new(duration, waveform, interruptions)
        scanner_constraints!.(waveform)
        return res
    end
end


GenericOverlapping(other_waveform::AbstractOverlapping) = GenericOverlapping(duration(other_waveform), waveform_sequence(other_waveform), interruptions(other_waveform))

make_generic(ao::AbstractOverlapping) = GenericOverlapping(
    duration(ao), 
    waveform_sequence(ao), 
    [(index=i.index, time=i.time, object=make_generic(i.object)) for i in interruptions(ao)]
)

waveform_sequence(go::GenericOverlapping) = go.waveform
interruptions(go::GenericOverlapping) = go.interruptions

make_generic(wait::WaitBlock) = GenericOverlapping(wait.duration, [], [])

end