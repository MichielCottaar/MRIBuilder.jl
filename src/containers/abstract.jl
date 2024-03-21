module Abstract
import ...Variables: AbstractBlock, duration, effective_time, gradient_strength, amplitude, phase, VariableType
import ...Components.Readouts: readout_times
import ...Components: BaseComponent, InstantPulse, InstantGradient, ReadoutComponent

"""
Parent type for `BuildingBlock` or `BaseSequence`, i.e., any building block that contains other MRI components/blocks.

Iterate over them to get the individual components.
"""
abstract type ContainerBlock <: AbstractBlock end


"""
    start_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

Also see [`duration`](@ref), [`end_time`](@ref), and [`effective_time`](@ref)
"""
function start_time(container::ContainerBlock, index1, index2, indices...)
    start_time(container, index1) + start_time(container[index1], index2, indices...)
end


"""
    end_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

Also see [`duration`](@ref), [`start_time`](@ref), and [`effective_time`](@ref)
"""
end_time(container::ContainerBlock, index, indices...) = start_time(container, index) + end_time(container[index], indices...)
end_time(block::AbstractBlock) = duration(block)
end_time(block::Tuple{<:VariableType, <:AbstractBlock}) = block[1] + end_time(block[2])

"""
    effective_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

This will crash if the component does not have an [`effective_time`](@ref) (e.g., if it is (part of) a gradient waveform).

Also see [`duration`](@ref), [`start_time`](@ref), and [`end_time`](@ref)
"""
effective_time(container::ContainerBlock, index, indices...) = start_time(container, index) + effective_time(container[index], indices...)
effective_time(block::Tuple{<:VariableType, <:AbstractBlock}) = block[1] + effective_time(block[2])


"""
    edge_times(container)

Returns all the edge times during a sequence in ms.

Edges are defined as any time, when:
- the edge of a building block
- the slope of the gradient profile changes suddenly
- an RF pulse starts or ends
"""
function edge_times end

"""
    gradient_strength(sequence, time)

Returns the gradient strength at a particular time within the sequence.
"""
function gradient_strength end


"""
    amplitude(sequence, time)

Returns the RF amplitude at a particular time within the sequence.
"""
function amplitude end


"""
    phase(sequence, time)

Returns the RF phase at a particular time within the sequence.
NaN if the ampltiude is zero.
"""
function phase end

"""
    iter(sequence, get_type)

Helper functions for any `iter_*` functions.
"""
function iter(container::ContainerBlock, get_type::Val) 
    [(start_time(container, key) + t, component) for key in keys(container) for (t, component) in iter(container[key], get_type)]
end
iter(container::ContainerBlock, get_type::Symbol) = iter(container, Val(get_type))
iter(component::BaseComponent, get_type::Val) = []
iter(component::Tuple{<:Number, <:BaseComponent}, get_type::Val) = []

"""
    iter_blocks(sequence)

Returns all the building blocks in the sequence with the time they will start
"""
iter_blocks(container::ContainerBlock) = iter(container, Val(:block))

"""
    iter_instant_pulses(sequence)

Returns all the [`InstantPulse`](@ref) within the sequence with their timings
"""
iter_instant_pulses(container::ContainerBlock) = iter(container, Val(:instantpulse))
iter(component::Tuple{<:Number, <:InstantPulse}, ::Val{:instantpulse}) = [component]

"""
    iter_instant_gradients(sequence)

Returns all the [`InstantGradient`](@ref) within the sequence with their timings
"""
iter_instant_gradients(container::ContainerBlock) = iter(container, Val(:instantgradient))
iter(component::Tuple{<:Number, <:InstantGradient}, ::Val{:instantgradient}) = [component]

"""
    readout_times(sequence)

Returns all the times that the sequence will readout.
"""
readout_times(container::ContainerBlock) = [time for (time, _) in iter(container, Val(:readout))]
iter(component::Tuple{<:Number, <:ReadoutComponent}, ::Val{:readout}) = [(time, nothing) for time in readout_times(component[2])]

end