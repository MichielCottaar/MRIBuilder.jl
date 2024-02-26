module Abstract
import ...Variables: AbstractBlock, duration, effective_time
import ...Components.Readouts: readout_times

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

"""
    effective_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

This will crash if the component does not have an [`effective_time`](@ref) (e.g., if it is (part of) a gradient waveform).

Also see [`duration`](@ref), [`start_time`](@ref), and [`end_time`](@ref)
"""
effective_time(container::ContainerBlock, index, indices...) = start_time(container, index) + effective_time(container[index], indices...)


"""
    readout_times(sequence/ADC)

Return the times of all of the readout events in ms.
"""
function readout_times end

"""
    edge_times(container)

Returns all the edge times during a sequence in ms.

Edges are defined as any time, when:
- the edge of a building block
- the slope of the gradient profile changes suddenly
- an RF pulse starts or ends
"""
function edge_times end



end