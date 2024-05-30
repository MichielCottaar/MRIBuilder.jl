module Abstract
import ...Variables: AbstractBlock, variables, VariableType, get_pulse, get_gradient, @defvar
import ...Components: BaseComponent, InstantPulse, InstantGradient, ReadoutComponent, NoGradient, RFPulseComponent, GradientWaveform

"""
Parent type for `BuildingBlock` or `BaseSequence`, i.e., any building block that contains other MRI components/blocks.

Iterate over them to get the individual components.
"""
abstract type ContainerBlock <: AbstractBlock end


"""
    start_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

Also see [`variables.duration`](@ref), [`end_time`](@ref), and [`effective_time`](@ref)
"""
function start_time(container::ContainerBlock, index1, index2, indices...)
    start_time(container, index1) + start_time(container[index1], index2, indices...)
end


"""
    end_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

Also see [`variables.duration`](@ref), [`start_time`](@ref), and [`effective_time`](@ref)
"""
end_time(container::ContainerBlock, index, indices...) = start_time(container, index) + end_time(container[index], indices...)
end_time(block::AbstractBlock) = variables.duration(block)
end_time(block::Tuple{<:VariableType, <:AbstractBlock}) = variables.duration(block[2])

@defvar begin
    effective_time(container::ContainerBlock, index, indices...) = start_time(container, index) + variables.effective_time(container[index], indices...)
    effective_time(block::Tuple{<:VariableType, <:AbstractBlock}) = block[1] + variables.effective_time(block[2])
end
"""
    effective_time(container, indices...)

Returns the start time of component with given `indices` with respect to the start of the [`ContainerBlock`](@ref).

This will crash if the component does not have an [`effective_time`](@ref) (e.g., if it is (part of) a gradient waveform).

Also see [`variables.duration`](@ref), [`start_time`](@ref), and [`end_time`](@ref)
"""
effective_time


"""
    gradient_strength(sequence, time)

Returns the gradient strength at a particular time within the sequence.
"""
gradient_strength = variables.gradient_strength


"""
    amplitude(sequence, time)

Returns the RF amplitude at a particular time within the sequence in kHz.
"""
amplitude = variables.amplitude


"""
    phase(sequence, time)

Returns the RF phase at a particular time within the sequence in degrees.

NaN is returned if there is no pulse activate at that `time`.
"""
phase = variables.phase

"""
    frequency(sequence, time)

Returns the RF frequency at a particular time within the sequence in kHz.

NaN is returned if there is no pulse activate at that `time`.
"""
frequency = variables.frequency

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
    get_pulse(container, time)

Gets the pulse running at a particular `time` (in ms) during a sequence of building block.

If there is a RF pulse, this function will return a tuple with 2 elements:
1. The [`RFPulseComponent`](@ref) itself
2. The time since the start of the pulse
If there is no active RF pulse, `nothing` is returned.
"""
function get_pulse end

"""
    get_gradient(container, time)

Gets the gradient running at a particular `time` (in ms) during a sequence of building block.

This function will return a tuple with 2 elements:
1. The [`GradientWaveform`](@ref) itself (which could be a [`NoGradient`](@ref) object).
2. The time since the start of the gradient
"""
function get_gradient end

@defvar readout_times(container::ContainerBlock) = [time for (time, _) in iter(container, Val(:readout))]
"""
    readout_times(sequence)

Returns all the times that the sequence will readout.
"""
readout_times

iter(component::Tuple{<:Number, <:ReadoutComponent}, ::Val{:readout}) = [(time, nothing) for time in readout_times(component[2])]

end