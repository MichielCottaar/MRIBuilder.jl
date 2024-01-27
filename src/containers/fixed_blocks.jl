
module FixedBlocks

import ...BuildingBlocks: ContainerBlock
import ...Gradients: FixedGradient
import ...Pulses: FixedPulse
import ...Variables: duration, variables
import ...WaitBlocks: WaitBlock

"""
    FixedBlock(duration=nothing; pulse_delay=0., pulse=nothing, gradient_delay=0., gradient=nothing, readout_times=[])

Create a [`BuildingBlock`](@ref) object that has been fully defined.
It can contain an RF pulse ([`FixedPulse`](@ref)), a gradient ([`FixedGradient`](@ref)), and any number of readouts.

The pulse and gradient will play after their respective delays in ms (relative to the start of this block).

If the `duration` is specified the block will last that long.
Otherwise, the block will last as long as the latest block.
"""
struct FixedBlock <: ContainerBlock
    duration :: Float64
    pulse_delay :: Float64
    pulse :: Union{Nothing, FixedPulse}
    gradient_delay :: Float64
    gradient :: Union{Nothing, FixedGradient}
    readout_times :: Vector{Float64}
    function FixedBlock(duration=nothing; pulse_delay=0., pulse=nothing, gradient_delay=0., gradient=nothing, readout_times=[])
        @assert pulse_delay >= 0.
        @assert gradient_delay >= 0.
        if pulse isa Tuple
            pulse = FixedPulse(pulse...)
        end
        if gradient isa Tuple
            gradient = FixedGradient(gradient...)
        end
        real_duration = max(pulse_delay + duration(pulse), gradient_delay + duration(gradient), length(readout_times) == 0 ? 0. : maximum(readout_times))
        if isnothing(duration)
            duration = real_duration
        else
            @assert duration >= real_duration
        end
        new(duration, pulse_delay, pulse, gradient_delay, gradient, Float64.(readout_times))
    end
end

duration(fb::FixedBlock) = fb.duration

variables(::Type{<:FixedBlock}) = []

fixed(wb::WaitBlock) = FixedBlock(value(duration(wb)))


end