module ConstantPulses
import JuMP: VariableRef, @constraint, @variable, value
import ....Sequences: RFPulse
import ...BuildingBlocks: BuildingBlock, properties, BuildingBlockPlaceholder, set_simple_constraints!, duration, to_mcmr_components
import ...SequenceBuilders: SequenceBuilder, owner_model, start_time, end_time
import ..Properties: flip_angle, phase, amplitude, frequency, bandwidth

"""
    ConstantPulse(; variables...)

Represents an radio-frequency pulse with a constant amplitude and frequency (i.e., a rectangular function).

## Properties
- [`flip_angle`](@ref): rotation expected for on-resonance spins in degrees.
- [`duration`](@ref): duration of the RF pulse in ms.
- [`amplitude`](@ref): amplitude of the RF pulse in kHz.
- [`phase`](@ref): phase at the start of the RF pulse in degrees.
- [`frequency`](@ref): frequency of the RF pulse relative to the Larmor frequency (in kHz).
"""
struct ConstantPulse <: BuildingBlock
    builder :: SequenceBuilder
    amplitude :: VariableRef
    duration :: VariableRef
    phase :: VariableRef
    frequency :: VariableRef
end

ConstantPulse(; kwargs...) = BuildingBlockPlaceholder{ConstantPulse}(; kwargs...)
function ConstantPulse(builder::SequenceBuilder; kwargs...) 
    model = owner_model(builder)
    res = ConstantPulse(
        builder,
        @variable(model),
        @variable(model),
        @variable(model),
        @variable(model)
    )
    @constraint model amplitude(res) >= 0
    set_simple_constraints!(model, res, kwargs)
    return res
end

amplitude(pulse::ConstantPulse) = pulse.amplitude
duration(pulse::ConstantPulse) = pulse.duration
phase(pulse::ConstantPulse) = pulse.phase
frequency(pulse::ConstantPulse) = pulse.frequency
flip_angle(pulse::ConstantPulse) = amplitude(pulse) * duration(pulse) * 360
bandwidth(pulse::ConstantPulse) = 3.79098854 / duration(pulse)

properties(::Type{<:ConstantPulse}) = [amplitude, duration, phase, frequency, flip_angle, bandwidth]

function to_mcmr_components(block::ConstantPulse)
    final_phase = phase(block) + duration(block) * frequency(block) * 360
    return RFPulse(
        value.([start_time(block), end_time(block)]),
        value.([amplitude(block), amplitude(block)]),
        value.([phase(block), final_phase])
    )
end


end