module ConstantPulses
import JuMP: VariableRef, @constraint, @variable, value, Model
import ...BuildingBlocks: BuildingBlock, set_simple_constraints!
import ...ConcreteBlocks: ConcreteBlock, to_concrete_block
import ...Variables: variables, get_free_variable, flip_angle, phase, amplitude, frequency, bandwidth, start_time, end_time, VariableType, duration
import ...GlobalModel: @global_model_constructor

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
    model :: Model
    amplitude :: VariableType
    duration :: VariableType
    phase :: VariableType
    frequency :: VariableType
end

@global_model_constructor ConstantPulse

function ConstantPulse(model::Model; amplitude=nothing, duration=nothing, phase=nothing, frequency=nothing, kwargs...) 
    res = ConstantPulse(
        model,
        [get_free_variable(model, value) for value in (amplitude, duration, phase, frequency)]...
    )
    @constraint model res.amplitude >= 0
    set_simple_constraints!(model, res, kwargs)
    return res
end

amplitude(pulse::ConstantPulse) = pulse.amplitude
duration(pulse::ConstantPulse) = pulse.duration
phase(pulse::ConstantPulse) = pulse.phase
frequency(pulse::ConstantPulse) = pulse.frequency
flip_angle(pulse::ConstantPulse) = amplitude(pulse) * duration(pulse) * 360
bandwidth(pulse::ConstantPulse) = 3.79098854 / duration(pulse)

variables(::Type{<:ConstantPulse}) = [amplitude, duration, phase, frequency, flip_angle, bandwidth]

function to_concrete_block(block::ConstantPulse)
    d = value(duration(block))
    final_phase = value(phase(block)) + d * value(frequency(block)) * 360
    return ConcreteBlock(value(duration(block)), pulse=[
        ([0., d]),
        value.([amplitude(block), amplitude(block)]),
        value.([phase(block), final_phase])
    ])
end


end