module ConstantPulses
import JuMP: VariableRef, @constraint, @variable, value
import ...BuildingBlocks: RFPulseBlock, set_simple_constraints!, fixed
import ...Variables: variables, get_free_variable, flip_angle, phase, amplitude, frequency, start_time, end_time, VariableType, duration, effective_time, inverse_bandwidth
import ...BuildSequences: global_model
import ..FixedPulses: FixedPulse

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
struct ConstantPulse <: RFPulseBlock
    amplitude :: VariableType
    duration :: VariableType
    phase :: VariableType
    frequency :: VariableType
    scale :: Union{Nothing, Symbol}
end

function ConstantPulse(amplitude=nothing, duration=nothing, phase=nothing, frequency=nothing, scale=nothing, kwargs...) 
    res = ConstantPulse(
        [get_free_variable(value) for value in (amplitude, duration, phase, frequency)]...,
        scale

    )
    @constraint global_model() res.amplitude >= 0
    set_simple_constraints!(res, kwargs)
    return res
end

amplitude(pulse::ConstantPulse) = pulse.amplitude
duration(pulse::ConstantPulse) = pulse.duration
phase(pulse::ConstantPulse) = pulse.phase
frequency(pulse::ConstantPulse) = pulse.frequency
flip_angle(pulse::ConstantPulse) = amplitude(pulse) * duration(pulse) * 360
inverse_bandwidth(pulse::ConstantPulse) = duration(pulse) * 4
effective_time(pulse::ConstantPulse) = duration(pulse) / 2

variables(::Type{<:ConstantPulse}) = [amplitude, duration, phase, frequency, flip_angle, inverse_bandwidth]

function fixed(block::ConstantPulse)
    d = value(duration(block))
    final_phase = phase(block) + d * frequency(block) * 360
    return FixedPulse(
        [0., d], 
        value.([amplitude(block), amplitude(block)]),
        value.([phase(block), final_phase]),
    )
end


end