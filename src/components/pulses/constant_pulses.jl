module ConstantPulses
import JuMP: @constraint
import ...AbstractTypes: RFPulseComponent
import ....BuildSequences: global_model
import ....Variables: duration, amplitude, effective_time, flip_angle, phase, inverse_bandwidth, VariableType, set_simple_constraints!, frequency
import ..GenericPulses: GenericPulse

"""
    ConstantPulse(; variables...)

Represents an radio-frequency pulse with a constant amplitude and frequency (i.e., a rectangular function).

## Parameters
- `group`: name of the group to which this pulse belongs. This is used for scaling or adding phases/off-resonance frequencies.

## Variables
- [`flip_angle`](@ref): rotation expected for on-resonance spins in degrees.
- [`duration`](@ref): duration of the RF pulse in ms.
- [`amplitude`](@ref): amplitude of the RF pulse in kHz.
- [`phase`](@ref): phase at the start of the RF pulse in degrees.
- [`frequency`](@ref): frequency of the RF pulse relative to the Larmor frequency (in kHz).
"""
struct ConstantPulse <: RFPulseComponent
    amplitude :: VariableType
    duration :: VariableType
    phase :: VariableType
    frequency :: VariableType
    group :: Union{Nothing, Symbol}
end

function ConstantPulse(; amplitude=nothing, duration=nothing, phase=nothing, frequency=nothing, group=nothing, kwargs...) 
    res = ConstantPulse(
        [get_free_variable(value) for value in (amplitude, duration, phase, frequency)]...,
        group
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

function make_generic(block::ConstantPulse)
    d = duration(block)
    final_phase = phase(block) + d * frequency(block) * 360
    return GenericPulse(
        [0., d], 
        [amplitude(block), amplitude(block)],
        [phase(block), final_phase],
        effective_time(block)
    )
end


end