module ConstantPulses
import JuMP: @constraint
import ...AbstractTypes: RFPulseComponent, split_timestep
import ....Variables: VariableType, set_simple_constraints!, make_generic, get_free_variable, adjust_internal, variables, @defvar, apply_simple_constraint!
import ..GenericPulses: GenericPulse

"""
    ConstantPulse(; variables...)

Represents an radio-frequency pulse with a constant amplitude and frequency (i.e., a rectangular function).

## Parameters
- `group`: name of the group to which this pulse belongs. This is used for scaling or adding phases/off-resonance frequencies.

## Variables
- [`variables.flip_angle`](@ref): rotation expected for on-resonance spins in degrees.
- [`variables.duration`](@ref): duration of the RF pulse in ms.
- [`variables.amplitude`](@ref): amplitude of the RF pulse in kHz.
- [`variables.phase`](@ref): phase at the start of the RF pulse in degrees.
- [`variables.frequency`](@ref): frequency of the RF pulse relative to the Larmor frequency (in kHz).
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
    apply_simple_constraint!(res.amplitude, :>=, 0)
    set_simple_constraints!(res, kwargs)
    return res
end

@defvar duration(pulse::ConstantPulse) = pulse.duration
@defvar effective_time(pulse::ConstantPulse) = variables.duration(pulse) / 2

@defvar pulse begin
    amplitude(pulse::ConstantPulse) = pulse.amplitude
    phase(pulse::ConstantPulse) = pulse.phase
    frequency(pulse::ConstantPulse) = pulse.frequency
    amplitude(pulse::ConstantPulse, time::Number) = variables.amplitude(pulse)
end

@defvar pulse begin
    flip_angle(pulse::ConstantPulse) = variables.amplitude(pulse) * variables.duration(pulse) * 360
    inverse_bandwidth(pulse::ConstantPulse) = variables.duration(pulse) * 4

    phase(pulse::ConstantPulse, time::Number) = variables.phase(pulse) + variables.frequency(pulse) * (time - variables.effective_time(pulse)) * 360.
    frequency(pulse::ConstantPulse, time::Number) = variables.frequency(pulse)
end

function make_generic(block::ConstantPulse)
    d = variables.duration(block)
    final_phase = variables.phase(block) + d * variables.frequency(block) * 360
    return GenericPulse(
        [0., d], 
        [variables.amplitude(block), variables.amplitude(block)],
        [variables.phase(block), final_phase],
        variables.effective_time(block)
    )
end


split_timestep(pulse::ConstantPulse, precision) = Inf

function adjust_internal(block::ConstantPulse; scale=1., frequency=0., stretch=1.)
    ConstantPulse(
        block.amplitude * scale,
        block.duration * stretch,
        block.phase,
        block.frequency + frequency,
        block.group,
    )
end

end