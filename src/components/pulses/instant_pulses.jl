module InstantPulses
import JuMP: @constraint
import ...AbstractTypes: RFPulseComponent
import ....Variables: VariableType, make_generic, get_free_variable, adjust_internal, variables, @defvar, apply_simple_constraint!, add_cost_function!

"""
    InstantPulse(; flip_angle=nothing, phase=nothing, group=nothing)

Return an instant RF pulse that rotates all spins by `flip_angle` around an axis that has an angle of `phase` with the X-Y plane.

## Parameters
- `group`: name of the group to which this pulse belongs. This is used for scaling or adding phases/off-resonance frequencies.

## Variables
- [`variables.flip_angle`](@ref): angle by which spins are rotated in degrees.
- [`variables.phase`](@ref): angle of axis around which spins are rotated in degrees.
"""
struct InstantPulse <: RFPulseComponent
    flip_angle :: VariableType
    phase :: VariableType
    group :: Union{Nothing, Symbol}
end

function InstantPulse(; flip_angle=nothing, phase=nothing, group=nothing) 
    res = InstantPulse(
        get_free_variable(flip_angle),
        get_free_variable(phase),
        group
    )
    apply_simple_constraint!(res.flip_angle, :>=, 0)
    add_cost_function!(res.phase^2)
    add_cost_function!((res.flip_angle - 90)^2)
    return res
end

@defvar begin
    flip_angle(instant::InstantPulse) = instant.flip_angle
    phase(instant::InstantPulse) = instant.phase
    duration(::InstantPulse) = 0.
    effective_time(::InstantPulse) = 0.
    inverse_bandwidth(::InstantPulse) = 0.
end


make_generic(block::InstantPulse) = block

function adjust_internal(block::InstantPulse; scale=1., frequency=0.)
    InstantPulse(block.flip_angle * scale, block.phase, block.group)
end

end