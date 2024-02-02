module InstantPulses
import JuMP: @constraint, @variable, VariableRef, value
import ...BuildingBlocks: RFPulseBlock
import ...Variables: flip_angle, phase, start_time, variables, duration, get_free_variable, VariableType, effective_time, inverse_bandwidth
import ...BuildSequences: global_model
import ..FixedPulses: FixedInstantPulse

struct InstantRFPulseBlock <: RFPulseBlock
    flip_angle :: VariableType
    phase :: VariableType
    scale :: Union{Nothing, Symbol}
end

function InstantRFPulseBlock(; flip_angle=nothing, phase=nothing, scale=nothing) 
    res = InstantRFPulseBlock(
        get_free_variable(flip_angle),
        get_free_variable(phase),
        scale
    )
    @constraint global_model() res.flip_angle >= 0
    return res
end

flip_angle(instant::InstantRFPulseBlock) = instant.flip_angle
phase(instant::InstantRFPulseBlock) = instant.phase
duration(::InstantRFPulseBlock) = 0.
variables(::Type{<:InstantRFPulseBlock}) = [flip_angle, phase]
effective_time(::InstantRFPulseBlock) = 0.
inverse_bandwidth(::InstantRFPulseBlock) = 0.


function fixed(block::InstantRFPulseBlock)
    return FixedInstantPulse(value(flip_angle(block)), value(phase(block)))
end

end