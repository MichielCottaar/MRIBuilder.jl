module InstantPulses
import JuMP: @constraint, @variable, VariableRef, value, Model
import ...BuildingBlocks: RFPulseBlock, fixed
import ...Variables: flip_angle, phase, start_time, variables, duration, get_free_variable, VariableType, effective_time, inverse_bandwidth
import ...BuildSequences: @global_model_constructor
import ..FixedPulses: FixedInstantPulse

struct InstantRFPulseBlock <: RFPulseBlock
    model :: Model
    flip_angle :: VariableType
    phase :: VariableType
end

@global_model_constructor InstantRFPulseBlock

function InstantRFPulseBlock(model::Model; flip_angle=nothing, phase=nothing) 
    res = InstantRFPulseBlock(
        model,
        get_free_variable(model, flip_angle),
        get_free_variable(model, phase)
    )
    @constraint model res.flip_angle >= 0
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