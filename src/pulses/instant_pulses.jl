module InstantPulses
import JuMP: @constraint, @variable, VariableRef, value, Model
import ...BuildingBlocks: BuildingBlock
import ...ConcreteBlocks: to_concrete_block, AbstractConcreteBlock
import ...Variables: flip_angle, phase, start_time, variables, duration, get_free_variable, VariableType
import ...GlobalModel: @global_model_constructor

struct InstantRFPulseBlock <: BuildingBlock
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
    set_simple_constraints!(model, res, kwargs)
    return res
end

flip_angle(instant::InstantRFPulseBlock) = instant.flip_angle
phase(instant::InstantRFPulseBlock) = instant.phase
duration(instant::InstantRFPulseBlock) = 0.
variables(::Type{<:InstantRFPulseBlock}) = [flip_angle, phase]


"""
    ConcreteInstantRFPulse(flip_angle, phase)

Instantaneous RF pulse with no free variables.

See [`InstantRFPulseBlock`](@ref) for a version where [`flip_angle`](@ref) and [`phase`](@ref) are variables.
"""
struct ConcreteInstantRFPulse <: AbstractConcreteBlock
    flip_angle :: Number
    phase :: Number
end

function to_concrete_block(block::InstantRFPulseBlock)
    return ConcreteInstantRFPulse(value(flip_angle(block)), value(phase(block)))
end

end