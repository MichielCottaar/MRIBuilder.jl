module InstantPulses
import JuMP: @constraint, @variable, VariableRef, value
import ...BuildingBlocks: BuildingBlock, properties, BuildingBlockPlaceholder, set_simple_constraints!, duration
import ...SequenceBuilders: SequenceBuilder, owner_model, start_time
import ...ConcreteBlocks: to_concrete_block, AbstractConcreteBlock, Sequence, AbstractSequence
import ..Properties: flip_angle, phase

struct InstantRFPulseBlock <: BuildingBlock
    builder :: SequenceBuilder
    flip_angle :: VariableRef
    phase :: VariableRef
end

InstantRFPulseBlock(; kwargs...) = BuildingBlockPlaceholder{InstantRFPulseBlock}(; kwargs...)
function InstantRFPulseBlock(builder::SequenceBuilder; kwargs...) 
    model = owner_model(builder)
    res = InstantRFPulseBlock(
        builder,
        @variable(model),
        @variable(model)
    )
    @constraint model flip_angle(res) >= 0
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
    builder :: AbstractSequence
    flip_angle :: Number
    phase :: Number
end

function to_concrete_block(builder::Sequence, block::InstantRFPulseBlock)
    return ConcreteInstantRFPulse(builder, value(flip_angle(block)), value(phase(block)))
end

end