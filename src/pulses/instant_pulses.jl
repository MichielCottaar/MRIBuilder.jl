module InstantPulses
import JuMP: @constraint, @variable, VariableRef, value
import ....Sequences: InstantRFPulse
import ...BuildingBlocks: BuildingBlock, properties, BuildingBlockPlaceholder, set_simple_constraints!, duration, to_mcmr_components
import ...SequenceBuilders: SequenceBuilder, owner_model, start_time
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
properties(::Type{<:InstantRFPulseBlock}) = [flip_angle, phase]


function to_mcmr_components(block::InstantRFPulseBlock)
    return InstantRFPulse(
        time=value(start_time(block)),
        flip_angle=value(flip_angle(block)),
        phase=value(phase(block)),
    )
end

end