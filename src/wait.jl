module Wait
import JuMP: Model, @constraint, @variable, VariableRef, owner_model
import ..BuildingBlocks: BuildingBlock, duration, properties, apply_simple_constraint!, BuildingBlockPlaceholder, to_mcmr_components
import ..SequenceBuilders: SequenceBuilder, to_block
import ...Scanners: Scanner

"""
    WaitBlock(duration)

An empty [`BuildingBlock`](@ref) of given `duration` (in ms).

Duration can be set to one of:
- numeric value to fix it
- `:min` to minimise its value given any external constraints
- `:max` to maximise its value given any external constraints
- `nothing` to make it fully determined by external constraints and objectives
"""
struct WaitBlock <: BuildingBlock
    builder :: SequenceBuilder
    duration :: VariableRef
end

function WaitBlock(builder::SequenceBuilder, duration_constraint=nothing)
    model = owner_model(builder)
    res = WaitBlock(builder, @variable(model))
    @constraint model duration(res) >= 0
    if !isnothing(duration_constraint)
        apply_simple_constraint!(model, duration(res), duration_constraint)
    end
    return res
end

WaitBlock(duration_constraint=nothing) = BuildingBlockPlaceholder{WaitBlock}(duration_constraint)

to_block(builder::SequenceBuilder, time::Union{Number, Symbol, Nothing, Val{:min}, Val{:max}}) = WaitBlock(builder, time)


properties(::Type{WaitBlock}) = [duration]

duration(wb::WaitBlock) = wb.duration

scanner_constraints!(::Model, ::WaitBlock, ::Scanner) = nothing

to_mcmr_components(::WaitBlock) = []

end