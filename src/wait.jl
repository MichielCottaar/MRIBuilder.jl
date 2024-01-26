module Wait
import JuMP: Model, @constraint, @variable, VariableRef, owner_model, value
import ..Variables: VariableType, variables, duration, get_free_variable
import ..BuildingBlocks: BuildingBlock
import ..ConcreteBlocks: to_concrete_block, ConcreteBlock
import ..GlobalModel: @global_model_constructor
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
    model :: Model
    duration :: VariableType
end

@global_model_constructor WaitBlock

function WaitBlock(model::Model, duration=nothing)
    res = WaitBlock(
        model, 
        get_free_variable(model, duration),
        VariableType[]
    )
    @constraint model duration(res) >= 0
    return res
end

to_block(time::Union{VariableType, Symbol, Nothing, Val{:min}, Val{:max}}) = WaitBlock(time)


variables(::Type{WaitBlock}) = [duration]

duration(wb::WaitBlock) = wb.duration

to_concrete_block(wb::WaitBlock) = ConcreteBlock(value(duration(wb)))

end