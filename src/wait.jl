module Wait
import JuMP: Model, @constraint, @variable, VariableRef, owner_model, value
import ..Variables: VariableType, variables, duration, get_free_variable
import ..BuildingBlocks: BuildingBlock, to_block
import ..ConcreteBlocks: to_concrete_block, ConcreteBlock
import ..BuildSequences: @global_model_constructor
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
    function WaitBlock(model::Model, duration=nothing)
        res = new(
            model, 
            get_free_variable(model, duration),
        )
        @constraint model res.duration >= 0
        return res
    end
end

@global_model_constructor WaitBlock

"""
    to_block(JuMP variable/:min or :max/nothing/number)

Converts object into a [`WaitBlock`](@ref).
- if a Number or a JuMP variable the wait time will match that value.
- if `:min`, the wait time will be minimised in the final sequence (limited by any external constraints).
- if `:max`, the wait time will be maximised in the final sequence (limited by any external constraints).
- if `nothing`, the wait time is only constrained by external factors.
"""
to_block(time::Union{VariableType, Symbol, Nothing, Val{:min}, Val{:max}}) = WaitBlock(time)


variables(::Type{WaitBlock}) = [duration]

duration(wb::WaitBlock) = wb.duration

to_concrete_block(wb::WaitBlock) = ConcreteBlock(value(duration(wb)))

end