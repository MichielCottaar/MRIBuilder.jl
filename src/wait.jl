module Wait
import JuMP: @constraint, @variable, VariableRef, value
import ..Variables: VariableType, variables, duration, get_free_variable
import ..BuildingBlocks: BuildingBlock, to_block
import ..BuildSequences: global_model
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
    duration :: VariableType
    function WaitBlock(duration=nothing)
        res = new(
            get_free_variable(duration),
        )
        @constraint model res.duration >= 0
        return res
    end
end

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

end