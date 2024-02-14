module WaitBlocks
import JuMP: @constraint
import ...BuildSequences: get_global_model
import ...Variables: get_free_variable, VariableType, duration
import ...Components: NoGradient
import ..BaseBuildingBlocks: BaseBuildingBlock

struct Wait <: BaseBuildingBlock
    duration :: VariableType
    function Wait(var)
        res = new(get_free_variable(var))
        if !(res.duration isa Number)
            @constraint get_global_model() res.duration >= 0
        end
    end
end

duration(wb::Wait) = wb.duration
Base.keys(::Wait) = (:empty)
Base.getindex(wb::Wait) = NoGradient{1}(wb.duration)

end