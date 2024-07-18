module Printing
import JuMP: value
import Printf: @sprintf
import ..Variables: VariableType, variables, AbstractBlock, Variable

function _robust_value(possible_number::VariableType)
    try
        return round(value(possible_number), sigdigits=3)
    catch
        return nothing
    end
end

function _robust_value(possible_vector::AbstractArray) 
    result = _robust_value.(possible_vector)
    if any(isnothing.(result))
        return nothing
    end
    return result
end

_robust_value(possible_tuple::Tuple) = _robust_value([possible_tuple...])
_robust_value(possible_tuple::NamedTuple) = NamedTuple(k => _robust_value(v) for (k, v) in pairs(possible_tuple))


function Base.show(io::IO, block::AbstractBlock)
    print(io, nameof(typeof(block)), "(")
    for name in propertynames(block)
        ft = fieldtype(typeof(block), name)
        if (
            ft == VariableType ||
            (ft <: AbstractVector && eltype(ft) == VariableType) ||
            string(name)[1] == '_'
        )
            continue
        end
        if (
            ft <: Union{Nothing, AbstractBlock} ||
            (ft <: AbstractVector && eltype(ft) <: AbstractBlock) ||
            (ft <: AbstractVector && eltype(ft) <: Pair)
        )
            continue
        end

        print(io, name, "=", repr(getproperty(block, name)), ", ")
    end

    for fn_name in names(variables, all=true)
        fn = getproperty(variables, fn_name)
        if !(fn isa Variable)
            continue
        end

        try
            numeric_value = _robust_value(fn(block))
            if isnothing(numeric_value)
                continue
            end
            as_string = string(numeric_value)
            if length(as_string) > 50
                continue
            end
            print(io, "$(fn_name)=$(as_string), ")
        catch e
            continue
        end
    end
    print(io, ")")
end

end