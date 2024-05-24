module Printing
import JuMP: value
import Printf: @sprintf
import ..Variables: VariableType, variables, AbstractBlock 

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

    for fn in values(variables)
        try
            numeric_value = _robust_value(fn(block))
            if isnothing(numeric_value)
                continue
            end
            print(io, "$(nameof(fn))=$(numeric_value), ")
        catch
            continue
        end
    end
    print(io, ")")
end

end