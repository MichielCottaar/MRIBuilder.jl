module Printing
import JuMP: value
import Printf: @sprintf
import ..Variables: VariableType, variables, AbstractBlock, VariableNotAvailable, alternative_variables


function _robust_value(possible_number::VariableType)
    try
        return value(possible_number)
    catch
        return nothing
    end
end

function _robust_value(possible_vector::AbstractVector) 
    result = _robust_value.(possible_vector)
    if any(isnothing.(result))
        return nothing
    end
    return result
end

_robust_value(possible_tuple::Tuple) = _robust_value([possible_tuple...])


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
            (ft <: AbstractVector && eltype(ft) <: AbstractBlock)
        )
            continue
        end

        print(io, name, "=", repr(getproperty(block, name)), ", ")
    end

    for fn in values(variables)
        if fn in [fn_alt for (fn_alt, _, _, _) in values(alternative_variables)]
            continue
        end
        try
            numeric_value = _robust_value(fn(block))
            if isnothing(numeric_value)
                continue
            end
            if numeric_value isa AbstractVector
                printed_value = "[" * join(map(v -> @sprintf("%.3g", v), numeric_value), ", ") * "]"
            else
                printed_value = @sprintf("%.3g", numeric_value)
            end
            print(io, "$(nameof(fn))=$(printed_value), ")
        catch e
            if e isa VariableNotAvailable
                continue
            end
            rethrow()
        end
    end
    print(io, ")")
end

end