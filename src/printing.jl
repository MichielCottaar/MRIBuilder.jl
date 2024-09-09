module Printing
import JuMP: value, AbstractJuMPScalar
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
_robust_value(other) = other


function Base.show(io::IO, block::AbstractBlock)
    show_block(io, block, 0)
end

function show_block(io::IO, block::AbstractBlock, nspaces::Int64)
    print(io, nameof(typeof(block)), "(\n")
    for name in propertynames(block)
        value = _robust_value(getproperty(block, name))
        if (
            value isa AbstractJuMPScalar ||
            (value isa AbstractVector{<:AbstractJuMPScalar}) ||
            string(name)[1] == '_' ||
            isnothing(value) ||
            (value isa AbstractVector{<:AbstractBlock}) ||
            (value isa AbstractVector{<:Pair})
        )
            continue
        end

        print(io, repeat(' ', nspaces + 2), name, "=")
        if value isa AbstractBlock
            show_block(io, value, nspaces + 2)
        else
            print(io, repr(value))
        end
        print(io, ",\n")
    end
    print(io, ")")
end

end