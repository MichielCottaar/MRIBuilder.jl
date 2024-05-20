module BasicParsers

"""
    parse_pulseq_dict(line, names, dtypes)

Parse a line of integers/floats with known names and dtypes.

This is useful to parse most of the columnar data in Pulseq, such as in BLOCKS, RF, GRADIENTS, etc.
"""
function parse_pulseq_dict(line, names, dtypes)
    parts = split(line)
    @assert length(parts) == length(names)
    values = parse.(dtypes, split(line))
    @assert names[1] == :id
    return Dict(Symbol(name) => value for (name, value) in zip(names, values))
end

"""
    parse_pulseq_properties(lines)

Parse any `pulseq` section formatted as:
```
<name> <value>
<name2> <value2>
...
```

This includes the VERSION, DEFINITIONS, and part of the SHAPES
"""
function parse_pulseq_properties(strings::Vector{<:AbstractString})
    result = Dict{String, Any}()
    for s in strings
        (name, value) = split(s, limit=2)
        result[name] = parse_value(value)
    end
    return result
end

"""
    parse_value(string)

Tries to value the `string` as a number of sequence of numbers.
If that does not work, the `string` is returned.
"""
function parse_value(value::AbstractString)
    for t in (Int, Float64)
        try
            return parse(t, value)
        catch
        end
    end
    for t in (Int, Float64)
        try
            return parse.(t, split(value))
        catch
        end
    end
    return value
end


end