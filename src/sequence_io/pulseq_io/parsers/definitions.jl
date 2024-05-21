function parse_section(section:: PulseqSection{:definitions}; version, kwargs...)
    props = parse_pulseq_properties(section.content)
    if version == v"1.3.1"
        for name in [
            "BlockDurationRaster", 
            "RadiofrequencyRasterTime", 
            "GradientRasterTime",
        ]
            props[name] = 1e-6
        end
        props["AdcRasterTime"] = 1e-9
    end
    return NamedTuple(Symbol(key) => value for (key, value) in props)
end

_to_string_value(value::AbstractString) = value
_to_string_value(value::Symbol) = string(value)
_to_string_value(value::Number) = string(value)
_to_string_value(value::Vector) = join(_to_string_value.(value), " ")

function gen_section(seq:: PulseqSequence, ::Val{:definitions})
    definitions = seq.definitions
    return PulseqSection{:definitions}(["$(key) $(_to_string_value(value))" for (key, value) in pairs(definitions)])
end