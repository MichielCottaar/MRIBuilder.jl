function parse_section(section:: PulseqSection{:definitions}; kwargs...)
    props = parse_pulseq_properties(section.content)
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