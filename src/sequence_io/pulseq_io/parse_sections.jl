"""
Translate between sets of [`PulseqSection`](@ref) objects and [`PulseqSequence`](@ref).
"""
module ParseSections

import ..Types: PulseqSequence, PulseqSection
import ..Components: PulseqComponents
import ..Parsers: parse_section, gen_section

"""
    parse_all_sections(sections)

Parses the sections read from a Pulseq file.

Returns a [`PulseqSequence`](@ref)
"""
function parse_all_sections(sections:: Dict{String, PulseqSection})
    sections = copy(sections)
    all_parts = Dict{Symbol, Any}()
    for name in ["version", "definitions", "delays", "shapes", "rf", "gradients", "trap", "adc", "extensions", "blocks"]
        if name in keys(sections)
            section = pop!(sections, name)
            all_parts[Symbol(name)] = parse_section(section; all_parts...)
        end
    end
    if length(sections) > 0
        @warn "Following sections in pulseq input file are not being used: $(keys(sections))"
    end
    return PulseqSequence(all_parts[:version], all_parts[:definitions], all_parts[:blocks])
end

function gen_all_sections(seq:: PulseqSequence)
    sections = Dict{Symbol, PulseqSection}()
    sections[:version] = gen_section(seq, Val(:version))
    sections[:definitions] = gen_section(seq, Val(:definitions))
    (section, extension_mapping) = gen_section(seq, Val(:extensions))
    sections[:extensions] = section

    comp = PulseqComponents()
    sections[:blocks] = gen_section(seq, comp, Val(:blocks), extension_mapping)
    for symbol in [:rf, :gradients, :trap, :adc, :shapes]
        sections[symbol] = gen_section(comp, Val(symbol))
    end
    return sections
end

end