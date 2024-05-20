module Parsers
import ..Types: PulseqSection, PulseqSequence, PulseqRFPulse, PulseqGradient, PulseqTrapezoid, PulseqExtensionDefinition, PulseqExtension, PulseqADC, PulseqBlock, PulseqShape
import ..Components: PulseqComponents, add_components!
import ..BasicParsers: parse_pulseq_dict, parse_pulseq_properties

"""
    parse_section(section)

Parses any [`PulseqSection`](@ref) and return the appropriate type.

The opposite is [`gen_section`](@ref).
"""
function parse_section end


"""
    gen_section(sequence, Val(:<title>))

Creates a specific [`PulseqSection`](@ref){<title>} from a part of the PulseqSequence.

This is the opposite of [`parse_section`](@ref)
"""
function gen_section end

_get_component(id::Int, shapes::Dict) = iszero(id) ? nothing : shapes[id]

include("version.jl")
include("definitions.jl")
include("delays.jl")
include("shapes.jl")
include("rf.jl")
include("gradients.jl")
include("trapezoids.jl")
include("adc.jl")
include("extensions.jl")
include("blocks.jl")
end