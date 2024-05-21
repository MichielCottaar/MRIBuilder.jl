"""Define IO for [`PulseqSection`](@ref)."""
module SectionsIO
import ..Types: PulseqSection

"""
    parse_pulseq_sections(io)

Reads a Pulseq file into a dictionary of [`PulseqSection`](@ref) objects.
"""
function parse_pulseq_sections(io::IO)
    sections = Dict{String, PulseqSection}()
    current_title = ""
    for line in readlines(io)
        line = strip(line)
        if length(line) == 0 || line[1] == '#'
            continue  # ignore comments
        end
        if line[1] == '[' && line[end] == ']'
            # new section starts
            current_title = lowercase(line[2:end-1])
            sections[current_title] = PulseqSection{Symbol(current_title)}(String[])
        elseif length(current_title) > 0
            push!(sections[current_title].content, line)
        else
            error("Content found in pulseq file before first section")
        end
    end
    return sections
end

"""
    write_pulseq_section(io, section::PulseqSection)

Writes a Pulseq `section` to the `IO`.
"""
function write_pulseq_section(io::IO, section::PulseqSection{T}) where {T}
    title = uppercase(string(T))
    write(io, "[$title]\n")
    for line in section.content
        if iszero(length(line)) || line[end] != '\n'
            line = line * '\n'
        end
        write(io, line)
    end
    write(io, "\n")
    write(io, "\n")
end

end