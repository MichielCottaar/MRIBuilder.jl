function parse_section(section::PulseqSection{:delays}; version, kwargs...)
    if version > v"1.3.1"
        error("Did not expect a [DELAYS] section in pulseq file with version $(version)")
    end
    delays = Dict{Int, Int}()
    for line in section.content
        (id, delay) = parse.(Int, split(line))
        delays[id] = delay
    end
    return delays
end