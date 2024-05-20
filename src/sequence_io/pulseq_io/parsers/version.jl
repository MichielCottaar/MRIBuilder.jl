function parse_section(section:: PulseqSection{:version}; kwargs...)
    props = parse_pulseq_properties(section.content)
    return VersionNumber(props["major"], props["minor"], props["revision"])
end

function gen_section(seq::PulseqSequence, ::Val{:version})
    version = seq.version
    return PulseqSection{:version}([
        "major $(version.major)",
        "minor $(version.minor)",
        "revision $(version.patch)",
    ])
end