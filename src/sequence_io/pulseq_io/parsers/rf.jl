function parse_section(section::PulseqSection{:rf}; shapes::Dict{Int, PulseqShape}, version, kwargs...)
    result = Dict{Int, PulseqRFPulse}()
    for line in section.content
        if version == v"1.3.1"
            props = parse_pulseq_dict(
                line,
                [:id, :amp, :mag_id, :phase_id, :delay, :freq, :phase], 
                [Int, Float64, Int, Int, Int, Float64, Float64],
            )
            props[:time_id] = 0
        else
            props = parse_pulseq_dict(
                line,
                [:id, :amp, :mag_id, :phase_id, :time_id, :delay, :freq, :phase], 
                [Int, Float64, Int, Int, Int, Int, Float64, Float64],
            )
        end
        result[props[:id]] = PulseqRFPulse(
            props[:amp],
            _get_component(props[:mag_id], shapes),
            _get_component(props[:phase_id], shapes),
            _get_component(props[:time_id], shapes),
            props[:delay],
            props[:freq],
            props[:phase],
        )
    end
    return result
end

function gen_section(comp:: PulseqComponents, ::Val{:rf})
    res = PulseqSection{:rf}(String[])
    for i in sort([keys(comp.pulses)...])
        pulse = comp.pulses[i]
        values = string.(Any[
            i,
            pulse.amplitude,
            add_components!(comp, pulse.magnitude),
            add_components!(comp, pulse.phase),
            add_components!(comp, pulse.time),
            pulse.delay,
            pulse.frequency,
            pulse.phase_offset
        ])
        push!(res.content, join(values, " "))
    end
    return res
end