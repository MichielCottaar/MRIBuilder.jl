function parse_section(section::PulseqSection{:adc}; kwargs...)
    result = Dict{Int, PulseqADC}()
    for line in section.content
        props = parse_pulseq_dict(
            line,
            [:id, :num, :dwell, :delay, :freq, :phase],
            [Int, Int, Float64, Int, Float64, Float64],
        )
        result[props[:id]] = PulseqADC(
            props[:num],
            props[:dwell],
            props[:delay],
            props[:freq],
            props[:phase],
        ) 
    end
    return result
end

function gen_section(comp:: PulseqComponents, ::Val{:adc})
    res = PulseqSection{:adc}(String[])
    for i in sort([keys(comp.adc)...])
        adc = comp.adc[i]
        values = string.(Any[
            i,
            adc.num,
            adc.dwell,
            adc.delay,
            adc.frequency,
            adc.phase,
        ])
        push!(res.content, join(values, " "))
    end
    return res
end