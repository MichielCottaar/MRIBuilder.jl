function parse_section(section::PulseqSection{:trap}; kwargs...)
    result = Dict{Int, PulseqTrapezoid}()
    for line in section.content
        props = parse_pulseq_dict(
            line,
            [:id, :amp, :rise, :flat, :fall, :delay],
            [Int, Float64, Int, Int, Int, Int],
        )
        result[props[:id]] = PulseqTrapezoid(
            props[:amp],
            props[:rise],
            props[:flat],
            props[:fall],
            props[:delay],
        ) 
    end
    return result
end

function gen_section(comp:: PulseqComponents, ::Val{:trap})
    res = PulseqSection{:trap}(String[])
    for i in sort([keys(comp.grads)...])
        grad = comp.grads[i]
        if !(grad isa PulseqTrapezoid)
            continue
        end
        values = string.(Any[
            i,
            grad.amplitude,
            grad.rise,
            grad.flat,
            grad.fall,
            grad.delay,
        ])
        push!(res.content, join(values, " "))
    end
    return res
end
