function parse_section(section::PulseqSection{:gradients}; shapes::Dict{Int, PulseqShape}, version::VersionNumber, kwargs...)
    result = Dict{Int, PulseqGradient}()
    for line in section.content
        if version == v"1.3.1"
            props = parse_pulseq_dict(
                line,
                [:id, :amp, :shape_id, :delay],
                [Int, Float64, Int, Int]
            )
            props[:time_id] = 0
        else
            props = parse_pulseq_dict(
                line,
                [:id, :amp, :shape_id, :time_id, :delay],
                [Int, Float64, Int, Int, Int]
            )
        end
        result[props[:id]] = PulseqGradient(
            props[:amp],
            _get_component(props[:shape_id], shapes),
            _get_component(props[:time_id], shapes),
            props[:delay]
        ) 
    end
    return result
end

function gen_section(comp:: PulseqComponents, ::Val{:gradients})
    res = PulseqSection{:gradients}(String[])
    for i in sort([keys(comp.grads)...])
        grad = comp.grads[i]
        if !(grad isa PulseqGradient)
            continue
        end
        values = string.(Any[
            i,
            grad.amplitude,
            add_components!(comp, grad.shape),
            add_components!(comp, grad.time),
            grad.delay,
        ])
        push!(res.content, join(values, " "))
    end
    return res
end