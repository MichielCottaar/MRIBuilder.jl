function parse_section(section::PulseqSection{:blocks}; version, rf=Dict(), gradients=Dict(), trap=Dict(), adc=Dict(), extensions=Dict(), delays=nothing, kwargs...)
    all_grad = merge(gradients, trap)
    
    res = Vector{PulseqBlock}()
    for line in section.content
        props = parse_pulseq_dict(
            line,
            [:id, :duration, :rf, :gx, :gy, :gz, :adc, :ext],
            fill(Int, 8),
        )
        if version == v"1.3.1"
            props[:duration] = iszero(props[:duration]) ? 0 : delays[props[:duration]]
        end
        @assert length(res) + 1 == props[:id]
        push!(res, PulseqBlock(
            props[:duration],
            _get_component(props[:rf], rf),
            _get_component(props[:gx], all_grad),
            _get_component(props[:gy], all_grad),
            _get_component(props[:gz], all_grad),
            _get_component(props[:adc], adc),
            iszero(props[:ext]) ? [] : _get_component(props[:ext], extensions),
        ))
    end
    return res
end

function gen_section(seq::PulseqSequence, comp:: PulseqComponents, ::Val{:blocks}, extension_mapping::Dict{Vector, Int})
    res = PulseqSection{:blocks}(String[])

    for (i, block) in enumerate(seq.blocks)
        values = Any[i, block.duration]
        for (search_vec, part) in [
            (comp.pulses, block.rf)
            (comp.grads, block.gx) 
            (comp.grads, block.gy) 
            (comp.grads, block.gz) 
            (comp.adc, block.adc)
        ]
            push!(values, add_components!(comp, search_vec, part))
        end
        if length(block.ext) > 0
            push!(values, extension_mapping[block.ext])
        else
            push!(values, 0)
        end
        push!(res.content, join(string.(values), " "))
    end
    return res
end