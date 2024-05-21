struct CompressedPulseqShape
    num :: Int
    samples :: Vector{Float64}
end


function parse_section(section:: PulseqSection{:shapes}; kwargs...)
    current_id = -1
    shapes = Dict{Int, CompressedPulseqShape}()
    for line in section.content
        if startswith(lowercase(line), "shape_id")
            current_id = parse(Int, line[9:end])
            continue
        end
        for text in ("num_uncompressed", "num_samples")
            if startswith(lowercase(line), text)
                @assert current_id != -1
                if current_id in keys(shapes)
                    error("Multiple shapes with the same ID detected.")
                end
                shapes[current_id] = CompressedPulseqShape(parse(Int, line[length(text)+1:end]), Float64[])
            end
        end
        if !startswith(lowercase(line), "num")
            push!(shapes[current_id].samples, parse(Float64, line))
        end
    end
    return Dict(key => uncompress(shape) for (key, shape) in shapes)
end

function uncompress(compressed::CompressedPulseqShape)
    if compressed.num == length(compressed.samples)
        # not actually compressed
        return PulseqShape(compressed.samples)
    end
    amplitudes = [compressed.samples[1]]
    repeating = false
    prev_sample = compressed.samples[1]
    for sample in compressed.samples[2:end]
        if repeating
            for _ in 1:Int(sample)
                push!(amplitudes, prev_sample + amplitudes[end])
            end
            repeating = false
        else
            push!(amplitudes, sample + amplitudes[end])
            if sample == prev_sample
                repeating = true
            else
                prev_sample = sample
            end
        end
    end
    if length(amplitudes) != compressed.num
        error("Uncompressing shape did not produce correct number of elements.")
    end
    return PulseqShape(amplitudes)
end

function gen_section(comp:: PulseqComponents, ::Val{:shapes})
    res = PulseqSection{:shapes}(String[])
    for index in sort([keys(comp.shapes)...])
        shape = comp.shapes[index]
        append!(res.content, [
            "",
            "shape_id $index",
            "num_samples $(length(shape.samples))"
        ])
        for sample in shape.samples
            push!(res.content, string(sample))
        end
    end
    return res
end