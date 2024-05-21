function parse_section(section::PulseqSection{:extensions}; kwargs...)
    current_extension = -1
    pre_amble = true
    linked_list = Dict{Int, NTuple{3, Int}}()
    extensions = Dict{Int, PulseqExtension}()
    for line in section.content
        if startswith(line, "extension ")
            pre_amble = false
            (_, name, str_id) = split(line)
            current_extension = int(str_id)
            extensions[current_extension] = PulseqExtension(name, String[])
        elseif pre_ample
            (id, type, ref, next) = int.(split(line))
            linked_list[id] = (type=type, ref=ref, next=next)
        else
            push!(extensions[current_extension], line)
        end
    end

    function get_extension_list(key::Int)
        if iszero(key)
            return Tuple{PulseqExtension, Int}[]
        else
            base = get_extension_list(linked_list[key].next)
            pushfirst!(base, (extensions[linked_list[key].type], linked_list[key].ref))
            return base
        end
    end

    return Dict(key => get_extension_list(key) for key in keys(linked_list))
end