import ..Extensions: parse_extension, get_extension_name, add_extension_definition!

function parse_section(section::PulseqSection{:extensions}; kwargs...)
    current_extension = -1
    pre_amble = true
    linked_list = Dict{Int, NamedTuple}()
    extensions = Dict{Int, PulseqExtensionDefinition}()
    for line in section.content
        if startswith(line, "extension ")
            pre_amble = false
            (_, name, str_id) = split(line)
            current_extension = parse(Int, str_id)
            extensions[current_extension] = PulseqExtensionDefinition{Symbol(name)}(String[])
        elseif pre_amble
            (id, type, ref, next) = parse.(Int, split(line))
            linked_list[id] = (type=type, ref=ref, next=next)
        else
            push!(extensions[current_extension].content, line)
        end
    end

    extension_mappers = Dict(key => parse_extension(ext) for (key, ext) in extensions)

    function get_extension_list(key::Int)
        if iszero(key)
            return []
        else
            base = get_extension_list(linked_list[key].next)
            pushfirst!(base, extension_mappers[linked_list[key].type][linked_list[key].ref])
            return base
        end
    end

    return Dict(key => get_extension_list(key) for key in keys(linked_list))
end

function gen_section(seq:: PulseqSequence, ::Val{:extensions})
    definitions = Dict{Symbol, PulseqExtensionDefinition}()
    extensions_ref_id = Dict{Any, Int}()

    all_ext_vec = [block.ext for block in seq.blocks if length(block.ext) > 0]
    for ext_vec in all_ext_vec
        for ext in ext_vec
            name = get_extension_name(ext)
            if !(name in keys(definitions))
                definitions[name] = PulseqExtensionDefinition{name}(String[])
            end
            extensions_ref_id[ext] = add_extension_definition!(definitions[name].content, ext)
        end
    end

    definitions_order = [keys(definitions)...]

    extension_vec_id = Dict{Vector, Int}()

    labels = Tuple{Int, Int, Int}[]
    function add_extension_vec!(ext_vec::Vector)
        if length(ext_vec) == 0
            return 0
        end
        next_id = add_extension_vec!(ext_vec[2:end])
        ext = ext_vec[1]
        type_id = findfirst(isequal(get_extension_name(ext)), definitions_order)
        ref_id = extensions_ref_id[ext]
        all_ids = (type_id, ref_id, next_id)
        if all_ids in labels
            extension_vec_id[ext_vec] = findfirst(isequal(all_ids), labels)
        else
            push!(labels, all_ids)
            extension_vec_id[ext_vec] = length(labels)
        end
    end
    add_extension_vec!.(all_ext_vec)

    content = String[]
    for (id0, (id1, id2, id3)) in enumerate(labels)
        push!(content, "$id0 $id1 $id2 $id3")
    end
    push!(content, "")

    for (index, name) in enumerate(definitions_order)
        def = definitions[name]
        push!(content, "extension $(string(name)) $(string(index))")
        append!(content, def.content)
        push!(content, "")
    end

    return PulseqSection{:extensions}(content), extension_vec_id
end