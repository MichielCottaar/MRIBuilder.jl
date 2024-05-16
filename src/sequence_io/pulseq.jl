module Pulseq


"""
    read_pulseq(IO)

Reads a sequence from a pulseq file (http://pulseq.github.io/).
Pulseq files can be produced using matlab (http://pulseq.github.io/) or python (https://pypulseq.readthedocs.io/en/master/).
"""
function read_pulseq(io::IO)
    sections = parse_pulseq_sections(io)
    return parse_all_sections(sections)
end

"""
    write_pulseq(IO, sequence)

Writes a sequence to an output IO file.
"""
function write_pulseq(io::IO, sequence::PulseqSequence)
    if sequence.version < v"1.4"
        error("Can only write to pulseq version 1.4 or later.")
    end
    sections = gen_all_sections(sequence)
    write_section.(io, sections)
end


"""
    parse_pulseq_dict(line, names, dtypes)

Parse a line of integers/floats with known names and dtypes.

This is useful to parse most of the columnar data in Pulseq, such as in BLOCKS, RF, GRADIENTS, etc.
"""
function parse_pulseq_dict(line, names, dtypes)
    parts = split(line)
    @assert length(parts) == length(names)
    values = parse.(dtypes, split(line))
    @assert names[1] == :id
    return Dict(Symbol(name) => value for (name, value) in zip(names, values))
end

"""
    parse_pulseq_properties(lines)

Parse any `pulseq` section formatted as:
```
<name> <value>
<name2> <value2>
...
```

This includes the VERSION, DEFINITIONS, and part of the SHAPES
"""
function parse_pulseq_properties(strings::Vector{<:AbstractString})
    result = Dict{String, Any}()
    for s in strings
        (name, value) = split(s, limit=2)
        result[name] = _parse_value(value)
    end
    return result
end

function _parse_value(value::AbstractString)
    for t in (Int, Float64)
        try
            return parse(t, value)
        catch
        end
    end
    for t in (Int, Float64)
        try
            return parse.(t, split(value))
        catch
        end
    end
    return value
end



"""
    PulseqSection(:<title>)(lines)

Represents a section in the pulseq file format.
"""
struct PulseqSection{T}
    content :: Vector{String}
end

"""
    parse_pulseq_sections(io)

Reads a Pulseq file into a dictionary of [`PulseqSection`](@ref) objects.
"""
function parse_pulseq_sections(io::IO)
    sections = Dict{String, PulseqSection}()
    current_title = ""
    for line in readlines(io)
        line = strip(line)
        if length(line) == 0 || line[1] == '#'
            continue  # ignore comments
        end
        if line[1] == '[' && line[end] == ']'
            # new section starts
            current_title = lowercase(line[2:end-1])
            sections[current_title] = PulseqSection{Symbol(current_title)}(String[])
        elseif length(current_title) > 0
            push!(sections[current_title].content, line)
        else
            error("Content found in pulseq file before first section")
        end
    end
    return sections
end

# Writing pulseq files

struct PulseqShape
    samples :: Vector{Float64}
end

struct PulseqRFPulse
    amplitude :: Float64
    magnitude :: PulseqShape
    phase :: PulseqShape
    time :: Union{Nothing, PulseqShape}
    delay :: Int
    frequency :: Float64
    phase_offset :: Float64
end

abstract type AnyPulseqGradient end
struct PulseqGradient <: AnyPulseqGradient
    amplitude :: Float64
    shape :: PulseqShape
    time :: Union{Nothing, PulseqShape}
    delay :: Int
end

struct PulseqTrapezoid <:AnyPulseqGradient
    amplitude :: Float64
    rise :: Int
    flat :: Int
    fall :: Int
    delay :: Int
end

struct PulseqADC
    nsamples :: Int
    dwell_time :: Float64
    delay :: Int
    frequency :: Float64
    phase :: Float64
end

struct PulseqExtension
    name :: String
    content :: Vector{String}
end

struct PulseqBlock
    duration :: Int
    rf :: Union{Nothing, PulseqRFPulse}
    gx :: Union{Nothing, AnyPulseqGradient}
    gy :: Union{Nothing, AnyPulseqGradient}
    gz :: Union{Nothing, AnyPulseqGradient}
    adc :: Union{Nothing, PulseqADC}
    ext :: Vector{<:Tuple{PulseqExtension, Int}}
end

struct PulseqSequence
    version:: VersionNumber
    definitions:: NamedTuple
    blocks:: Vector{PulseqBlock}
end



# I/O all sections
function parse_all_sections(sections:: Dict{String, PulseqSection})
    sections = copy(sections)
    all_parts = Dict{Symbol, Any}()
    for name in ["version", "definitions", "delays", "shapes", "rf", "gradients", "trap", "adc", "extensions", "blocks"]
        if name in keys(sections)
            section = pop!(sections, name)
            all_parts[Symbol(name)] = parse_section(section; all_parts...)
        end
    end
    if length(sections) > 0
        @warn "Following sections in pulseq input file are not being used: $(keys(sections))"
    end
    return PulseqSequence(all_parts[:version], all_parts[:definitions], all_parts[:blocks])
end

function gen_all_sections(seq:: PulseqSequence)
    sections = [gen_section(seq, Val(symbol)) for symbol in [:version, :definitions, :shapes, :rf, :gradients, :trap, :adc, :extensions, :blocks]]
    return [section for section in sections if length(section.content) > 0]
end

# Version I/O

"""
    parse_section(section)

Parses any [`PulseqSection`](@ref) and return the appropriate type.

The opposite is [`gen_section`](@ref).
"""
function parse_section(section:: PulseqSection{:version}; kwargs...)
    props = parse_pulseq_properties(section.content)
    return VersionNumber(props["major"], props["minor"], props["revision"])
end

"""
    gen_section(sequence, Val(:<title>))

Creates a specific [`PulseqSection`](@ref){<title>} from a part of the PulseqSequence.

This is the opposite of [`parse_section`](@ref)
"""
function gen_section(seq::PulseqSequence, ::Val{:version})
    version = seq.version
    return PulseqSection{:version}([
        "major $(version.major)",
        "minor $(version.minor)",
        "revision $(version.patch)",
    ])
end

# Definitions I/O

function parse_section(section:: PulseqSection{:definitions}; kwargs...)
    props = parse_pulseq_properties(section.content)
    return NamedTuple(Symbol(key) => value for (key, value) in props)
end

_to_string_value(value::AbstractString) = value
_to_string_value(value::Number) = string(value)
_to_string_value(value::Vector) = join(_to_string_value.(value), " ")

function gen_section(seq:: PulseqSequence, ::Val{:definitions})
    definitions = seq.definitions
    return PulseqSection{:definitions}(["$(key) $(_to_string_value(value))" for (key, value) in pairs(definitions)])
end


# Shapes I/O
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
                push!(amplitudes, prev_sample)
            end
            repeating = false
        else
            push!(amplitudes, sample)
            if sample == prev_sample
                repeating = true
            end
        end
    end
    return PulseqShape(amplitudes)
end

same_shape(shape1::PulseqShape, shape2::PulseqShape) = length(shape1.samples) == length(shape2.samples) && all(shape1.samples .≈ shape2.samples)

"""
    shapes(sequence::PulseqSequence)

Return all the unique shapes used in the MR sequence.
"""
function shapes(sequence:: PulseqSequence)
    result = PulseqShape[]
    for block in sequence.blocks
        for new_shape in _all_shapes(block)
            if !any(s -> same_shape(s, new_shape), result)
                push!(result, new_shape)
            end
        end
    end
    return result
end

_all_shapes(block::PulseqBlock) = vcat(_all_shapes.([block.gx, block.gy, block.gz, block.rf])...)
_all_shapes(block::PulseqRFPulse) = vcat(_all_shapes.([block.magnitude, block.phase, block.time])...)
_all_shapes(block::PulseqGradient) = vcat(_all_shapes.([block.shape, block.time])...)
_all_shapes(shape::PulseqShape) = PulseqShape[shape]
_all_shapes(::PulseqTrapezoid) = PulseqShape[]
_all_shapes(::Nothing) = PulseqShape[]


function gen_section(seq:: PulseqSequence, ::Val{:shapes})
    res = PulseqSection{:shapes}(String[])
    for (index, shape) in enumerate(shapes(seq))
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

# Components I/O
_get_component(id::Int, shapes::Dict) = iszero(id) ? nothing : shapes[id]

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

function _get_shape_id(shapes, shape)
    for (i, s) in enumerate(shapes)
        if same_shape(s, shape)
            return i
        end
    end
    error("Shape not found.")
end


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
            iszero(props[:ext]) ? Tuple{PulseqExtension, Int}[] : _get_component(props[:ext], extensions),
        ))
    end
    return res
end

end