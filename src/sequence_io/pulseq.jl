module Pulseq
import ...Variables: duration
import ...Scanners: Scanner
import ...Components: GenericPulse, GradientWaveform, ADC
import ...Containers: BuildingBlock, Sequence, Wait
import DataStructures: OrderedDict
import Interpolations: linear_interpolation, Flat
import StaticArrays: SVector

struct PulseqSection
    title :: String
    content :: Vector{String}
end

struct CompressedPulseqShape
    id :: Int
    num :: Int
    samples :: Vector{Float64}
end


"""
    control_points(shape::CompressedPulseqShape)

Returns a tuple with:
- vector of times
- vector of amplitudes
"""
function control_points(pulseq::CompressedPulseqShape)
    compressed = length(pulseq.samples) != pulseq.num
    if !compressed
        times = range(0, 1, length=pulseq.num)
        return (times, pulseq.samples)
    end
    times = [zero(Float64)]
    amplitudes = [pulseq.samples[1]]
    repeating = false
    time_norm = 1 / (pulseq.num - 1)
    prev_sample = pulseq.samples[1]
    prev_applied = true
    for sample in pulseq.samples[2:end]
        if repeating
            nrepeats = Int(sample) + 2
            if prev_applied
                nrepeats -= 1
            end
            push!(times, times[end] + nrepeats * time_norm)
            push!(amplitudes, amplitudes[end] + nrepeats * prev_sample)
            repeating = false
            prev_applied = true
        elseif sample == prev_sample
            repeating = true
        else
            if !prev_applied
                push!(times, times[end] + time_norm)
                push!(amplitudes, amplitudes[end] + prev_sample)
            end
            prev_sample = sample
            prev_applied = false
        end
    end
    if !prev_applied
        push!(times, times[end] + time_norm)
        push!(amplitudes, amplitudes[end] + prev_sample)
    end
    return (times, amplitudes)
end

"""
    read_pulseq(filename; scanner=Scanner(B0=B0), B0=3., TR=<sequence duration>)

Reads a sequence from a pulseq file (http://pulseq.github.io/).
Pulseq files can be produced using matlab (http://pulseq.github.io/) or python (https://pypulseq.readthedocs.io/en/master/).
"""
function read_pulseq(filename; kwargs...)
    keywords = open(read_pulseq_sections, filename)
    name = basename(filename)
    if endswith(name, ".seq")
        name = name[1:end-4]
    end
    build_sequence(; name=Symbol(name), kwargs..., keywords...)
end

function read_pulseq_sections(io::IO)
    sections = PulseqSection[]
    version = nothing
    title = ""
    for line in readlines(io)
        line = strip(line)
        if length(line) == 0 || line[1] == '#'
            continue  # ignore comments
        end
        if line[1] == '[' && line[end] == ']'
            if title == "VERSION"
                version = parse_pulseq_section(sections[end]).second
            end
            # new section starts
            title = line[2:end-1]
            push!(sections, PulseqSection(title, String[]))
        elseif length(sections) > 0
            push!(sections[end].content, line)
        else
            error("Content found in pulseq file before first section")
        end
    end
    Dict(filter(pair -> !isnothing(pair), parse_pulseq_section.(sections, version))...)
end

function parse_pulseq_ordered_dict(strings::Vector{<:AbstractString}, names, dtypes)
    as_dict = OrderedDict{Int, NamedTuple{Tuple(names), Tuple{dtypes...}}}()
    for line in strings
        parts = split(line)
        @assert length(parts) == length(names)
        values = parse.(dtypes, split(line))
        @assert names[1] == :id
        as_dict[values[1]] = (; zip(names, values)...)
    end
    return as_dict
end

function parse_pulseq_properties(strings::Vector{<:AbstractString})
    result = Dict{String, Any}()
    for s in strings
        (name, value) = split(s, limit=2)
        result[name] = value
    end
    return result
end

section_headers = Dict(
    "BLOCKS" => ([:id, :duration, :rf, :gx, :gy, :gz, :adc, :ext], fill(Int, 8)),
    ("BLOCKS", v"1.3.1") => ([:id, :delay, :rf, :gx, :gy, :gz, :adc, :ext], fill(Int, 8)),
    ("DELAYS", v"1.3.1") => ([:id, :delay], [Int, Int]),
    ("RF", v"1.3.1") => ([:id, :amp, :mag_id, :phase_id, :delay, :freq, :phase], [Int, Float64, Int, Int, Int, Float64, Float64]),
    "RF" => ([:id, :amp, :mag_id, :phase_id, :time_id, :delay, :freq, :phase], [Int, Float64, Int, Int, Int, Int, Float64, Float64]),
    ("GRADIENTS", v"1.3.1") => ([:id, :amp, :shape_id, :delay], [Int, Float64, Int, Int]),
    "GRADIENTS" => ([:id, :amp, :shape_id, :time_id, :delay], [Int, Float64, Int, Int, Int]),
    "TRAP" => ([:id, :amp, :rise, :flat, :fall, :delay], [Int, Float64, Int, Int, Int, Int]),
    "ADC" => ([:id, :num, :dwell, :delay, :freq, :phase], [Int, Int, Float64, Int, Float64, Float64]),
)

function parse_pulseq_section(section::PulseqSection, version=nothing)
    if section.title == "VERSION"
        props = parse_pulseq_properties(section.content)
        result = VersionNumber(
            parse(Int, props["major"]),
            parse(Int, props["minor"]),
            parse(Int, props["revision"]),
        )
    elseif section.title == "DEFINITIONS"
        result = parse_pulseq_properties(section.content)
    elseif (section.title, version) in keys(section_headers)
        result = parse_pulseq_ordered_dict(section.content, section_headers[(section.title, version)]...)
    elseif section.title in keys(section_headers)
        result = parse_pulseq_ordered_dict(section.content, section_headers[section.title]...)
    elseif section.title == "EXTENSION"
        println("Ignoring all extensions in pulseq")
    elseif section.title == "SHAPES"
        current_id = -1
        shapes = CompressedPulseqShape[]
        for line in section.content
            if length(line) > 8 && lowercase(line[1:8]) == "shape_id"
                current_id = parse(Int, line[9:end])
                continue
            end
            for text in ("num_uncompressed", "num_samples")
                if startswith(lowercase(line), text)
                    @assert current_id != -1
                    push!(shapes, CompressedPulseqShape(current_id, parse(Int, line[length(text)+1:end]), Float64[]))
                    current_id = -1
                    break
                end
            end
            if !startswith(lowercase(line), "num")
                @assert current_id == -1
                push!(shapes[end].samples, parse(Float64, line))
            end
        end
        result = Dict(
            [s.id => (s.num, control_points(s)...) for s in shapes]...
        )
    elseif section.title in ["SIGNATURE"]
        # silently ignore these sections
        return nothing
    else
        error("Unrecognised pulseq section: $(section.title)")
    end
    return Symbol(lowercase(section.title)) => result
end

function align_in_time(pairs...)
    interps = [linear_interpolation(time, ampl; extrapolation_bc=Flat()) for (time, ampl) in pairs]
    all_times = sort(unique(vcat([time for (time, _) in pairs]...)))
    return (all_times, [interp.(all_times) for interp in interps]...)
end


function build_sequence(; scanner=nothing, B0=3., TR=nothing, definitions, version, blocks, rf=nothing, gradients=nothing, trap=nothing, delays=nothing, shapes=nothing, adc=nothing, name=:Pulseq)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    if version == v"1.4.0"
        # load raster times (converting seconds to milliseconds)
        convert = key -> parse(Float64, definitions[key]) * Float64(1e3)
        gradient_raster = convert("GradientRasterTime")
        rf_raster = convert("RadiofrequencyRasterTime")
        adc_raster = convert("AdcRasterTime")
        block_duration_raster = convert("BlockDurationRaster")
    elseif version == v"1.3.1"
        gradient_raster = rf_raster = block_duration_raster = Float64(1e-3) # 1 microsecond as default raster
        adc_raster = Float64(1e-6) # ADC dwell time is in ns by default
    else
        error("Can only load pulseq files with versions v1.3.1 and v1.4.0, not $(version)")
    end

    full_blocks = BuildingBlock[]

    for block in values(blocks)
        block_duration = 0.
        events = []
        if !iszero(block.rf)
            proc = rf[block.rf]
            (num, times_shape, amplitudes_shape) = shapes[proc.mag_id]
            block_duration = max(num * rf_raster + proc.delay * 1e-3, block_duration)
            ampl_times = times_shape .* (num * rf_raster)
            ampl_size = amplitudes_shape .* (proc.amp * 1e-3)
            if iszero(proc.phase_id)
                phase_times = ampl_times
                phase_size = ampl_times .* 0
            else
                (num, times_shape, phase_shape) = shapes[proc.phase_id]
                block_duration = max(num * rf_raster + proc.delay * 1e-3, block_duration)
                phase_times = times_shape .* (num * rf_raster)
                phase_size = rad2deg.(phase_shape) .+ phase_times .* (proc.freq * 1e-3 * 360) .+ rad2deg(proc.phase)
            end
            if version != v"1.3.1" && !iszero(proc.time_id)
                (num, time_shape) = shapes[proc.time_id]
                times = time_shape.amplitudes .* rf_raster
                ampl = ampl_size
                phase = phase_size
            else
                (times, ampl, phase) = align_in_time((ampl_times, ampl_size), (phase_times, phase_size))
            end
            push!(events, (proc.delay * 1e-3, GenericPulse(times, ampl, phase)))
        end
        if !iszero(block.adc)
            proc = adc[block.adc]
            push!(events, (proc.delay, ADC(proc.num, proc.dwell * 1e-6, proc.dwell * proc.num * 1e-6 / 2, 0.)))
            block_duration = max(proc.delay * 1e-3 + proc.dwell * proc.num * 1e-6, block_duration)
        end
        grad_shapes = []
        for symbol_grad in [:gx, :gy, :gz]
            grad_id = getfield(block, symbol_grad)
            if iszero(grad_id)
                push!(grad_shapes, (Float64[], Float64[]))
                continue
            end
            if !isnothing(gradients) && grad_id in keys(gradients)
                proc = gradients[grad_id]
                start_time = proc.delay * 1e-3
                (num, grad_shape) = shapes[proc.shape_id]

                if version != v"1.3.1" && !iszero(proc.time_id)
                    (num, time_shape) = shapes[proc.time_id]
                    times = time_shape.amplitudes .* gradient_raster .+ start_time
                else
                    times = grad_shape.times .* (num * gradient_raster) .+ start_time
                end
                push!(grad_shapes, (times, grad_shape.amplitudes .* (proc.amp * 1e-9)))
                block_duration = max(proc.delay * 1e-3 + num * gradient_raster, block_duration)
            elseif !isnothing(trap) && grad_id in keys(trap)
                proc = trap[grad_id]
                start_time = proc.delay * 1e-3
                times = (cumsum([0, proc.rise, proc.flat, proc.fall]) .* 1e-3) .+ start_time
                push!(grad_shapes, (times, [0, proc.amp * 1e-9, proc.amp * 1e-9, 0]))
                block_duration = max((proc.delay + proc.rise + proc.flat + proc.fall) * 1e-3, block_duration)
            else
                error("Gradient ID $grad_id not found in either of [GRADIENTS] or [TRAP] sections")
            end
        end
        if version == v"1.3.1"
            if !iszero(block.delay)
                actual_duration = max(block_duration, delays[block.delay].delay * 1e-3)
            else
                actual_duration = block_duration
            end
        else
            actual_duration = block.duration * block_duration_raster
            @assert actual_duration >= block_duration
        end
        function extend_grad!(pair)
            (times, amplitudes) = pair
            if iszero(length(times)) || !iszero(times[1])
                pushfirst!(times, 0.)
                pushfirst!(amplitudes, 0.)
            end
            if times[end] != actual_duration
                push!(times, actual_duration)
                push!(amplitudes, 0.)
            end
        end
        extend_grad!.(grad_shapes)
        arrs = align_in_time(grad_shapes...)
        waveform = [(t, SVector{3, Float64}(gx, gy, gz)) for (t, gx, gy, gz) in zip(arrs...)]
        push!(full_blocks, BuildingBlock(waveform, events))
    end
    full_block_duration = sum(duration.(full_blocks))
    if isnothing(TR)
        total_duration = parse(Float64, get(definitions, "TotalDuration", "0")) * 1000
        TR = iszero(total_duration) ? full_block_duration : total_duration
    end
    if TR > full_block_duration
        push!(full_blocks, Wait(TR - full_block_duration))
    elseif TR < full_block_duration
        @warn "Given TR ($TR) is shorter than the total duration of all the building blocks ($full_block_duration). Ignoring TR..."
    end
    Sequence(full_blocks; scanner=scanner, name=name)
end

end