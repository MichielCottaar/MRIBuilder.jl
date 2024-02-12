module BaseBuildingBlocks
import ...ContainerBlocks: ContainerBlock, GradientWaveform, EventComponent, NoGradientBlock, ChangingGradientBlock, ConstantGradientBlock, split_gradient
import ...Components: BaseComponent
import ...Variables: qvec, bmat_gradient

"""
Basic BuildingBlock, which can consist of a gradient waveforms with any number of RF pulses/readouts overlaid

Main interface:
- iteration will give the gradient waveforms interspersed by RF pulses/readouts.
    - Indiviual indices can be accessed using `keys(building_block)`
- [`waveform_sequence`](@ref) returns just the gradient waveform as a sequence of [`GradientWaveform`](@ref) objects.
- [`waveform`](@ref) returns just the gradient waveform as a sequence of (time, gradient_strength) tuples.
- [`events`](@ref) returns the RF pulses and readouts.
- [`qvec`](@ref) returns area under curve for (part of) the gradient waveform.

Sub-types need to implement:
- `Base.keys`: returns sequence of keys to all the components
- `Base.getindex`: returns the actual component for each key
"""
abstract type BaseBuildingBlock <: ContainerBlock end

# Iterator interface
Base.length(c::BaseBuildingBlock) = length(keys(c))
Base.eltype(::Type{<:BaseBuildingBlock}) = BaseComponent
Base.iter(c::BaseBuildingBlock) = iter(c, 1)
Base.iter(c::BaseBuildingBlock, index::Integer) = length(c) < index ? nothing : (c[keys(c)[index]], index + 1)


"""
    events(building_block)

Returns just the non-gradient (i.e., RF pulses/readouts) events as a sequence of [`EventComponent`](@ref) objects (with their keys).
"""
function waveform_sequence(bb::BaseBuildingBlock)
    return [(key, bb[key]) for key in keys(bb) if bb[key] isa EventComponent]
end

"""
    waveform_sequence(building_block)

Returns just the gradient waveform as a sequence of [`GradientWaveform`](@ref) objects (with their keys).
"""
function waveform_sequence(bb::BaseBuildingBlock)
    return [(key, bb[key]) for key in keys(bb) if bb[key] isa GradientWaveform]
end

function ndim_grad(bb::BaseBuildingBlock)
    g = waveform_sequence(bb)
    if izero(length(g))
        return 0
    end
    for N in (1, 3)
        if all(isa.(g, GradientWaveform{N}))
            return N
        end
    end
    error("$(typeof(bb)) contains both 1D and 3D gradient waveforms.")
end


"""
    waveform(building_block)

Returns the gradient waveform of any [`BaseBuildingBlock`](@ref) as a sequence of control points.

Each control point is stored as a tuple with the time in ms and the gradient as a length-3 vector.
The gradient is linearly interpolated between these points (see [`waveform_sequence`](@ref)).
"""
function waveform(bb::BaseBuildingBlock)
    ndim = ndim_grad(bb)
    if ndim == 3
        result = Tuple{VariableType, SVector{3, VariableType}}[(0., zero(SVector{3, Float64}))]
    elseif ndim == 1
        result = Tuple{VariableType, SVector{3, VariableType}}[(0., 0.)]
    else
        return []
    end
    for (_, block) in waveform_sequence(bb)
        new_time = result[end][1] + max(duration(block), 0)
        prev_grad = result[end][2]
        if block isa NoGradientBlock
            @assert all(abs.(prev_grad) <= 1e-12) "$(typeof(bb)) inserts NoGradientBlock before the gradient is zero. This is probably caused by an improper implementation of this BuildingBlock."
            push!(result, (new_time, prev_grad))
        elseif block isa ConstantGradientBlock
            @assert all(gradient_strength(block) .≈ prev_grad) "$(typeof(bb)) inserts ConstantGradientBlock that does not match previous gradient strength. This is probably caused by an improper implementation of this BuildingBlock."
            push!(result, (new_time, prev_grad))
        elseif block isa ChangingGradientBlock
            @assert all(block.gradient_strength_start .≈ prev_grad) "$(typeof(bb)) inserts ChangingGradientBlock that does not match previous gradient strength. This is probably caused by an improper implementation of this BuildingBlock."
            push!(result, (new_time, prev_grad .+ slew_rate(block) .* duration(block)))
        else
            error("Unrecognised block type in BuildingBlock: $(typeof(bb)).")
        end
    end
    @assert all(abs.(result[end][2]) <= 1e-12) "$(typeof(bb)) does not end up with a gradient of zero. This is probably caused by an improper implementation of this BuildingBlock."
    return result
end


# Pathway support
"""
    waveform_sequence(building_block, first, last)

Gets the sequence of [`GradientWaveform`](@ref) from the event with key `first` till the event with key `last`.

Setting `first` to nothing indicates to start from the beginning of the `building_block`.
Similarly, setting `last` to nothing indicates to continue till the end of the `building_block`.
"""
function waveform_seqeuence(bb::BaseBuildingBlock, first, last)
    started = isnothing(first)
    current_grad = current_start = nothing
    parts = GradientWaveform[]
    for key in bb
        if bb[key] isa GradientWaveform
            if started
                push!(parts, isnothing(current_start) ? current_grad : split_gradient(current_grad, current_start)[2])
            end
            current_grad = bb[key]
            current_start = nothing
        end
        if key == first
            @assert !started
            started = true
            current_started = effective_time(bb[key])
        end
        if key == last
            @assert started
            if isnothing(current_started)
                push!(parts, split_gradient(current_grad, effective_time(bb[key]))[1])
            else
                push!(parts, split_gradient(current_grad, current_started, effective_time(bb[key]))[2])
            end
        end
    end
    return parts
end


"""
    qvec(overlapping[, first_event, last_event])

Computes the area under the curve for the gradient waveform in [`BaseBuildingBlock`](@ref).

If `first_event` is set to something else than `nothing`, only the gradient waveform after this RF pulse/Readout will be considered.
Similarly, if `last_event` is set to something else than `nothing`, only the gradient waveform up to this RF pulse/Readout will be considered.
"""
function qvec(bb::BaseBuildingBlock, index1::Union{Nothing, Integer}, index2::Union{Nothing, Integer})
    @assert isnothing(index1) || isnothing(index2) || index2 >= index1
    if (index1 isa Number) && (index1 == index2)
        return zeros(3)
    end
    sum(qvec.(get_parts(bb, index1, index2)); init=zero(SVector{3, Float64}))
end
qvec(bb::BaseBuildingBlock) = qvec(bb, nothing, nothing)

function bmat_gradient(bb::BaseBuildingBlock, qstart, index1::Union{Nothing, Integer}, index2::Union{Nothing, Integer})
    @assert isnothing(index1) || isnothing(index2) || index2 >= index1
    if (index1 isa Number) && (index1 == index2)
        return zeros(3, 3)
    end
    result = Matrix{VariableType}(zeros(3, 3))
    qcurrent = Vector{VariableType}(qstart)

    for part in get_parts(bb, index1, index2)
        result = result .+ bmat_gradient(part, qcurrent)
        qcurrent = qcurrent .+ qvec(part, qcurrent)
    end
    return result
end
bmat_gradient(bb::BaseBuildingBlock, qstart) = bmat_gradient(bb, qstart, nothing, nothing)

end

end