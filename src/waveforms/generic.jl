module Generic
import ...BuildingBlocks: BuildingBlock, ContainerBlock, RFPulseBlock
import ...Wait: WaitBlock
import ...Readouts: InstantReadout
import ...Gradients: GradientBlock, split_gradient
import ...Variables: duration, start_time, qvec, bmat_gradient, gradient_strength, slew_rate, VariableType
import ...Variables: flip_angle, amplitude, phase, frequency, bandwidth, inverse_bandwidth, N_left, N_right, slice_thickness, all_variables_symbols

"""
Parent type for all objects, where gradients, RF pulses, and/or readouts might overlap with each other.

Do not sub-type from this type directly, use [`SpecificWaveform`](@ref) instead.
"""
abstract type AbstractWaveform <: ContainerBlock end

"""
Parent type for any specific gradient waveform that might overlap with RF pulses and/or readouts.

New waveforms should sub-type from this one.
At the very least they should implement a conversion to [`GenericWaveform`](@ref).
"""
abstract type SpecificWaveform <: AbstractWaveform end


"""
    GenericWaveform(duration, waveform, interruptions)

Generic description that can capture any overlapping gradients, RF pulses, and/or readouts.

Interruptions are stores as tuples with 3 fields:
- `index`: which part of the `waveform` is being interrupted (cannot be a free variable).
- `time`: time of the interruption relative to the start of the `waveform` part (between 0 and the length of this part). This can be a variable.
- `block`: [`RFPulseBlock`](@ref) or [`InstantReadout`](@ref)
"""
struct GenericWaveform <: AbstractWaveform
    duration :: VariableType
    waveform :: Vector{Union{WaitBlock, GradientBlock}}
    interruptions :: Vector{NamedTuple{(:index, :time, :block), (Int64, <:VariableType, <:Union{RFPulseBlock, InstantReadout})}}
end


"""
    SingleInterrupted(grad_like_block, interruptions)

Represents a single part of the waveform within [`GenericWaveform`](@ref).
"""
struct SingleInterrupted{T<:Union{WaitBlock, GradientBlock}} <: AbstractWaveform
    block :: T
    interruptions :: Vector{NamedTuple{(:index, :time, :block), (Int64, Float64, <:Union{RFPulseBlock, InstantReadout})}}
end


"""
    GenericWaveform(other_overlapping_block)

Converts any sub-type of [`SpecificWaveform`](@ref) into a [`GenericWaveform`](@ref).
This needs to be defined for every sub-type.
"""
function GenericWaveform end

waveform(go::GenericWaveform) = go.waveform
interruptions(go::GenericWaveform) = go.interruptions

duration(go::GenericWaveform) = go.duration
duration(ao::AbstractWaveform) = sum(duration.(waveform(ao)))
duration(single::SingleInterrupted) = duration(single.block)

pulses(ao::AbstractWaveform) = [pulse for (_, _, pulse) in ao.interruptions if pulse isa RFPulseBlock]
readouts(ao::AbstractWaveform) = [readout for (_, _, readout) in ao.interruptions if readout isa InstantReadout]

# For any overlapping building blocks with a single RF pulse, get information about that RF pulse.
for (func_symbol, _) in Dict(all_variables_symbols)[:pulse]
    @eval function $func_symbol(go::AbstractWaveform, args...; kwargs...)
        single_pulse = pulses(go)
        if iszero(length(single_pulse))
            error("Building block does not contain any RF pulse, so cannot compute $func_symbol.")
        elseif lenth(single_pulse) > 1
            error("Building block has multipl pulses, so cannot compute $func_symbol.")
        end
        $func_symbol(single_pulse[1], args...; kwargs...)
    end
end

# Acting as a valid container
get_children_indices(go::Union{SpecificWaveform, GenericWaveform}) = 1:length(waveform(go))
function Base.getindex(go::GenericWaveform, index::Integer)
    grad_like_block = waveform(go)[index]
    first_interrupt = findfirst(int => int[2] == index, interruptions(go))
    if isnothing(first_interrupt)
        return grad_like_block
    else
        last_interrupt = findlast(int => int[2] == index, go.interruptions)
        return SingleInterrupted(grad_like_block, go.interruptions[first_interrupt:last_interrupt])
    end
end

function start_time(go::Union{SpecificWaveform, GenericWaveform}, index::Integer)
    sum(duration.(waveform(go)[1:index]))
end

function get_part(so::SpecificWaveform, first::Union{Nothing, Number}, last::Union{Nothing, Number})
    if isnothing(first)
        if isnothing(last)
            part = so.block
        else
            (part, _) = split_gradient(so, so.interruptions[last][1])
        end
    else
        tfirst = so.interruptions[first][1]
        if isnothing(last)
            (_, part) = split_gradient(so, tfirst)
        else
            (_, part, _) = split_gradient(so, tfirst, so.interruptions[last][1])
        end
    end
    return part
end

function get_parts(go::AbstractWaveform, first::Union{Nothing, Integer}, last::Union{Nothing, Integer})
    inter = interruptions(go)
    form = waveform(go)
    whole_start_index = isnothing(first) ? 0 : inter[first].index
    whole_final_index = isnothing(last) ? length(form) + 1 : inter[first].index
    if whole_start_index == whole_final_index
        return [split_gradient(form[whole_start_index].block, inter[first].time, inter[last].time)]
    end

    parts = form[whole_start_index+1:whole_final_index-1]
    if !isnothing(first)
        pushfirst!(parts, split_gradient(form[whole_start_index].block, inter[first].time)[2])
    end
    if !isnothing(last)
        push!(parts, split_gradient(form[whole_final_index].block, inter[last].time)[1])
    end
    return parts
end


# Computing gradient properties
for func in (:qvec, :bmat_gradient)
    @eval $func(ao::AbstractWaveform, args..., kwargs...) = $func(GenericWaveform(ao), args...; kwargs...)
end

"""
    qvec(overlapping[, first_interruption, last_interruption])

Computes the area under the curve for the gradient waveform in [`AbstractOverlapping`](@ref).

If `first_interruption` is set to something else than `nothing`, only the gradient waveform after this RF pulse/Readout will be considered.
Similarly, if `last_interruption` is set to something else than `nothing`, only the gradient waveform up to this RF pulse/Readout will be considered.
"""
function qvec(go::GenericOverlapping, index1::Union{Nothing, Integer}, index2::Union{Nothing, Integer})
    @assert index2 >= index1
    if (index1 isa Number) && (index1 == index2)
        return zeros(3)
    end
    sum(qvec.(get_parts(go, index1, index2)))
end
qvec(ao::AbstractOverlapping) = qvec(ao, nothing, nothing)

function bmat_gradient(go::GenericOverlapping, qstart, index1::Union{Nothing, Integer}, index2::Union{Nothing, Integer})
    @assert index2 >= index1
    if (index1 isa Number) && (index1 == index2)
        return zeros(3, 3)
    end
    result = Matrix{VariableType}(zeros(3, 3))
    qcurrent = Vector{VariableType}(qstart)

    for part in get_parts(go, index1, index2)
        result = result .+ bmat_gradient(part, qcurrent)
        qcurrent = qcurrent .+ qvec(part, qcurrent)
    end
    return result
end
bmat_gradient(ao::AbstractOverlapping, qstart) = bmat_gradient(ao, qstart, nothing, nothing)


end