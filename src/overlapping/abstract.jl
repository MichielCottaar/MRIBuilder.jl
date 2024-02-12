module Abstract

import StaticArrays: SVector
import ...BuildingBlocks: BuildingBlock, ContainerBlock, RFPulseBlock, get_children_indices, VariableNotAvailable
import ...Wait: WaitBlock
import ...Readouts: InstantReadout
import ...Gradients: GradientBlock, split_gradient, ConstantGradientBlock, ChangingGradientBlock
import ...Variables: duration, start_time, qvec, bmat_gradient, gradient_strength, slew_rate, VariableType, effective_time, gradient_strength, slew_rate
import ...Variables: flip_angle, amplitude, phase, frequency, bandwidth, inverse_bandwidth, N_left, N_right, slice_thickness, all_variables_symbols

"""
Parent type for any gradient waveform that might overlap with RF pulses and/or readouts.

New waveforms should sub-type from this one.
At the very least they should implement:
- [`waveform_sequence`](@ref): return the subsequent [`GradientBlock`](@ref) and [`WaitBlock`](@ref) objects describing the gradient waveform. To get the wavform in a more useful format use [`waveform`](@ref).
- [`interruptions`](@ref): return a vector of all [`RFPulseBlock`](@ref) and [`InstantReadout`](@ref) with information on which part of the waveform they interrupt.
"""
abstract type AbstractOverlapping <: ContainerBlock end

"""
    SingleInterrupted(grad_like_block, interruptions)

Represents a single part of the waveform within [`GenericOverlapping`](@ref).
"""
struct SingleInterrupted{T<:Union{WaitBlock, GradientBlock}}
    block :: T
    interruptions :: Vector{NamedTuple{(:index, :time, :block), Tuple{Int64, Float64, <:Union{RFPulseBlock, InstantReadout}}}}
end

"""
    waveform_sequence(overlapping)

Returns the gradient waveform as a sequence of [`WaitBlock`](@ref), [`ChangingGradientBlock`](@ref), and [`ConstantGradientBlock`](@ref) objects.

This should be implemented for any sub-type of [`AbstractOverlapping`](@ref).
"""
function waveform_sequence end


"""
    waveform(overlapping)

Returns the gradient waveform of a [`AbstractOverlapping`](@ref) as a sequence of control points.
Each control point is stored as a tuple with the time in ms and the gradient as a length-3 vector.
The gradient is linearly interpolated between these points (see [`waveform_sequence`](@ref)).
"""
function waveform(ao::AbstractOverlapping)
    result = Tuple{VariableType, SVector{3, VariableType}}[(0, zero(SVector{3, Float64}))]
    for block in waveform_sequence(ao)
        new_time = result[end][1] + max(duration(block), 0)
        prev_grad = result[end][2]
        if block isa WaitBlock
            @assert all(abs.(gradient_strength(block)) <= 1e-12) "$(typeof(ao)) inserts WaitBlock before the gradient is zero. This is probably caused by an improper implementation of this gradient waveform."
            push!(result, (new_time, prev_grad))
        elseif block isa ConstantGradientBlock
            @assert all(gradient_strength(block) .≈ prev_grad) "$(typeof(ao)) inserts ConstantGradientBlock that does not match previous gradient strength. This is probably caused by an improper implementation of this gradient waveform."
            push!(result, (new_time, prev_grad))
        elseif block isa ChangingGradientBlock
            @assert all(block.gradient_strength_start .≈ prev_grad) "$(typeof(ao)) inserts ChangingGradientBlock that does not match previous gradient strength. This is probably caused by an improper implementation of this gradient waveform."
            push!(result, (new_time, prev_grad .+ slew_rate(block) .* duration(block)))
        else
            error("Unrecognised block type in gradient waveform: $(typeof(ao)).")
        end
    end
    if duration(ao) > result[end][1]
        push!(result, (duration(ao), zero(SVector{3, Float64})))
    end
    return result
end


"""
    interruptions(overlapping)

Returns the interruptions as a sequence named tuples with 3 elements:
- `index`: which object in the [`waveform`](@ref) should be interrupted.
- `time`: time since the start of the [`waveform`](@ref) object that the interruptions takes place ([`effective_time`](@ref) for RF pulses).
- `object`: [`RFPulseBlock`](@ref) or [`InstantReadout`](@ref).

This should be implemented for any sub-type of [`AbstractOverlapping`](@ref).
"""
function interruptions end

duration(ao::AbstractOverlapping) = sum(duration.(waveform_sequence(ao)))
duration(single::SingleInterrupted) = duration(single.block)

pulses(ao::AbstractOverlapping) = [pulse for (_, _, pulse) in interruptions(ao) if pulse isa RFPulseBlock]
readouts(ao::AbstractOverlapping) = [readout for (_, _, readout) in interruptions(ao) if readout isa InstantReadout]

# For any overlapping building blocks with a single RF pulse, get information about that RF pulse.
for (func_symbol, _) in [Dict(all_variables_symbols)[:pulse]..., :effective_time => ""]
    if func_symbol == :slice_thickness
        continue
    end
    @eval function $func_symbol(ao::AbstractOverlapping, args...; kwargs...)
        single_pulse = pulses(ao)
        if iszero(length(single_pulse))
            throw(VariableNotAvailable(typeof(ao), $func_symbol))
        elseif length(single_pulse) > 1
            error("Building block has multipl pulses, so cannot compute ", $func_symbol, ".")
        end
        $func_symbol(single_pulse[1], args...; kwargs...)
    end
end

# Acting as a valid container
get_children_indices(ao::AbstractOverlapping) = 1:length(waveform_sequence(ao))
function Base.getindex(ao::AbstractOverlapping, index::Integer)
    grad_like_block = waveform_sequence(ao)[index]
    inter = interruptions(ao)
    first_interrupt = findfirst(int -> int[2] == index, inter)
    if isnothing(first_interrupt)
        return grad_like_block
    else
        last_interrupt = findlast(int -> int[2] == index, inter)
        return SingleInterrupted(grad_like_block, inter[first_interrupt:last_interrupt])
    end
end

start_time(ao::AbstractOverlapping, index::Integer) = sum(duration.(waveform_sequence(ao)[1:index]))

function get_part(si::SingleInterrupted, first::Union{Nothing, Number}, last::Union{Nothing, Number})
    if isnothing(first)
        if isnothing(last)
            part = si.block
        else
            (part, _) = split_gradient(si, si.interruptions[last][1])
        end
    else
        tfirst = si.interruptions[first][1]
        if isnothing(last)
            (_, part) = split_gradient(si, tfirst)
        else
            (_, part, _) = split_gradient(si, tfirst, so.interruptions[last][1])
        end
    end
    return part
end

function get_parts(ao::AbstractOverlapping, first::Union{Nothing, Integer}, last::Union{Nothing, Integer})
    inter = interruptions(ao)
    form = waveform_sequence(ao)
    whole_start_index = isnothing(first) ? 0 : inter[first].index
    whole_final_index = isnothing(last) ? length(form) + 1 : inter[last].index
    if whole_start_index == whole_final_index
        return [split_gradient(form[whole_start_index], inter[first].time, inter[last].time)]
    end

    parts = form[whole_start_index+1:whole_final_index-1]
    if !isnothing(first)
        pushfirst!(parts, split_gradient(form[whole_start_index], inter[first].time)[2])
    end
    if !isnothing(last)
        push!(parts, split_gradient(form[whole_final_index], inter[last].time)[1])
    end
    return parts
end


"""
    qvec(overlapping[, first_interruption, last_interruption])

Computes the area under the curve for the gradient waveform in [`AbstractOverlapping`](@ref).

If `first_interruption` is set to something else than `nothing`, only the gradient waveform after this RF pulse/Readout will be considered.
Similarly, if `last_interruption` is set to something else than `nothing`, only the gradient waveform up to this RF pulse/Readout will be considered.
"""
function qvec(ao::AbstractOverlapping, index1::Union{Nothing, Integer}, index2::Union{Nothing, Integer})
    @assert isnothing(index1) || isnothing(index2) || index2 >= index1
    if (index1 isa Number) && (index1 == index2)
        return zeros(3)
    end
    sum(qvec.(get_parts(ao, index1, index2)); init=zero(SVector{3, Float64}))
end
qvec(ao::AbstractOverlapping) = qvec(ao, nothing, nothing)

function bmat_gradient(ao::AbstractOverlapping, qstart, index1::Union{Nothing, Integer}, index2::Union{Nothing, Integer})
    @assert isnothing(index1) || isnothing(index2) || index2 >= index1
    if (index1 isa Number) && (index1 == index2)
        return zeros(3, 3)
    end
    result = Matrix{VariableType}(zeros(3, 3))
    qcurrent = Vector{VariableType}(qstart)

    for part in get_parts(ao, index1, index2)
        result = result .+ bmat_gradient(part, qcurrent)
        qcurrent = qcurrent .+ qvec(part, qcurrent)
    end
    return result
end
bmat_gradient(ao::AbstractOverlapping, qstart) = bmat_gradient(ao, qstart, nothing, nothing)

end