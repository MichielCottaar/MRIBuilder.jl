"""
Defines [`BaseBuildingBlock`](@ref), [`BuildingBlock`](@ref) and [`Wait`](@ref).
"""
module BuildingBlocks
import LinearAlgebra: norm
import JuMP: @constraint
import StaticArrays: SVector
import ..Abstract: ContainerBlock, start_time
import ...BuildSequences: global_model
import ...Components: BaseComponent, GradientWaveform, EventComponent, NoGradient, ChangingGradient, ConstantGradient, split_gradient, DelayedEvent, RFPulseComponent, ReadoutComponent
import ...Variables: qval, bmat_gradient, effective_time, get_free_variable, qval3, slew_rate, gradient_strength
import ...Variables: VariableType, duration, make_generic, get_pulse, get_readout, scanner_constraints!

"""
Basic BuildingBlock, which can consist of a gradient waveforms with any number of RF pulses/readouts overlaid

Main interface:
- iteration will give the gradient waveforms interspersed by RF pulses/readouts.
    - Indiviual indices can be accessed using `keys(building_block)`
- [`waveform_sequence`](@ref) returns just the gradient waveform as a sequence of [`GradientWaveform`](@ref) objects.
- [`waveform`](@ref) returns just the gradient waveform as a sequence of (time, gradient_strength) tuples.
- [`events`](@ref) returns the RF pulses and readouts.
- [`qval`](@ref) returns area under curve for (part of) the gradient waveform.

Sub-types need to implement:
- `Base.keys`: returns sequence of keys to all the components.
- `Base.getindex`: returns the actual component for each key.
"""
abstract type BaseBuildingBlock <: ContainerBlock end

# Iterator interface
Base.length(c::BaseBuildingBlock) = length(keys(c))
Base.eltype(::Type{<:BaseBuildingBlock}) = BaseComponent
Base.iterate(c::BaseBuildingBlock) = Base.iterate(c, 1)
Base.iterate(c::BaseBuildingBlock, index::Integer) = length(c) < index ? nothing : (c[keys(c)[index]], index + 1)
Base.getindex(bb::BaseBuildingBlock, s::Symbol) = bb[Val(s)]


"""
    events(building_block)

Returns just the non-gradient (i.e., RF pulses/readouts) events as a sequence of [`EventComponent`](@ref) objects (with their keys).
"""
function events(bb::BaseBuildingBlock)
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
    g = [ws[2] for ws in waveform_sequence(bb)]
    if iszero(length(g))
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
        result = Tuple{VariableType, VariableType}[(0., 0.)]
    else
        return []
    end
    for (_, block) in waveform_sequence(bb)
        new_time = result[end][1] + max(duration(block), 0)
        prev_grad = result[end][2]
        if block isa NoGradient
            @assert all(abs.(prev_grad) <= 1e-12) "$(typeof(bb)) inserts NoGradient before the gradient is zero. This is probably caused by an improper implementation of this BuildingBlock."
            push!(result, (new_time, prev_grad))
        elseif block isa ConstantGradient
            @assert all(gradient_strength(block) .≈ prev_grad) "$(typeof(bb)) inserts ConstantGradient that does not match previous gradient strength. This is probably caused by an improper implementation of this BuildingBlock."
            push!(result, (new_time, prev_grad))
        elseif block isa ChangingGradient
            @assert all(block.gradient_strength_start .≈ prev_grad) "$(typeof(bb)) inserts ChangingGradient that does not match previous gradient strength. This is probably caused by an improper implementation of this BuildingBlock."
            push!(result, (new_time, prev_grad .+ slew_rate(block) .* duration(block)))
        else
            error("Unrecognised block type in BuildingBlock: $(typeof(bb)).")
        end
    end
    @assert all(abs.(result[end][2]) .<= 1e-12) "$(typeof(bb)) does not end up with a gradient of zero. This is probably caused by an improper implementation of this BuildingBlock."
    return result
end

equal_key(i1::Val, i2) = i1 == Val(i2)
equal_key(i1, i2::Val) = Val(i1) == i2
equal_key(i1::Val, i2::Val) = i1 == i2
equal_key(i1, i2) = i1 == i2

function start_time(building_block::BaseBuildingBlock, index)
    time = 0.
    prev_time = 0.
    for key in keys(building_block)
        if equal_key(key, index)
            return prev_time
        end
        component = building_block[key]
        if component isa GradientWaveform
            prev_time = time
            time += duration(component)
        end
    end
    error("Building block with index '$index' not found")
end

duration(bb::BaseBuildingBlock) = sum([duration(wv) for (_, wv) in waveform_sequence(bb)])

# Pathway support
"""
    waveform_sequence(building_block, first, last)

Gets the sequence of [`GradientWaveform`](@ref) from the event with key `first` till the event with key `last`.

Setting `first` to nothing indicates to start from the beginning of the `building_block`.
Similarly, setting `last` to nothing indicates to continue till the end of the `building_block`.
"""
function waveform_sequence(bb::BaseBuildingBlock, first, last)
    started = isnothing(first)
    current_grad_key = current_start = nothing
    parts = Tuple{Any, GradientWaveform}[]
    for key in keys(bb)
        if bb[key] isa GradientWaveform
            if started && !isnothing(current_grad_key)
                push!(parts, (current_grad_key, isnothing(current_start) ? bb[current_grad_key] : split_gradient(bb[current_grad_key], current_start)[2]))
            end
            current_grad_key = key
            current_start = nothing
        end
        if equal_key(key, first)
            @assert !started
            started = true
            current_start = effective_time(bb[key])
        end
        if equal_key(key, last)
            @assert started
            if isnothing(current_start)
                push!(parts, (current_grad_key, split_gradient(bb[current_grad_key], effective_time(bb[key]))[1]))
            else
                push!(parts, (current_grad_key, split_gradient(bb[current_grad_key], current_start, effective_time(bb[key]))[2]))
            end
            return parts
        end
    end
    if !started
        error("Starting index of $first not recognised.")
    end
    if !isnothing(last)
        error("Final index of $last not recognised.")
    end
    push!(parts, (current_grad_key, isnothing(current_start) ? bb[current_grad_key] : split_gradient(bb[current_grad_key], current_start)[2]))
    return parts
end


"""
    qval(overlapping[, first_event, last_event])

Computes the area under the curve for the gradient waveform in [`BaseBuildingBlock`](@ref).

If `first_event` is set to something else than `nothing`, only the gradient waveform after this RF pulse/Readout will be considered.
Similarly, if `last_event` is set to something else than `nothing`, only the gradient waveform up to this RF pulse/Readout will be considered.
"""
function qval(bb::BaseBuildingBlock, index1, index2)
    if (!isnothing(index1)) && (index1 == index2)
        return 0.
    end
    sum([qval(wv) for (_, wv) in waveform_sequence(bb, index1, index2)]; init=0.)
end
qval(bb::BaseBuildingBlock) = qval(bb, nothing, nothing)

function bmat_gradient(bb::BaseBuildingBlock, qstart, index1, index2)
    if (!isnothing(index1)) && (index1 == index2)
        return zeros(3, 3)
    end
    result = Matrix{VariableType}(zeros(3, 3))
    qcurrent = Vector{VariableType}(qstart)

    for (_, part) in waveform_sequence(bb, index1, index2)
        result = result .+ bmat_gradient(part, qcurrent)
        qcurrent = qcurrent .+ qval3(part, qcurrent)
    end
    return result
end
bmat_gradient(bb::BaseBuildingBlock, qstart) = bmat_gradient(bb, qstart, nothing, nothing)

"""
    BuildingBlock(waveform, events; duration=nothing, orientation=nothing, group)

Generic [`BaseBuildingBlock`](@ref) that can capture any overlapping gradients, RF pulses, and/or readouts.
The gradients cannot contain any free variables.

Scanner constraints are automatically applied.

## Arguments
- `waveform`: Sequence of 2-element tuples with (time, (Gx, Gy, Gz)). If `orientation` is set then the tuple is expected to look like (time, G). This cannot contain any free variables.
- `events`: Sequence of 2-element tuples with (index, pulse/readout). The start time of the pulse/readout at the start of the gradient waveform element with index `index` (use [`DelayedEvent`](@ref) to make this earlier or later).
- `duration`: duration of this `BuildingBlock`. If not set then it will be assumed to be the time of the last element in `waveform`.
- `orientation`: orientation of the gradients in the waveform. If not set, then the full gradient vector should be given explicitly.
- `group`: group of the gradient waveform
"""
struct BuildingBlock <: BaseBuildingBlock
    parts :: Vector{<:BaseComponent}
end

function BuildingBlock(waveform::AbstractVector, events::AbstractVector; orientation=nothing, group=nothing)
    events = Any[events...]
    waveform = Any[waveform...]
    ndim = isnothing(orientation) ? 1 : 3
    zero_grad = isnothing(orientation) ? zeros(3) : 0.
    if length(waveform) == 0 || waveform[1][1] > 0.
        pushfirst!(waveform, (0., zero_grad))
        events = [(i+1, e) for (i, e) in events]
    end

    components = BaseComponent[]
    for (index_grad, ((prev_time, prev_grad), (time, grad))) in enumerate(zip(waveform[1:end-1], waveform[2:end]))
        duration = time - prev_time
        if norm(prev_grad) <= 1e-12 && norm(grad) <= 1e-12
            push!(components, NoGradient{ndim}(duration))
        elseif norm(prev_grad) ≈ norm(grad)
            push!(components, ConstantGradient(prev_grad, orientation, duration, group))
        else
            push!(components, ChangingGradient(prev_grad, (grad .- prev_grad) ./ duration, orientation, duration, group))
        end
        while length(events) > 0 && index_grad == events[1][1]
            (_, event) = popfirst!(events)
            push!(components, event)
        end
    end
    #for comp in components
    #    scanner_constraints!(comp)
    #end
    return BuildingBlock(components)
end

make_generic(other_block::BaseBuildingBlock) = BuildingBlock([other_block...])
Base.keys(bb::BuildingBlock) = 1:length(bb.parts)
Base.getindex(bb::BuildingBlock, i::Integer) = bb.parts[i]

function get_pulse(bb::BuildingBlock)
    pulses = [p for (_, p) in events(bb) if p isa RFPulseComponent]
    if length(pulses) == 0
        error("BuildingBlock does not contain any pulses.")
    end
    if length(pulses) == 1
        return pulses[1]
    end
    error("BuildingBlock contains more than one pulse. Not sure which one to return.")
end

function get_readout(bb::BuildingBlock)
    readouts = [r for (_, r) in events(bb) if r isa ReadoutComponent]
    if length(readouts) == 0
        error("BuildingBlock does not contain any readouts.")
    end
    if length(readouts) == 1
        return readouts[1]
    end
    error("BuildingBlock contains more than one readout. Not sure which one to return.")
end

function effective_time(bb::BuildingBlock)
    index = [i for (i, r) in events(bb) if r isa Union{RFPulseComponent, ReadoutComponent}]
    if length(index) == 0 
        error("BuildingBlock does not contain any RF pulse or readout events, so `effective_time` is not defined.")
    elseif length(index) > 1
        error("BuildingBlock contains multiple RF pulse or readout events, so `effective_time` is not defined.")
    end
    index = index[1]
    return effective_time(bb, index)
end

struct Wait <: BaseBuildingBlock
    duration :: VariableType
    function Wait(var)
        res = new(get_free_variable(var))
        if !(res.duration isa Number)
            @constraint global_model() res.duration >= 0
        end
        return res
    end
end

duration(wb::Wait) = wb.duration
Base.keys(::Wait) = (Val(:empty),)
Base.getindex(wb::Wait, ::Val{:empty}) = NoGradient{1}(wb.duration)

end