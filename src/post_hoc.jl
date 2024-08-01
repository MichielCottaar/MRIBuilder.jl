"""
Define post-fitting adjustments of the sequences
"""
module PostHoc

import ..Variables: AbstractBlock, adjust_internal, adjustable, adjust
import ..Components: GradientWaveform, RFPulseComponent, BaseComponent, NoGradient
import ..Containers: ContainerBlock, Sequence, Wait

"""
    adjust(block; kwargs...)

Generate one or more new sequences/building_blocks/components with some post-fitting adjustments.

The following adjustments are allowed:
- for MR gradients
    - `orientation`: set the orientation to a given vector.
    - `rotation`: rotate the gradient orientations using a rotations from [`Rotations.jl`](https://juliageometry.github.io/Rotations.jl/stable/).
    - `scale`: multiply the gradient strength by the given value. Note that if you use a value not between -1 and 1 you might break the scanner's maximum gradient or slew rate.
- for RF pulses:
    - `frequency`: shift the off-resonance frequency by the given value (in kHz).
    - `scale`: multiply the RF pulse amplitude by the given value (used to model the B1 transmit field).

A vector of multiple values can be passed on to any of these in order to create multiple sequences with different adjustments.
The will usually be merged together. You can get the individual sequences by passing on `merge=false`.
The time between these repeated sequences can be adjusted using the keywords described in [`merge_sequences`](@ref) passed on to the merge keyword:
e.g., `merge=(wait_time=10, )` adds a wait time of 10 ms between each repeated sequence.

Specific sequence components that can be adjusted are identified by their `group` name.
For example, `adjust(sequence, diffusion=(orientation=[0, 1, 0], ))` will set any gradient in the group `:diffusion` to point in the y-direction.

To affect all gradients or pulses, use `gradient=` or `pulse`, e.g.
`adjust(sequence, pulse=(scale=0.5, ))`
will divide the amplitude of all RV pulses by two.
"""
function adjust(block::AbstractBlock; merge=true, kwargs...) 
    used_names = Set{Symbol}()
    n_adjust, kwargs_list = adjust_kwargs_list(; kwargs...)
    if isnothing(n_adjust)
        res = adjust_helper(block, used_names; kwargs_list[1]...)
    else
        res = [adjust_helper(block, used_names; kw...) for kw in kwargs_list]
        if merge !== false
            if merge === true
                merge = NamedTuple()
            end
            res = merge_sequences(res...; merge...)
        end
    end
    unused_names = filter(keys(kwargs)) do key
        !(key in used_names)
    end
    if length(unused_names) > 0
        @warn "Some group/type names were not used in call to `MRIBuilder.adjust`, namely: $(unused_names)."
    end
    res
end

function adjust_kwargs_list(; kwargs...)
    n_adjust = nothing
    for (_, named_tuple) in kwargs
        for key in keys(named_tuple)
            value = named_tuple[key]
            if key == :orientation && value isa AbstractVector{<:Number}
                continue
            end
            if value isa AbstractVector
                if isnothing(n_adjust)
                    n_adjust = length(value)
                else
                    @assert length(value) == n_adjust
                end
            end
        end
    end
    use_n_adjust = isnothing(n_adjust) ? 1 : n_adjust
    kwargs_list = [Dict{Symbol, Any}([field=>Dict{Symbol, Any}() for field in keys(kwargs)]...) for _ in 1:use_n_adjust]
    for (field, named_tuple) in kwargs
        for key in keys(named_tuple)
            value = named_tuple[key]
            for index in 1:use_n_adjust
                if (key == :orientation && value isa AbstractVector{<:Number}) || !(value isa AbstractVector)
                    kwargs_list[index][field][key] = value
                else
                    kwargs_list[index][field][key] = value[index]
                end
            end
        end
    end
    return (n_adjust, kwargs_list)
end


function adjust_helper(block::AbstractBlock, used_names::Set{Symbol}; gradient=(), pulse=(), kwargs...)
    params = []
    adjust_type = adjustable(block)
    if adjust_type != :false && (!isnothing(block.group) && (block.group in keys(kwargs)))
        push!(used_names, block.group)
        new_block = adjust_internal(block; kwargs[block.group]...)
    elseif adjust_type == :gradient
        push!(used_names, :gradient)
        new_block = adjust_internal(block; gradient...)
    elseif adjust_type == :pulse
        push!(used_names, :pulse)
        new_block = adjust_internal(block; pulse...)
    else
        new_block = block
    end
    for prop_name in propertynames(new_block)
        push!(params, adjust_helper(getproperty(new_block, prop_name), used_names; gradient=gradient, pulse=pulse, kwargs...))
    end
    return typeof(block)(params...)
end

adjust_helper(some_value, used_names::Set{Symbol}; kwargs...) = some_value
adjust_helper(array_variable::AbstractArray, used_names::Set{Symbol}; kwargs...) = map(array_variable) do v adjust_helper(v, used_names; kwargs...) end
adjust_helper(dict_variable::AbstractDict, used_names::Set{Symbol}; kwargs...) = typeof(dict_variable)(k => adjust_helper(v, used_names; kwargs...) for (k, v) in pairs(dict_variable))
adjust_helper(tuple_variable::Tuple, used_names::Set{Symbol}; kwargs...) = map(tuple_variable) do v adjust_helper(v, used_names; kwargs...) end
adjust_helper(pair:: Pair, used_names::Set{Symbol}; kwargs...) = adjust_helper(pair[1], used_names; kwargs...) => adjust_helper(pair[2], used_names; kwargs...)


"""
    merge_sequences(sequences...; wait_time=0.)

Merge multiple sequences together.

Sequences will be run one after each other with `wait_time` in between.
"""
merge_sequences(sequences::Sequence{S}...; kwargs...) where {S} = merge_internal(sequences...; name=S, kwargs...)

merge_sequences(sequences::Sequence...; kwargs...) = merge_internal(sequences...; kwargs...)

function merge_internal(sequences...; name=:Sequence, wait_time=0.)
    wb = Wait(wait_time)
    new_blocks = ContainerBlock[sequences[1]]
    for seq in sequences[2:end]
        if !iszero(wait_time)
            push!(new_blocks, wb)
        end
        push!(new_blocks, seq)
    end

    return Sequence(new_blocks; scanner=sequences[1].scanner, name=name)
end

end