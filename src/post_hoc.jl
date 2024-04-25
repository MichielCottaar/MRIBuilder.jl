"""
Define post-fitting adjustments of the sequences
"""
module PostHoc

import ..Variables: AbstractBlock, adjust_internal, adjustable
import ..Components: GradientWaveform, RFPulseComponent, BaseComponent, NoGradient
import ..Containers: ContainerBlock, Sequence, Wait

"""
    adjust(block; kwargs...)

Generate a new sequence/building_block/component with some post-fitting adjustments.

The following adjustments are allowed:
- for MR gradients
    - `orientation`: set the orientation to a given vector.
    - `rotation`: rotate the gradient orientations using a rotations from [`Rotations.jl`](https://juliageometry.github.io/Rotations.jl/stable/).
    - `scale`: multiply the gradient strength by the given value. Note that if you use a value not between -1 and 1 you might break the scanner's maximum gradient or slew rate.
- for RF pulses:
    - `frequency`: shift the off-resonance frequency by the given value (in kHz).
    - `scale`: multiply the RF pulse amplitude by the given value (used to model the B1 transmit field).

Specific sequence components that can be adjusted are identified by their `group` name.
For example, `adjust(sequence, dwi=(orientation=[0, 1, 0], ))` will set any gradient in the group `dwi` to point in the y-direction.

To affect all gradients or pulses, use `gradient=` or `pulse`, e.g.
`adjust(sequence, pulse=(scale=0.5, ))`
will divide the amplitude of all RV pulses by two.
"""
function adjust(block::AbstractBlock; kwargs...) 
    used_names = Set{Symbol}()
    res = adjust_helper(block, used_names; kwargs...)
    unused_names = filter(keys(kwargs)) do key
        !(key in used_names)
    end
    if length(unused_names) > 0
        @warn "Some group/type names were not used in call to `MRIBuilder.adjust`, namely: $(unused_names)."
    end
    res
end

function adjust_helper(block::AbstractBlock, used_names::Set{Symbol}; gradient=(), pulse=(), kwargs...)
    params = []
    adjust_type = adjustable(block)
    if adjust_type == :false
        for prop_name in propertynames(block)
            push!(params, adjust_helper(getproperty(block, prop_name), used_names; gradient=gradient, pulse=pulse, kwargs...))
        end
        return typeof(block)(params...)
    else
        if !isnothing(block.group) && (block.group in keys(kwargs))
            push!(used_names, block.group)
            return adjust_internal(block; kwargs[block.group]...)
        elseif adjust_type == :gradient
            push!(used_names, :gradient)
            return adjust_internal(block; gradient...)
        elseif adjust_type == :pulse
            push!(used_names, :pulse)
            return adjust_internal(block; pulse...)
        end
    end
end

adjust_helper(some_value, used_names::Set{Symbol}; kwargs...) = some_value
adjust_helper(array_variable::AbstractArray, used_names::Set{Symbol}; kwargs...) = map(array_variable) do v adjust_helper(v, used_names; kwargs...) end
adjust_helper(dict_variable::AbstractDict, used_names::Set{Symbol}; kwargs...) = typeof(dict_variable)(k => adjust_helper(v, used_names; kwargs...) for (k, v) in pairs(dict_variable))
adjust_helper(tuple_variable::Tuple, used_names::Set{Symbol}; kwargs...) = map(tuple_variable) do v adjust_helper(v, used_names; kwargs...) end
adjust_helper(pair:: Pair, used_names::Set{Symbol}; kwargs...) = adjust_helper(pair[1], used_names; kwargs...) => adjust_helper(pair[2], used_names; kwargs...)


"""
    merge(sequences...; wait_time=0.)

Merge multiple sequences together.

Sequences will be run one after each other with `wait_time` in between.
"""
merge(sequences::Sequence{S}...; kwargs...) where {S} = merge_internal(sequences...; name=S, kwargs...)

merge(sequences::Sequence...; kwargs...) = merge_internal(sequences...; kwargs...)

function merge_internal(sequences...; name=:Sequence, wait_time=0.)
    wb = Wait(wait_time)
    new_blocks = ContainerBlock[]
    for seq in sequences
        append!(new_blocks, map(b -> b[2], seq.blocks))
        push!(new_blocks, wb)
    end
    new_blocks = new_blocks[1:end-1]

    return Sequence(new_blocks; scanner=sequences[1].scanner, name=name)
end

end