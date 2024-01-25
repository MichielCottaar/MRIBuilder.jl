module IntegrateGradients
import ...BuildingBlocks: BuildingBlock
import ...SequenceBuilders: SequenceBuilder, TR, duration, builder, index


"""
    qval(blocks)
    qval(builder::SequenceBuilder, indices)

Computes the total q-value summed over multiple gradient [`BuildingBlock`](@ref) objects.

The `blocks` should contain the actual [`BuildingBlock`](@ref) objects, possibly interspersed with one of:
- `:TR` wait for one TR (has no effect, but included for consistency with [`bval`](@ref)).
- `:flip` flips the current `qval` value (e.g., if a refocus pulse has happened).

When a `builder` is provided, the `indices` are expected to be the integer indices of the individual blocks rather than the actual [`BuildingBlock`](@ref) objects.
`:TR` and `:flip` can also still be provided.
In addition, in this interface one can provide negative indices to indicate that a specific gradient should have a negative contribution to the q-value (this is an alternative for using :flip).

The integral can occur over multiple repetition times by including [`BuildingBlock`](@ref) objects out of order (or by using :TR).
"""
qval(builder::SequenceBuilder, indices::AbstractVector) = full_integral(builder, indices)[1]
qval(indices::AbstractVector) = full_integral(indices)[1]

"""
    bval(blocks)
    bval(builder::SequenceBuilder, indices)

Computes the total b-value combined over multiple gradient [`BuildingBlock`](@ref) objects.

The `blocks` should contain the actual [`BuildingBlock`](@ref) objects, possibly interspersed with one of:
- `:TR` wait for one TR (has no effect, but included for consistency with [`bval`](@ref)).
- `:flip` flips the current `qval` value (e.g., if a refocus pulse has happened).

When a `builder` is provided, the `indices` are expected to be the integer indices of the individual blocks rather than the actual [`BuildingBlock`](@ref) objects.
`:TR` and `:flip` can also still be provided.
In addition, in this interface one can provide negative indices to indicate that a specific gradient should have a negative contribution to the q-value (this is an alternative for using :flip).

The integral can occur over multiple repetition times by including [`BuildingBlock`](@ref) objects out of order (or by using :TR).
"""
bval(builder::SequenceBuilder, indices::AbstractVector) = full_integral(builder, indices)[2]
bval(indices::AbstractVector) = full_integral(indices)[2]

function full_integral(blocks::AbstractVector)
    actual_blocks = filter(b -> b isa BuildingBlock)
    if length(actual_blocks) == 0
        return (0., 0.)
    end
    sb = builder(blocks[1])

    return full_integral(sb, map(b -> b isa BuildingBlock ? index(sb, b) : b))
end

function full_integral(builder::SequenceBuilder, indices::AbstractVector)
    qval_current = 0.
    current_index = 0
    bval_current = 0.
    for index in indices
        if index == :flip
            qval_current = -qval_current
        elseif index == :TR
            bval_current = bval_current + qval_current^2 * TR(builder)
        elseif index isa Integer
            next_gradient = index < 0 ? -index : index
            wait_time = duration(builder, current_index + 1, next_gradient - 1)

            bval_current = (
                bval_current + 
                qval_current^2 * wait_time +
                bval(builder[next_gradient], index < 0 ? -qval_current : qval_current)
            )
            qval_current = (
                qval_current +
                sign(index) * qval(builder[next_gradient])
            )
        else
            error("qval/bval indices should refer to gradient blocks or :flip/:TR, not $index")
        end
    end
    return (qval_current, bval_current)
end

end