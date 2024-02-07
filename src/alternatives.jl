module Alternatives
import JuMP: @constraint
import ..BuildingBlocks: BuildingBlock, match_blocks!
import ..BuildSequences: global_model
import ..Variables: duration

"""
    AlternativeBlocks(name, blocks)

Represents a part of the sequence where there are multiple possible alternatives.

Variables can be matched across these alternatives using [`match_blocks!`](@ref).

The `name` is a symbol that is used to identify this `AlternativeBlocks` in the broader sequence.
"""
struct AlternativeBlocks <: BuildingBlock
    name :: Symbol
    options :: Vector{<:BuildingBlock}
end

Base.getindex(alt::AlternativeBlocks, index::Int) = alt.options[index]

duration(alt::AlternativeBlocks) = maximum(duration.(alt.options))

"""
    match_blocks!(alternatives, function)

Matches the outcome of given `function` on each of the building blocks in [`AlternativeBlocks`](@ref).

For example, `match_blocks!(alternatives, duration)` will ensure that all the alternative building blocks have the same duration.
"""
function match_blocks!(alternatives::AlternativeBlocks, func)
    baseline = func(alternatives[1])
    for other_block in alternatives.options[2:end]
        if baseline isa AbstractVector
            @constraint global_model() baseline == func(other_block)
        else
            @constraint global_model() baseline .== func(other_block)
        end
    end
end



end