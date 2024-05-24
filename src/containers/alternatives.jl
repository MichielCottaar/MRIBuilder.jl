module Alternatives
import JuMP: @constraint
import ..Abstract: ContainerBlock
import ...BuildSequences: global_model, fixed
import ...Variables: @defvar, make_generic

"""
    AlternativeBlocks(name, blocks)

Represents a part of the sequence where there are multiple possible alternatives.

Variables can be matched across these alternatives using [`match_blocks!`](@ref).

The `name` is a symbol that is used to identify this `AlternativeBlocks` in the broader sequence.
"""
struct AlternativeBlocks <: ContainerBlock
    name :: Symbol
    options :: Dict{Any, <:ContainerBlock}
end

AlternativeBlocks(name::Symbol, options_vector::AbstractVector) = AlternativeBlocks(name, Dict(index => value for (index, value) in enumerate(options_vector)))

Base.getindex(alt::AlternativeBlocks, index) = alt.options[index]
Base.length(alt::AlternativeBlocks) = length(alt.options)

@defvar duration(alt::AlternativeBlocks) = maximum(duration.(values(alt.options)))

"""
    match_blocks!(alternatives, function)

Matches the outcome of given `function` on each of the building blocks in [`AlternativeBlocks`](@ref).

For example, `match_blocks!(alternatives, duration)` will ensure that all the alternative building blocks have the same duration.
"""
function match_blocks!(alternatives::AlternativeBlocks, func)
    options = [values(alternatives.options)...]
    baseline = func(options[1])
    for other_block in options[2:end]
        if baseline isa AbstractVector
            @constraint global_model() baseline == func(other_block)
        else
            @constraint global_model() baseline .== func(other_block)
        end
    end
end

fixed(alt::AlternativeBlocks) = AlternativeBlocks(alt.name, Dict(key=>fixed(value) for (key, value) in alt.options))
make_generic(alt::AlternativeBlocks) = AlternativeBlocks(alt.name, Dict(key=>make_generic(value) for (key, value) in alt.options))


end