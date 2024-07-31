module Alternatives
import JuMP: @constraint
import ..Abstract: ContainerBlock
import ...BuildSequences: fixed
import ...Variables: @defvar, make_generic, apply_simple_constraint!

"""
Parent type for all blocks that can take different MR sequence components between multiple repetitions of the sequence.

They can be extended into their individual components using `adjust(<name>=:all)`.

Each subtype of [`AbstractAlternativeBlock`](@ref) needs to implement two methods:
- [`get_alternatives_name`](@ref): returns the `name` used to identify this block in `adjust`
- [`get_alternatives_options`](@ref): returns a dictionary mapping the name of the different options to a [`ContainerBlock`](@ref) with the actual sequence building block.
"""
abstract type AbstractAlternativeBlocks <: ContainerBlock end

"""
    AlternativeBlocks(name, blocks)

Represents a part of the sequence where there are multiple possible alternatives.

Variables can be matched across these alternatives using [`match_blocks!`](@ref).

The `name` is a symbol that is used to identify this `AlternativeBlocks` in the broader sequence (as in `adjust`).
"""
struct AlternativeBlocks <: AbstractAlternativeBlocks
    name :: Symbol
    options :: Dict{Any, <:ContainerBlock}
end

AlternativeBlocks(name::Symbol, options_vector::AbstractVector) = AlternativeBlocks(name, Dict(index => value for (index, value) in enumerate(options_vector)))

"""
    get_alternatives_name(alternative_block)

Get the name with which any [`AbstractAlternativeBlocks`](@ref) will be identified in a call to `adjust`.
"""
get_alternatives_name(alt::AlternativeBlocks) = alt.name

"""
    get_alternatives_options(alternative_block)

Get the options available for a [`AbstractAlternativeBlocks`](@ref).
"""
get_alternatives_options(alt::AlternativeBlocks) = alt.options

Base.getindex(alt::AbstractAlternativeBlocks, index) = get_alternatives_options(alt)[index]
Base.length(alt::AbstractAlternativeBlocks) = length(get_alternatives_options(alt))
Base.keys(alt::AbstractAlternativeBlocks) = keys(get_alternatives_options(alt))
Base.values(alt::AbstractAlternativeBlocks) = values(get_alternatives_options(alt))

@defvar duration(alt::AbstractAlternativeBlocks) = maximum(variables.duration.(values(alt.options)))

"""
    match_blocks!(alternatives, function)

Matches the outcome of given `function` on each of the building blocks in [`AlternativeBlocks`](@ref).

For example, `match_blocks!(alternatives, duration)` will ensure that all the alternative building blocks have the same duration.
"""
function match_blocks!(alternatives::AbstractAlternativeBlocks, func)
    options = [values(get_alternatives_options(alternatives.options))...]
    if length(options) <= 1
        return
    end
    baseline = func(options[1])
    for other_block in options[2:end]
        apply_simple_constraint!(func(other_block), baseline)
    end
end

fixed(alt::AlternativeBlocks) = AlternativeBlocks(alt.name, Dict(key=>fixed(value) for (key, value) in alt.options))
make_generic(alt::AbstractAlternativeBlocks) = AlternativeBlocks(get_alternatives_name(alt), get_alternatives_options(alt))


end