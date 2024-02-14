module BaseSequences
import ...ContainerBlocks: ContainerBlock
import ...AllBuildingBlocks: BaseBuildingBlock
import ...Variables: TR, duration

"""
Super-type of any sequence of non-overlapping building blocks that should be played after each other.

It contains `N` [`ContainerBlock`](@ref) objects (e.g., building blocks or other sequences).

Main interface:
- Acts as an iterable containing the blocks and sequences.
    - Indiviual blocks/sequences can be obtained using indexing.

Sub-types need to implement:
- `Base.getindex`: returns the actual component for each integer index in 1:N.
"""
abstract type BaseSequence{N} <: ContainerBlock end


Base.iterate(bs::BaseSequence) = Base.iterate(bs, 1)
Base.iterate(bs::BaseSequence{N}, index::Integer) where {N} = index > N ? nothing : (bs[index], index + 1)
Base.length(::BaseSequence{N}) where {N} = N
Base.eltype(::Type{<:BaseSequence}) = ContainerBlock


"""
    nrepeat(sequence)

How often sequence should be repeated. 

Set to 0 to repeat indefinetely (default value if not overriden for a particular sequence type).
"""
nrepeat(bs::BaseSequence) = 0

"""
    TR(sequence)

On what timescale the sequence should repeat itself in ms (will have no effect if [`nrepeat`](@ref) is one).

By default this is set to the total duration of the sequence.
"""
TR(bs::BaseSequence) = duration(bs)


duration(bs::BaseSequence{0}) = 0.
duration(bs::BaseSequence) = sum(duration.(bs))


end
