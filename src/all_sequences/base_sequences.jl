module BaseSequences
import ...ContainerBlocks: ContainerBlock, start_time
import ...AllBuildingBlocks: BaseBuildingBlock
import ...Variables: TR, duration

"""
Super-type of any sequence of non-overlapping building blocks that should be played after each other.

It contains `N` [`ContainerBlock`](@ref) objects (e.g., building blocks or other sequences).

Main interface:
- Acts as an iterable containing the blocks and sequences.
    - Indiviual blocks/sequences can be obtained using indexing.
    - If there is a finite number of repeats, the iteration will continue over all repeats.

Sub-types need to implement:
- `get_index_single_TR`: return the index assuming it is between 1 and N
- [`nrepeat`](@ref): how often the sequence should repeat (if not implemented, this will be 1).
- [`TR`](@ref): time scale on which to repeat (if not implemented, this will be the sum of the individual block durations).
"""
abstract type BaseSequence{N} <: ContainerBlock end

function Base.getindex(bs::BaseSequence{N}, index::Integer) where {N}
    if nrepeat(bs) > 0 && (index < 1 || index > length(bs))
        throw(BoundsError(bs, index))
    end
    base_index = ((index - 1) % N) + 1
    return get_index_single_TR(bs, base_index)
end
Base.iterate(bs::BaseSequence) = Base.iterate(bs, 1)
Base.iterate(bs::BaseSequence{N}, index::Integer) where {N} = index > length(bs) ? nothing : (bs[index], index + 1)
Base.length(bs::BaseSequence{N}) where {N} = iszero(nrepeat(bs)) ? N : (nrepeat(bs) * N)
Base.eltype(::Type{<:BaseSequence}) = ContainerBlock

function start_time(bs::BaseSequence{N}, index::Integer) where {N}
    nTR = div(index-1, N, RoundDown)
    if 0 < nrepeat(bs) <= nTR
        throw(BoundsError(bs, index))
    end
    base_index = ((index - 1) % N) + 1
    base_time = sum(i -> duration(bs[i]), 1:base_index-1; init=0.)
    if iszero(nTR)
        return base_time
    else
        return nTR * TR(bs) + base_time
    end
end


"""
    get_index_single_TR(sequence, index)

Used internally by any [`BaseSequence`](@ref) to get a specific block.
The `index` should be between 1 and N.
It should be implemented for any sub-classes of [`BaseSequence`](@ref).
"""
get_index_single_TR(seq::BaseSequence, index::Integer) = get_index_single_TR(seq, Val(index))


"""
    nrepeat(sequence)

How often sequence should be repeated. 
"""
nrepeat(bs::BaseSequence) = 1

"""
    TR(sequence)

On what timescale the sequence should repeat itself in ms (will have no effect if [`nrepeat`](@ref) is one).

By default this is set to the total duration of the sequence.
"""
TR(bs::BaseSequence) = duration(bs)


duration(bs::BaseSequence{0}) = 0.
duration(bs::BaseSequence) = sum(duration.(bs); init=0.)


end
