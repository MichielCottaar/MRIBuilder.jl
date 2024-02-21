"""
Defines [`BaseSequence`](@ref) and [`Sequence`](@ref)
"""
module Sequences
import StaticArrays: SVector
import JuMP: @constraint
import ...Variables: get_free_variable, TR, VariableType, duration, variables, VariableNotAvailable, Variables, set_simple_constraints!
import ...BuildSequences: global_model
import ...Components: EventComponent
import ..Abstract: ContainerBlock, start_time
import ..BuildingBlocks: Wait, BuildingBlock, BaseBuildingBlock

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

function start_time(bs::BaseSequence{N}, s::Symbol) where {N}
    index = findfirst(p -> p[1] == s, bs.blocks)
    if isnothing(index)
        error("Sequence does not contain a block with name $s")
    end
    return start_time(bs, index)
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

"""
    Sequence(blocks; TR=:min, nrepeat=0)
    Sequence(blocks...; TR=:min, nrepeat=0)

Defines an MRI sequence from a vector of building blocks.

## Arguments
- [`nrepeat`](@ref): how often the sequence will repeat itself (keep at default of 0 to repeat indefinetely).
- [`TR`](@ref): how long between repeats in ms (defaults to the duration of the sequence). Can be set to `nothing` to be a free variable.
"""
struct Sequence{S, N} <: BaseSequence{N}
    blocks :: SVector{N, Pair{<:Union{Symbol, Nothing}, <:ContainerBlock}}
end

function Sequence(blocks::AbstractVector; name=:Sequence, variables...)
    blocks = to_block_pair.(blocks)
    res = Sequence{name, length(blocks)}(SVector{length(blocks)}(blocks))
    set_simple_constraints!(res, variables)
    return res
end

Base.show(io::IO, ::Type{<:Sequence{S, N}}) where {S, N} = print(io, S, "{$N}")

Sequence(blocks...; kwargs...) = Sequence([blocks...]; kwargs...)

get_index_single_TR(s::Sequence, i::Integer) = s.blocks[i][2]
Base.getindex(seq::Sequence, sym::Symbol) = seq[findfirst(p -> p[1] == sym, seq.blocks)]
nrepeat(::Sequence) = 0

to_block_pair(pair::Pair) = pair[1] => to_block(pair[2])
to_block_pair(other) = nothing => to_block(other)

"""
    to_block(block_like)

Converst object into something that can be included in the sequence:
- :min/:max/number/variable/nothing => [`Wait`](@ref).
- `building_block` or `sequence` => no change.
- RF pulse/readout => will be embedded within a [`BuildingBlock`](@ref).
"""
to_block(cb::ContainerBlock) = cb
to_block(s::Symbol) = to_block(Val(s))
to_block(s::Union{VariableType, Nothing, Val{:min}, Val{:max}}) = Wait(s)
to_block(ec::EventComponent) = BuildingBlock([(0., 0.), (duration(ec), 0.)], [(1, ec)])

end
