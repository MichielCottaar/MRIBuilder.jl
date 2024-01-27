module Sequences
import JuMP: Model
import ...BuildSequences: @global_model_constructor
import ...Variables: variables, start_time, duration, VariableType
import ...BuildingBlocks: BuildingBlock, to_block

"""
    Sequence(building_blocks...)
    Sequence([building_blocks])

Represents a series of [`BuildingBlock`](@ref) objects run in turn.
"""
struct Sequence <: BuildingBlock
    model :: Model
    blocks :: Vector{<:BuildingBlock}
end

@global_model_constructor Sequence

Sequence(model::Model, blocks...) = return Sequence(model, to_block.([blocks...]))

Base.length(seq::Sequence) = length(seq)
Base.getindex(seq::Sequence, index) = seq[index]

"""
    start_time(sequence::Sequence, index::Integer, args...)

Returns the starting time of the [`BuildingBlock`](@ref) with index `index`.
Additional `args` can be used to select a sub-block of that [`BuildingBlock`](@ref).
The starting time is returned with respect to the start of this sequence.
"""
start_time(seq::Sequence) = 0.
start_time(seq::Sequence, index::Integer) = iszero(index) ? start_time(seq) : (start_time(seq, index-1) + duration(seq[index]))
start_time(seq::Sequence, index::Integer, args...) = start_time(seq, index) + start_time(seq[index], args...)

"""
    end_time(sequence::Sequence, index::Integer, args...)

Returns the end time of the [`BuildingBlock`](@ref) with index `index`.
Additional `args` can be used to select a sub-block of that [`BuildingBlock`](@ref).
The end time is returned with respect to the start of this sequence.
"""
end_time(seq::Sequence, index::Integer) = start_time(seq, index) + duration(seq[index])
end_time(seq::Sequence, index::Integer, args...) = start_time(seq, index) + end_time(seq[index], args...)

duration(seq::Sequence) = end_time(seq, length(seq))



end




