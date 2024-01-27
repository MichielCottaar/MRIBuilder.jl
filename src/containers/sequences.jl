module Sequences
import JuMP: Model, @constraint
import ...BuildSequences: @global_model_constructor
import ...Variables: variables, start_time, duration, VariableType, get_free_variable, TR
import ...BuildingBlocks: BuildingBlock, ContainerBlock, to_block, get_children_blocks, scanner_constraints!

"""
    Sequence(building_blocks...; TR=nothing, scanner=nothing)
    Sequence([building_blocks]; TR=nothing, scanner=nothing)

Represents a series of [`BuildingBlock`](@ref) objects run in turn.

Providing a [`Scanner`](@ref) will lead to [`scanner_constraints!`](@ref) to be called on all building blocks.

## Variables
- [`TR`](@ref): repetition time of sequence in ms.
"""
struct Sequence <: ContainerBlock
    model :: Model
    blocks :: Vector{<:BuildingBlock}
    TR :: VariableType
    function Sequence(model::Model, blocks::AbstractVector; TR=nothing, scanner=nothing)
        seq = new(
            model,
            to_block.(blocks),
            get_free_variable(model, TR),
        )
        @constraint model seq.TR >= duration(seq)
        if !isnothing(scanner)
            scanner_constraints!(model, seq, scanner)
        end
        return seq
    end
end

@global_model_constructor Sequence

Sequence(model::Model, blocks...; TR=nothing, scanner=nothing) = Sequence(model, [blocks...]; TR=TR, scanner=scanner)

Base.length(seq::Sequence) = length(seq.blocks)
Base.getindex(seq::Sequence, index) = seq.blocks[index]
get_children_blocks(seq::Sequence) = enumerate(seq.blocks)

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

TR(seq::Sequence) = seq.TR
variables(::Type{<:Sequence}) = [TR]


end




