"""
Define the [`Sequence`](@ref) building block.
"""
module Sequences
import JuMP: @constraint
import ...BuildSequences: global_model
import ...Variables: variables, start_time, duration, VariableType, get_free_variable, TR, end_time
import ...BuildingBlocks: BuildingBlock, ContainerBlock, to_block, get_children_indices, fixed, fixed, make_generic

"""
    Sequence(building_blocks...; TR=nothing)
    Sequence([building_blocks]; TR=nothing)

Represents a series of [`BuildingBlock`](@ref) objects run in turn.

This can be used as a top-level NMR/MRI sequence (in which case the [`TR`](@ref) variable is relevant)
or be embedded as a [`BuildingBlock`](@ref) into higher-order `Sequence` or other [`ContainerBlock`](@ref) objects.

## Variables
- [`TR`](@ref): repetition time of sequence in ms.
"""
struct Sequence <: ContainerBlock
    _blocks :: Vector{<:BuildingBlock}
    TR :: VariableType
    function Sequence(blocks::AbstractVector; TR=nothing)
        seq = new(
            to_block.(blocks),
            get_free_variable(TR),
        )
        if !(TR isa Number && (isinf(TR) || duration(seq) isa Number))
            @constraint global_model() seq.TR >= duration(seq)
        end
        return seq
    end
end

Sequence(blocks...; TR=nothing) = Sequence([blocks...]; TR=TR)
fixed(seq::Sequence) = Sequence(fixed.(seq._blocks), TR=fixed(seq.TR))
make_generic(seq::Sequence) = Sequence(make_generic.(seq._blocks), TR=seq.TR)

Base.length(seq::Sequence) = length(seq._blocks)
Base.getindex(seq::Sequence, index) = seq._blocks[index]
get_children_indices(seq::Sequence) = eachindex(seq._blocks)

"""
    start_time(sequence::Sequence, index::Integer, args...)

Returns the starting time of the [`BuildingBlock`](@ref) with index `index`.
Additional `args` can be used to select a sub-block of that [`BuildingBlock`](@ref).
The starting time is returned with respect to the start of this sequence.
"""
start_time(seq::Sequence, index::Integer) = isone(index) ? start_time(seq) : (start_time(seq, index-1) + duration(seq[index-1]))

duration(seq::Sequence) = end_time(seq, length(seq))

TR(seq::Sequence) = seq.TR


end


