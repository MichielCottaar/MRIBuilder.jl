"""
Define the [`Sequence`](@ref) building block.
"""
module Sequences
import JuMP: Model, @constraint
import ...BuildSequences: @global_model_constructor
import ...Variables: variables, start_time, duration, VariableType, get_free_variable, TR, end_time
import ...BuildingBlocks: BuildingBlock, ContainerBlock, to_block, get_children_indices, scanner_constraints!, fixed, BuildingBlockPrinter
import ..FixedBlocks: FixedBlock

abstract type AbstractSequence <: ContainerBlock end

"""
    Sequence(building_blocks...; TR=nothing, scanner=nothing)
    Sequence([building_blocks]; TR=nothing, scanner=nothing)

Represents a series of [`BuildingBlock`](@ref) objects run in turn.

Providing a [`Scanner`](@ref) will lead to [`scanner_constraints!`](@ref) to be called on all building blocks.

This can be used as a top-level NMR/MRI sequence (in which case the [`TR`](@ref) variable is relevant)
or be embedded as a [`BuildingBlock`](@ref) into higher-order `Sequence` or other [`ContainerBlock`](@ref) objects.

## Variables
- [`TR`](@ref): repetition time of sequence in ms.
"""
struct Sequence <: AbstractSequence
    model :: Model
    _blocks :: Vector{<:BuildingBlock}
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

Base.length(seq::AbstractSequence) = length(seq._blocks)
Base.getindex(seq::AbstractSequence, index) = seq._blocks[index]
get_children_indices(seq::AbstractSequence) = eachindex(seq._blocks)

"""
    start_time(sequence::Sequence, index::Integer, args...)

Returns the starting time of the [`BuildingBlock`](@ref) with index `index`.
Additional `args` can be used to select a sub-block of that [`BuildingBlock`](@ref).
The starting time is returned with respect to the start of this sequence.
"""
start_time(seq::AbstractSequence, index::Integer) = isone(index) ? start_time(seq) : (start_time(seq, index-1) + duration(seq[index-1]))

duration(seq::AbstractSequence) = end_time(seq, length(seq))

TR(seq::AbstractSequence) = seq.TR
variables(::Type{<:Sequence}) = [TR]

# print timings when printing sequences
Base.show(io::IO, seq::AbstractSequence) = print(io, BuildingBlockPrinter(seq, nothing, 0., 0))

"""
Represents a series of [`BuildingBlock`](@ref) run in sequence.

Each [`BuildingBlock`](@ref) is fixed, i.e, it does not contain any free parameters.

A fixed version of any [`Sequence`](@ref) can be created by calling [`fixed`](@ref) on that sequence.
"""
struct FixedSequence <: AbstractSequence
    _blocks :: Vector{<:BuildingBlock}
    TR :: Number
end

function fixed(seq::Sequence)
    FixedSequence(
        fixed.(seq._blocks),
        value(TR(seq))
    )
end

variables(::Type{<:FixedSequence}) = []

end




