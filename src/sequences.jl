"""
Define the [`Sequence`](@ref) building block.
"""
module Sequences
import Printf: @sprintf
import JuMP: @constraint
import ...BuildSequences: global_model
import ...Variables: variables, start_time, duration, VariableType, get_free_variable, TR, end_time
import ...BuildingBlocks: BuildingBlock, ContainerBlock, to_block, get_children_indices, fixed, BuildingBlockPrinter, _robust_value, fixed, make_generic

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
        if TR isa Number
            if !isinf(TR)
                @assert seq.TR >= duration(seq)
            end
        else
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

# print timings when printing sequences
Base.show(io::IO, seq::Sequence) = print(io, BuildingBlockPrinter(seq, 0., 0))

function Base.show(io::IO, printer::BuildingBlockPrinter{<:Sequence})
    seq = printer.bb
    print(io, "Sequence(")
    d = _robust_value(duration(seq))
    if !isnothing(d)
        if !isnothing(printer.start_time)
            print(io, "t=", @sprintf("%.3g", printer.start_time), "-", @sprintf("%.3g", printer.start_time + d), ", ")
        else
            print(io, "duration=", @sprintf("%.3g", d), ",")
        end
    end
    TR = _robust_value(seq.TR)
    if !isnothing(TR) && isfinite(TR)
        print(io, "TR=", Int(round(TR)))
    end
    print(io, "):")

    for child_index in get_children_indices(seq)
        child_block = seq[child_index]
        child_printer = BuildingBlockPrinter(
            child_block,
            isnothing(printer.start_time) ? nothing : _robust_value(start_time(seq, child_index) + printer.start_time),
            printer.spaces + 2
        )
        print(io, "\n", repeat(' ', printer.spaces + 2), "- ", child_index, ": ", child_printer)
    end

end

end


