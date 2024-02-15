module Sequences
import StaticArrays: SVector
import JuMP: @constraint
import ...Variables: get_free_variable, TR, VariableType, duration
import ...AllBuildingBlocks: Wait, BuildingBlock
import ...BuildSequences: global_model
import ...Components: EventComponent
import ..BaseSequences: ContainerBlock, BaseSequence, nrepeat, get_index_single_TR

"""
    Sequence(blocks; TR=:min, nrepeat=0)
    Sequence(blocks...; TR=:min, nrepeat=0)

Defines an MRI sequence from a vector of building blocks.

## Arguments
- [`nrepeat`](@ref): how often the sequence will repeat itself (keep at default of 0 to repeat indefinetely).
- [`TR`](@ref): how long between repeats in ms (defaults to the duration of the sequence). Can be set to `nothing` to be a free variable.
"""
struct Sequence{N} <: BaseSequence{N}
    blocks :: SVector{N, <:ContainerBlock}
    TR :: VariableType
end

function Sequence(blocks::AbstractVector; TR=:min)
    blocks = to_block.(blocks)
    actual_duration = sum(duration, blocks; init=0.)
    if TR == :min
        TR = actual_duration
    end
    res = Sequence{length(blocks)}(SVector{length(blocks)}(blocks), get_free_variable(TR))
    if !(res.TR isa Number) || !(duration(res) isa Number)
        @constraint global_model() res.TR >= actual_duration
    end
end

Sequence(blocks...; kwargs...) = Sequence([blocks...]; kwargs...)

get_index_single_TR(s::Sequence, i::Integer) = s.blocks[i]
nrepeat(::Sequence) = 0

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
to_block(ec::EventComponent) = BuildingBlock([], [(0, ec)]; duration=duration(ec))


end