module InstantReadouts
import ...BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, duration, properties
import ...SequenceBuilders: SequenceBuilder, start_time, to_block, AbstractSequence
import ...ConcreteBlocks: AbstractConcreteBlock, to_concrete_block

"""
    InstantReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct InstantReadout <: AbstractConcreteBlock
    builder::AbstractSequence
end

InstantReadout() = BuildingBlockPlaceholder{InstantReadout}()

properties(::Type{<:InstantReadout}) = []
to_block(builder::SequenceBuilder, cls::Type{<:InstantReadout}) = cls(builder)
to_concrete_block(builder::AbstractSequence, ::InstantReadout) = InstantReadout(builder)

end