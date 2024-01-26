module InstantReadouts
import ...BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, duration, properties
import ...SequenceBuilders: SequenceBuilder, start_time, to_block
import ...ConcreteBlocks: AbstractConcreteBlock

"""
    InstantReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct InstantReadout <: AbstractConcreteBlock
    builder::SequenceBuilder
end

InstantReadout() = BuildingBlockPlaceholder{InstantReadout}()

properties(::Type{<:InstantReadout}) = []
to_block(builder::SequenceBuilder, cls::Type{<:InstantReadout}) = cls(builder)

end