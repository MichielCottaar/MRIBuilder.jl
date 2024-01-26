module InstantReadouts
import ...BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, duration, properties
import ...SequenceBuilders: SequenceBuilder, start_time, to_block
import ...ConcreteBlocks: ConcreteBlock, to_concrete_block

"""
    InstantReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct InstantReadout <: BuildingBlock
    builder::SequenceBuilder
end

InstantReadout() = BuildingBlockPlaceholder{InstantReadout}()

duration(::InstantReadout) = 0.
properties(::Type{<:InstantReadout}) = []
to_block(builder::SequenceBuilder, cls::Type{<:InstantReadout}) = cls(builder)

to_mcmr_components(readout::InstantReadout) = Readout(value(start_time(readout)))
to_concrete_block(builder::SequenceBuilder, ::InstantReadout) = ConcreteBlock(builder, 0., readout_times=[0.])
end