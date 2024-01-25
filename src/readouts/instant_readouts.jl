module InstantReadouts
import ...BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, duration, properties, to_mcmr_components
import ...SequenceBuilders: SequenceBuilder, start_time, to_block

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
end