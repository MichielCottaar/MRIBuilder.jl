module InstantReadouts
import ...BuildingBlocks: BuildingBlock, to_block, fixed
import ...Variables: variables, duration

"""
    InstantReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct InstantReadout <: BuildingBlock
end

variables(::Type{<:InstantReadout}) = []
to_block(::Type{<:InstantReadout}) = InstantReadout()
duration(::InstantReadout) =0.
fixed(i::InstantReadout) = i
end