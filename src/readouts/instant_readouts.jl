module InstantReadouts
import ...BuildingBlocks: BuildingBlock
import ...ConcreteBlocks: AbstractConcreteBlock, to_concrete_block
import ...Variables: variables

"""
    InstantReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct InstantReadout <: AbstractConcreteBlock
end

variables(::Type{<:InstantReadout}) = []
to_concrete_block(::InstantReadout) = InstantReadout()
end