module InstantReadouts
import ...BuildingBlocks: BuildingBlock, to_block, effective_time, make_generic
import ...Variables: variables, duration

"""
    InstantReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct InstantReadout <: BuildingBlock
end

to_block(::Type{<:InstantReadout}) = InstantReadout()
duration(::InstantReadout) = 0.
effective_time(::InstantReadout) = 0.

make_generic(ir::InstantReadout) = ir
end