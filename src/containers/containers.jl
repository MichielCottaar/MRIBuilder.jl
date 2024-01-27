"""
Defines [`ContainerBlock`](@ref) sub-types, i.e., [`BuildingBlocks`](@ref) objects that contain other BuildingBlocks.

The most important of these is [`Sequence`](@ref).

To get the children blocks use [`get_children_blocks`](@ref).
"""
module Containers
include("fixed_blocks.jl")
include("sequences.jl")
import ..BuildingBlocks: ContainerBlock, get_children_blocks

import .FixedBlocks: FixedBlock
import .Sequences: Sequence, FixedSequence
end