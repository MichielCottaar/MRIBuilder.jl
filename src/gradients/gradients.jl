"""
Module defining sub-types of the [`GradientBlock`](@ref), i.e., any [`BuildingBlock`](@ref) that only defines a gradient profile.

There are only three types of [`GradientBlock`] objects:
- [`ChangingGradientBlock`](@ref): any gradient changing linearly in strength.
- [`ConstantGradientBlock`](@ref): any gradient staying constant in strength. These can overlap with a pulse (`SliceSelectPulse`).
- [`InstantGradientBlock`](@ref): an infinitely short gradient pulse.

Combinations of these should sub-type from [`SpecificWaveform`](@ref).
"""
module Gradients

include("changing_gradient_blocks.jl")
include("constant_gradient_blocks.jl")
include("instant_gradients.jl")


import ..BuildingBlocks: GradientBlock
import .ChangingGradientBlocks: ChangingGradientBlock, split_gradient
import .ConstantGradientBlocks: ConstantGradientBlock
import .InstantGradients: InstantGradientBlock
end