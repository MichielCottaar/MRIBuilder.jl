"""
Module defining sub-types of the [`GradientBlock`](@ref), i.e., any [`BuildingBlock`](@ref) that only defines a gradient profile.

In principle, there are only three types of [`GradientBlock`]:
- [`ChangingGradientBlock`](@ref): any gradient changing linearly in strength.
- [`ConstantGradientBlock`](@ref): any gradient staying constant in strength. These can overlap with a pulse (`SliceSelectPulse`).
- [`InstantGradientBlock`](@ref): an infinitely short gradient pulse.

All other gradient profiles are made up of [`ChangingGradientBlock`](@ref) and [`ConstantGradientBlock`](@ref) and hence are technically sequences.
They are still included here for clarity:
- [`PulsedGradient`](@ref)
- [`GradientWaveform`](@ref)
- [`FixedGradient`](@ref)
"""
module Gradients

include("changing_gradient_blocks.jl")
include("constant_gradient_blocks.jl")
include("slice_selects.jl")

#include("integrate_gradients.jl")
include("fixed_gradients.jl")
include("pulsed_gradients.jl")
include("instant_gradients.jl")

import ..BuildingBlocks: GradientBlock

import .PulsedGradients: PulsedGradient
import .InstantGradients: InstantGradientBlock
import .FixedGradients: FixedGradient
end