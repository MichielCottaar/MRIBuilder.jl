"""
Module defining sub-types of the [`GradientBlock`](@ref), i.e., any [`BuildingBlock`](@ref) that only defines a gradient profile.

- [`PulsedGradient`](@ref)
- [`InstantGradientBlock`](@ref)
- [`FixedGradient`](@ref)
"""
module Gradients

#include("integrate_gradients.jl")
include("fixed_gradients.jl")
include("pulsed_gradients.jl")
include("instant_gradients.jl")

import ..BuildingBlocks: GradientBlock

import .PulsedGradients: PulsedGradient
import .InstantGradients: InstantGradientBlock
import .FixedGradients: FixedGradient
end