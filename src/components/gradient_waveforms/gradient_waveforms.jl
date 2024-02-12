"""
Module defining sub-types of the [`GradientWaveform`](@ref).

There are only three types of [`GradientBlock`] objects:
- [`ChangingGradientBlock`](@ref): any gradient changing linearly in strength.
- [`ConstantGradientBlock`](@ref): any gradient staying constant in strength. These can overlap with a pulse (`SliceSelectPulse`).
- [`NoGradientBlock`](@ref): any part of the gradient waveform when no gradient is active.

These parts are combined into a full gradient waveform in a `BuildingBlock`.

Each part of this gradient waveform can compute:
- [`gradient_strength`]: maximum gradient strength in each dimension.
- [`slew_rate`]: maximum slew rate in each dimension.
- [`qvec`]: area under curve in each dimension
- [`bmat_gradient`]: diffusion weighting (scalar in 1D or matrix in 3D).
"""
module GradientWaveforms

include("changing_gradient_blocks.jl")
include("constant_gradient_blocks.jl")
include("no_gradient_blocks.jl")


import ..AbstractTypes: GradientWaveform
import .NoGradientBlocks: NoGradient
import .ChangingGradientBlocks: ChangingGradient, split_gradient
import .ConstantGradientBlocks: ConstantGradient

end