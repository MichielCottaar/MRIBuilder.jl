"""
Builds and optimises NMR/MRI sequences.
"""
module MRIBuilder

include("build_sequences.jl")
include("scanners.jl")
include("variables.jl")
include("building_blocks.jl")
include("concrete_blocks.jl")
include("wait.jl")
include("containers/containers.jl")
include("gradients/gradients.jl")
include("pulses/pulses.jl")
include("readouts/readouts.jl")

import .BuildSequences: build_sequence
export build_sequence

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

import .Variables: variables, duration, start_time, end_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, area_under_curve, δ, rise_time, flat_time, slew_rate, gradient_strength
export variables, duration, start_time, end_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, area_under_curve, δ, rise_time, flat_time, slew_rate, gradient_strength

import .BuildingBlocks: BuildingBlock
export BuildingBlocks

import .ConcreteBlocks: ConcreteBlock, AbstractConcreteBlock
export ConcreteBlocks, AbstractConcreteBlock

import .Wait: WaitBlock
export WaitBlock

import .Containers: Sequence, FixedBlock
export Sequence, FixedBlock

import .Gradients: PulsedGradient, InstantGradientBlock, FixedGradient
export PulsedGradient, InstantGradientBlock, FixedGradient

import .Pulses: InstantRFPulseBlock, ConstantPulse, SincPulse, FixedPulse
export InstantRFPulseBlock, ConstantPulse, SincPulse, FixedPulse

import .Readouts: InstantReadout
export InstantReadout

import JuMP: @constraint, @objective, objective_function, optimize!, has_values, value
export @constraint, @objective, objective_function, optimize!, has_values, value

end
