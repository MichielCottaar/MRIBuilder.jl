"""
Builds and optimises NMR/MRI sequences.
"""
module MRIBuilder

include("scanners.jl")
include("building_blocks.jl")
include("sequence_builders.jl")
include("concrete_blocks.jl")
include("wait.jl")
include("gradients/gradients.jl")
include("pulses/pulses.jl")
include("readouts/readouts.jl")

import .BuildingBlocks: BuildingBlock, scanner_constraints!
export BuildingBlock, scanner_constraints!

import .SequenceBuilders: SequenceBuilder, start_time, end_time, duration, TR
export SequenceBuilder, start_time, end_time, duration, TR

import .ConcreteBlocks: ConcreteBlock, Sequence
export ConcreteBlock, Sequence

import .Wait: WaitBlock
export WaitBlock

import .Gradients: PulsedGradient, InstantGradientBlock, qval, rise_time, flat_time, slew_rate, gradient_strength, bval
export PulsedGradient, InstantGradientBlock, qval, rise_time, flat_time, slew_rate, gradient_strength, bval

import .Pulses: InstantRFPulseBlock, ConstantPulse, SincPulse, flip_angle, phase, frequency, bandwidth, N_left, N_right
export InstantRFPulseBlock, ConstantPulse, SincPulse, flip_angle, phase, frequency, bandwidth, N_left, N_right

import .Readouts: InstantReadout
export InstantReadout

import .Scanners: Scanner, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

using JuMP
export @constraint, @objective, objective_function, optimize!, has_values, value

end
