"""
Builds and optimises NMR/MRI sequences.
"""
module MRIBuilder

include("build_sequences.jl")
include("scanners.jl")
include("variables.jl")
include("building_blocks.jl")
include("wait.jl")
include("gradients/gradients.jl")
include("pulses/pulses.jl")
include("containers/containers.jl")
include("readouts/readouts.jl")
include("pathways.jl")

import .BuildSequences: build_sequence
export build_sequence

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

import .Variables: variables, duration, start_time, end_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength
export variables, duration, start_time, end_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength

import .BuildingBlocks: BuildingBlock, fixed
export BuildingBlocks, fixed

import .Wait: WaitBlock
export WaitBlock

import .Pulses: InstantRFPulseBlock, ConstantPulse, SincPulse, FixedPulse
export InstantRFPulseBlock, ConstantPulse, SincPulse, FixedPulse

import .Gradients: PulsedGradient, InstantGradientBlock, FixedGradient
export PulsedGradient, InstantGradientBlock, FixedGradient

import .Readouts: InstantReadout
export InstantReadout

import .Containers: Sequence, FixedBlock, FixedSequence
export Sequence, FixedBlock, FixedSequence

import .Pathways: Pathway, duration_transverse, duration_dephase, bval, bmat
export duration_transverse, duration_dephase, bval, bmat

import JuMP: @constraint, @objective, objective_function, optimize!, has_values, value, owner_model, Model
export @constraint, @objective, objective_function, optimize!, has_values, value, owner_model, Model

end
