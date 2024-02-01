"""
Builds and optimises NMR/MRI sequences.
"""
module MRIBuilder

include("variables.jl")
include("building_blocks.jl")
include("scanners.jl")
include("build_sequences.jl")
include("wait.jl")
include("gradients/gradients.jl")
include("pulses/pulses.jl")
include("readouts/readouts.jl")
include("overlapping/overlapping.jl")
include("sequences.jl")
include("pathways.jl")
include("helper_functions.jl")

import .BuildSequences: build_sequence
export build_sequence

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

import .Variables: variables, duration, start_time, end_time, effective_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength, qvec, qval_square
export variables, duration, start_time, end_time, effective_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength, qvec, qval_square

import .BuildingBlocks: BuildingBlock, fixed
export BuildingBlocks, fixed

import .Wait: WaitBlock
export WaitBlock

import .Pulses: InstantRFPulseBlock, ConstantPulse, SincPulse, FixedPulse
export InstantRFPulseBlock, ConstantPulse, SincPulse, FixedPulse

import .Gradients: InstantGradientBlock
export InstantGradientBlock

import .Readouts: InstantReadout
export InstantReadout

import .Overlapping: TrapezoidGradient, SpoiltSliceSelect, interruptions, waveform
export TrapezoidGradient, SpoiltSliceSelect, interruptions, waveform

import .Sequences: Sequence
export Sequence

import .Pathways: Pathway, duration_transverse, duration_dephase, bval, bmat
export Pathway, duration_transverse, duration_dephase, bval, bmat

import .HelperFunctions: excitation_pulse, refocus_pulse
export excitation_pulse, refocus_pulse

import JuMP: @constraint, @objective, objective_function, optimize!, has_values, value, owner_model, Model
export @constraint, @objective, objective_function, optimize!, has_values, value, owner_model, Model

end
