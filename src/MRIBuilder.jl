"""
Builds and optimises NMR/MRI sequences.
"""
module MRIBuilder

include("scanners.jl")
include("build_sequences.jl")
include("variables.jl")
include("building_blocks.jl")
include("wait.jl")
include("gradients/gradients.jl")
include("pulses/pulses.jl")
include("readouts/readouts.jl")
include("overlapping/overlapping.jl")
include("sequences.jl")
include("alternatives.jl")
include("pathways.jl")
include("helper_functions.jl")
include("printing.jl")

import .BuildSequences: build_sequence, global_model, global_scanner
export build_sequence, global_model, global_scanner

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

import .Variables: variables, duration, start_time, end_time, effective_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength, qvec, qval_square, slice_thickness, inverse_slice_thickness, fov, inverse_fov, voxel_size, inverse_voxel_size, resolution
export variables, duration, start_time, end_time, effective_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength, qvec, qval_square, slice_thickness, inversne_slice_thickness, fov, inverse_fov, voxel_size, inverse_voxel_size, resolution

import .BuildingBlocks: BuildingBlock, fixed, get_children_blocks, get_children_indices, make_generic
export BuildingBlocks, fixed, get_children_blocks, get_children_indices, make_generic

import .Wait: WaitBlock
export WaitBlock

import .Pulses: InstantRFPulseBlock, ConstantPulse, SincPulse, GenericPulse
export InstantRFPulseBlock, ConstantPulse, SincPulse, GenericPulse

import .Gradients: InstantGradientBlock
export InstantGradientBlock

import .Readouts: InstantReadout, ADC
export InstantReadout, ADC

import .Overlapping: TrapezoidGradient, SpoiltSliceSelect, interruptions, waveform, SingleLine, opposite_kspace_lines, waveform_sequence
export TrapezoidGradient, SpoiltSliceSelect, interruptions, waveform, SingleLine, opposite_kspace_lines, waveform_sequence

import .Sequences: Sequence
export Sequence

import .Alternatives: AlternativeBlocks
export AlternativeBlocks

import .Pathways: Pathway, duration_transverse, duration_dephase, bval, bmat
export Pathway, duration_transverse, duration_dephase, bval, bmat

import .HelperFunctions: excitation_pulse, refocus_pulse, epi_readout, single_line_readout
export excitation_pulse, refocus_pulse, epi_readout, single_line_readout

import JuMP: @constraint, @objective, objective_function, value, Model
export @constraint, @objective, objective_function, value, Model

end
