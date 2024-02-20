"""
Builds and optimises NMR/MRI sequences.
"""
module MRIBuilder

include("scanners.jl")
include("build_sequences.jl")
include("variables.jl")
include("components/components.jl")
include("containers/containers.jl")
include("pathways.jl")
include("parts/parts.jl")
#include("printing.jl")

import .BuildSequences: build_sequence, global_model, global_scanner
export build_sequence, global_model, global_scanner

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

import .Variables: variables, duration, effective_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength, qvec, qval_square, slice_thickness, inverse_slice_thickness, fov, inverse_fov, voxel_size, inverse_voxel_size, resolution, nsamples, oversample, dwell_time, ramp_overlap, spoiler_scale, TR, Δ
export variables, duration, effective_time, flip_angle, amplitude, phase, frequency, bandwidth, N_left, N_right, qval, δ, rise_time, flat_time, slew_rate, gradient_strength, qvec, qval_square, slice_thickness, inversne_slice_thickness, fov, inverse_fov, voxel_size, inverse_voxel_size, resolution, nsamples, oversample, dwell_time, ramp_overlap, spoiler_scale, TR, Δ

import .Components: InstantPulse, ConstantPulse, SincPulse, GenericPulse, InstantGradient, SingleReadout, ADC
export InstantPulse, ConstantPulse, SincPulse, GenericPulse, InstantGradient, SingleReadout, ADC

import .Containers: ContainerBlock, start_time, end_time, waveform, waveform_sequence, events, BaseBuildingBlock, BuildingBlock, Wait, BaseSequence, nrepeat, Sequence, AlternativeBlocks, match_blocks!, get_index_single_TR
export ContainerBlock, start_time, end_time, waveform, waveform_sequence, events, BaseBuildingBlock, BuildingBlock, Wait, BaseSequence, nrepeat, Sequence, AlternativeBlocks, match_blocks!, get_index_single_TR

import .Pathways: Pathway, duration_transverse, duration_dephase, bval, bmat
export Pathway, duration_transverse, duration_dephase, bval, bmat

import .Parts: excitation_pulse, refocus_pulse, epi_readout, single_line_readout, Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines, SpoiltSliceSelect, SliceSelectRephase, EPIReadout
export excitation_pulse, refocus_pulse, epi_readout, single_line_readout, Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines, SpoiltSliceSelect, SliceSelectRephase, EPIReadout

import JuMP: @constraint, @objective, objective_function, value, Model
export @constraint, @objective, objective_function, value, Model

end
