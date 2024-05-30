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
include("post_hoc.jl")
include("sequences/sequences.jl")
include("printing.jl")
include("sequence_io/sequence_io.jl")
include("plot.jl")

import .BuildSequences: build_sequence, global_model, global_scanner, fixed
export build_sequence, global_model, global_scanner, fixed

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra, Default_Scanner
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra, Default_Scanner

import .Variables: variables, make_generic, @defvar, get_pulse, get_readout, get_pathway, get_gradient
export variables, make_generic, @defvar, get_pulse, get_readout, get_pathway, get_gradient

import .Components: InstantPulse, ConstantPulse, SincPulse, GenericPulse, InstantGradient, SingleReadout, ADC, CompositePulse, edge_times
export InstantPulse, ConstantPulse, SincPulse, GenericPulse, InstantGradient, SingleReadout, ADC, CompositePulse, edge_times

import .Containers: ContainerBlock, start_time, end_time, waveform, waveform_sequence, events, BaseBuildingBlock, BuildingBlock, Wait, BaseSequence, nrepeat, Sequence, AlternativeBlocks, match_blocks!, get_index_single_TR, iter_blocks, iter_instant_gradients, iter_instant_pulses
export ContainerBlock, start_time, end_time, waveform, waveform_sequence, events, BaseBuildingBlock, BuildingBlock, Wait, BaseSequence, nrepeat, Sequence, AlternativeBlocks, match_blocks!, get_index_single_TR, iter_blocks, iter_instant_gradients, iter_instant_pulses

import .Pathways: Pathway
export Pathway

import .Parts: dwi_gradients, readout_event, excitation_pulse, refocus_pulse, Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines, SpoiltSliceSelect, SliceSelectRephase, EPIReadout, interpret_image_size
export dwi_gradients, readout_event, excitation_pulse, refocus_pulse, Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines, SpoiltSliceSelect, SliceSelectRephase, EPIReadout, interpret_image_size

import .PostHoc: adjust, merge_sequences
export adjust, merge_sequences

import .Sequences: GradientEcho, SpinEcho, DiffusionSpinEcho, DW_SE, DWI
export GradientEcho, SpinEcho, DiffusionSpinEcho, DW_SE, DWI

import .SequenceIO: read_sequence, write_sequence
export read_sequence, write_sequence

import .Plot: plot_sequence
export plot_sequence

import JuMP: @constraint, @objective, objective_function, value, Model
export @constraint, @objective, objective_function, value, Model

end
