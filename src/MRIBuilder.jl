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

import .BuildSequences: build_sequence, global_scanner, fixed
export build_sequence, global_scanner, fixed

import .Scanners: Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra, Default_Scanner
export Scanner, B0, Siemens_Connectom, Siemens_Prisma, Siemens_Terra, Default_Scanner

import .Variables: variables, make_generic, @defvar, get_pulse, get_readout, get_pathway, get_gradient, add_cost_function!, apply_simple_constraint!, set_simple_constraints!
export variables, make_generic, @defvar, get_pulse, get_readout, get_pathway, get_gradient, add_cost_function!, apply_simple_constraint!, set_simple_constraints!

import .Components: InstantPulse, ConstantPulse, SincPulse, GenericPulse, InstantGradient, SingleReadout, ADC, CompositePulse, edge_times, BinomialPulse
export InstantPulse, ConstantPulse, SincPulse, GenericPulse, InstantGradient, SingleReadout, ADC, CompositePulse, edge_times, BinomialPulse

import .Containers: ContainerBlock, start_time, end_time, waveform, waveform_sequence, events, BaseBuildingBlock, BuildingBlock, Wait, BaseSequence, nrepeat, Sequence, AbstractAlternativeBlocks, get_alternatives_name, get_alternatives_options, AlternativeBlocks, match_blocks!, get_index_single_TR, iter_blocks, iter_instant_gradients, iter_instant_pulses
export ContainerBlock, start_time, end_time, waveform, waveform_sequence, events, BaseBuildingBlock, BuildingBlock, Wait, BaseSequence, nrepeat, Sequence, AbstractAlternativeBlocks, get_alternatives_name, get_alternatives_options, AlternativeBlocks, match_blocks!, get_index_single_TR, iter_blocks, iter_instant_gradients, iter_instant_pulses

import .Pathways: Pathway
export Pathway

import .Parts: dwi_gradients, readout_event, excitation_pulse, refocus_pulse, Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines, SpoiltSliceSelect, SliceSelectRephase, EPIReadout, interpret_image_size, saturation_pulse
export dwi_gradients, readout_event, excitation_pulse, refocus_pulse, Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines, SpoiltSliceSelect, SliceSelectRephase, EPIReadout, interpret_image_size, saturation_pulse

import .PostHoc: adjust, merge_sequences
export adjust, merge_sequences

import .Sequences: GradientEcho, SpinEcho, DiffusionSpinEcho, DW_SE, DWI, MagnetisationTransfer, MultiSpinEcho
export GradientEcho, SpinEcho, DiffusionSpinEcho, DW_SE, DWI, MagnetisationTransfer, MultiSpinEcho

import .SequenceIO: read_sequence, write_sequence
export read_sequence, write_sequence

import .Plot: plot_sequence
export plot_sequence

end
