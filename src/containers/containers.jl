module Containers
include("abstract.jl")
include("building_blocks.jl")
include("base_sequences.jl")
include("alternatives.jl")

import .Abstract: ContainerBlock, start_time, end_time, iter_blocks, iter_instant_gradients, iter_instant_pulses
import .BuildingBlocks: BaseBuildingBlock, BuildingBlock, Wait, waveform, waveform_sequence, events, ndim_grad
import .BaseSequences: BaseSequence, Sequence, nrepeat, get_index_single_TR
import .Alternatives: AbstractAlternativeBlocks, get_alternatives_name, get_alternatives_options, AlternativeBlocks, match_blocks!

end