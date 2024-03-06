module Containers
include("abstract.jl")
include("building_blocks.jl")
include("base_sequences.jl")
include("alternatives.jl")
include("linearise.jl")

import .Abstract: ContainerBlock, start_time, end_time, effective_time, readout_times, edge_times
import .BuildingBlocks: BaseBuildingBlock, BuildingBlock, Wait, waveform, waveform_sequence, events, ndim_grad
import .BaseSequences: BaseSequence, Sequence, nrepeat, get_index_single_TR
import .Alternatives: AlternativeBlocks, match_blocks!
import .Linearise: SequencePart, split_times, split_timestep, linearise

end