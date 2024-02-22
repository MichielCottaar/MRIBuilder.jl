module Containers
include("abstract.jl")
include("building_blocks.jl")
include("base_sequences.jl")
include("alternatives.jl")

import .Abstract: ContainerBlock, start_time, end_time, effective_time
import .BuildingBlocks: BaseBuildingBlock, BuildingBlock, Wait, waveform, waveform_sequence, events
import .BaseSequences: BaseSequence, Sequence, nrepeat, get_index_single_TR
import .Alternatives: AlternativeBlocks, match_blocks!

end