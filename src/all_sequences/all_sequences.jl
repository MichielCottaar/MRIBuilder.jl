module AllSequences
include("base_sequences.jl")
include("sequences.jl")
include("slice_select_rephases.jl")
include("epi_readouts.jl")

import .BaseSequences: BaseSequence, nrepeat
import .Sequences: Sequence
import .SliceSelectRephases: SliceSelectRephase
import .EPIReadouts: EPIReadout


end