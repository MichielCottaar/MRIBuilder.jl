module AllSequences
include("base_sequences.jl")
include("sequences.jl")
include("slice_select_rephases.jl")

import .BaseSequences: BaseSequence, nrepeat
import .Sequences: Sequence
import .SliceSelectRephases: SliceSelectRephase

end