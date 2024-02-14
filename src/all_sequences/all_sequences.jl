module AllSequences
include("base_sequences.jl")
include("sequences.jl")

import .BaseSequences: BaseSequence, nrepeat
import .Sequences: Sequence

end