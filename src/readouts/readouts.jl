module Readouts
include("instant_readouts.jl")
include("ADCs.jl")

import .InstantReadouts: InstantReadout
import .ADCs: ADC
end