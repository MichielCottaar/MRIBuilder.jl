module Readouts
include("ADCs.jl")
include("single_readouts.jl")

import ..AbstractTypes: ReadoutComponent
import .ADCs: ADC, readout_times
import .SingleReadouts: SingleReadout
end