module Readouts
include("ADCs.jl")
include("single_readouts.jl")

import ..AbstractTypes: ReadoutComponent, split_timestep
import .ADCs: ADC, readout_times
import .SingleReadouts: SingleReadout

split_times(rc::ReadoutComponent, precision) = Inf
end