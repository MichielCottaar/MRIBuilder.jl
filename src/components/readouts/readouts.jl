module Readouts
include("ADCs.jl")
include("single_readouts.jl")

import ..AbstractTypes: ReadoutComponent, split_times
import .ADCs: ADC, readout_times
import .SingleReadouts: SingleReadout

split_times(rc::ReadoutComponent, t1, t2, precision) = [t1, t2]
end