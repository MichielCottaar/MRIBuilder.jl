module Readouts
include("ADCs.jl")
include("single_readouts.jl")

import ..AbstractTypes: ReadoutComponent, split_timestep
import .ADCs: ADC
import .SingleReadouts: SingleReadout

split_timestep(rc::ReadoutComponent, precision) = Inf
end