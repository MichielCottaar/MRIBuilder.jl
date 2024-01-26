module Pulses
include("instant_pulses.jl")
include("constant_pulses.jl")
include("sinc_pulses.jl")

import .InstantPulses: InstantRFPulseBlock
import .ConstantPulses: ConstantPulse
import .SincPulses: SincPulse

end