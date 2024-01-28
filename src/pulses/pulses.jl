module Pulses
include("fixed_pulses.jl")
include("instant_pulses.jl")
include("constant_pulses.jl")
include("sinc_pulses.jl")

import ..Variables: effective_time
import .FixedPulses: FixedPulse
import .InstantPulses: InstantRFPulseBlock
import .ConstantPulses: ConstantPulse
import .SincPulses: SincPulse

end