module Pulses
include("properties.jl")
include("instant_pulses.jl")
include("constant_pulses.jl")
include("sinc_pulses.jl")

import .Properties: flip_angle, phase, amplitude, frequency, bandwidth
import .InstantPulses: InstantRFPulseBlock
import .ConstantPulses: ConstantPulse
import .SincPulses: SincPulse, N_left, N_right

end