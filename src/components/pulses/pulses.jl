module Pulses
include("generic_pulses.jl")
include("instant_pulses.jl")
include("constant_pulses.jl")
include("sinc_pulses.jl")

import ..AbstractTypes: RFPulseComponent
import .GenericPulses: GenericPulse
import .InstantPulses: InstantPulse
import .ConstantPulses: ConstantPulse
import .SincPulses: SincPulse

end