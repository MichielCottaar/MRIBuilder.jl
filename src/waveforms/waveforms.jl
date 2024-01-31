module Waveforms
include("generic.jl")
include("trapezoid_gradients.jl")

import .Generic: GenericWaveform, AbstractWaveform, SpecificWaveform
import .TrapezoidGradients: TrapezoidGradient

end