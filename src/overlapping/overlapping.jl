module Overlapping
include("abstract.jl")
include("generic.jl")
include("gradient_pulses//gradient_pulses.jl")

import .Abstract: AbstractOverlapping, interruptions, waveform
import .Generic: GenericOverlapping
import .GradientPulses: TrapezoidGradient, SpoiltSliceSelect

end