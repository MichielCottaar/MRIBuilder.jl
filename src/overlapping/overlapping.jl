module Waveforms
include("abstract.jl")
include("generic.jl")
include("trapezoid_gradients.jl")

import .Abstract: AbstractOverlapping, interruptions, waveform
import .Generic: GenericOverlapping
import .TrapezoidGradients: TrapezoidGradient

end