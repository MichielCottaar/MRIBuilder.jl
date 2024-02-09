module Overlapping
include("abstract.jl")
include("generic.jl")
include("gradient_pulses/gradient_pulses.jl")
include("gradient_readouts/gradient_readouts.jl")

import .Abstract: AbstractOverlapping, interruptions, waveform, get_parts, waveform_sequence
import .Generic: GenericOverlapping
import .GradientPulses: TrapezoidGradient, SpoiltSliceSelect
import .GradientReadouts: SingleLine, opposite_kspace_lines

end