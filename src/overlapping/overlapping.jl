module Overlapping
include("abstract.jl")
include("generic.jl")
include("trapezoid_gradients.jl")
include("spoilt_slice_selects.jl")

import .Abstract: AbstractOverlapping, interruptions, waveform
import .Generic: GenericOverlapping
import .TrapezoidGradients: TrapezoidGradient
import .SpoiltSliceSelects: SpoiltSliceSelect

end