module GradientPulses
include("trapezoid_gradients.jl")
include("spoilt_slice_selects.jl")

import .TrapezoidGradients: TrapezoidGradient
import .SpoiltSliceSelects: SpoiltSliceSelect

end