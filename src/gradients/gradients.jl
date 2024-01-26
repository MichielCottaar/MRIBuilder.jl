module Gradients
include("integrate_gradients.jl")
include("pulsed_gradients.jl")
include("instant_gradients.jl")

import .PulsedGradients: PulsedGradient
import .InstantGradients: InstantGradientBlock
end