module Sequences
include("gradient_echoes.jl")
include("spin_echoes.jl")

import .GradientEchoes: GradientEcho
import .SpinEchoes: SpinEcho

end