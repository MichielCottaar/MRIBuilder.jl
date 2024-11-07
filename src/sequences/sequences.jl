module Sequences
include("gradient_echoes.jl")
include("spin_echoes.jl")
include("diffusion_spin_echoes.jl")
include("magnetisation_transfers.jl")
include("multi_spin_echoes.jl")

import .GradientEchoes: GradientEcho
import .SpinEchoes: SpinEcho
import .DiffusionSpinEchoes: DiffusionSpinEcho, DW_SE, DWI
import .MagnetisationTransfers: MagnetisationTransfer
import .MultiSpinEchoes: MultiSpinEcho

end