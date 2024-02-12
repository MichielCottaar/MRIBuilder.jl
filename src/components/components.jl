module Components
include("abstract_types.jl")
include("gradient_waveforms/gradient_waveforms.jl")
include("instant_gradients.jl")
include("pulses/pulses.jl")
include("readouts/readouts.jl")

import .AbstractTypes: BaseComponent, GradientWaveform, EventComponent, RFPulseComponent, ReadoutComponent
import .GradientWaveforms: ConstantGradientBlock, ChangingGradientBlock, NoGradientBlock
import .InstantGradients: InstantGradient1D, InstantGradient3D
import .Pulses: GenericPulse, InstantRFPulse, SincPulse, ConstantPulse
import .Readouts: ADC, SingleReadout

end