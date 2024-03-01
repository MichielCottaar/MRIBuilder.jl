module Components
include("abstract_types.jl")
include("gradient_waveforms/gradient_waveforms.jl")
include("instant_gradients.jl")
include("pulses/pulses.jl")
include("readouts/readouts.jl")
include("delayed_events.jl")

import .AbstractTypes: BaseComponent, GradientWaveform, EventComponent, RFPulseComponent, ReadoutComponent, split_times
import .GradientWaveforms: ConstantGradient, ChangingGradient, NoGradient, split_gradient
import .InstantGradients: InstantGradient
import .Pulses: GenericPulse, InstantPulse, SincPulse, ConstantPulse
import .Readouts: ADC, SingleReadout
import .DelayedEvents: DelayedEvent

end