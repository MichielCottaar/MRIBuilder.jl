module InstantReadouts
import ...AbstractTypes: ReadoutComponent
import ...Variables: variables, duration, effective_time
import ..ADCs: ADC

"""
    SingleReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or properties to set.
"""
struct SingleReadout <: ReadoutComponent
end

duration(::SingleReadout) = 0.
effective_time(::SingleReadout) = 0.

make_generic(::SingleReadout) = ADC(1, 0., 0., 0.)
end