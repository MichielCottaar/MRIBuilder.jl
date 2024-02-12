module AbstracTypes
import ..Variables: AbstractBlock, duration, variables, effective_time

"""
Super-type for all individual components that form an MRI sequence (i.e., RF pulse, gradient waveform, or readout event).

RF pulses, instant gradients, and readouts are grouped together into [`EventComponent`](@ref).

These all should have a [`duration`](@ref) in addition to any other relevant [`variables`](@ref).
"""
abstract type BaseComponent <: AbstractBlock end

"""
Super-type for all parts of a gradient waveform.

N should be 1 for a 1D gradient waveform or 3 for a 3D one.
"""
abstract type GradientWaveform{N} <: BaseComponent end


"""
Super-type for all RF pulses, instant gradients and readouts that might play out during a gradient waveform.

These all have an [`effective_time`](@ref), which should quantify at what single time one can approximate the RF pulse or readout to have taken place.
"""
abstract type EventComponent <: BaseComponent end

"""
Super type for all RF pulses.
"""
abstract type RFPulseComponent <: EventComponent end

"""
Super type for all readout events.
"""
abstract type ReadoutComponent <: EventComponent end

end