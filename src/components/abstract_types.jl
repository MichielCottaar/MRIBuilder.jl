module AbstractTypes
import ...Variables: AbstractBlock, variables, adjustable, gradient_orientation, @defvar

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

@defvar begin
    function slew_rate end
    function gradient_strength end
end

"""
    slew_rate(gradient)

Maximum 3D slew rate of the gradient in kHz/um/ms.
"""
variables.slew_rate

"""
    gradient_strength(gradient)

Maximum 3D gradient strength of the gradient in kHz/um.
"""
variables.gradient_strength

"""
Super-type for all RF pulses, instant gradients and readouts that might play out during a gradient waveform.

These all have an [`effective_time`](@ref), which should quantify at what single time one can approximate the RF pulse or readout to have taken place.
"""
abstract type EventComponent <: BaseComponent end

"""
Super type for all RF pulses.
"""
abstract type RFPulseComponent <: EventComponent end

@defvar pulse begin
    function phase end
    function amplitude end
    function flip_angle end
    function frequency end
end

"""
    phase(pulse)

Return the phase of an [`RFPulseComponent`](@ref) in degrees.
"""
variables.phase

"""
    amplitude(pulse)

Return the amplitude of an [`RFPulseComponent`](@ref) in kHz.
"""
variables.amplitude

"""
    frequency(pulse)

Return the off-resonance frequency of an [`RFPulseComponent`](@ref) in kHz.
"""
variables.frequency

"""
    flip_angle(pulse)

Return the flip angle of an [`RFPulseComponent`](@ref) in degrees.
"""
variables.flip_angle

"""
    bandwidth(pulse)

Return the bandwidth of an [`RFPulseComponent`](@ref) in kHz.
"""
variables.bandwidth

"""
Super type for all readout events.
"""
abstract type ReadoutComponent <: EventComponent end


"""
    split_timestep(component, precision)

Indicates the maximum timestep that a component can be linearised with and still achieve the required `precision`.

Typically, this will be determined by the maximum second derivative:

``\\sqrt{\\frac{2 \\epsilon}{max(|d^2y/dx^2|)}}``

It should be infinite if the component is linear.
"""
split_timestep(comp_tuple::Tuple{<:Number, <:EventComponent}, precision::Number) = split_timestep(comp_tuple[2], precision)


adjustable(::RFPulseComponent) = :pulse
adjustable(::GradientWaveform) = :gradient
gradient_orientation(gw::GradientWaveform{1}) = gw.orientation


"""
    edge_times(container/component; tol=1e-6)

Returns all the edge times during a sequence in ms.

Edges are defined as any time, when:
- the edge of a building block
- the slope of the gradient profile changes suddenly
- an RF pulse starts or ends

Edges that are within `tol` ms of each other are considered to be one edge (default: 1 ns).
"""
edge_times(comp::BaseComponent) = [0., variables.duration(comp)]

end