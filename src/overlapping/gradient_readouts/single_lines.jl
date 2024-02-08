module SingleLines
import JuMP: @constraint
import ....BuildingBlocks: set_simple_constraints!
import ....BuildSequences: global_model
import ....Readouts: ADC
import ....Variables: dwell_time, inverse_fov, resolution, nsamples, effective_time, gradient_strength, slew_rate, inverse_voxel_size, rise_time, duration, variables, fov, voxel_size, flat_time, ramp_overlap, VariableType, get_free_variable
import ...Abstract: AbstractOverlapping, waveform, interruptions
import ...GradientPulses: TrapezoidGradient

struct SingleLine <: AbstractOverlapping
    adc :: ADC
    grad :: TrapezoidGradient
    ramp_overlap :: VariableType
end

"""
    SingleLine(rotate=:FOV, orientation=[1, 0, 0], ramp_overlap=1., variables...)

Overlays a [`TrapezoidGradient`](@ref) and [`ADC`](@ref) event.

## Parameters
- `rotate`: with which user-set parameter will this gradient be rotated (default: :FOV).
- `orientation`: gradient orientation before rotation (default: [1, 0, 0]).
- `ramp_overlap`: how much the gradient ramp should overlap with the ADC. 0 for no overlap, 1 for full overlap (default: 1). Can be set to `nothing` to become a free variable.

## Variables
New variables added at this level:
- [`fov`](@ref): FOV of the output image along this single k-space line in mm.
- [`voxel_size`](@ref): size of each voxel along this single k-space line in mm.

ADC variables:
- [`resolution`](@ref): number of voxels in the readout direction. This can be a non-integer value during optimisation.
- [`nsamples`](@ref): number of samples in the readout. This can be a non-integer value during optimisation.
- [`dwell_time`](@ref): Time between each readout sample in ms.
- [`effective_time`](@ref): time till the center of k-space from the start of the gradient (might not be the same as the start of the ADC).

Gradient variables:
- [`gradient_strength`](@ref): Maximum gradient strength achieved during the pulse in kHz/um
"""
function SingleLine(; rotate=:FOV, orientation=[1, 0, 0], ramp_overlap=1., resolution=nothing, oversample=nothing, nsamples=nothing, dwell_time=nothing, kwargs...)
    res = SingleLine(
        ADC(; resolution=resolution, oversample=oversample, nsamples=nsamples, dwell_time=dwell_time),
        TrapezoidGradient(orientation=orientation, duration=:min, rotate=rotate),
        get_free_variable(ramp_overlap),
    )
    @constraint global_model() (res.ramp_overlap * rise_time(res.grad) * 2 + flat_time(res.grad)) == duration(res.adc)
    set_simple_constraints!(res, kwargs)
    return res
end

waveform(sl::SingleLine) = waveform(sl.grad)
interruptions(sl::SingleLine) = [(index=2, time=effective_time(sl.adc), object=sl.adc)]

ramp_overlap(sl::SingleLine) = sl.ramp_overlap
effective_time(sl::SingleLine) = (1. - ramp_overlap(sl)) * rise_time(sl.grad) + effective_time(sl.adc)
dwell_time(sl::SingleLine) = dwell_time(sl.adc)
slew_rate(sl::SingleLine) = slew_rate(sl.grad)[1]
gradient_strenth(sl::SingleLine) = slew_rate(sl) * rise_time(sl.grad)
inverse_fov(sl::SingleLine) = 1e3 * dwell_time(sl) * gradient_strenth(sl)
inverse_voxel_size(sl::SingleLine) = 1e3 * duration(sl.adc) * gradient_strenth(sl)
nsamples(sl::SingleLine) = nsamples(sl.adc)
resolution(sl::SingleLine) = resolution(sl.adc)
duration(sl::SingleLine) = duration(sl.grad)


"""
    opposite_kspace_lines(rotate=:FOV, orientation=[1, 0, 0], ramp_overlap=1., variables...)

Create two [`SingleLine`](@ref) objects that run in the opposite direction through k-space.

The first will run in the `+orientation` direction, the other in the `-orienation` direction.

Parameters and variables are the same as for [`SingleLine`](@ref).
"""
function opposite_kspace_lines(; kwargs...)
    positive_line = SingleLine(; kwargs...)
    grad = positive_line.grad
    negative_grad = TrapezoidGradient(
        grad.rise_time,
        grad.flat_time,
        -grad.slew_rate,
        grad.rotate,
        grad.scale,
        grad.time_before_pulse,
        grad.pulse,
        grad.time_after_pulse,
        grad.slew_rate_1d
    )
    negative_line = SingleLine(
        positive_line.adc,
        negative_grad,
        positive_line.ramp_overlap
    )
    return (positive_line, negative_line)
end

end