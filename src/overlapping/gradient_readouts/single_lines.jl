module SingleLines
import JuMP: @constraint
import ....BuildingBlocks: set_simple_constraints!
import ....BuildSequences: global_model
import ....Readouts: ADC
import ....Variables: dwell_time, inverse_fov, resolution, nsamples, effective_time, gradient_strength, slew_rate, inverse_voxel_size, rise_time, duration, variables, fov, voxel_size, flat_time, ramp_overlap
import ...Abstract: AbstractOverlapping, waveform, interruptions
import ...GradientPulses: TrapezoidGradient

struct SingleLine <: AbstractOverlapping
    adc :: ADC
    grad :: TrapezoidGradient
    ramp_overlap :: VariableType
end

"""
    SingleLine(adc; rotate=:FOV, orientation=[1, 0, 0], ramp_overlap=1., variables...)

Adds a [`TrapezoidGradient`](@ref) to underly the [`ADC`](@ref) event.

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
function SingleLine(adc; rotate=:FOV, orientation=[1, 0, 0], ramp_overlap=1., kwargs...)
    res = SingleLine(
        adc,
        TrapezoidGradient(orientation=orientation, duration=:min, rotate=rotate),
        get_free_variable(ramp_overlap),
    )
    @constraint global_model() (res.ramp_overlap * rise_time(res.grad) + flat_time(res.grad)) == duration(adc)
    set_simple_constraints!(res, kwargs)
    return res
end

waveform(sl::SingleLine) = waveform(sl.grad)
interruptions(sl::SingleLine) = [(index=2, time=effective_time(sl.adc), object=sl.adc)]

ramp_overlap(sl::SingleLine) = sl.ramp_overlap
effective_time(sl::SingleLine) = rise_time(sl.grad) + effective_time(sl.adc)
dwell_time(sl::SingleLine) = dwell_time(sl.adc)
slew_rate(sl::SingleLine) = slew_rate(sl.grad)[1]
gradient_strenth(sl::SingleLine) = slew_rate(sl) * rise_time(sl.grad)
inverse_fov(sl::SingleLine) = 1e3 * dwell_time(sl) * gradient_strenth(sl)
inverse_voxel_size(sl::SingleLine) = 1e3 * duration(sl.adc) * gradient_strenth(sl)
nsamples(sl::SingleLine) = nsamples(sl.adc)
resolution(sl::SingleLine) = resolution(sl.adc)
duration(sl::SingleLine) = duration(sl.grad)

end