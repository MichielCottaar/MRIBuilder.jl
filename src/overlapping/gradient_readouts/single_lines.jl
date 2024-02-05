module SingleLines
import ....BuildingBlocks: set_simple_constraints!
import ....Readouts: ADC
import ....Variables: dwell_time, fov, resolution, nsamples, effective_time, gradient_strength, slew_rate, voxel_size, rise_time, duration, variables
import ...Abstract: AbstractOverlapping, waveform, interruptions
import ...GradientPulses: TrapezoidGradient

struct SingleLine <: AbstractOverlapping
    adc :: ADC
    grad :: TrapezoidGradient
end

function SingleLine(; nsamples=nothing, dwell_time=nothing, time_to_center=nothing, adc_duration=nothing, center_halfway=true, rotate=:FOV, kwargs...)
    readout = ADC(nsamples=nsamples, dwell_time=dwell_time, time_to_center=time_to_center, duration=adc_duration, center_halfway=center_halfway)
    res = SingleLine(
        readout,
        TrapezoidGradient(orientation=[1, 0, 0], flat_time=duration(readout), duration=:min, rotate=rotate)
    )
    set_simple_constraints!(res, kwargs)
    return res
end

waveform(sl::SingleLine) = waveform(sl.grad)
interruptions(sl::SingleLine) = [(index=2, time=effective_time(sl.adc), object=sl.adc)]

effective_time(sl::SingleLine) = rise_time(sl.grad) + effective_time(sl.adc)
dwell_time(sl::SingleLine) = dwell_time(sl.adc)
slew_rate(sl::SingleLine) = slew_rate(sl.grad)[1]
gradient_strenth(sl::SingleLine) = slew_rate(sl) * rise_time(sl.grad)
fov_inverse(sl::SingleLine) = 1e3 * dwell_time(sl) * gradient_strenth(sl)
voxel_size_inverse(sl::SingleLine) = 1e3 * duration(sl.adc) * gradient_strenth(sl)
nsamples(sl::SingleLine) = nsamples(sl.adc)
resolution(sl::SingleLine) = resolution(sl.adc)
duration(sl::SingleLine) = duration(sl.grad)

end