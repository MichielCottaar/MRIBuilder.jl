module HelperFunctions
import JuMP: @constraint
import ..BuildSequences: global_model
import ..Sequences: Sequence
import ..Pulses: SincPulse, ConstantPulse, InstantRFPulseBlock
import ..Overlapping: TrapezoidGradient, SpoiltSliceSelect
import ..Variables: qvec


function _get_pulse(shape, flip_angle, phase, frequency, Nzero, scale, bandwidth, duration)
    if shape == :sinc
        pulse = SincPulse(flip_angle=flip_angle, phase=phase, frequency=frequency, N_lobes=Nzero, scale=scale, bandwidth=bandwidth, duration=duration)
    elseif shape in (:constant, :hard)
        pulse = ConstantPulse(flip_angle=flip_angle, phase=phase, frequency=frequency, scale=scale, bandwidth=bandwidth, duration=duration)
    elseif shape == :instant
        pulse = InstantRFPulseBlock(flip_angle=flip_angle, phase=phase, scale=scale)
    end
    return pulse
end

"""
    excitation_pulse(; parameters...)

Create an excitation RF pulse.

By default there is no slice-selective gradient. 
To enable slice selection `min_slice_thickness` has to be set to a number or to :min.
If `min_slice_thickness` is not set or is set to `:min`, then either `bandwidth` or `duration` should be set, otherwise the optimisation might be unconstrained (ignore this for `shape=:instant`).

## Parameters
For an [`InstantRFPulseBlock`](@ref) (i.e., `shape=:instant`), only the `flip_angle`, `phase`, and `scale` will be used. All other parameters are ignored.
### Pulse parameters
- `shape`: The shape of the RF pulse. One of `:sinc` (for [`SincPulse`](@ref)), `:constant`/`:hard` (for [`ConstantPulse`](@ref)), or `:instant` (for [`InstantRFPulseBlock`](@ref)).
- `flip_angle`: size of the flip due to the RF pulse in degrees (default: 90).
- `phase`: angle of the RF pulse in the x-y plane in degrees (default: 0).
- `frequency`: frequency of the RF pulse relative to the Larmor frequency in kHz (default: 0).
- `bandwidth`: width of the RF pulse in Fourier space in kHz (default: free variable).
- `duration`: duration of the RF pulse in ms (default: free variable).
- `Nzero`: sets the number of zero crossings for a [`SincPulse`](@ref) (default: 3). Can be set to a tuple of two numbers to set a different number of zero crossings before and after the pulse maximum.
- `scale`: name of the parameter with which the RF pulse amplitude can be modulated after sequence optimisation (default: `:transmit_B1`).

### Slice selection
- `slice_thickness`: minimum slice thickness that should be possible without adjusting the sequence timings in um (not mm!) (default: no slice selection). Can be set to `:min` to indicate that this should be minimised given the scanner constraints and user values for `bandwidth` or `duration`.
- `rephase`: set to false to disable the spin rephasing after the RF pulse.
- `rotate_grad`: name of the parameter with which the slice selection gradient will be rotated after sequence optimisation (default: `:FOV`).
"""
function excitation_pulse(; flip_angle=90, phase=0., frequency=0., shape=:sinc, slice_thickness=Inf, rephase=true, Nzero=3, scale=:transmit_B1, rotate_grad=:FOV, bandwidth=nothing, duration=nothing)
    pulse = _get_pulse(shape, flip_angle, phase, frequency, Nzero, scale, bandwidth, duration)
    if pulse isa InstantRFPulseBlock
        if !isinf(slice_thickness)
            error("An instant RF pulse always affects all spins equally, so using `shape=:instant` is incompatible with setting `slice_thickness`.")
        end
        return pulse
    end
    if isinf(slice_thickness)
        return pulse
    end

    grad = TrapezoidGradient(pulse=pulse, duration=:min, slice_thickness=[Inf, Inf, slice_thickness], orientation=[0, 0, 1.], rotate=rotate_grad)
    if !rephase
        return grad
    end
    seq = Sequence(
        grad,
        TrapezoidGradient(orientation=[0, 0, -1.], rotate=rotate_grad, duration=:min);
        TR=Inf
    )
    @constraint global_model() qvec(grad, 1, nothing)[3] == -qvec(seq[2])[3]
    return seq
end



"""
    refocus_pulse(; parameters...)

Create an excitation RF pulse.

By default there is no slice-selective gradient. 
To enable slice selection `slice_thickness` has to be set to a number or to :min.
If `slice_thickness` is not set or is set to `:min`, then either `bandwidth` or `duration` should be set, otherwise the optimisation might be unconstrained (ignore this for `shape=:instant`).

## Parameters
For an [`InstantRFPulseBlock`](@ref) (i.e., `shape=:instant`), only the `flip_angle`, `phase`, and `scale` will be used. All other parameters are ignored.
### Pulse parameters
- `shape`: The shape of the RF pulse. One of `:sinc` (for [`SincPulse`](@ref)), `:constant`/`:hard` (for [`ConstantPulse`](@ref)), or `:instant` (for [`InstantRFPulseBlock`](@ref)).
- `flip_angle`: size of the flip due to the RF pulse in degrees (default: 180).
- `phase`: angle of the RF pulse in the x-y plane in degrees (default: 0).
- `frequency`: frequency of the RF pulse relative to the Larmor frequency in kHz (default: 0).
- `bandwidth`: width of the RF pulse in Fourier space in kHz (default: free variable).
- `duration`: duration of the RF pulse in ms (default: free variable).
- `Nzero`: sets the number of zero crossings for a [`SincPulse`](@ref) (default: 3). Can be set to a tuple of two numbers to set a different number of zero crossings before and after the pulse maximum.
- `scale`: name of the parameter with which the RF pulse amplitude can be modulated after sequence optimisation (default: `:transmit_B1`).

### Slice selection and spoilers
- `slice_thickness`: minimum slice thickness that should be possible without adjusting the sequence timings in um (not mm!) (default: no slice selection). Can be set to `:min` to indicate that this should be minimised given the scanner constraints and user values for `bandwidth` or `duration`.
- `spoiler`: set to the spatial scale on which the spins should be dephased in mm. For rotating spoilers, this does include the contribution from the slice select gradient as well.
- `rotate_grad`: name of the parameter with which the slice selection and spoiler gradient will be rotated after sequence optimisation (default: `:FOV`).
"""
function refocus_pulse(; flip_angle=180, phase=0., frequency=0., shape=:sinc, slice_thickness=Inf, Nzero=3, scale=:transmit_B1, rotate_grad=:FOV, bandwidth=nothing, duration=nothing, spoiler=Inf)
    pulse = _get_pulse(shape, flip_angle, phase, frequency, Nzero, scale, bandwidth, duration)
    if pulse isa InstantRFPulseBlock && !isinf(slice_thickness)
        error("An instant RF pulse always affects all spins equally, so using `shape=:instant` is incompatible with setting `slice_thickness`.")
    end
    if isinf(spoiler)
        if isinf(slice_thickness)
            return pulse
        else
            return TrapezoidGradient(pulse=pulse, duration=:min, slice_thickness=[Inf, Inf, slice_thickness], orientation=[0, 0, 1.], rotate=rotate_grad)
        end
    else
        orientation=isnothing(rotate_grad) ? [1, 1, 1] : [0, 0, 1]
        if isinf(slice_thickness)
            grad = TrapezoidGradient(orientation=orientation, duration=:min, rotate=rotate_grad)
            @constraint global_model() qvec(grad)[3] == 2π * 1e-3 / spoiler
            return Sequence(grad, pulse, grad; TR=Inf)
        else
            res = SpoiltSliceSelect(pulse; orientation=orientation, duration=:min, rotate=rotate_grad, slice_thickness=slice_thickness, spoiler_scale=spoiler)
            return res
        end
    end
end


end