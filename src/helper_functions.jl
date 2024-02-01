module HelperFunctions
import JuMP: @constraint
import ..BuildSequences: get_global_model
import ..Sequences: Sequence
import ..Pulses: SincPulse, ConstantPulse, InstantRFPulseBlock
import ..Overlapping: TrapezoidGradient
import ..Variables: qvec

"""
    excitation_pulse(; parameters...)

Create an excitation RF pulse.

By default there is no slice-selective gradient. 
To enable slice selection `min_slice_thickness` has to be set to a number or to :min.
if `min_slice_thickness` is not set or is set to `:min`, then either `bandwidth` or `duration` should be set, otherwise the problem might be unconstrained (ignore this for `shape=:instant`).

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
- `min_slice_thickness`: minimum slice thickness that should be possible without adjusting the sequence timings in um (default: no slice selection). Can be set to `:min` to indicate that this should be minimised given the scanner constraints and user values for `bandwidth` or `duration`.
- `rephase`: set to false to disable the spin rephasing after the RF pulse.
- `rotate_grad`: name of the parameter with which the slice selection gradient will be rotated after sequence optimisation (default: `:FOV`).
"""
function excitation_pulse(; flip_angle=90, phase=0., frequency=0., shape=:sinc, min_slice_thickness=Inf, rephase=true, Nzero=3, scale=:transmit_B1, rotate_grad=:FOV, bandwidth=nothing, duration=nothing)
    if shape == :sinc
        pulse = SincPulse(flip_angle=flip_angle, phase=phase, frequency=frequency, N_lobes=Nzero, scale=scale, bandwidth=bandwidth, duration=duration)
    elseif shape in (:constant, :hard)
        pulse = ConstantPulse(flip_angle=flip_angle, phase=phase, frequency=frequency, scale=scale, bandwidth=bandwidth, duration=duration)
    elseif shape == :instant
        pulse = InstantRFPulseBlock(flip_angle=flip_angle, phase=phase, scale=scale)
        if !isinf(min_slice_thickness)
            error("An instant RF pulse always affects all spins equally, so using `shape=:instant` is incompatible with setting `min_slice_thickness`.")
        end
        return pulse
    end
    if isinf(min_slice_thickness)
        return pulse
    end

    grad = TrapezoidGradient(pulse=pulse, rise_time=:min, slice_thickness=[Inf, Inf, min_slice_thickness], orientation=[0, 0, 1.], rotate=rotate_grad)
    if !rephase
        return grad
    end
    seq = Sequence(
        grad,
        TrapezoidGradient(orientation=[0, 0, -1.], rotate=rotate_grad, duration=:min);
        TR=Inf
    )
    @constraint get_global_model() qvec(grad, 1, nothing)[2] == -qvec(seq[2])[2]
    return seq
end


end