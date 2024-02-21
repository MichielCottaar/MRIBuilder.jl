module HelperFunctions
import JuMP: @constraint
import ...Containers: AlternativeBlocks, match_blocks!, BuildingBlock
import ..Trapezoids: Trapezoid, opposite_kspace_lines, SliceSelect
import ..SpoiltSliceSelects: SpoiltSliceSelect
import ..SliceSelectRephases: SliceSelectRephase
import ..EPIReadouts: EPIReadout
import ...BuildSequences: global_model, build_sequence
import ...Containers: Sequence
import ...Components: SincPulse, ConstantPulse, InstantPulse, SingleReadout
import ...Variables: qvec, flat_time, rise_time


function _get_pulse(shape, flip_angle, phase, frequency, Nzeros, group, bandwidth, duration)
    if shape == :sinc
        pulse = SincPulse(flip_angle=flip_angle, phase=phase, frequency=frequency, Nzeros=Nzeros, group=group, bandwidth=bandwidth, duration=duration)
    elseif shape in (:constant, :hard)
        pulse = ConstantPulse(flip_angle=flip_angle, phase=phase, frequency=frequency, group=group, bandwidth=bandwidth, duration=duration)
    elseif shape == :instant
        pulse = InstantPulse(flip_angle=flip_angle, phase=phase, group=group)
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
For an [`InstantPulse`](@ref) (i.e., `shape=:instant`), only the `flip_angle`, `phase`, and `group` will be used. All other parameters are ignored.
- `optimise`: set to true to optimise this RF pulse separately from the embedding sequence.

### Pulse parameters
- `shape`: The shape of the RF pulse. One of `:sinc` (for [`SincPulse`](@ref)), `:constant`/`:hard` (for [`ConstantPulse`](@ref)), or `:instant` (for [`InstantPulse`](@ref)).
- `flip_angle`: size of the flip due to the RF pulse in degrees (default: 90).
- `phase`: angle of the RF pulse in the x-y plane in degrees (default: 0).
- `frequency`: frequency of the RF pulse relative to the Larmor frequency in kHz (default: 0).
- `bandwidth`: width of the RF pulse in Fourier space in kHz (default: free variable).
- `duration`: duration of the RF pulse in ms (default: free variable).
- `Nzeros`: sets the number of zero crossings for a [`SincPulse`](@ref) (default: 3). Can be set to a tuple of two numbers to set a different number of zero crossings before and after the pulse maximum.
- `group`: name of the group of which the RF pulse is part. This is used to add transformations after the sequence is optimised.

### Slice selection
- `slice_thickness`: minimum slice thickness that should be possible without adjusting the sequence timings in um (not mm!) (default: no slice selection). Can be set to `:min` to indicate that this should be minimised given the scanner constraints and user values for `bandwidth` or `duration`.
- `rephase`: set to false to disable the spin rephasing after the RF pulse.
- `rotate_grad`: name of the parameter with which the slice selection gradient will be rotated after sequence optimisation (default: `:FOV`).
- `scanner`: overrides the [`global_scanner`](@ref) for this part of the sequence. Recommended to set only if not part of a larger sequence.
"""
function excitation_pulse(; flip_angle=90, phase=0., frequency=0., shape=:sinc, slice_thickness=Inf, rephase=true, Nzeros=3, group=nothing, rotate_grad=:FOV, bandwidth=nothing, duration=nothing, scanner=nothing, optimise=false, other_kwargs...)
    bad_keys = [key for (key, value) in pairs(other_kwargs) if !isnothing(value)]
    if length(bad_keys) > 0
        error("Unrecognised keyword arguments in call of `excitation_pulse`: $bad_keys")
    end
    build_sequence(scanner; optimise=optimise) do
        pulse = _get_pulse(shape, flip_angle, phase, frequency, Nzeros, group, bandwidth, duration)
        if pulse isa InstantPulse
            if !isinf(slice_thickness)
                error("An instant RF pulse always affects all spins equally, so using `shape=:instant` is incompatible with setting `slice_thickness`.")
            end
            return pulse
        end
        if isinf(slice_thickness)
            return pulse
        end

        cls = rephase ? SliceSelectRephase : SliceSelect
        return cls(pulse; duration=:min, slice_thickness=slice_thickness, orientation=[0, 0, 1.], group=:FOV)
    end
end



"""
    refocus_pulse(; parameters...)

Create an excitation RF pulse.

By default there is no slice-selective gradient. 
To enable slice selection `slice_thickness` has to be set to a number or to :min.
If `slice_thickness` is not set or is set to `:min`, then either `bandwidth` or `duration` should be set, otherwise the optimisation might be unconstrained (ignore this for `shape=:instant`).

## Parameters
- `optimise`: set to true to optimise this RF pulse separately from the embedding sequence.

### Pulse parameters
For an [`InstantPulse`](@ref) (i.e., `shape=:instant`), only the `flip_angle`, `phase`, and `group` will be used. All other parameters are ignored.
- `shape`: The shape of the RF pulse. One of `:sinc` (for [`SincPulse`](@ref)), `:constant`/`:hard` (for [`ConstantPulse`](@ref)), or `:instant` (for [`InstantPulse`](@ref)).
- `flip_angle`: size of the flip due to the RF pulse in degrees (default: 180).
- `phase`: angle of the RF pulse in the x-y plane in degrees (default: 0).
- `frequency`: frequency of the RF pulse relative to the Larmor frequency in kHz (default: 0).
- `bandwidth`: width of the RF pulse in Fourier space in kHz (default: free variable).
- `duration`: duration of the RF pulse in ms (default: free variable).
- `Nzeros`: sets the number of zero crossings for a [`SincPulse`](@ref) (default: 3). Can be set to a tuple of two numbers to set a different number of zero crossings before and after the pulse maximum.
- `group`: name of the group of which the RF pulse is part. This is used to add transformations after the sequence is optimised.

### Slice selection and spoilers
- `slice_thickness`: minimum slice thickness that should be possible without adjusting the sequence timings in um (not mm!) (default: no slice selection). Can be set to `:min` to indicate that this should be minimised given the scanner constraints and user values for `bandwidth` or `duration`.
- `spoiler`: set to the spatial scale on which the spins should be dephased in mm. For rotating spoilers, this does include the contribution from the slice select gradient as well.
- `rotate_grad`: name of the parameter with which the slice selection and spoiler gradient will be rotated after sequence optimisation (default: `:FOV`).
- `scanner`: overrides the [`global_scanner`](@ref) for this part of the sequence. Recommended to set only if not part of a larger sequence.
- `orientation`: vector with orientation of slice select gradient and pulses (defaults: z-direction).
"""
function refocus_pulse(; flip_angle=180, phase=0., frequency=0., shape=:sinc, slice_thickness=Inf, Nzeros=3, group=nothing, bandwidth=nothing, duration=nothing, spoiler=Inf, scanner=nothing, optimise=false, orientation=[0, 0, 1], other_kwargs...)
    bad_keys = [key for (key, value) in pairs(other_kwargs) if !isnothing(value)]
    if length(bad_keys) > 0
        error("Unrecognised keyword arguments in call of `refocus_pulse`: $bad_keys")
    end
    build_sequence(scanner; optimise=optimise) do
        pulse = _get_pulse(shape, flip_angle, phase, frequency, Nzeros, group, bandwidth, duration)
        if pulse isa InstantPulse && !isinf(slice_thickness)
            error("An instant RF pulse always affects all spins equally, so using `shape=:instant` is incompatible with setting `slice_thickness`.")
        end
        if isinf(spoiler)
            if isinf(slice_thickness)
                return pulse
            else
                return SliceSelect(pulse; duration=:min, slice_thickness=slice_thickness, orientation=orientation, group=:FOV)
            end
        else
            res = SpoiltSliceSelect(pulse; orientation=orientation, duration=:min, group=:FOV, slice_thickness=slice_thickness, spoiler_scale=spoiler)
            return res
        end
    end
end


"""
    readout_event(; type, optimise=false, variables...)

Adds a readout event to the sequence.

## Parameters
- `type`: A symbol describing the type of readout. One of the following:
    - `:epi`: EPI readout. See [`EPIReadout`](@ref) for the relevant `variables`.
    - `:instant`: single isolated readout event [`SingleReadout`](@ref) (e.g., for NMR). Does not expect any `variables`.
- `optimise`: Whether to optimise this readout event in isolation from the rest of the sequence. Use this with caution. It can speed up the optimisation (and for very complicated sequences make it more robust), however the resulting parameters might not represent the optimal solution of any external constraints (which are ignored if the readout is optimised in isolation).
- `scanner`: Used for testing. Do not set this parameter at this level (instead set it for the total sequence using [`build_sequence`](@ref)).
"""
function readout_event(; type, optimise=false, scanner=nothing, variables...)
    if type == :instant
        optimise = false # there is nothing to optimise
    end
    real_variables = Dict(key => value for (key, value) in pairs(variables) if !isnothing(variables))
    build_sequence(scanner; optimise=optimise) do 
        func_dict = Dict(
            :epi => EPIReadout,
            :instant => SingleReadout,
        )
        if !(type in keys(func_dict))
            error("Readout event type `$type` has not been implemented. Please use one of $(keys(func_dict)).")
        end
        return func_dict[type](; real_variables...)
    end
end



end