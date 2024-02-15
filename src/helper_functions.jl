module HelperFunctions
import JuMP: @constraint
import ..AllBuildingBlocks: BuildingBlock, Trapezoid, SpoiltSliceSelect, opposite_kspace_lines
import ..BuildSequences: global_model, build_sequence
import ..AllSequences: Sequence
import ..Components: SincPulse, ConstantPulse, InstantPulse
import ..Variables: qvec, flat_time, rise_time
import ..Alternatives: AlternativeBlocks, match_blocks!


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
- `optimise`: set to true to optimise this RF pulse separately from the embedding sequence.

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
- `scanner`: overrides the [`global_scanner`](@ref) for this part of the sequence. Recommended to set only if not part of a larger sequence.
"""
function excitation_pulse(; flip_angle=90, phase=0., frequency=0., shape=:sinc, slice_thickness=Inf, rephase=true, Nzero=3, scale=:transmit_B1, rotate_grad=:FOV, bandwidth=nothing, duration=nothing, scanner=nothing, optimise=false)
    build_sequence(scanner; optimise=optimise) do
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
For an [`InstantRFPulseBlock`](@ref) (i.e., `shape=:instant`), only the `flip_angle`, `phase`, and `scale` will be used. All other parameters are ignored.
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
- `scanner`: overrides the [`global_scanner`](@ref) for this part of the sequence. Recommended to set only if not part of a larger sequence.
"""
function refocus_pulse(; flip_angle=180, phase=0., frequency=0., shape=:sinc, slice_thickness=Inf, Nzero=3, scale=:transmit_B1, rotate_grad=:FOV, bandwidth=nothing, duration=nothing, spoiler=Inf, scanner=nothing, optimise=false)
    build_sequence(scanner; optimise=optimise) do
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

"""
    epi_readout(; resolution=, fov/voxel_size=; scanner=nothing, optimise=false, variables...)

Creates an EPI readout with given `resolution` and `fov` or `voxel_size`.

## Parameters
- `optimise`: set to true to optimise this readout separately from the embedding sequence.
- `scanner`: overrides the [`global_scanner`](@ref) for this part of the sequence. Recommended to set only if not part of a larger sequence.
- `resolution`: 2-element vector with number of voxels of the output image in the x- and y-direction. This parameter is required.
- `fov`: 2-element vector with size of the output image in the x- and y-direction in mm. Either this parameter or `voxel_size` should be set.
- `voxel_size`: 2-element vector with size of the voxels in the output image in the x- and y-direction in mm. Either this parameter or `fov` should be set.
- `oversample`: by how much to oversample in the x-direction.
"""
function epi_readout(; resolution, fov=nothing, voxel_size=nothing, kwargs...)
    start_lines = [-resolution[2]]
    readout_lines = 0:(2 * resolution[2])
    fov = _get_fov(fov, voxel_size, resolution)
    return cartesian_readout(start_lines, readout_lines, fov, resolution[1]; kwargs...)
end

"""
    single_line_readout(; resolution=, fov/voxel_size=; scanner=nothing, optimise=false, variables...)

Creates a readout that scans a single line in k-space at a time with given `resolution` and `fov` or `voxel_size`.

## Parameters
- `optimise`: set to true to optimise this readout separately from the embedding sequence.
- `scanner`: overrides the [`global_scanner`](@ref) for this part of the sequence. Recommended to set only if not part of a larger sequence.
- `resolution`: 2-element vector with number of voxels of the output image in the x- and y-direction. This parameter is required.
- `fov`: 2-element vector with size of the output image in the x- and y-direction in mm. Either this parameter or `voxel_size` should be set.
- `voxel_size`: 2-element vector with size of the voxels in the output image in the x- and y-direction in mm. Either this parameter or `fov` should be set.
- `oversample`: by how much to oversample in the x-direction.
"""
function single_line_readout(; resolution, fov=nothing, voxel_size=nothing, kwargs...)
    start_lines = -resolution[2]:resolution[2]
    readout_lines = [0]
    fov = _get_fov(fov, voxel_size, resolution)
    return cartesian_readout(start_lines, readout_lines, fov, resolution[1]; kwargs...)
end

"""
Helper function that detects consistency between `fov`, `voxel_size`, and `resolution`.
It computes the `fov` if needed.
"""
function _get_fov(fov, voxel_size, resolution)
    if isnothing(fov)
        if isnothing(voxel_size)
            error("Set either FOV or voxel size when creating an single line readout.")
        end
        return voxel_size .* resolution
    elseif !isnothing(voxel_size)
        @assert all(fov .≈ voxel_size .* resolution) "FOV, resolution, and voxel_size have been set to inconsistent values."
    end
    return fov
end

"""
    cartesian_readout(start_lines, readout_lines, fov, resolution_x; scanner=nothing, optimise=false, kwargs...)

Helper function used to build the readout for any Cartesian readout, i.e.:
- [`epi_readout`](@ref)
- [`single_line_readout`](@ref)
"""
function cartesian_readout(start_lines, readout_lines, fov, resolution_x; scanner=nothing, optimise=false, kwargs...)
    @assert iszero(readout_lines[1])
    build_sequence(scanner; optimise=optimise) do
        (pos_line, neg_line) = opposite_kspace_lines(; rotate=:FOV, fov=fov[1], resolution=resolution_x, kwargs...)
        ky = @. start_lines * 1e-3 / fov[2]
        if length(start_lines) == 1
            prepare_kspace = TrapezoidGradient(rotate=:FOV, duration=:min, qvec=[-qvec(pos_line, nothing, 1)[1], ky[1], 0.])
        else
            prepare_kspace = AlternativeBlocks(
                :readout_segment,
                [TrapezoidGradient(rotate=:FOV, duration=:min, qvec=[-qvec(pos_line, nothing, 1)[1], ky_prep, 0.]) for ky_prep in ky]
            )
            match_blocks!(prepare_kspace, flat_time)
            match_blocks!(prepare_kspace, rise_time)
        end

        steps = (readout_lines[2:end] - readout_lines[1:end-1])
        blips = Dict(
            s => TrapezoidGradient(orientation=[0, 1, 0], rotate=:FOV, duration=:min)
            for s in steps
        )
        for (s, trap) in blips
            @constraint global_model() qvec(trap)[2] == 1e-3 * s / fov[2]
        end

        result = BuildingBlock[prepare_kspace, pos_line]
        next_line = neg_line
        for s in steps
            append!(result, [blips[s], next_line])
            next_line = next_line == pos_line ? neg_line : pos_line
        end
        return Sequence(result...; TR=Inf)
    end
end



end