module HelperFunctions
import JuMP: @constraint
import ...Containers: AlternativeBlocks, match_blocks!, BuildingBlock
import ..Trapezoids: Trapezoid, opposite_kspace_lines, SliceSelect
import ..SpoiltSliceSelects: SpoiltSliceSelect
import ..SliceSelectRephases: SliceSelectRephase
import ..EPIReadouts: EPIReadout
import ...BuildSequences: global_model, build_sequence
import ...Containers: Sequence
import ...Components: SincPulse, ConstantPulse, InstantPulse, SingleReadout, InstantGradient
import ...Variables: qvec, flat_time, rise_time, qval, apply_simple_constraint!, variables


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
- `shape`: The shape of the RF pulse. One of `:sinc` (for [`SincPulse`](@ref)), `:constant`/`:hard` (for [`ConstantPulse`](@ref)), or `:instant` (for [`InstantPulse`](@ref)). Default is :sinc for slice-selective pulses or :instant otherwise.
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
function excitation_pulse(; flip_angle=90, phase=0., frequency=0., shape=nothing, slice_thickness=Inf, rephase=true, Nzeros=3, group=nothing, rotate_grad=:FOV, bandwidth=nothing, duration=nothing, scanner=nothing, optimise=false, other_kwargs...)
    bad_keys = [key for (key, value) in pairs(other_kwargs) if !isnothing(value)]
    if length(bad_keys) > 0
        error("Unrecognised keyword arguments in call of `excitation_pulse`: $bad_keys")
    end
    if isnothing(slice_thickness) || slice_thickness isa Symbol
        error("slice thickness cannot be a free parameter")
    end
    if isnothing(shape)
        shape = isinf(slice_thickness) ? :instant : :sinc
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
- `shape`: The shape of the RF pulse. One of `:sinc` (for [`SincPulse`](@ref)), `:constant`/`:hard` (for [`ConstantPulse`](@ref)), or `:instant` (for [`InstantPulse`](@ref)). Default is :sinc for slice-selective pulses or :instant otherwise.
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
function refocus_pulse(; flip_angle=180, phase=0., frequency=0., shape=nothing, slice_thickness=Inf, Nzeros=3, group=nothing, bandwidth=nothing, duration=nothing, spoiler=Inf, scanner=nothing, optimise=false, orientation=[0, 0, 1], other_kwargs...)
    bad_keys = [key for (key, value) in pairs(other_kwargs) if !isnothing(value)]
    if length(bad_keys) > 0
        error("Unrecognised keyword arguments in call of `refocus_pulse`: $bad_keys")
    end
    if isnothing(slice_thickness) || slice_thickness isa Symbol
        error("slice thickness cannot be a free parameter")
    end
    if isnothing(shape)
        shape = isinf(slice_thickness) ? :instant : :sinc
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
- `type`: A symbol describing the type of readout. It will default to `:epi` if a resolution has been set and `:instant` otherwise. Can be set to one of the following:
    - `:epi`: EPI readout. See [`EPIReadout`](@ref) for the relevant `variables`.
    - `:instant`: single isolated readout event [`SingleReadout`](@ref) (e.g., for NMR). Does not expect any `variables`.
- `optimise`: Whether to optimise this readout event in isolation from the rest of the sequence. Use this with caution. It can speed up the optimisation (and for very complicated sequences make it more robust), however the resulting parameters might not represent the optimal solution of any external constraints (which are ignored if the readout is optimised in isolation).
- `scanner`: Used for testing. Do not set this parameter at this level (instead set it for the total sequence using [`build_sequence`](@ref)).
"""
function readout_event(; type=nothing, optimise=false, scanner=nothing, oversample=2, duration=:min, all_variables...)
    real_variables = Dict(key => value for (key, value) in pairs(all_variables) if !(isnothing(value) || (value isa AbstractVector && all(isnothing.(value)))))
    if isnothing(type)
        resolution = get(real_variables, :resolution, nothing)
        type = (isnothing(resolution) || (resolution isa Union{Tuple, AbstractVector} && all(isnothing.(resolution)))) ? :instant : :epi
    end
    if type == :instant
        optimise = false # there is nothing to optimise
        oversample = nothing
    end
    build_sequence(scanner; optimise=optimise) do 
        func_dict = Dict(
            :epi => EPIReadout,
            :instant => SingleReadout,
        )
        if !(type in keys(func_dict))
            error("Readout event type `$type` has not been implemented. Please use one of $(keys(func_dict)).")
        end
        return func_dict[type](; duration=duration, oversample=oversample, real_variables...)
    end
end


"""
    dwi_gradients(; type, optimise=false, refocus=true, orientation=[1, 0, 0], group=:DWI, variables...)

Returns two DWI gradients that are guaranteed to cancel each other out.

## Parameters
- `type`: A symbol describing the type of gradient. One of:
    - `:trapezoid`: Pulsed trapezoidal gradient profile. See [`Trapezoid`](@ref) for the relevant `variables`.
    - `:instant`: instantaneous gradient (e.g., to test short-pulse approximations). See [`InstantGradient`](@ref) for the relevant `variables`.
- `optimise`: Whether to optimise this readout event in isolation from the rest of the sequence. Use this with caution. It can speed up the optimisation (and for very complicated sequences make it more robust), however the resulting parameters might not represent the optimal solution of any external constraints (which are ignored if the readout is optimised in isolation).
- `scanner`: Used for testing. Do not set this parameter at this level (instead set it for the total sequence using [`build_sequence`](@ref)).
"""
function dwi_gradients(; type=:trapezoid, optimise=false, scanner=nothing, refocus=true, orientation=[1, 0, 0], group=:DWI, match=nothing, all_variables...)
    real_variables = Dict(key => value for (key, value) in pairs(all_variables) if !(isnothing(value) || (value isa AbstractVector && all(isnothing.(value)))))

    func_dict = Dict(
        :trapezoid => Trapezoid,
        :instant => InstantGradient,
    )
    if !(type in keys(func_dict))
        error("DWI gradients type `$type` has not been implemented. Please use one of $(keys(func_dict)).")
    end

    if isnothing(match)
        match = get(Dict(
            :trapezoid => [:rise_time, :flat_time, :slew_rate],
        ), type, [])
    end

    get_index(var::Union{Tuple, AbstractVector}, i) = length(var) == 2 ? var[i] : var
    get_index(var::NamedTuple, i) = var
    get_index(var, i) = var

    build_sequence(scanner; optimise=optimise) do 
        other_orientation = isnothing(orientation) ? nothing : (refocus ? orientation : -orientation)
        (g1, g2) = [func_dict[type](;
            group=group, orientation=o, Dict(key => get_index(value, i) for (key, value) in pairs(real_variables))...
            ) for (i, o) in enumerate((orientation, other_orientation))]
        if !isnothing(orientation) || refocus
            apply_simple_constraint!(qval(g1), qval(g2))
        else
            apply_simple_constraint!(qval(g1), -qval(g2))
        end
        for var_func in match
            if var_func isa Symbol
                var_func = variables[var_func]
            end
            apply_simple_constraint!(var_func(g1), var_func(g2))
        end
        return (g1, g2)
    end
end


"""
    interpret_image_size(fov, resolution, voxel_size, slice_thickness)

Combines the user-provided information of the image size to produce a tuple with:
- `slice_thickness`: if not set explicitly, will be set to the third element of `voxel_size`
- `resolution_z`: number of voxels in the slice-select direction.
- keywords parameters for the [`readout_event`](@ref) functions with the `fov`, `resolution`, and `voxel_size` in the x- and y- dimensions.
"""
function interpret_image_size(fov, resolution, voxel_size, slice_thickness)
    _real_value(n::Number) = true
    _real_value(n::NTuple{N, <:Number}) where {N} = true
    _real_value(n::AbstractVector{<:Number}) = true
    _real_value(other) = false

    if all(isnothing.((fov, resolution, voxel_size)))
        return (isnothing(slice_thickness) ? Inf : slice_thickness, nothing, (fov=nothing, resolution=nothing, voxel_size=nothing))
    end
    if !_real_value(resolution)
        if !(_real_value(fov) && _real_value(voxel_size))
            error("`resolution` needs to be fully defined; either by setting it directly or by setting `fov` and `voxel_size` to concrete values.")
        end
        resolution = Int.(div.(fov, voxel_size, RoundUp))
        fov = resolution .* voxel_size
    end

    getz(v::Union{Tuple, AbstractVector}) = length(v) == 2 ? nothing : v[3]
    getz(v) = v
    if isnothing(slice_thickness)
        slice_thickness = getz(voxel_size)
        if isnothing(slice_thickness)
            fovz = getz(fov)
            if isnothing(fovz)
                slice_thickness = Inf
            else
                slice_thickness = fovz / getz(resolution)
            end
        end
    end

    getxy(v::Union{Tuple, AbstractVector}) = [v[1], v[2]]
    getxy(v) = [v, v]

    fov=getxy(fov)
    resolution=getxy(resolution)
    voxel_size=getxy(voxel_size)
    if all(isnothing.(voxel_size))
        return (slice_thickness, getz(resolution), (resolution=resolution, fov=fov))
    else
        return (slice_thickness, getz(resolution), (resolution=resolution, voxel_size=voxel_size))
    end
end

end