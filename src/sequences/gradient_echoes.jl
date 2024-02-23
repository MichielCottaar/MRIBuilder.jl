module GradientEchoes
import ...Containers: Sequence
import ...Parts: excitation_pulse, readout_event, interpret_image_size, Trapezoid, gradient_spoiler
import ...Variables: get_pulse, get_readout, echo_time, duration_transverse
import ...Pathways: Pathway, get_pathway
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner

const GradientEcho = Sequence{:GradientEcho}

"""
    GradientEcho(; echo_time, excitation=(), readout=(), resolution/fov/voxel_size/slice_thickness)

Defines a gradient echo sequence with a single readout event.

By default, an instant excitation pulse and readout event are used.
If image parameters are provided, this will switch to a sinc pulse and EPI readout.

## Parameters
- [`excitation`](@ref): properties of the excitation pulse as described in [`excitation_pulse`](@ref).
- [`readout`](@ref): properties of the readout as described in [`readout_event`](@ref).
- [`spoiler`](@ref): if set adds a spoiler [`gradient_spoiler`](@ref) gradient after the readout (e.g., `spoiler=()` to add a gradient in the z-direction of the `FOV` coordinate system that fully dephases spins over 1 mm).
- Image parameters ([`resolution`](@ref)/[`fov`](@ref)/[`voxel_size`](@ref)/[`slice_thickness`](@ref)): describe the properties of the resulting image. See [`interpret_image_size`](@ref) for details.
- [`optim`](@ref): parameters to pass on to the Ipopt optimiser (see https://coin-or.github.io/Ipopt/OPTIONS.html).

## Variables
- [`TE`](@ref)/[`echo_time`](@ref): echo time between excitation pulse and readout in ms (required).
- [`TR`](@ref)/[`repetition_time`](@ref)/[`duration`](@ref): total duration of the sequence from start of excitation pulse to end of readout or spoiler in ms.
"""
function GradientEcho(; excitation=(), readout=(), optim=(), spoiler=nothing, resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, scanner=Default_Scanner, variables...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        parts = Any[
            :excitation => excitation_pulse(; slice_thickness=slice_thickness, excitation...),
            nothing,
            :readout => readout_event(; extra_readout_params..., readout...),
        ]
        if !isnothing(spoiler)
            push!(parts, gradient_spoiler(; spoiler...))
        end
        return Sequence(parts; name=:GradientEcho, variables...)
    end
end


get_pulse(ge::GradientEcho) = ge.excitation
get_readout(ge::GradientEcho) = ge.readout
get_pathway(ge::GradientEcho) = Pathway(ge, [90], 1)
echo_time(ge::GradientEcho) = duration_transverse(ge)



end