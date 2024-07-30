module GradientEchoes
import ...Containers: Sequence
import ...Parts: excitation_pulse, readout_event, interpret_image_size, Trapezoid, gradient_spoiler
import ...Variables: get_pulse, get_readout, variables, @defvar
import ...Pathways: Pathway, get_pathway
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner

const GradientEcho = Sequence{:GradientEcho}

"""
    GradientEcho(; echo_time, excitation=(), readout=(), optim=() resolution/fov/voxel_size/slice_thickness, scanner)

Defines a gradient echo sequence with a single readout event.

By default, an instant excitation pulse and readout event are used.
If image parameters are provided, this will switch to a sinc pulse and EPI readout.

## Parameters
- `excitation`: properties of the excitation pulse as described in [`excitation_pulse`](@ref).
- `readout`: properties of the readout as described in [`readout_event`](@ref).
- Image parameters ([`variables.resolution`](@ref)/[`variables.fov`](@ref)/[`variables.voxel_size`](@ref)/[`variables.slice_thickness`](@ref)): describe the properties of the resulting image. See [`interpret_image_size`](@ref) for details.
- `optim`: parameters to pass on to the Ipopt optimiser (see https://coin-or.github.io/Ipopt/OPTIONS.html).
- `scanner`: Sets the [`Scanner`](@ref) used to constraint the gradient parameters. If not set, the [`Default_Scanner`](@ref) will be used.

## Variables
- [`variables.TE`](@ref)/[`variables.echo_time`](@ref): echo time between excitation pulse and readout in ms (required).
- [`variables.duration`](@ref): total duration of the sequence from start of excitation pulse to end of readout in ms.
"""
function GradientEcho(; excitation=(), readout=(), optim=(), resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, scanner=Default_Scanner, variables...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        return Sequence([
            :excitation => excitation_pulse(; slice_thickness=slice_thickness, excitation...),
            nothing,
            :readout => readout_event(; extra_readout_params..., readout...),
            nothing
        ]; name=:GradientEcho, variables...)
    end
end


get_pulse(ge::GradientEcho) = ge.excitation
get_readout(ge::GradientEcho) = ge.readout
get_pathway(ge::GradientEcho) = Pathway(ge, [90], 1)
@defvar echo_time(ge::GradientEcho) = variables.duration_transverse(ge)



end