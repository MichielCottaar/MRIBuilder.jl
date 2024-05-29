module SpinEchoes
import ...Containers: Sequence
import ...Parts: excitation_pulse, readout_event, interpret_image_size, Trapezoid, gradient_spoiler, refocus_pulse
import ...Variables: get_pulse, get_readout, variables, @defvar
import ...Pathways: Pathway, get_pathway
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner

const SpinEcho = Sequence{:SpinEcho}

"""
    SpinEcho(; echo_time, delay=0., excitation=(), refocus=(), readout=(), optim=(), resolution/fov/voxel_size/slice_thickness, scanner)

Defines a gradient echo sequence with a single readout event.

By default, an instant excitation pulse and readout event are used.
If image parameters are provided, this will switch to a sinc pulse and EPI readout.

## Parameters
- [`excitation`](@ref): properties of the excitation pulse as described in [`excitation_pulse`](@ref).
- [`refocus`](@ref): properties of the refocus pulse as described in [`refocus_pulse`](@ref).
- [`readout`](@ref): properties of the readout as described in [`readout_event`](@ref).
- [`spoiler`](@ref): if set adds a spoiler [`gradient_spoiler`](@ref) gradient after the readout (e.g., `spoiler=()` to add a gradient in the z-direction of the `FOV` coordinate system that fully dephases spins over 1 mm).
- Image parameters ([`resolution`](@ref)/[`fov`](@ref)/[`voxel_size`](@ref)/[`slice_thickness`](@ref)): describe the properties of the resulting image. See [`interpret_image_size`](@ref) for details.
- [`optim`](@ref): parameters to pass on to the Ipopt optimiser (see https://coin-or.github.io/Ipopt/OPTIONS.html).
- [`scanner`](@ref): Sets the [`Scanner`](@ref) used to constraint the gradient parameters. If not set, the [`Default_Scanner`](@ref) will be used.

## Variables
- [`TE`](@ref)/[`echo_time`](@ref): echo time between excitation pulse and spin echo in ms (required).
- [`delay`](@ref): delay between the readout and the peak of the spin echo in ms (positive number indicates that readout is after the spin echo). Defaults to zero.
- [`duration`](@ref): total duration of the sequence from start of excitation pulse to end of readout or spoiler in ms.
"""
function SpinEcho(; delay=0., excitation=(), refocus=(), readout=(), optim=(), spoiler=nothing, resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, scanner=Default_Scanner, variables...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        parts = Any[
            :excitation => excitation_pulse(; slice_thickness=slice_thickness, excitation...),
            nothing,
            :refocus => refocus_pulse(; slice_thickness=slice_thickness, refocus...),
            nothing,
            :readout => readout_event(; extra_readout_params..., readout...),
        ]
        if !isnothing(spoiler)
            push!(parts, gradient_spoiler(; spoiler...))
        end
        return Sequence(parts; name=:SpinEcho, delay=delay, variables...)
    end
end


get_pulse(ge::SpinEcho) = (excitation=ge.excitation, refocus=ge.refocus)
get_readout(ge::SpinEcho) = ge.readout
get_pathway(ge::SpinEcho) = Pathway(ge, [90, 180], 1)

@defvar begin
    echo_time(ge::SpinEcho) = 2 * (variables.effective_time(ge, :refocus) - variables.effective_time(ge, :excitation))
    delay(ge::SpinEcho) = variables.duration_transverse(ge) - variables.echo_time(ge)
end


end