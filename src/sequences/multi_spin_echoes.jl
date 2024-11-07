module MultiSpinEchoes
import ...Containers: Sequence, start_time
import ...Parts: excitation_pulse, readout_event, interpret_image_size, Trapezoid, gradient_spoiler, refocus_pulse, Repeat
import ...Variables: get_pulse, get_readout, variables, @defvar, apply_simple_constraint!
import ...Pathways: Pathway, get_pathway
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner

const MultiSpinEcho = Sequence{:MultiSpinEcho}

"""
    SpinEcho(; echo_time, nrepeat, delay=0., excitation=(), refocus=(), readout=(), optim=(), resolution/fov/voxel_size/slice_thickness, scanner)

Defines a gradient echo sequence with a single readout event.

By default, an instant excitation and refocus pulse and readout event are used.
If image parameters are provided, this will switch to a sinc pulse and EPI readout.

## Parameters
- [`excitation`](@ref): properties of the excitation pulse as described in [`excitation_pulse`](@ref).
- [`refocus`](@ref): properties of the refocus pulse as described in [`refocus_pulse`](@ref).
- [`readout`](@ref): properties of the readout as described in [`readout_event`](@ref).
- Image parameters ([`variables.resolution`](@ref)/[`variables.fov`](@ref)/[`variables.voxel_size`](@ref)/[`variables.slice_thickness`](@ref)): describe the properties of the resulting image. See [`interpret_image_size`](@ref) for details.
- [`optim`](@ref): parameters to pass on to the Ipopt optimiser (see https://coin-or.github.io/Ipopt/OPTIONS.html).
- [`scanner`](@ref): Sets the [`Scanner`](@ref) used to constraint the gradient parameters. If not set, the [`Default_Scanner`](@ref) will be used.

## Variables
- [`variables.TE`](@ref)/[`variables.echo_time`](@ref): echo time between excitation pulse and spin echo in ms (required).
- [`variables.delay`](@ref): delay between the readout and the peak of the spin echo in ms (positive number indicates that readout is after the spin echo). Defaults to zero.
- [`variables.duration`](@ref): total duration of the sequence from start of excitation pulse to end of the final readout.
- [`variables.nrepeat`](@ref): How often to repeat the readout.
"""
function MultiSpinEcho(; excitation=(), refocus=(), readout=(), optim=(), spoiler=nothing, resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, delay=0., scanner=Default_Scanner, vars...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        repeating_part = Sequence(Any[
            :refocus => refocus_pulse(; slice_thickness=slice_thickness, refocus...),
            nothing,
            :readout => readout_event(; extra_readout_params..., readout...),
            nothing,
        ])
        res = Repeat(repeating_part)
        seq = Sequence(Any[
            :excitation => excitation_pulse(; slice_thickness=slice_thickness, excitation...),
            nothing,
            :repeat => res
        ]; name=:MultiSpinEcho, delay=delay, vars...)
        apply_simple_constraint!(
            (variables.effective_time(repeating_part, :refocus) + start_time(seq, :repeat) - variables.effective_time(seq, :excitation)) * 2,
            variables.echo_time(seq)
        )
        return seq
    end
end


get_pulse(se::MultiSpinEcho) = (excitation=se.excitation, refocus=se.repeat.block.refocus)
get_readout(se::MultiSpinEcho) = se.repeat.block.readout

@defvar begin
    echo_time(se::MultiSpinEcho) = variables.duration(se.repeat.block)
    delay(se::MultiSpinEcho) = variables.effective_time(se.repeat.block, :readout) - variables.echo_time(se)/2
    nrepeat(se::MultiSpinEcho) = variables.nrepeat(se.repeat)
end

end
