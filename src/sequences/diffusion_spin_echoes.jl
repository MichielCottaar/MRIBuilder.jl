module DiffusionSpinEchoes
import ...Containers: Sequence
import ...Components: InstantGradient
import ...Parts: excitation_pulse, readout_event, interpret_image_size, Trapezoid, gradient_spoiler, refocus_pulse, dwi_gradients
import ...Containers: start_time, end_time
import ...Variables: get_pulse, get_readout, get_gradient, variables, @defvar, add_cost_function!
import ...Pathways: Pathway, get_pathway
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner, Scanner

const DiffusionSpinEcho = Sequence{:DiffusionSpinEcho}
const DW_SE = DiffusionSpinEcho
const DWI = DiffusionSpinEcho

"""
    DiffusionSpinEcho(; echo_time, delay=0., excitation=(), gradient=(), refocus=(), readout=(), optim=(), resolution/fov/voxel_size/slice_thickness, scanner)

Defines a diffusion-weighted spin echo (Stejskal-Tanner) sequence.

`DWI`, `DW_SE`, and `DiffusionSpinEcho` are all synonyms.

By default, an instant excitation pulse and readout event are used.
If image parameters are provided, this will switch to a sinc pulse and EPI readout.

## Parameters
- `excitation`: properties of the excitation pulse as described in [`excitation_pulse`](@ref).
- `gradient`: properties of the diffusion-weighting gradients as described in [`dwi_gradients`](@ref).
- `refocus`: properties of the refocus pulse as described in [`refocus_pulse`](@ref).
- `readout`: properties of the readout as described in [`readout_event`](@ref).
- Image parameters ([`variables.resolution`](@ref)/[`variables.fov`](@ref)/[`variables.voxel_size`](@ref)/[`variables.slice_thickness`](@ref)): describe the properties of the resulting image. See [`interpret_image_size`](@ref) for details.
- `optim`: parameters to pass on to the Ipopt optimiser (see https://coin-or.github.io/Ipopt/OPTIONS.html).
- `scanner`: Sets the [`Scanner`](@ref) used to constraint the gradient parameters. If not set, the [`Default_Scanner`](@ref) will be used.

## Variables
- [`variables.TE`](@ref)/[`variables.echo_time`](@ref): echo time between excitation pulse and spin echo in ms.
- [`variables.delay`](@ref): delay between the readout and the peak of the spin echo in ms (positive number indicates that readout is after the spin echo). Defaults to zero.
- [`variables.duration`](@ref): total duration of the sequence from start of excitation pulse to end of readout or spoiler in ms.
- [`variables.Δ`](@ref)/[`variables.diffusion_time`](@ref): Time from the start of one diffusion-weighted gradient till the other in ms.
"""
function DiffusionSpinEcho(; delay=0., excitation=(), gradient=(), refocus=(), readout=(), optim=(), spoiler=nothing, resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, scanner=Default_Scanner, vars...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        (g1, g2) = dwi_gradients(; gradient...)
        parts = Any[
            :excitation => excitation_pulse(; slice_thickness=slice_thickness, excitation...),
            nothing,
            :gradient => g1,
            nothing,
            :refocus => refocus_pulse(; slice_thickness=slice_thickness, refocus...),
            nothing,
            :gradient2 => g2,
            nothing,
            :readout => readout_event(; extra_readout_params..., readout...),
            nothing,
        ]
        if !isnothing(spoiler)
            push!(parts, gradient_spoiler(; spoiler...))
        end
        seq = Sequence(parts; name=:DiffusionSpinEcho, delay=delay, vars...)
        
        if g1 isa InstantGradient
            add_cost_function!((variables.duration(seq[2]) - variables.duration(seq[4]))^2)
            add_cost_function!((variables.duration(seq[6]) - variables.duration(seq[8]))^2)
        else
            # Minimise sequence duration with less emphasis on waiting block between second gradient and readout
            add_cost_function!(variables.duration(seq) - 0.5 * variables.duration(seq[8]), 2)
        end
        return seq
    end
end


get_pulse(ge::DiffusionSpinEcho) = (excitation=ge[:excitation], refocus=ge[:refocus])
get_gradient(ge::DiffusionSpinEcho) = (gradient=ge[:gradient], gradient2=ge[:gradient2])
get_readout(ge::DiffusionSpinEcho) = ge[:readout]
get_pathway(ge::DiffusionSpinEcho) = Pathway(ge, [90, 180], 1, group=:diffusion)

@defvar begin
    echo_time(ge::DiffusionSpinEcho) = 2 * (variables.effective_time(ge, :refocus) - variables.effective_time(ge, :excitation))
    delay(ge::DiffusionSpinEcho) = variables.duration_transverse(ge) - variables.echo_time(ge)
    diffusion_time(ge::DiffusionSpinEcho) = start_time(ge, :gradient2) - start_time(ge, :gradient)
end

"""
    diffusion_time(diffusion_sequence)

Returns the diffusion time of a [`DiffusionSpinEcho`](@ref) in ms.
"""
variables.diffusion_time

"""
    delay(sequence)

Returns the offset beetween the readout and the spin echo in ms.
"""
variables.delay

end