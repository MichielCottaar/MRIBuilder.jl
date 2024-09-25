module MagnetisationTransfers
import ...Containers: Sequence
import ...Parts: excitation_pulse, readout_event, interpret_image_size, saturation_pulse
import ...Variables: get_pulse, get_readout, @defvar, variables
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner

const MagnetisationTransfer = Sequence{:MagnetisationTransfer}


"""
    MagnetisationTransfer(; saturation=(), excitation=(), readout=(), optim=(), resolution/fov/voxel_size, slice_thickness, scanner, variables...)

Creates a magnetisation transfer sequence of a long saturation pulse followed by an excitation pulse and then readout.

# Parameters
- `saturation`: properties of the saturation pulse as described in [`saturation_pulse`](@ref).
- `excitation`: properties of the excitation pulse as described in [`excitation_pulse`](@ref).
- `readout`: properties of the readout as described in [`readout_event`](@ref).
- Image parameters ([`variables.resolution`](@ref)/[`variables.fov`](@ref)/[`variables.voxel_size`](@ref)/[`variables.slice_thickness`](@ref)): describe the properties of the resulting image. See [`interpret_image_size`](@ref) for details.
- `scanner`: Sets the [`Scanner`](@ref) used to constraint the gradient parameters. If not set, the [`Default_Scanner`](@ref) will be used.

# Variables
- [`variables.duration`](@ref): total duration of the sequence in ms.
- [`variables.echo_time`](@ref): time between the excitation pulse and the readout.
"""
function MagnetisationTransfer(; saturation=(), excitation=(), readout=(), optim=(), resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, scanner=Default_Scanner, vars...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        parts = Any[
            :saturation => saturation_pulse(; saturation...),
            nothing,
            :excitation => excitation_pulse(; slice_thickness=slice_thickness, excitation...),
            nothing,
            :readout => readout_event(; extra_readout_params..., readout...),
            nothing,
        ]
        seq = Sequence(parts; name=:MagnetistionTransfer, vars...)
        return seq
    end
end


get_pulse(mt::MagnetisationTransfer) = (saturation=mt[:saturation].base_pulse, excitation=mt[:excitation])
get_readout(mt::MagnetisationTransfer) = mt[:readout]

@defvar begin
    echo_time(mt::MagnetisationTransfer) = effective_time(mt, :readout) - effective_time(mt, :excitation)
end


end