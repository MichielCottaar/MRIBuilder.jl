module MagnetisationTransfers
import ...Containers: Sequence
import ...Parts: excitation_pulse, readout_event, interpret_image_size, saturation_pulse_binomial
import ...Variables: get_pulse, get_readout
import ...BuildSequences: build_sequence
import ...Scanners: Default_Scanner

const MagnetisationTransfer = Sequence{:MagnetisationTransfer}


"""
    MagnetisationTransfer(; saturation=(), excitation=(), readout=(), optim=())
"""
function MagnetisationTransfer(; saturation=(), excitation=(), readout=(), optim=(), resolution=nothing, fov=nothing, voxel_size=nothing, slice_thickness=nothing, scanner=Default_Scanner, vars...)
    build_sequence(scanner; optim...) do
        (slice_thickness, _, extra_readout_params) = interpret_image_size(fov, resolution, voxel_size, slice_thickness)
        parts = Any[
            :saturation => saturation_pulse_binomial(; saturation...),
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


get_pulse(mt::MagnetisationTransfer) = (saturation=mt[:excitation].base_pulse, refocus=ge[:refocus])
get_readout(mt::MagnetisationTransfer) = ge[:readout]


end