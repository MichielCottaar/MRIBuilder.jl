module EPIReadouts
import ...AllBuildingBlocks: LineReadout, Trapezoid, opposite_kspace_lines
import ...Components: ADC
import ...Variables: get_free_variable, VariableType, qval, qvec, set_simple_constraints!, resolution, inverse_voxel_size, inverse_fov, resolution, get_readout, apply_simple_constraint!
import ..BaseSequences: BaseSequence

"""
    EPIReadout(resolution; ky_lines=-resolution[2]:resolution[2], recenter=false, group=:FOV, variables...)

Defines an (accelerated) EPI readout.

## Parameters
- [`resolution`](@ref): Resolution of the final image in the frequency- and phase-encode directions.
- `recenter`: if true, the signal will be recentred in k-space after the EPI readout.
- `group`: name of the group used to rotate the readout gradients (default: :FOV).

## Variables:
- [`voxel_size`](@ref): size of the voxel in the frequency- and phase-encode directions.
- [`fov`](@ref): size of the FOV in the frequency- and phase-encode directions.
- [`ramp_overlap`](@ref): what fraction of the gradient ramp should overlap with the readout.
- [`oversample`](@ref): by how much to oversample in the frequency-encode direcion.
- [`dwell_time`](@ref): dwell time in the frequency-encode direction
"""
struct EPIReadout <: BaseSequence{0}
    start_gradient :: Trapezoid
    positive_line :: LineReadout
    negative_line :: LineReadout
    blips :: Dict{Integer, Trapezoid}
    recenter_gradient :: Union{Nothing, Trapezoid}
    ky_step :: VariableType
    ky_lines :: AbstractVector{<:Integer}
end

function EPIReadout(resolution::AbstractVector{<:Integer}; recenter=false, group=:FOV, ky_lines=nothing, variables...)
    if length(resolution) != 2
        error("EPIReadout expects the image resolution in the x- and y-direction.")
    end
    if isnothing(ky_lines)
        ky_lines = -resolution[2]:resolution[2]
    end
    (pos, neg) = opposite_kspace_lines(; resolution=resolution[1], orientation=[1, 0, 0], group=group)
    res = EPIReadout(
        Trapezoid(group=group),
        pos,
        neg,
        Dict{Integer, Trapezoid}(),
        recenter ? Trapezoid(group=group) : nothing,
        get_free_variable(nothing),
        ky_lines
    )
    apply_simple_constraint!(qvec(res.start_gradient), VariableType[-qval(pos)/2, ky_lines[1] * res.ky_step, 0.])
    if recenter
        sign = isodd(length(ky_lines)) ? -1 : 1
        apply_simple_constraint!(qvec(res.recenter_gradient), VariableType[sign * qval(pos)/2, -ky_lines[end] * res.ky_step, 0.])
    end
    set_simple_constraints!(res, variables)
    return res
end

inverse_fov(epi::EPIReadout) = [inverse_fov(epi.positive_line), 1e3 * epi.ky_step]
inverse_voxel_size(epi::EPIReadout) = [inverse_voxel_size(epi.positive_line), 1e3 * epi.ky_step * maximum(abs.(epi.ky_lines))]
resolution(epi::EPIReadout) = [resolution(epi.positive_line), maximum(abs.(epi.ky_lines))]
get_readout(epi::EPIReadout) = epi.positive_line


end