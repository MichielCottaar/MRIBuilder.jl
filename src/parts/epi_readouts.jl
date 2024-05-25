module EPIReadouts
import ...Containers: BaseSequence, get_index_single_TR
import ..Trapezoids: Trapezoid, opposite_kspace_lines, LineReadout
import ...Components: ADC
import ...Variables: get_free_variable, VariableType, set_simple_constraints!, get_readout, apply_simple_constraint!, variables, @defvar
import ...Pathways: PathwayWalker, update_walker_till_time!, walk_pathway!

"""
    EPIReadout(; resolution, ky_lines=-resolution[2]:resolution[2], recenter=false, group=:FOV, variables...)

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
struct EPIReadout{N} <: BaseSequence{N}
    start_gradient :: Trapezoid
    positive_line :: LineReadout
    negative_line :: LineReadout
    blips :: Dict{Integer, Trapezoid}
    recenter_gradient :: Union{Nothing, Trapezoid}
    ky_step :: VariableType
    ky_lines :: AbstractVector{<:Integer}
end

function EPIReadout(; resolution::AbstractVector{<:Integer}, recenter=false, group=:FOV, ky_lines=nothing, vars...)
    if length(resolution) != 2
        error("EPIReadout expects the image resolution in the x- and y-direction.")
    end
    if isnothing(ky_lines)
        ky_lines = -resolution[2]:resolution[2]
    end
    (pos, neg) = opposite_kspace_lines(; resolution=resolution[1], orientation=[1, 0, 0], group=group)
    res = EPIReadout{2 * length(ky_lines) + recenter}(
        Trapezoid(group=group),
        pos,
        neg,
        Dict{Integer, Trapezoid}(),
        recenter ? Trapezoid(group=group) : nothing,
        get_free_variable(nothing),
        ky_lines
    )
    apply_simple_constraint!(variables.qvec(res.start_gradient), VariableType[-variables.qvec(pos)[1]/2, ky_lines[1] * res.ky_step, 0.])
    if recenter
        sign = isodd(length(ky_lines)) ? -1 : 1
        apply_simple_constraint!(variables.qvec(res.recenter_gradient), VariableType[sign * variables.qvec(pos)[1]/2, -ky_lines[end] * res.ky_step, 0.])
    end
    for shift in unique(ky_lines[2:end] - ky_lines[1:end-1])
        res.blips[shift] = Trapezoid(orientation=[0, shift > 0 ? 1 : -1, 0], group=group, qval=abs(shift) * res.ky_step)
    end
    set_simple_constraints!(res, vars)
    return res
end

@defvar readout begin
    inverse_fov(epi::EPIReadout) = [variables.inverse_fov(epi.positive_line), 1e3 * epi.ky_step]
    inverse_voxel_size(epi::EPIReadout) = [variables.inverse_voxel_size(epi.positive_line), 1e3 * epi.ky_step * maximum(abs.(epi.ky_lines))]
    resolution(epi::EPIReadout) = [variables.resolution(epi.positive_line), maximum(abs.(epi.ky_lines))]
end
get_readout(epi::EPIReadout) = epi.positive_line
@defvar function effective_time(epi::EPIReadout)
    index = findfirst(iszero, epi.ky_lines)
    if isnothing(index)
        error("EPI readout does not pass through the centre of k-space")
    end
    return variables.effective_time(epi, index * 2)
end

function get_index_single_TR(epi::EPIReadout, index::Integer)
    if index == 1
        return epi.start_gradient
    elseif !isnothing(epi.recenter_gradient) && index == length(epi)
        return epi.recenter_gradient
    end
    current_line = div(index, 2, RoundDown)
    if isodd(index)
        line_shift = epi.ky_lines[current_line + 1] - epi.ky_lines[current_line]
        return epi.blips[line_shift]
    else
        if isodd(current_line)
            return epi.positive_line
        else
            return epi.negative_line
        end
    end
end

# overrides default pathway behaviour, so that EPI is counted as a single readout
function walk_pathway!(epi::EPIReadout, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType)
    if length(pulse_effects) > 0
        return false
    end
    nreadout[] -= 1
    if nreadout[] > 0
        return false
    end
    update_walker_till_time!(walker, block_start_time + variables.effective_time(epi))
    return true
end

end