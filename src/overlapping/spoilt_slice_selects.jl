module SpoiltSliceSelects

import LinearAlgebra: norm
import StaticArrays: SVector
import JuMP: @variable, @constraint, @objective, objective_function
import ...BuildingBlocks: RFPulseBlock, set_simple_constraints!
import ...BuildSequences: global_model, global_scanner
import ...Variables: VariableType, variables, duration, rise_time, flat_time, effective_time, qvec, gradient_strength, slew_rate
import ...Gradients: ChangingGradientBlock, ConstantGradientBlock
import ..Abstract: interruptions, waveform, AbstractOverlapping


"""
    SpoiltSliceSelect(pulse; parameters, variables...)

Adds slice selection to the `pulse` and surrounds it with spoiler gradients.

## Parameters
- `orientation`: vector with orientation of the slice selection and the spoilers (default: [0, 0, 1])
- `rotate`: with which user-set parameter will this gradient be rotated (e.g., :FOV). Default: no rotation.

## Variables
- `duration`: total duration of the block in ms.
- `slice_thickness`: slice thickness in mm.
- `spoiler_scale`: length scale on which the spoilers achieve 2π dephasing in mm. This sets the minimum spoiling. If this spoiling level is not achieved by the slice-select gradient alone, then there will be additional gradients added.
"""
struct SpoiltSliceSelect <: AbstractOverlapping
    orientation :: SVector{3, Float64}
    rise_time1 :: VariableType
    flat_time1 :: VariableType
    diff_time :: VariableType
    pulse :: RFPulseBlock
    flat_time2 :: VariableType
    fall_time2 :: VariableType
    rotate :: Union{Nothing, Symbol}
    slew_rate :: Number
end

function SpoiltSliceSelect(pulse; orientation=[0, 0, 1], rotate=nothing, spoiler_scale=nothing, kwargs...)
    model = global_model()
    res = SpoiltSliceSelect(
        orientation ./ norm(orientation),
        @variable(model, start=0.1),
        @variable(model, start=0.1),
        @variable(model, start=0.05),
        pulse,
        @variable(model, start=0.1),
        @variable(model, start=0.1),
        rotate,
        slew_rate(global_scanner())
    )
    for time_var in (res.rise_time1, res.flat_time1, res.diff_time, res.flat_time2, res.fall_time2)
        @constraint model time_var >= 0
    end
    @constraint model res.diff_time <= res.rise_time1
    @constraint model res.diff_time <= res.fall_time2

    dim = findfirst(v -> !iszero(v), res.orientation)
    @constraint model qvec(res, nothing, 1)[dim] == qvec(res, 1, nothing)[dim]
    if !isnothing(spoiler_scale)
        if spoiler_scale == :min # note that spoiler_scale is inverse of qvec
            @objective model Min objective_function(model) - qvec(res, nothing, 1)[dim]
        elseif spoiler_scale == :max
            @objective model Min objective_function(model) + qvec(res, nothing, 1)[dim]
        else
            @constraint model (qvec(res, nothing, 1)[dim] ./ res.orientation[dim]) >= 2π * 1e-3 / spoiler_scale
        end
    end
    set_simple_constraints!(res, kwargs)

    max_time = gradient_strength(global_scanner()) / res.slew_rate
    @constraint model rise_time(res)[1] <= max_time
    @constraint model fall_time(res)[2] <= max_time
    return res
end

duration_trap1(spoilt::SpoiltSliceSelect) = 2 * spoilt.rise_time1 + spoilt.flat_time1 - spoilt.diff_time
duration_trap2(spoilt::SpoiltSliceSelect) = 2 * spoilt.fall_time2 + spoilt.flat_time2 - spoilt.diff_time

function waveform(spoilt::SpoiltSliceSelect)
    slew = slew_rate(spoilt)
    (grad1, grad2, grad3) = [g .* spoilt.orientation for g in all_gradient_strengths(spoilt)]
    return [
        ChangingGradientBlock(zeros(3), slew, rise_time(spoilt)[1], spoilt.rotate, nothing),
        ConstantGradientBlock(grad1, flat_time(spoilt)[1], spoilt.rotate, nothing),
        ChangingGradientBlock(grad1, -slew, fall_time(spoilt)[1], spoilt.rotate, nothing),
        ConstantGradientBlock(grad2, duration(spoilt.pulse), spoilt.rotate, nothing),
        ChangingGradientBlock(grad2, slew, rise_time(spoilt)[2], spoilt.rotate, nothing),
        ConstantGradientBlock(grad3, flat_time(spoilt)[2], spoilt.rotate, nothing),
        ChangingGradientBlock(grad3, -slew, fall_time(spoilt)[2], spoilt.rotate, nothing),
    ]
end
interruptions(spoilt::SpoiltSliceSelect) = [(index=4, time=effective_time(spoilt.pulse), object=spoilt.pulse)]

rise_time(spoilt::SpoiltSliceSelect) = (spoilt.rise_time1, spoilt.fall_time2 - spoilt.diff_time)
flat_time(spoilt::SpoiltSliceSelect) = (spoilt.flat_time1, spoilt.flat_time2)
fall_time(spoilt::SpoiltSliceSelect) = (spoilt.rise_time1 - spoilt.diff_time, spoilt.fall_time2)
duration(spoilt::SpoiltSliceSelect) = sum(rise_time(spoilt)) + sum(flat_time(spoilt)) + sum(flat_time(spoilt)) + duration(spoilt.pulse)
slew_rate(spoilt::SpoiltSliceSelect) = spoilt.orientation .* spoilt.slew_rate
inverse_slice_thickness(spoilt::SpoiltSliceSelect) = spoilt.slew_rate * spoilt.diff_time * duration(spoilt.pulse) * 1e3
function all_gradient_strengths(spoilt::SpoiltSliceSelect)
    grad1 = spoilt.slew_rate .* rise_time(spoilt)[1]
    grad2 = grad1 .- spoilt.slew_rate .* flat_time(spoilt)[1]
    grad3 = spoilt.slew_rate .* fall_time(spoilt)[2]
    return [grad1, grad2, grad3]
end

variables(::Type{<:SpoiltSliceSelect}) = [duration, qvec, slew_rate, inverse_slice_thickness, all_gradient_strengths, rise_time, flat_time, fall_time]


end