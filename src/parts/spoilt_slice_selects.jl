module SpoiltSliceSelects

import LinearAlgebra: norm
import StaticArrays: SVector
import JuMP: @constraint, @objective, objective_function
import ...BuildSequences: global_model, global_scanner
import ...Variables: VariableType, duration, rise_time, flat_time, effective_time, qval, gradient_strength, slew_rate, inverse_slice_thickness, get_free_variable, get_pulse, set_simple_constraints!, gradient_orientation
import ...Components: ChangingGradient, ConstantGradient, RFPulseComponent
import ...Containers: BaseBuildingBlock


"""
    SpoiltSliceSelect(pulse; parameters, variables...)

Adds slice selection to the `pulse` and surrounds it with spoiler gradients.

## Parameters
- `orientation`: vector with orientation of the slice selection and the spoilers (default: [0, 0, 1])
- `group`: name of the group of the gradient. This will be used to scale and rotate the gradients after optimisation. Scaling is not recommended as this might ruin the spoiling.

## Variables
- `duration`: total duration of the block in ms.
- `slice_thickness`: slice thickness in mm.
- `spoiler_scale`: length scale on which the spoilers achieve 2π dephasing in mm. This sets the minimum spoiling. If this spoiling level is not achieved by the slice-select gradient alone, then there will be additional gradients added.
"""
struct SpoiltSliceSelect <: BaseBuildingBlock
    orientation :: SVector{3, Float64}
    rise_time1 :: VariableType
    flat_time1 :: VariableType
    diff_time :: VariableType
    pulse :: RFPulseComponent
    flat_time2 :: VariableType
    fall_time2 :: VariableType
    group :: Union{Nothing, Symbol}
    slew_rate :: Number
end

function SpoiltSliceSelect(pulse::RFPulseComponent; orientation=[0, 0, 1], group=nothing, spoiler_scale=nothing, slice_thickness=nothing, kwargs...)
    model = global_model()
    res = nothing
    if slice_thickness isa Number && isinf(slice_thickness)
        rise_time_var = get_free_variable(nothing)
        flat_time_var = get_free_variable(nothing)
        res = SpoiltSliceSelect(
            orientation ./ norm(orientation),
            rise_time_var,
            flat_time_var,
            0.,
            pulse,
            flat_time_var,
            rise_time_var,
            group,
            slew_rate(global_scanner()),
        )
        for time_var in (rise_time_var, flat_time_var)
            @constraint model time_var >= 0
        end
    else
        res = SpoiltSliceSelect(
            orientation ./ norm(orientation),
            get_free_variable(nothing; start=0.1),
            get_free_variable(nothing; start=0.1),
            get_free_variable(nothing; start=0.05),
            pulse,
            get_free_variable(nothing; start=0.1),
            get_free_variable(nothing; start=0.1),
            group,
            slew_rate(global_scanner()),
        )
        for time_var in (res.rise_time1, res.flat_time1, res.diff_time, res.flat_time2, res.fall_time2)
            @constraint model time_var >= 0
        end
        @constraint model res.diff_time <= res.rise_time1
        @constraint model res.diff_time <= res.fall_time2
        @constraint model qval(res, nothing, :pulse) == qval(res, :pulse, nothing)
    end

    if !isnothing(spoiler_scale)
        if spoiler_scale == :min # note that spoiler_scale is inverse of qval
            @objective model Min objective_function(model) - qval(res, nothing, 1)
        elseif spoiler_scale == :max
            @objective model Min objective_function(model) + qval(res, nothing, 1)
        else
            @constraint model qval(res, nothing, :pulse) >= 2π * 1e-3 / spoiler_scale
        end
    end
    set_simple_constraints!(res, kwargs)

    max_time = gradient_strength(global_scanner()) / res.slew_rate
    @constraint model rise_time(res)[1] <= max_time
    @constraint model fall_time(res)[2] <= max_time
    return res
end

gradient_orientation(spoilt::SpoiltSliceSelect) = spoilt.orientation
duration_trap1(spoilt::SpoiltSliceSelect) = 2 * spoilt.rise_time1 + spoilt.flat_time1 - spoilt.diff_time
duration_trap2(spoilt::SpoiltSliceSelect) = 2 * spoilt.fall_time2 + spoilt.flat_time2 - spoilt.diff_time

Base.keys(::SpoiltSliceSelect) = Val.((:rise1, :flat1, :fall1, :flat_pulse, :pulse, :rise2, :flat2, :fall2))
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:rise1}) = ChangingGradient(0., slew_rate(spoilt), gradient_orientation(spoilt), rise_time(spoilt)[1], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:flat1}) = ConstantGradient(slew_rate(spoilt) * rise_time(spoilt)[1], gradient_orientation(spoilt), flat_time(spoilt)[1], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:fall1}) = ChangingGradient(slew_rate(spoilt) * rise_time(spoilt)[1], -slew_rate(spoilt), gradient_orientation(spoilt), fall_time(spoilt)[1], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:flat_pulse}) = ConstantGradient(slew_rate(spoilt) * spoilt.diff_time, gradient_orientation(spoilt), duration(spoilt.pulse), spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:pulse}) = spoilt.pulse
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:rise2}) = ChangingGradient(slew_rate(spoilt) * spoilt.diff_time, slew_rate(spoilt), gradient_orientation(spoilt), rise_time(spoilt)[2], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:flat2}) = ConstantGradient(slew_rate(spoilt) * fall_time(spoilt)[2], gradient_orientation(spoilt), flat_time(spoilt)[2], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:fall2}) = ChangingGradient(slew_rate(spoilt) * fall_time(spoilt)[2], -slew_rate(spoilt), gradient_orientation(spoilt), fall_time(spoilt)[2], spoilt.group)

rise_time(spoilt::SpoiltSliceSelect) = (spoilt.rise_time1, spoilt.fall_time2 - spoilt.diff_time)
flat_time(spoilt::SpoiltSliceSelect) = (spoilt.flat_time1, spoilt.flat_time2)
fall_time(spoilt::SpoiltSliceSelect) = (spoilt.rise_time1 - spoilt.diff_time, spoilt.fall_time2)
duration(spoilt::SpoiltSliceSelect) = sum(rise_time(spoilt)) + sum(flat_time(spoilt)) + sum(flat_time(spoilt)) + duration(spoilt.pulse)
slew_rate(spoilt::SpoiltSliceSelect) = spoilt.slew_rate
inverse_slice_thickness(spoilt::SpoiltSliceSelect) = spoilt.slew_rate * spoilt.diff_time * duration(spoilt.pulse) * 1e3
gradient_strength(spoilt::SpoiltSliceSelect) = slew_rate(spoilt) * max(spoilt.rise_time1, spoilt.fall_time2)
get_pulse(spoilt::SpoiltSliceSelect) = spoilt.pulse

function all_gradient_strengths(spoilt::SpoiltSliceSelect)
    grad1 = spoilt.slew_rate * rise_time(spoilt)[1]
    grad2 = grad1 - spoilt.slew_rate * flat_time(spoilt)[1]
    grad3 = spoilt.slew_rate * fall_time(spoilt)[2]
    return [grad1, grad2, grad3]
end


end