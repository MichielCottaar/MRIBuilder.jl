module SpoiltSliceSelects

import LinearAlgebra: norm
import StaticArrays: SVector
import JuMP: @constraint, @objective, objective_function
import ...BuildSequences: global_model, global_scanner
import ...Variables: VariableType, get_pulse, set_simple_constraints!, variables, @defvar, gradient_orientation
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

function SpoiltSliceSelect(pulse::RFPulseComponent; orientation=[0, 0, 1], group=nothing, slice_thickness=nothing, kwargs...)
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
        @constraint model qvec(res, nothing, :pulse) == qvec(res, :pulse, nothing)
    end

    set_simple_constraints!(res, kwargs)

    max_time = gradient_strength(global_scanner()) / res.slew_rate
    @constraint model rise_time(res)[1] <= max_time
    @constraint model fall_time(res)[2] <= max_time
    return res
end

gradient_orientation(spoilt::SpoiltSliceSelect) = spoilt.orientation

@defvar begin
    duration_trap1(spoilt::SpoiltSliceSelect) = 2 * spoilt.rise_time1 + spoilt.flat_time1 - spoilt.diff_time
    duration_trap2(spoilt::SpoiltSliceSelect) = 2 * spoilt.fall_time2 + spoilt.flat_time2 - spoilt.diff_time
end

Base.keys(::SpoiltSliceSelect) = Val.((:rise1, :flat1, :fall1, :flat_pulse, :pulse, :rise2, :flat2, :fall2))
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:rise1}) = ChangingGradient(0., variables.slew_rate(spoilt), gradient_orientation(spoilt), variables.rise_time(spoilt)[1], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:flat1}) = ConstantGradient(variables.slew_rate(spoilt) * variables.rise_time(spoilt)[1], gradient_orientation(spoilt), variables.flat_time(spoilt)[1], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:fall1}) = ChangingGradient(variables.slew_rate(spoilt) * variables.rise_time(spoilt)[1], -variables.slew_rate(spoilt), gradient_orientation(spoilt), variables.fall_time(spoilt)[1], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:flat_pulse}) = ConstantGradient(variables.slew_rate(spoilt) * spoilt.diff_time, gradient_orientation(spoilt), variables.duration(spoilt.pulse), spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:pulse}) = spoilt.pulse
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:rise2}) = ChangingGradient(variables.slew_rate(spoilt) * spoilt.diff_time, variables.slew_rate(spoilt), gradient_orientation(spoilt), variables.rise_time(spoilt)[2], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:flat2}) = ConstantGradient(variables.slew_rate(spoilt) * variables.fall_time(spoilt)[2], gradient_orientation(spoilt), variables.flat_time(spoilt)[2], spoilt.group)
Base.getindex(spoilt::SpoiltSliceSelect, ::Val{:fall2}) = ChangingGradient(variables.slew_rate(spoilt) * variables.fall_time(spoilt)[2], -variables.slew_rate(spoilt), gradient_orientation(spoilt), variables.fall_time(spoilt)[2], spoilt.group)

@defvar begin
    rise_time(spoilt::SpoiltSliceSelect) = (spoilt.rise_time1, spoilt.fall_time2 - spoilt.diff_time)
    flat_time(spoilt::SpoiltSliceSelect) = (spoilt.flat_time1, spoilt.flat_time2)
    fall_time(spoilt::SpoiltSliceSelect) = (spoilt.rise_time1 - spoilt.diff_time, spoilt.fall_time2)
    slew_rate(spoilt::SpoiltSliceSelect) = spoilt.slew_rate
end

@defvar begin
    duration(spoilt::SpoiltSliceSelect) = sum(variables.rise_time(spoilt)) + sum(variables.flat_time(spoilt)) + sum(variables.flat_time(spoilt)) + variables.duration(spoilt.pulse)
    inverse_slice_thickness(spoilt::SpoiltSliceSelect) = variables.slew_rate(spoilt) * spoilt.diff_time * variables.duration(spoilt.pulse) * 1e3
    gradient_strength(spoilt::SpoiltSliceSelect) = variables.slew_rate(spoilt) * max(spoilt.rise_time1, spoilt.fall_time2)
end

get_pulse(spoilt::SpoiltSliceSelect) = spoilt.pulse

@defvar function all_gradient_strengths(spoilt::SpoiltSliceSelect)
    grad1 = spoilt.slew_rate * variables.rise_time(spoilt)[1]
    grad2 = grad1 - spoilt.slew_rate * variables.flat_time(spoilt)[1]
    grad3 = spoilt.slew_rate * variables.fall_time(spoilt)[2]
    return [grad1, grad2, grad3]
end

"""
    all_gradient_strengths(spoilt_slice_select)

Returns the gradient strength before, during, and after the pulse in [`SpoiltSliceSelect`](@ref).
"""
variables.all_gradient_strengths


"""
    fall_time(spoilt_slice_select)

Returns the time of the [`SpoiltSliceSelect`](@ref) to return to zero.
"""
variables.all_gradient_strengths

end