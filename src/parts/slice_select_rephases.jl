module SliceSelectRephases
import ...Containers: BaseSequence, get_index_single_TR
import ..Trapezoids: SliceSelect, Trapezoid
import ...Variables: get_pulse, qval, apply_simple_constraint!, effective_time
import ...Components: RFPulseComponent

"""
    SliceSelectRephase(pulse; kwargs...)

Creates an excitatory RF pulse with slice selection and a rephasing gradient.

Parameters are the same as for [`SliceSelect`](@ref).
"""
struct SliceSelectRephase{N} <: BaseSequence{2}
    slice_select :: SliceSelect{N}
    rephase :: Trapezoid{N}
end

function SliceSelectRephase(pulse::RFPulseComponent; orientation=nothing, group=nothing, slice_thickness=nothing, kwargs...)
    N = isnothing(orientation) ? 3 : 1
    rephase_orientation = isnothing(orientation) ? nothing : -orientation 
    res = SliceSelectRephase{N}(
        SliceSelect(pulse; orientation=orientation, group=group, slice_thickness=slice_thickness, kwargs...),
        Trapezoid(; orientation=rephase_orientation, group=group, kwargs...)
    )
    if N == 1
        apply_simple_constraint!(qval(res.slice_select, :pulse, nothing), qval(res.rephase))
    else
        apply_simple_constraint!(qval(res.slice_select, :pulse, nothing), -qval(res.rephase))
    end
    return res
end

get_index_single_TR(ss::SliceSelectRephase, ::Val{1}) = ss.slice_select
get_index_single_TR(ss::SliceSelectRephase, ::Val{2}) = ss.rephase
get_pulse(ssr::SliceSelectRephase) = ssr.slice_select
effective_time(ssr::SliceSelectRephase) = effective_time(ssr, 1)

end