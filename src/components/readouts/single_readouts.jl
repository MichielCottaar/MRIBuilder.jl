module SingleReadouts
import JuMP: @constraint
import ...AbstractTypes: ReadoutComponent
import ....Variables: duration, effective_time, make_generic
import ..ADCs: ADC, readout_times

"""
    SingleReadout()

Represents an instantaneous `Readout` of the signal.

It has no parameters or variables to set.
"""
struct SingleReadout <: ReadoutComponent
    SingleReadout(s::Symbol) = new()
end

function SingleReadout(; duration=nothing, kwargs...)
    bad_keys = [key for (key, value) in pairs(kwargs) if !(isnothing(value) || (value isa Union{AbstractArray, Tuple} && all(isnothing.(value))))]
    if length(bad_keys) > 0
        error("SingleReadout does not expect any variables, yet it received $bad_keys.")
    end
    if !(duration in (nothing, :min, :max, 0, 0.))
        error("SingleReadout did not expect a value for duration, yet received $duration.")
    end
    return SingleReadout(:prevent_loop)
end

duration(::SingleReadout) = 0.
effective_time(::SingleReadout) = 0.
readout_times(::SingleReadout) = [0.]

make_generic(::SingleReadout) = ADC(1, 0., 0., 0.)
end