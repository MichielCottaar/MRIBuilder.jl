module ADCs
import ...Variables: variables, dwell_time, duration, effective_time, get_free_variable, VariableType
import ...BuildingBlocks: BuildingBlock, apply_simple_constraint!, set_simple_constraints!

"""
    ADC(nsamples; dwell_time=nothing, duration=nothing, time_to_center=duration/2)

Adds a readout event with `nsamples` readouts.

## Parameters
- `nsamples`: number of samples in the readout.

## Variables
- `dwell_time`: Time between each readout sample in ms.
- `duration`: Total duration of the ADC event in ms.
- `time_to_center`: time till the center of k-space (assumed to be `duration/2` if not set explicitly).
"""
struct ADC <: BuildingBlock
    nsamples :: Integer
    dwell_time :: VariableType
    time_to_center :: Union{Nothing, VariableType}
end

function ADC(nsamples; dwell_time=nothing, time_to_center=nothing, kwargs...)
    res = ADC(
        nsamples,
        get_free_variable(dwell_time),
        time_to_center isa VariableType ? time_to_center : nothing
    )
    if !(time_to_center isa VariableType)
        apply_simple_constraint!(effective_time(res), time_to_center)
    end
    set_simple_constraints!(res, kwargs)
    return res
end

dwell_time(adc::ADC) = adc.dwell_time
duration(adc::ADC) = adc.nsamples * dwell_time(adc)
time_to_center(adc::ADC) = isnothing(adc.time_to_center) ? duration(adc) / 2 : adc.time_to_center
effective_time(adc::ADC) = time_to_center(adc)

variables(::Type{<:ADC}) = [dwell_time, duration, effective_time]

end