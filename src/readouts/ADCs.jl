module ADCs
import JuMP: @constraint
import ...Variables: variables, dwell_time, duration, effective_time, get_free_variable, VariableType
import ...BuildingBlocks: BuildingBlock, apply_simple_constraint!, set_simple_constraints!
import ...BuildSequences: global_model

"""
    ADC(nsamples; dwell_time=nothing, duration=nothing, time_to_center=duration/2)

Adds a readout event with `nsamples` readouts.

## Parameters
- `nsamples`: number of samples in the readout.
- `center_halfway`: by default the `time_to_center` is assumed to be half of the `duration`. Set this to false to disable this assumption.

## Variables
- `dwell_time`: Time between each readout sample in ms.
- `duration`: Total duration of the ADC event in ms.
- `time_to_center`: time till the center of k-space.
"""
struct ADC <: BuildingBlock
    nsamples :: Integer
    dwell_time :: VariableType
    time_to_center :: Union{Nothing, VariableType}
end

function ADC(nsamples; dwell_time=nothing, time_to_center=nothing, center_halfway=true, kwargs...)
    res = ADC(
        nsamples,
        get_free_variable(dwell_time),
        get_free_variable(time_to_center),
    )
    @constraint global_model() res.dwell_time >= 0
    if center_halfway
        @constraint global_model() 2 * res.time_to_center == duration(res)
    else
        @constraint global_model() res.time_to_center >= 0
        @constraint global_model() res.time_to_center <= duration(res)
    end
    set_simple_constraints!(res, kwargs)
    return res
end

dwell_time(adc::ADC) = adc.dwell_time
duration(adc::ADC) = adc.nsamples * dwell_time(adc)
time_to_center(adc::ADC) = adc.time_to_center
effective_time(adc::ADC) = time_to_center(adc)

variables(::Type{<:ADC}) = [dwell_time, duration, time_to_center]

end