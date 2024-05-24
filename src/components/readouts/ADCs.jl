module ADCs
import JuMP: @constraint, value
import ...AbstractTypes: ReadoutComponent
import ....BuildSequences: global_model, fixed
import ....Variables: VariableType, apply_simple_constraint!, set_simple_constraints!, get_free_variable, make_generic, variables, @defvar


"""
    ADC(; center_halfway=true, oversample=1, variables...)

Adds a readout event.

## Parameters
- `center_halfway`: by default the `time_to_center` is assumed to be half of the `duration`. Set this to false to disable this assumption.
- `oversample`: by how much the ADC should oversample (minimum of 1).

## Variables
- `resolution`: number of voxels in the readout direction. This can be a non-integer value during optimisation.
- `nsamples`: number of samples in the readout. This can be a non-integer value during optimisation.
- `dwell_time`: Time between each readout sample in ms.
- `duration`: Total duration of the ADC event in ms.
- `time_to_center`: time till the center of k-space from start of ADC in ms.
- `effective_time`: same as `time_to_center`.
"""
struct ADC <: ReadoutComponent
    resolution :: VariableType
    dwell_time :: VariableType
    time_to_center :: VariableType
    oversample :: VariableType
end

function ADC(; resolution=nothing, dwell_time=nothing, time_to_center=nothing, center_halfway=true, oversample=nothing, kwargs...)
    res = ADC(
        get_free_variable(resolution),
        get_free_variable(dwell_time),
        get_free_variable(time_to_center),
        get_free_variable(oversample),
    )
    @constraint global_model() res.dwell_time >= 0
    @constraint global_model() res.oversample >= 1
    if center_halfway
        apply_simple_constraint!(duration(res), 2 * res.time_to_center)
    else
        @constraint global_model() res.time_to_center >= 0
        @constraint global_model() res.time_to_center <= duration(res)
    end
    set_simple_constraints!(res, kwargs)
    return res
end

@defvar begin
    readout_times(adc::ADC) = ((1:Int(nsamples(adc))) .- 0.5) .* dwell_time(adc)
    oversample(adc::ADC) = adc.oversample
    nsamples(adc::ADC) = resolution(adc) * oversample(adc)
    dwell_time(adc::ADC) = adc.dwell_time
    duration(adc::ADC) = nsamples(adc) * dwell_time(adc)
    time_to_center(adc::ADC) = adc.time_to_center
    effective_time(adc::ADC) = time_to_center(adc)
    resolution(adc::ADC) = adc.resolution
end

function fixed(adc::ADC)
    # round nsamples during fixing
    r = Int(round(value(resolution(adc)), RoundNearest))
    n = Int(round(value(nsamples(adc)), RoundNearest))
    if iszero(n)
        return ADC(0, NaN, NaN, NaN)
    end
    oversample = n // r
    return ADC(r, value(adc.dwell_time), value(adc.time_to_center), oversample)
end

make_generic(adc::ADC) = adc

end