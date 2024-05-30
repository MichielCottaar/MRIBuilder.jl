module CompositePulses
import JuMP: @constraint
import ...AbstractTypes: RFPulseComponent, split_timestep, edge_times
import ....BuildSequences: global_model
import ....Variables: VariableType, set_simple_constraints!, make_generic, get_free_variable, adjust_internal, variables, @defvar
import ..GenericPulses: GenericPulse

"""
    CompositePulse(; variables...)

A composite RF pulse formed by repeating a base RF pulse.

# Parameters
- `base_pulse`: The base RF pulse that will be repeated.
- `nweights`: The number of repeated pulses. This will be ignored if a vector of `weights` is explicitly provided.

# Variables
- `weights`: The weight of each of the base RF pulses.
- `pulse_time`: Time between the center of the RF pulses.
- `scale_amplitude`: How strongly one should scale the amplitude versus the duration to achieve the desired weights. If set to 1 only the RF pulse amplitude will be scaled. If set to 0 only the RF pulse duration will be scaled.
"""
struct CompositePulse <: RFPulseComponent
    pulses :: Vector{RFPulseComponent}
    pulse_time :: VariableType
    scale_amplitude :: VariableType
end

function CompositePulse(; base_pulse::RFPulseComponent, nweights=nothing, weights=nothing, pulse_time=nothing, scale_amplitude=nothing, variables...)
    if isnothing(weights)
        if isnothing(nweights)
            error("Either `nweights` or `weights` should be set when constructing a CompositePulse.")
        end
        weights = [get_free_variable(nothing) for _ in 1:nweights]
    elseif !isnothing(nweights)
        @assert length(weights) == nweights
    end
    scale_amplitude = get_free_variable(scale_amplitude)
    pulses = [
        adjust_internal(
            base_pulse; 
            scale=w * scale_amplitude + sign(w) * (1 - scale_amplitude),
            stretch=abs(w) * (1 - scale_amplitude) + scale_amplitude,
        )
        for w in weights
    ]
    res = CompositePulse(
        pulses,
        get_free_variable(pulse_time),
        scale_amplitude
    )
    return res
    @constraint global_model() minimum(wait_times(res)) >= 0.
    if !(res.scale_amplitude isa Number)
        @constraint global_model() res.scale_amplitude >= 0.
        @constraint global_model() res.scale_amplitude <= 1.
    end
    return res
end

Base.length(comp::CompositePulse) = length(comp.pulses)
function wait_times(comp::CompositePulse)
    d = duration.(comp.pulses)
    comp.pulse_time - (d[1:end-1] + d[2:end]) / 2
end

@defvar duration(pulse::CompositePulse) = (
        0.5 * variables.duration(pulse.pulses[1]) + 
        0.5 * variables.duration(pulse.pulses[end]) + 
        pulse.pulse_time * (length(pulse) - 1)
    )

@defvar begin
    flip_angle(pulse::CompositePulse) = sum(variables.flip_angle.(pulse.pulses))
    effective_time(pulse::CompositePulse) = variables.duration(pulse) / 2
end


function get_pulse_index(comp::CompositePulse, time::Number)
    t_first_center = 0.5 * variables.duration(comp.pulses[1])
    index = Int(div(time - t_first_center, comp.pulse_time, RoundNearest)) + 1
    if index < 1 || index > length(comp)
        return (index, 0.)
    else
        relative_to_center = time - t_first_center - (index - 1) * comp.pulse_time
        relative_to_start = relative_to_center + 0.5 * variables.duration(comp.pulses[index])
        return (index, relative_to_start)
    end
end

for (fn, default) in [
    (:amplitude, 0.),
    (:phase, NaN),
    (:frequency, NaN),
]
    @eval function variables.$fn.f(pulse::CompositePulse, time::Number)
        (index, rtime) = get_pulse_index(pulse, time)
        if (
            index < 1 || index > length(pulse) ||
            rtime < 0 || rtime > variables.duration(pulse.pulses[index])
        )
            return $default
        end
        return $fn(pulse.pulses[index], rtime)
    end
end

function edge_times(comp::CompositePulse)
    res = Float64[]
    t1 = variables.duration(comp.pulses[1]) * 0.5
    for (index, pulse) in enumerate(comp.pulses)
        d = variables.duration(pulse)
        push!(res, t1 - d/2 + comp.pulse_time * (index - 1))
        push!(res, t1 + d/2 + comp.pulse_time * (index - 1))
    end
    return res
end


function make_generic(comp::CompositePulse)
    sub_generics = make_generic.(comp.pulses)
    t1 = variables.duration(comp.pulses[1]) * 0.5
    times = Float64[]
    amplitude = Float64[]
    phase = Float64[]
    for (index, generic) in enumerate(sub_generics)
        start_time = t1 + (index - 1) * comp.pulse_time - variables.duration(comp.pulses[index]) / 2
        push!(times, start_time)
        push!(amplitude, 0.)
        push!(phase, generic.phase[1])
        append!(times, start_time .+ generic.time)
        append!(amplitude, generic.amplitude)
        append!(phase, generic.phase)
        push!(times, times[end])
        push!(amplitude, 0.)
        push!(phase, generic.phase[end])
    end
    return GenericPulse(times, amplitude, phase)
end

split_timestep(comp::CompositePulse, precision) = minimum(split_timestep.(comp.pulses, precision))

function adjust_internal(comp::CompositePulse; stretch=1., kwargs...)
    return CompositePulse(
        adjust_internal.(comp.pulses; stretch=stretch, kwargs...),
        comp.pulse_time * stretch,
        comp.scale_amplitude,
    )
end

end