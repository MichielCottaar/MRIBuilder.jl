module SincPulses
import JuMP: @constraint
import QuadGK: quadgk
import ....BuildSequences: global_model
import ....Variables: VariableType, set_simple_constraints!, make_generic, get_free_variable, adjust_internal, variables, @defvar
import ...AbstractTypes: RFPulseComponent, split_timestep
import ..GenericPulses: GenericPulse

"""
    SincPulse(; Nzeros=3, apodise=true, variables...)

Represents a radio-frequency pulse with a sinc-like amplitude and constant frequency.

## Parameters
- `Nzeros`: Number of zero-crossings on each side of the sinc pulse. Can be set to a tuple with two values to have a different number of zero crossings on the left and the right of the sinc pulse.
- `apodise`: if true (default) applies a Hanning apodising window to the sinc pulse.
- `group`: name of the group to which this pulse belongs. This is used for scaling or adding phases/off-resonance frequencies.

## Variables
- [`flip_angle`](@ref): rotation expected for on-resonance spins in degrees.
- [`duration`](@ref): duration of the RF pulse in ms.
- [`amplitude`](@ref): amplitude of the RF pulse in kHz.
- [`phase`](@ref): phase at the start of the RF pulse in degrees.
- [`frequency`](@ref): frequency of the RF pulse relative to the Larmor frequency (in kHz).
- [`bandwidth`](@ref): width of the rectangular function in frequency space (in kHz). If the `duration` is short (compared with 1/`bandwidth`), this bandwidth will only be approximate.
"""
struct SincPulse <: RFPulseComponent
    apodise :: Bool
    Nzeros :: Tuple{Integer, Integer}
    norm_flip_angle :: Tuple{Float64, Float64}
    amplitude :: VariableType
    phase :: VariableType
    frequency :: VariableType
    lobe_duration :: VariableType
    group :: Union{Nothing, Symbol}
end

function SincPulse(; 
    Nzeros=3, apodise=true,
    amplitude=nothing, phase=nothing, frequency=nothing, lobe_duration=nothing, group=nothing, kwargs...
) 
    if Nzeros isa Number
        Nzeros = (Nzeros, Nzeros)
    end
    res = SincPulse(
        apodise,
        Nzeros,
        integral_nzero.(Nzeros, apodise),
        [get_free_variable(value) for value in (amplitude, phase, frequency, lobe_duration)]...,
        group
    )
    if !(res.amplitude isa Number)
        @constraint global_model() res.amplitude >= 0
    end
    if !(res.lobe_duration isa Number)
        @constraint global_model() res.lobe_duration >= 0
    end
    set_simple_constraints!(res, kwargs)
    return res
end

function normalised_function(x, Nleft, Nright; apodise=false)
    if iszero(x)
        return 1.
    end
    if apodise
        if x < 0
            return (0.54 + 0.46 * cos(π * x / Nleft)) * sin(π * x) / (π * x)
        else
            return (0.54 + 0.46 * cos(π * x / Nright)) * sin(π * x) / (π * x)
        end
    else
        return sin(π * x) / (π * x)
    end
end

function integral_nzero(Nzeros, apodise)
    f = x -> normalised_function(x, Nzeros, Nzeros; apodise=apodise)
    return quadgk(f, 0, Nzeros)[1]
end

@defvar begin
    amplitude(pulse::SincPulse) = pulse.amplitude
    N_left(pulse::SincPulse) = pulse.Nzeros[1]
    N_right(pulse::SincPulse) = pulse.Nzeros[2]
    duration(pulse::SincPulse) = (N_left(pulse) + N_right(pulse)) * lobe_duration(pulse)
    phase(pulse::SincPulse) = pulse.phase
    frequency(pulse::SincPulse) = pulse.frequency
    flip_angle(pulse::SincPulse) = (pulse.norm_flip_angle[1] + pulse.norm_flip_angle[2]) * amplitude(pulse) * lobe_duration(pulse) * 360
    lobe_duration(pulse::SincPulse) = pulse.lobe_duration
    inverse_bandwidth(pulse::SincPulse) = lobe_duration(pulse)
    effective_time(pulse::SincPulse) = N_left(pulse) * lobe_duration(pulse)

    amplitude(pulse::SincPulse, time::Number) = amplitude(pulse) * normalised_function(abs((time - effective_time(pulse))) / lobe_duration(pulse), N_left(pulse), N_right(pulse); apodise=pulse.apodise)
    phase(pulse::SincPulse, time::Number) = phase(pulse) + frequency(pulse) * (time - effective_time(pulse)) * 360.
    frequency(pulse::SincPulse, time::Number) = frequency(pulse)
end

function make_generic(block::SincPulse)
    normed_times = -N_left(block):0.1:N_right(block) + 1e-5
    times = max.(0., (normed_times .+ N_left(block))) .* lobe_duration(block)
    amplitudes = amplitude(block) .* (normalised_function.(normed_times, N_left(block), N_right(block); apodise=block.apodise))
    phases = [frequency(block) .* lobe_duration(block)] .* normed_times .* 360
    return GenericPulse(times, amplitudes, phases, effective_time(block))
end

function split_timestep(block::SincPulse, precision)
    max_second_derivative = π^2/3 * (amplitude(block) / lobe_duration(block))^2
    return sqrt(2 * precision / max_second_derivative)
end

function adjust_internal(block::SincPulse; scale=1., frequency=0., stretch=1.)
    SincPulse(
        block.apodise,
        block.Nzeros,
        block.norm_flip_angle,
        block.amplitude * scale,
        block.phase,
        block.frequency + frequency,
        block.lobe_duration * stretch,
        block.group
    )
end

end