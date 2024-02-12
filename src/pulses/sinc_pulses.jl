module SincPulses

import JuMP: VariableRef, @constraint, @variable
import QuadGK: quadgk
import Polynomials: fit, Polynomial
import ...BuildingBlocks: RFPulseBlock, set_simple_constraints!
import ...Variables: flip_angle, phase, amplitude, frequency, VariableType, variables, get_free_variable, duration, effective_time, inverse_bandwidth
import ...BuildSequences: global_model
import ..GenericPulses: GenericPulse, make_generic

"""
    SincPulse(; symmetric=true, max_Nlobes=nothing, apodise=true, variables...)

Represents an radio-frequency pulse with a constant amplitude and frequency.

## Parameters
- `symmetric`: by default the sinc pulse will be symmetric (i.e., `N_left` == `N_right`). Set `symmetric=false` to independently control `N_left` and `N_right`.
- `max_Nlobes`: by default the computed [`flip_angle`](@ref) is only approximated as `amplitude` * `lobe_duration`. By setting `max_Nlobes` the flip_angle will be given by the actual integral of the sinc function, which will be more accurate for small number of lobes. However, the number of lobes will not be able to exceed this `max_Nlobes`.
- `apodise`: if true (default) applies a Hanning apodising window to the sinc pulse.

## Variables
- [`N_left`](@ref): number of zero-crossings of the sinc function before the main peak (minimum of 1).
- [`N_right`](@ref): number of zero-crossings of the sinc function after the main peak (minimum of 1).
- [`flip_angle`](@ref): rotation expected for on-resonance spins in degrees (only approximate if `max_Nlobes` is not set).
- [`duration`](@ref): duration of the RF pulse in ms.
- [`amplitude`](@ref): amplitude of the RF pulse in kHz.
- [`phase`](@ref): phase at the start of the RF pulse in degrees.
- [`frequency`](@ref): frequency of the RF pulse relative to the Larmor frequency (in kHz).
- [`bandwidth`](@ref): width of the rectangular function in frequency space (in kHz). If the `duration` is short (compared with 1/`bandwidth`), this bandwidth will only be approximate.
"""
struct SincPulse <: RFPulseBlock
    symmetric :: Bool
    apodise :: Bool
    nlobe_integral :: Polynomial
    N_left :: VariableType
    N_right :: VariableType
    amplitude :: VariableType
    phase :: VariableType
    frequency :: VariableType
    lobe_duration :: VariableType
    scale :: Union{Nothing, Symbol}
end

function SincPulse(; 
    symmetric=true, max_Nlobes=nothing, apodise=true, N_lobes=nothing, N_left=nothing, N_right=nothing, 
    amplitude=nothing, phase=nothing, frequency=nothing, lobe_duration=nothing, scale=nothing, kwargs...
) 
    if symmetric
        N_lobes = get_free_variable(N_lobes)
        @assert isnothing(N_left) && isnothing(N_right) "N_left and N_right cannot be set if symmetric=true (default)"
        N_left_var = N_right_var = N_lobes
    else
        @assert isnothing(N_lobes) "N_lobes cannot be set if symmetric=true (default)"
        N_left_var = get_free_variable(N_left)
        N_right_var = get_free_variable(N_right)
    end
    res = SincPulse(
        symmetric,
        apodise,
        nlobe_integral_params(max_Nlobes, apodise),
        N_left_var,
        N_right_var,
        [get_free_variable(value) for value in (amplitude, phase, frequency, lobe_duration)]...,
        scale
    )
    model = global_model()
    @constraint model res.amplitude >= 0
    @constraint model res.N_left >= 1
    if !symmetric
        @constraint model res.N_right >= 1
    end
    set_simple_constraints!(res, kwargs)
    return res
end

function normalised_function(x; apodise=false)
    if iszero(x)
        return 1.
    end
    if apodise
        return (0.54 + 0.46 * cos(π * x)) * sin(π * x) / (π * x)
    else
        return sin(π * x) / (π * x)
    end
end

function nlobe_integral_params(Nlobe_max, apodise=false)
    if isnothing(Nlobe_max)
        return fit([1], [0.5], 0)
    end
    f = x -> normalised_function(x; apodise=apodise)
    integral_values = [quadgk(f, 0, i)[1] for i in 1:Nlobe_max]
    return fit(1:Nlobe_max, integral_values, Nlobe_max - 1)
end

amplitude(pulse::SincPulse) = pulse.amplitude
N_left(pulse::SincPulse) = pulse.N_left
N_right(pulse::SincPulse) = pulse.N_right
duration(pulse::SincPulse) = (N_left(pulse) + N_right(pulse)) * lobe_duration(pulse)
phase(pulse::SincPulse) = pulse.phase
frequency(pulse::SincPulse) = pulse.frequency
flip_angle(pulse::SincPulse) = (pulse.nlobe_integral(N_left(pulse)) + pulse.nlobe_integral(N_right(pulse))) * amplitude(pulse) * lobe_duration(pulse) * 360
lobe_duration(pulse::SincPulse) = pulse.lobe_duration
inverse_bandwidth(pulse::SincPulse) = lobe_duration(pulse)
effective_time(pulse::SincPulse) = N_left(pulse) * lobe_duration(pulse)

function make_generic(block::SincPulse)
    normed_times = -N_left(block):0.1:N_right(block) + 1e-5
    times = (normed_times .+ N_left(block)) .* lobe_duration(block)
    amplitudes = amplitude(block) .* (normalised_function.(normed_times; apodise=block.apodise))
    phases = [frequency(block) .* lobe_duration(block)] .* normed_times .* 360
    return GenericPulse(times, amplitudes, phases, effective_time(block))
end

end