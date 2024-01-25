module SincPulses

import JuMP: VariableRef, @constraint, @variable, value
import QuadGK: quadgk
import Polynomials: fit, Polynomial
import ....Sequences: RFPulse
import ...BuildingBlocks: BuildingBlock, properties, BuildingBlockPlaceholder, set_simple_constraints!, duration, to_mcmr_components
import ...SequenceBuilders: SequenceBuilder, owner_model, start_time, end_time
import ..Properties: flip_angle, phase, amplitude, frequency, bandwidth

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
struct SincPulse <: BuildingBlock
    builder :: SequenceBuilder
    symmetric :: Bool
    apodise :: Bool
    nlobe_integral :: Polynomial
    N_left :: VariableRef
    N_right :: VariableRef
    amplitude :: VariableRef
    phase :: VariableRef
    frequency :: VariableRef
    lobe_duration :: VariableRef
end

SincPulse(; kwargs...) = BuildingBlockPlaceholder{SincPulse}(; kwargs...)
function SincPulse(builder::SequenceBuilder; symmetric=true, max_Nlobes=nothing, apodise=true, kwargs...) 
    model = owner_model(builder)
    if symmetric
        N_lobes = @variable(model, integer=true)
        N_left_var = N_right_var = N_lobes
    else
        N_left_var = @variable(model, integer=true)
        N_right_var = @variable(model, integer=true)
    end
    res = SincPulse(
        builder,
        symmetric,
        apodise,
        nlobe_integral_params(max_Nlobes, apodise),
        N_left_var,
        N_right_var,
        @variable(model),
        @variable(model),
        @variable(model),
        @variable(model)
    )
    @constraint model amplitude(res) >= 0
    @constraint model N_left(res) >= 1
    if !symmetric
        @constraint model N_right(res) >= 1
    end
    set_simple_constraints!(model, res, kwargs)
    return res
end

function normalised_function(x; apodise=false)
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
bandwidth(pulse::SincPulse) = 1 / lobe_duration(pulse)
properties(::Type{<:SincPulse}) = [amplitude, N_left, N_right, duration, phase, frequency, flip_angle, lobe_duration, bandwidth]

function to_mcmr_components(block::SincPulse)
    normed_times = -value(N_left(block)):0.1:value(N_right(block)) + 1e-5
    times = ((normed_times .+ value(N_left(block))) .* value(lobe_duration(block))) .+ value(start_time(block))
    amplitudes = value(amplitude(block)) .* (normalised_function.(normed_times; apodise=block.apodise))
    phases = (value(frequency(block)) .* value(lobe_duration(block))) .* normed_times * 360
    return RFPulse(times, amplitudes, phases)
end

end