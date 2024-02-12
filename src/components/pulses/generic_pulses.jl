module GenericPulses

import ...AbstractTypes: RFPulseComponent
import ....Variables: duration, amplitude, effective_time, flip_angle


"""
    GenericPulse(time, amplitude, phase, effective_time=<halfway>)
    GenericPulse(time, amplitude; phase=0., frequency=0., effective_time=<halfway>)

Create a Pulse profile that has been fully defined by N control point.

All arguments should be arrays of the same length N defining these control points.

This pulse has no free variables.

- `time`: time since the start of this [`BuildingBlock`](@ref) in ms.
- `amplitude`: amplitude of the RF pulse at every timepoint in kHz.
- `phase`: phase of the RF pulse at every timpoint in degrees. If not set explicitly it will be determined by the provided starting `phase` (degrees) and the `frequency` (kHz).
- `effective_time`: the time that the RF pulse should be considered to have taken place when computing a `Pathway` (defaults: whenever half of the final flip angle has been achieved for on-resonance spins).
"""
struct GenericPulse <: RFPulseComponent
    time :: Vector{Float64}
    amplitude :: Vector{Float64}
    phase :: Vector{Float64}
    effective_time :: Float64
    function GenericPulse(time::AbstractVector{<:Number}, amplitude::AbstractVector{<:Number}, phase::AbstractVector{<:Number}, effective_time=nothing)
        @assert length(time) == length(amplitude)
        @assert length(time) == length(phase)
        @assert all(time .>= 0)
        if isnothing(effective_time)
            test_res = new(Float64.(time), Float64.(amplitude), Float64.(phase), 0.)
            effective_time = time_halfway_flip(test_res)
        end
        return new(Float64.(time), Float64.(amplitude), Float64.(phase), Float64(effective_time))
    end
end

function GenericPulse(time::AbstractVector{<:Number}, amplitude::AbstractVector{<:Number}; phase=0., frequency=0., effective_time=nothing)
    return GenericPulse(time, amplitude, (time .- time[1]) .* (frequency * 360) .+ phase, effective_time)
end

duration(fp::GenericPulse) = maximum(fp.time)
amplitude(fp::GenericPulse) = maximum(abs.(fp.amplitude))
effective_time(pulse::GenericPulse) = pulse.time[argmax(abs.(pulse.amplitude))]
flip_angle(pulse::GenericPulse) = sum(get_weights(pulse) .* pulse.amplitude) * 360
function time_halfway_flip(pulse::GenericPulse)
    w = get_weights(pulse)
    flip_so_far = cumsum(w .* pulse.amplitude)
    return pulse.times[findfirst(f -> f >= flip_so_far[end] / 2, flip_so_far)]
end

function get_weights(pulse::GenericPulse)
    weights_double = pulse.time[2:end] - pulse.time[1:end-1]
    return [weights_double[1] / 2, ((weights_double[1:end-1] + weights_double[2:end]) / 2)..., weights_double[end]/2]
end



make_generic(gp::GenericPulse) = gp


end
