module GenericPulses

import Polynomials: fit
import ...AbstractTypes: RFPulseComponent, split_timestep
import ....Variables: duration, amplitude, effective_time, flip_angle, make_generic, phase


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

"""
    GenericPulse(pulse, t1, t2)

Creates a new [`GenericPulse`](@ref) by slicing another pulse between `t1` and `t2`
"""
function GenericPulse(pulse::GenericPulse, t1::Number, t2::Number)
    if t1 < pulse.time[1] || t2 > pulse.time[end]
        error("Cannot extrapolate GenericPulse")
    end
    use = t1 .<= pulse.time .<= t2
    tnew = pulse.time[use]
    anew = pulse.amplitude[use]
    pnew = pulse.phase[use]
    @show tnew[end], anew[end]
    if !(t1 ≈ tnew[1])
        pushfirst!(tnew, t1)
        pushfirst!(anew, amplitude(pulse, t1))
        pushfirst!(pnew, phase(pulse, t1))
    elseif !(t2 ≈ tnew[end])
        push!(tnew, t2)
        push!(anew, amplitude(pulse, t2))
        push!(pnew, phase(pulse, t2))
    end
    @show tnew[end], anew[end]
    return GenericPulse(tnew .- t1, anew, pnew, pulse.effective_time - t1)
end

GenericPulse(pulse::RFPulseComponent, t1::Number, t2::Number) = GenericPulse(make_generic(pulse), t1, t2)

duration(fp::GenericPulse) = maximum(fp.time)
amplitude(fp::GenericPulse) = maximum(abs.(fp.amplitude))
effective_time(pulse::GenericPulse) = pulse.time[findmax(abs.(pulse.amplitude))]
phase(pulse::GenericPulse) = pulse.phase[findmax(abs.(pulse.amplitude))[2]]
flip_angle(pulse::GenericPulse) = sum(get_weights(pulse) .* pulse.amplitude) * 360
function time_halfway_flip(pulse::GenericPulse)
    w = get_weights(pulse)
    flip_so_far = cumsum(w .* pulse.amplitude)
    return pulse.time[findfirst(f -> f >= flip_so_far[end] / 2, flip_so_far)]
end

for fn in (:amplitude, :phase)
    @eval function $fn(fp::GenericPulse, time::Number)
        i2 = findfirst(t->t > time, fp.time)
        if isnothing(i2)
            @assert time ≈ fp.time[end]
            return fp.$(fn)[end]
        end
        i1 = i2 - 1
        if time ≈ fp.time[i1]
            return fp.$(fn)[i1]
        end
        w2 = (time - fp.time[i1]) / (fp.time[i2] - fp.time[i1])
        w1 = (fp.time[i2] - time) / (fp.time[i2] - fp.time[i1])
        return w1 * fp.$(fn)[i1] + w2 * fp.$(fn)[i2]
    end
end


function get_weights(pulse::GenericPulse)
    weights_double = pulse.time[2:end] - pulse.time[1:end-1]
    return [weights_double[1] / 2, ((weights_double[1:end-1] + weights_double[2:end]) / 2)..., weights_double[end]/2]
end



make_generic(gp::GenericPulse) = gp

function split_timestep(gp::GenericPulse, precision)
    function second_derivative(arr)
        max_second_der = 0.
        for index in 2:length(arr)-1
            poly = fit(gp.times[index-1:index+1], arr[index-1:index+1])
            second_der = abs(poly.coeffs[end])
            if second_der > max_second_der
                max_second_der = second_der
            end
        end
        return max_second_der
    end
    max_second_der = max(
        second_derivative(gp.amplitude ./ maximum(gp.amplitude)),
        second_derivative(gp.phase ./ 360),
    )
    if iszero(max_second_der)
        return Inf
    end
    return sqrt(2 * precision / max_second_der)
end


end
