module FixedPulses

import Printf: @sprintf
import ...BuildingBlocks: RFPulseBlock, fixed, BuildingBlockPrinter
import ...Variables: variables, duration, amplitude, effective_time, flip_angle


"""
    FixedPulse(time, amplitude, phase)
    FixedPulse(time, amplitude; phase=0., frequency=0.)

Create a Pulse profile that has been fully defined by N control point.
All arguments should be arrays of the same length N defining these control points.

- `time`: time since the start of this [`BuildingBlock`](@ref) in ms.
- `amplitude`: amplitude of the RF pulse at every timepoint in kHz.
- `phase`: phase of the RF pulse at every timpoint in degrees. If not set explicitly it will be determined by the provided starting `phase` (degrees) and the `frequency` (kHz).
"""
struct FixedPulse <: RFPulseBlock
    time :: Vector{Float64}
    amplitude :: Vector{Float64}
    phase :: Vector{Float64}
    function FixedPulse(time::AbstractVector{<:Number}, amplitude::AbstractVector{<:Number}, phase::AbstractVector{<:Number})
        @assert length(time) == length(amplitude)
        @assert length(time) == length(phase)
        @assert all(time .>= 0)
        new(Float64.(time), Float64.(amplitude), Float64.(phase))
    end
end

function FixedPulse(time::AbstractVector{<:Number}, amplitude::AbstractVector{<:Number}; phase=0., frequency=0.)
    return FixedPulse(time, amplitude, (time .- time[1]) .* (frequency * 360) .+ phase)
end

duration(fp::FixedPulse) = maximum(fp.time)
amplitude(fp::FixedPulse) = maximum(abs.(fp.amplitude))
effective_time(pulse::FixedPulse) = pulse.time[argmax(abs.(pulse.amplitude))]
function flip_angle(pulse::FixedPulse)
    weights_double = pulse.time[2:end] - pulse.time[1:end-1]
    weights = [weights_double[1] / 2, ((weights_double[1:end-1] + weights_double[2:end]) / 2)..., weights_double[end]/2]
    return sum(weights .* pulse.amplitude) * 360
end

function Base.show(io::IO, printer::BuildingBlockPrinter{<:FixedPulse})
    fp = printer.bb
    pp(value::Number) = @sprintf("%.3g", value)

    print(io, "FixedPulse(")
    if isnothing(printer.start_time)
        print(io, "duration=", pp(duration(fp)))
    else
        print(io, "t=", pp(printer.start_time), "-", pp(printer.start_time + duration(fp)))
    end
    print(io, ", flip_angle=", pp(flip_angle(fp)), ", effective_time=", pp(effective_time(fp)), ")")
end


"""
    FixedInstantPulse(flip_angle, phase)

Fixed version of `InstantRFPulseBlock`, where both `flip_angle` and `phase` are actual numbers.
"""
struct FixedInstantPulse <: RFPulseBlock
    flip_angle :: Number
    phase :: Number
end

duration(fip::FixedInstantPulse) = 0.

fixed(fp::Union{FixedPulse, FixedInstantPulse}) = fp


end
