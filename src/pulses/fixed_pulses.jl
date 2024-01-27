module FixedPulses

import ...BuildingBlocks: PulseBlock, fixed
import ...Variables: variables, duration


"""
    FixedPulse(time, amplitude, phase)
    FixedPulse(time, amplitude; phase=0., frequency=0.)

Create a Pulse profile that has been fully defined by N control point.
All arguments should be arrays of the same length N defining these control points.

- `time`: time since the start of this [`BuildingBlock`](@ref) in ms.
- `amplitude`: amplitude of the RF pulse at every timepoint in kHz.
- `phase`: phase of the RF pulse at every timpoint in degrees. If not set explicitly it will be determined by the provided starting `phase` (degrees) and the `frequency` (kHz).
"""
struct FixedPulse <: PulseBlock
    time :: Vector{Float64}
    amplitude :: Vector{Float64}
    phase :: Vector{Float64}
    function FixedPulse(time::AbstractVector{<:Number}, amplitude::AbstractVector{<:Number}, phase::AbstractVector{<:Number})
        @assert length(time) == length(amplitude)
        @assert length(time) == length(phase)
        new(Float64.(time), Float64.(amplitude), Float64.(phase))
    end
end

function FixedPulse(time::AbstractVector{<:Number}, amplitude::AbstractVector{<:Number}; phase=0., frequency=0.)
    return FixedPulse(time, amplitude, (time .- time[1]) .* (frequency * 360) .+ phase)
end

variables(::Type{<:FixedPulse}) = []

duration(fp::FixedPulse) = maximum(fp.time)

Base.show(io::IO, fp::FixedPulse) = print(io, "FixedPulse for $(duration(fp)) ms")


"""
    FixedInstantPulse(flip_angle, phase)

Fixed version of `InstantRFPulseBlock`, where both `flip_angle` and `phase` are actual numbers.
"""
struct FixedInstantPulse <: PulseBlock
    flip_angle :: Number
    phase :: Number
end

duration(fip::FixedInstantPulse) = 0.

fixed(fp::Union{FixedPulse, FixedInstantPulse}) = fp


end
