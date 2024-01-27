module FixedGradients

import ...BuildingBlock: GradientBlock, fixed
import ...Variables: variables, duration, qval


"""
    FixedGradient(time, Gx, Gy, Gz)
    FixedGradient(time, Gx)
    FixedGradient(time, gradients)

Create a Gradient profile that has been fully defined by N control point.
All arguments should be arrays of the same length N defining these control points.

- `time`: time since the start of this [`BuildingBlock`](@ref) in ms.
- `Gx`: gradient in the x-direction at every timepoint.
- `Gy`: gradient in the y-direction at every timepoint (assumed zero if not provided).
- `Gz`: gradient in the z-direction at every timepoint (assumed zero if not provided).
- `gradients`: length-3 array at every timepoint defining the x-, y-, and z- direction of the gradient.
"""
struct FixedGradient <: GradientBlock
    time :: Vector{Float64}
    Gx :: Vector{Float64}
    Gy :: Vector{Float64}
    Gz :: Vector{Float64}
    rotate :: Bool
    function FixedGradient(time::AbstractVector{<:Number}, Gx::AbstractVector{<:Number}, Gy::AbstractVector{<:Number}, Gz::AbstractVector{<:Number}; rotate=false)
        @assert length(time) == length(Gx)
        @assert length(time) == length(Gy)
        @assert length(time) == length(Gz)
        new(Float64.(time), Float64.(Gx), Float64.(Gy), Float64.(Gz), rotate)
    end
end

function FixedGradient(time::AbstractVector{<:Number}, arr::AbstractVector{<:AbstractVector{<:Number}}; kwargs...) 
    @assert all(length.(arr) .== 3)
    FixedGradient(
        time,
        [a[1] for a in arr],
        [a[2] for a in arr],
        [a[3] for a in arr];
        kwargs...
    )
end

FixedGradient(time::AbstractVector{<:Number}, Gx::AbstractVector{<:Number}; kwargs...) = FixedGradient(time, Gx, zeros(length(time)), zeros(length(time)); kwargs...)

variables(::Type{<:FixedGradient}) = []

duration(fg::FixedGradient) = maximum(fg.time)

Base.show(io::IO, fb::FixedGradient) = print(io, "FixedGradient for $(duration(fb)) ms")

"""
    FixedInstantGradient(orientation, qval)

Instantaneous MR gradient with no free variables.
"""
struct FixedInstantGradient <: GradientBlock
    orientation :: Any
    qval :: Number
end

duraction(instant::FixedInstantGradient) = 0.
qval(instant::FixedInstantGradient) = instant.qval

fixed(f::Union{FixedGradient, FixedInstantGradient}) = f


end