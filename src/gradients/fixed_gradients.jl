module FixedGradients

import Printf: @sprintf
import LinearAlgebra: norm
import ...BuildingBlocks: GradientBlock, fixed, BuildingBlock, BuildingBlockPrinter
import ...Variables: variables, duration, qval, gradient_strength, slew_rate


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

gradient_strength(fg::FixedGradient) = max(maximum(abs.(fg.Gx)), maximum(abs.(fg.Gy)), maximum(abs.(fg.Gz)))

slew_rate(fg::FixedGradient) = maximum(map((fg.Gx, fg.Gy, fg.Gz)) do gradient
    return maximum(abs.(gradient[2:end] - gradient[1:end-1]) ./ (fg.time[2:end] - fg.time[1:end-1]))
end)

qvec(fg::FixedGradient) = map((fg.Gx, fg.Gy, fg.Gz)) do gradient
    weights_double = fg.time[2:end] - fg.time[1:end-1]
    weights = [weights_double[1] / 2, ((weights_double[1:end-1] + weights_double[2:end]) / 2)..., weights_double[end]/2]
    return sum(weights .* gradient)
end

qval(fg::FixedGradient) = norm(qvec(fg))


function Base.show(io::IO, printer::BuildingBlockPrinter{<:FixedGradient})
    fp = printer.bb
    pp(value::Number) = @sprintf("%.3g", value)

    print(io, "FixedGradient(")
    if isnothing(printer.start_time)
        print(io, "duration=", pp(duration(fp)))
    else
        print(io, "t=", pp(printer.start_time), "-", pp(printer.start_time + duration(fp)))
    end
    for fn in (qval, gradient_strength, slew_rate)
        print(io, ", $(nameof(fn))=", pp(fn(fp)))
    end
    print(io, ")")
end

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