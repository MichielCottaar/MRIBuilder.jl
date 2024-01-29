module FixedGradients

import Printf: @sprintf
import LinearAlgebra: norm
import StaticArrays: SVector
import ...BuildingBlocks: ContainerBlock, fixed, BuildingBlock, BuildingBlockPrinter, get_children_blocks
import ...Variables: variables, duration, qvec, gradient_strength, slew_rate, start_time
import ..ChangingGradientBlocks: FixedChangingGradientBlock


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
    gradient_strength :: Vector{SVector{3, Float64}}
    rotate :: Union{Symbol, Nothing}
    scale :: Union{Symbol, Nothing}
    function FixedGradient(time::AbstractVector{<:Number}, Gx::AbstractVector{<:Number}, Gy::AbstractVector{<:Number}, Gz::AbstractVector{<:Number}; rotate=nothing, scale=nothing)
        @assert length(time) == length(Gx)
        @assert length(time) == length(Gy)
        @assert length(time) == length(Gz)
        grad = [SVector{3, Float64}(Gx[i], Gy[i], Gz[i]) for i in eachindex(Gx)]
        new(Float64.(time), grad, rotate, scale)
    end
end

get_children_indices(fg::FixedGradient) = 1:(length(fg.time)-1)
Base.getindex(fg::FixedGradient, index::Int) = FixedChangingGradientBlock(
    fg.gradient_strength[index],
    fg.gradient_strength[index + 1],
    fg.times[index + 1] - fg.times[index],
    rotate,
    scale,
)
start_time(fg::FixedGradient, index::Int) = fg.times[index]

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

function gradient_strength(fg::FixedGradient)
    if isnothing(fg.rotate)
        return maximum(map(g -> max(abs.(g)...), fg.gradient_strength))
    else
        return maximum(map(norm, fg.gradient_strength))
    end
end

function slew_rate(fg::FixedGradient)
    diff = (fg.gradient_strength[2:end] .- fg.gradient_strength[1:end-1]) ./ (fg.time[2:end] - fg.time[1:end-1])
    if isnothing(fg.rotate)
        return maximum(map(d -> max(abs.(d)...), diff))
    else
        return maximum(map(norm, diff))
    end
end


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
    qvec :: SVector{3, Number}
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end

duraction(instant::FixedInstantGradient) = 0.
qvec(instant::FixedInstantGradient) = instant.qvec

fixed(f::Union{FixedGradient, FixedInstantGradient}) = f


end