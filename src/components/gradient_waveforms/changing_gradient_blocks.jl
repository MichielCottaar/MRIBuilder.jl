module ChangingGradientBlocks
import StaticArrays: SVector
import Rotations: RotMatrix3
import LinearAlgebra: I
import ....Variables: VariableType, duration, qval, bmat_gradient, gradient_strength, slew_rate, get_free_variable, adjust_internal
import ...AbstractTypes: GradientWaveform


"""
    ChangingGradient(grad1_scalar, slew_rate_scalar, orientation, duration, group=nothing)
    ChangingGradient(grad1_vec, slew_rate_vec, duration, group=nothing)

Underlying type for any linearly changing part in a 1D (first constructor) or 3D (second constructor) gradient waveform.

Usually, you do not want to create this object directly, use a `BuildingBlock` instead.
"""
abstract type ChangingGradient{N} <: GradientWaveform{N} end
(::Type{ChangingGradient})(grad1::VariableType, slew_rate::VariableType, orientation::AbstractVector, duration::VariableType, group=nothing) = ChangingGradient1D(grad1, slew_rate, orientation, duration, group)
(::Type{ChangingGradient})(grad1::AbstractVector, slew_rate::AbstractVector, ::Nothing, duration::VariableType, group=nothing) = ChangingGradient3D(grad1, slew_rate, duration, group)
(::Type{ChangingGradient})(grad1::AbstractVector, slew_rate::AbstractVector, duration::VariableType, group=nothing) = ChangingGradient3D(grad1, slew_rate, duration, group)

struct ChangingGradient1D <: ChangingGradient{1}
    gradient_strength_start :: VariableType
    slew_rate :: VariableType
    orientation :: SVector{3, Float64}
    duration :: VariableType
    group :: Union{Nothing, Symbol}
end

struct ChangingGradient3D <: ChangingGradient{3}
    gradient_strength_start :: SVector{3, <:VariableType}
    slew_rate :: SVector{3, <:VariableType}
    duration :: VariableType
    group :: Union{Nothing, Symbol}
end


duration(cgb::ChangingGradient) = cgb.duration

grad_start(cgb::ChangingGradient) = cgb.gradient_strength_start
slew_rate(cgb::ChangingGradient) = cgb.slew_rate
grad_end(cgb::ChangingGradient) = grad_start(cgb) .+ slew_rate(cgb) .* duration(cgb)
gradient_strength(cgb::ChangingGradient) = max.(grad_start(cgb), grad_end(cgb))
qval(cgb::ChangingGradient) = (grad_start(cgb) .+ grad_end(cgb)) .* (duration(cgb) * π)

gradient_strength(cgb::ChangingGradient, time::Number) = slew_rate(cgb) .* time .+ grad_start(cgb)

_mult(g1::VariableType, g2::VariableType) = g1 * g2
_mult(g1::AbstractVector, g2::AbstractVector) = g1 .* permutedims(g2)

to_vec(cgb::ChangingGradient1D, g::VariableType) = cgb.orientation .* g
to_vec(::ChangingGradient3D, g::AbstractVector) = g

function bmat_gradient(cgb::ChangingGradient, qstart::AbstractVector)
    # grad = (g1 * (duration - t) + g2 * t) / duration
    #      = g1 + (g2 - g1) * t / duration
    # q = qstart + g1 * t + (g2 - g1) * t^2 / (2 * duration)
    # \int dt (qstart + t * grad)^2 = 
    #   qstart^2 * duration +
    #   qstart * g1 * duration^2 +
    #   qstart * (g2 - g1) * duration^2 / 3 +
    #   g1^2 * duration^3 / 3 +
    #   g1 * (g2 - g1) * duration^3 / 4 +
    #   (g2 - g1)^2 * duration^3 / 10
    grad_aver = to_vec(cgb, 2 .* grad_start(cgb) .+ grad_end(cgb))
    return (
        _mult(qstart, qstart) .* duration(cgb) .+
        duration(cgb)^2 .* _mult(qstart, grad_aver) .* 2π ./ 3 .+
        bmat_gradient(cgb)
    )
end

function bmat_gradient(cgb::ChangingGradient)
    gs = to_vec(cgb, grad_start(cgb))
    diff = to_vec(cgb, slew_rate(cgb) .* duration(cgb))
    return (2π)^2 .* (
        _mult(gs, gs) ./ 3 .+
        _mult(gs, diff) ./ 4 .+
        _mult(diff, diff) ./ 10
    ) .* duration(cgb)^3
end


"""
    split_gradient(constant/changing_gradient_block, times...)

Split a single gradient at a given times.

All times are relative to the start of the gradient block (in ms).
Times are assumed to be in increasing order and between 0 and the duration of the gradient block.

For N times this returns a vector with the N+1 replacement [`ConstantGradientBlock`](@ref) or [`ChangingGradientBlock`](@ref) objects.
"""
function split_gradient(cgb::ChangingGradient, times::VariableType...)
    all_times = [0., times...]
    durations = [times[1], [t[2] - t[1] for t in zip(times[1:end-1], times[2:end])]..., duration(cgb) - times[end]]
    if cgb isa ChangingGradient1D
        return [ChangingGradient1D(cgb.gradient_strength .+ cgb.slew_rate .* t, cgb.slew_rate, cgb.orientation, d, cgb.group) for (t, d) in zip(all_times, durations)]
    else
        return [ChangingGradient3D(cgb.gradient_strength .+ cgb.slew_rate .* t, cgb.slew_rate, d, cgb.group) for (t, d) in zip(all_times, durations)]
    end
end

function adjust_internal(cgb::ChangingGradient1D; orientation=nothing, scale=1., rotation=nothing)
    if !isnothing(orientation) && !isnothing(rotation)
        error("Cannot set both the gradient orientation and rotation.")
    end
    new_orientation = isnothing(orientation) ? (isnothing(rotation) ? cgb.orientation : rotation * cgb.orientation) : orientation
    return ChangingGradient1D(
        cgb.gradient_strength_start * scale,
        cgb.slew_rate * scale,
        new_orientation,
        cgb.duration,
        cgb.group
    )
end

function adjust_internal(cgb::ChangingGradient3D; scale=1., rotation=RotMatrix3(I(3)))
    return ChangingGradient3D(
        rotation * (cgb.gradient_strength_start .* scale),
        rotation * (cgb.slew_rate .* scale),
        cgb.duration,
        cgb.group
    )
end

end
