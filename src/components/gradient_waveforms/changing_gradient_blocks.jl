module ChangingGradientBlocks
import StaticArrays: SVector
import Rotations: RotMatrix3
import LinearAlgebra: I, norm
import ....Variables: @defvar, VariableType, variables, get_free_variable, adjust_internal
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

@defvar duration(cgb::ChangingGradient) = cgb.duration

@defvar gradient begin
    grad_start(cgb::ChangingGradient1D) = cgb.gradient_strength_start .* cgb.orientation
    grad_start(cgb::ChangingGradient3D) = cgb.gradient_strength_start
    slew_rate(cgb::ChangingGradient1D) = cgb.slew_rate .* cgb.orientation
    slew_rate(cgb::ChangingGradient3D) = cgb.slew_rate
    grad_end(cgb::ChangingGradient) = variables.grad_start(cgb) .+ variables.slew_rate(cgb) .* variables.duration(cgb)
    gradient_strength_norm(cgb::ChangingGradient1D) = max(
        abs(cgb.gradient_strength_start),
        abs(cgb.gradient_strength_start + cgb.slew_rate * cgb.duration)
    )
    gradient_strength_norm(cgb::ChangingGradient3D) = max(
        norm(cgb.grad_start),
        norm(cgb.grad_end)
    )
    qvec(cgb::ChangingGradient) = (variables.grad_start(cgb) .+ variables.grad_end(cgb)) .* (variables.duration(cgb) * π)

    gradient_strength(cgb::ChangingGradient, time::Number) = variables.slew_rate(cgb) .* time .+ variables.grad_start(cgb)
end

_mult(g1::VariableType, g2::VariableType) = g1 * g2
_mult(g1::AbstractVector, g2::AbstractVector) = g1 .* permutedims(g2)

@defvar gradient begin
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
        grad_aver = 2 .* variables.grad_start(cgb) .+ variables.grad_end(cgb)
        return (
            _mult(qstart, qstart) .* variables.duration(cgb) .+
            variables.duration(cgb)^2 .* _mult(qstart, grad_aver) .* 2π ./ 3 .+
            variables.bmat_gradient(cgb)
        )
    end

    function bmat_gradient(cgb::ChangingGradient)
        gs = variables.grad_start(cgb)
        diff = variables.slew_rate(cgb) .* variables.duration(cgb)
        return (2π)^2 .* (
            _mult(gs, gs) ./ 3 .+
            _mult(gs, diff) ./ 4 .+
            _mult(diff, diff) ./ 10
        ) .* variables.duration(cgb)^3
    end
end


"""
    split_gradient(constant/changing_gradient_block, times...)

Split a single gradient at a given times.

All times are relative to the start of the gradient block (in ms).
Times are assumed to be in increasing order and between 0 and the duration of the gradient block.

For N times this returns a vector with the N+1 replacement [`ConstantGradient`](@ref) or [`ChangingGradient`](@ref) objects.
"""
function split_gradient(cgb::ChangingGradient, times::VariableType...)
    all_times = [0., times...]
    durations = [times[1], [t[2] - t[1] for t in zip(times[1:end-1], times[2:end])]..., variables.duration(cgb) - times[end]]
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
