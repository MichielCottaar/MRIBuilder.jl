module ChangingGradientBlocks
import StaticArrays: SVector
import ....Variables: VariableType, duration, qvec, bmat_gradient, gradient_strength, slew_rate, get_free_variable
import ...AbstractTypes: GradientWaveform


"""
    ChangingGradient(grad1, slew_rate, duration)

Underlying type for any linearly changing part in a 1D or 3D gradient waveform (depending on whether `grad1` and `slew_rate` are a scalar or a vector)

Usually, you do not want to create this object directly, use a `BuildingBlock` instead.
"""
abstract type ChangingGradient{N} <: GradientWaveform{N} end
(::Type{ChangingGradient})(grad1::VariableType, slew_rate::VariableType, duration::VariableType) = ChangingGradient1D(grad1, slew_rate, duration)
(::Type{ChangingGradient})(grad1::AbstractVector, slew_rate::AbstractVector, duration::VariableType) = ChangingGradient3D(grad1, slew_rate, duration)

struct ChangingGradient1D <: ChangingGradient{1}
    gradient_strength_start :: VariableType
    slew_rate :: VariableType
    duration :: VariableType
end

struct ChangingGradient3D <: ChangingGradient{3}
    gradient_strength_start :: SVector{3, <:VariableType}
    slew_rate :: SVector{3, <:VariableType}
    duration :: VariableType
end


duration(cgb::ChangingGradient) = cgb.duration

grad_start(cgb::ChangingGradient) = cgb.gradient_strength_start
slew_rate(cgb::ChangingGradient) = cgb.slew_rate
grad_end(cgb::ChangingGradient) = grad_start(cgb) .+ slew_rate(cgb) .* duration(cgb)
gradient_strength(cgb::ChangingGradient) = max.(grad_start(cgb), grad_end(cgb))
qvec(cgb::ChangingGradient) = (grad_start(cgb) .+ grad_end(cgb)) .* (duration(cgb) * π)

_mult(g1::VariableType, g2::VariableType) = g1 * g2
_mult(g1::AbstractVector, g2::AbstractVector) = g1 .* permutedims(g2)

function bmat_gradient(cgb::ChangingGradient, qstart)
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
    return (
        _mult(qstart, qstart) .* duration(cgb) .+
        duration(cgb)^2 .* _mult(qstart, 2 .* grad_start(cgb) .+ grad_end(cgb)) .* 2π ./ 3 .+
        bmat_gradient(cgb)
    )
end

function bmat_gradient(cgb::ChangingGradient)
    diff = slew_rate(cgb) .* duration(cgb)
    return (2π)^2 .* (
        _mult(grad_start(cgb), grad_start(cgb)) ./ 3 .+
        _mult(grad_start(cgb), diff) ./ 4 .+
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
    return [cls(cgb)(cgb.gradient_strength .+ cgb.slew_rate .* t, cgb.slew_rate, d) for (t, d) in zip(all_times, durations)]
end


end
