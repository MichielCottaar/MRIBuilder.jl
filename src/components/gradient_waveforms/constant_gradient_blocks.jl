module ConstantGradientBlocks
import StaticArrays: SVector
import ....Variables: VariableType, duration, qvec, bmat_gradient, gradient_strength, slew_rate, get_free_variable
import ...AbstractTypes: GradientWaveform

"""
    ConstantGradient(gradient_strength, duration)

Underlying type for any flat part in a 1D or 3D gradient waveform (depending on whether `gradient_strength` is a scalar or a vector).

Usually, you do not want to create this object directly, use a `BuildingBlock` instead.
"""
abstract type ConstantGradient{N} <: GradientWaveform{N} end
(::Type{ConstantGradient})(grad1::VariableType, duration::VariableType) = ConstantGradient1D(grad1, duration)
(::Type{ConstantGradient})(grad1::AbstractVector, duration::VariableType) = ConstantGradient3D(grad1, duration)

struct ConstantGradient1D <: ConstantGradient{1}
    gradient_strength :: VariableType
    duration :: VariableType
end

struct ConstantGradient3D <: ConstanGradient{3}
    gradient_strength :: SVector{3, <:VariableType}
    duration :: VariableType
end

duration(cgb::ConstantGradient) = cgb.duration
gradient_strength(cgb::ConstantGradient) = cgb.gradient_strength
slew_rate(::ConstantGradient1D) = 0.
slew_rate(::ConstantGradient3D) = zero(SVector{3, Float64})
qvec(cgb::ConstantGradient1D) = duration(cgb) * gradient_strength(cgb) * 2π
qvec(cgb::ConstantGradient3D) = @. duration(cgb) * gradient_strength(cgb) * 2π

_mult(g1::VariableType, g2::VariableType) = g1 * g2
_mult(g1::AbstractVector, g2::AbstractVector) = g1 .* permutedims(g2)

function bmat_gradient(cgb::ConstantGradient)
    grad = 2π .* gradient_strength(cgb)
    return _mult(grad, grad) .* duration(cgb)^3 ./3
end

function bmat_gradient(cgb::ConstantGradient, qstart)
    # \int dt (qstart + t * grad)^2 = 
    #   qstart^2 * duration +
    #   qstart * grad * duration^2 +
    #   grad * grad * duration^3 / 3 +
    grad = 2π .* gradient_strength(cgb)
    return (
        _mult(qstart, qstart) .* duration(cgb) .+
        _mult(qstart, grad) .* duration(cgb)^2 .+
        bmat_gradient(cgb)
    )
end

function split_gradient(cgb::ConstantGradient, times::VariableType...)
    durations = [times[1], [t[2] - t[1] for t in zip(times[1:end-1], times[2:end])]..., duration(cgb) - times[end]]
    @assert all(durations >= 0.)
    return [typeof(cgb)(cgb.gradient_strength, d) for d in durations]
end

end
