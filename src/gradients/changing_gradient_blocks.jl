module ChangingGradientBlocks
import StaticArrays: SVector
import ...Variables: VariableType, variables
import ...BuildingBlocks: GradientBlock
import JuMP: @constraint, @variable, Model, owner_model
import ...Variables: qvec, bmat_gradient, gradient_strength, slew_rate, duration, variables, VariableType
import ...BuildingBlocks: GradientBlock, fixed, RFPulseBlock

"""
    ChangingGradientBlock(grad1, slew_rate, duration, rotate, scale)

Underlying type for any linearly changing part in a gradient waveform.

Usually, you do not want to create this object directly, use a gradient waveform instead.

## Parameters
- `rotate`: with which user-set parameter will this gradient be rotated (e.g., :bvec). Default is no rotation.
- `scale`: with which user-set parameter will this gradient be scaled (e.g., :bval). Default is no scaling.
"""
struct ChangingGradientBlock <: GradientBlock
    gradient_strength_start :: SVector{3, <:VariableType}
    slew_rate :: SVector{3, <:VariableType}
    duration :: VariableType
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end


duration(cgb::ChangingGradientBlock) = cgb.duration

grad_start(cgb::ChangingGradientBlock) = cgb.gradient_strength_start
slew_rate(cgb::ChangingGradientBlock) = cgb.slew_rate
grad_end(cgb::ChangingGradientBlock) = grad_start(cgb) .+ slew_rate(cgb) .* duration(cgb)
gradient_strength(cgb::ChangingGradientBlock) = max.(grad_start(cgb), grad_end(cgb))
qvec(cgb::ChangingGradientBlock) = (grad_start(cgb) .+ grad_end(cgb)) .* (duration(cgb) * π)

function bmat_gradient(cgb::ChangingGradientBlock, qstart)
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
        qstart .* permutedims(qstart) .* duration(cgb) .+
        duration(cgb)^2 .* qstart .* permutedims(
            2 .* grad_start(cgb) .+
            grad_end(cgb)) .* 2π ./ 3 .+
        bmat_gradient(cgb)
    )
end

function bmat_gradient(cgb::ChangingGradientBlock)
    diff = slew_rate(cgb) .* duration(cgb)
    return (2π)^2 .* (
        grad_start(cgb) .* permutedims(grad_start(cgb)) ./ 3 .+
        grad_start(cgb) .* permutedims(diff) ./ 4 .+
        diff .* permutedims(diff) ./ 10
    ) .* duration(cgb)^3
end


variables(::Type{<:ChangingGradientBlock}) = [duration, slew_rate, gradient_strength, qvec]

end
