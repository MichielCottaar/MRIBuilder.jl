module ConstantGradientBlocks
import StaticArrays: SVector
import ...Variables: VariableType, variables
import ...BuildingBlocks: GradientBlock
import JuMP: @constraint, @variable, Model, owner_model
import ...Variables: qvec, bmat_gradient, gradient_strength, slew_rate, duration, variables, VariableType
import ...BuildingBlocks: GradientBlock, fixed, RFPulseBlock

"""
    ConstantGradientBlock(gradient_strength, duration, rotate, scale)

Underlying type for any flat part in a gradient waveform.
It can overlap with an [`RFPulseBlock`], in which case it will be a `SliceSelectPulse`

Usually, you do not want to create this object directly, use a gradient waveform instead.

## Parameters
- `rotate`: with which user-set parameter will this gradient be rotated (e.g., :bvec). Default is no rotation.
- `scale`: with which user-set parameter will this gradient be scaled (e.g., :bval). Default is no scaling.
"""
struct ConstantGradientBlock{T<:VariableType} <: GradientBlock
    gradient_strength :: SVector{3, T}
    duration :: T
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end

const ConstantGradientBlock = ConstantGradientBlock{VariableType}

"""
    FixedConstantGradientBlock(gradient_strength, duration, rotate, scale)

Fixed equivalent of [`ConstantGradientBlock`](@ref).
"""
const FixedConstantGradientBlock = ConstantGradientBlock{Float64}

duration(cgb::ConstantGradientBlock) = cgb.duration
gradient_strength(cgb::ConstantGradientBlock) = cgb.gradient_strength
slew_rate(::ConstantGradientBlock) = zero(SVector{3, Float64})
qvec(cgb::ConstantGradientBlock) = duration(cgb) .* gradient_strength(cgb) .* 2π

function bmat_gradient(cgb::ConstantGradientBlock)
    grad = 2π .* gradient_strength(cgb)
    return (grad .* permutedims(grad)) .* duration(cgb)^3 ./3
end

function bmat_gradient(cgb::ConstantGradientBlock, qstart)
    # \int dt (qstart + t * grad)^2 = 
    #   qstart^2 * duration +
    #   qstart * grad * duration^2 +
    #   grad * grad * duration^3 / 3 +
    grad = 2π .* gradient_strength(cgb)
    return (
        qstart .* permutedims(qstart) .* duration(cgb) .+
        qstart .* permutedims(grad) .* duration(cgb)^2 .+
        bmat_gradient(cgb)
    )
end

variables(::Type{<:ConstantGradientBlock}) = [duration, gradient_strength, qvec]

end
