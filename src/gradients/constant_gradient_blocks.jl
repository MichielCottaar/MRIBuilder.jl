module ConstantGradientBlocks
import StaticArrays: SVector
import ...Variables: VariableType, variables
import ...BuildingBlocks: GradientBlock
import JuMP: @constraint, @variable, Model, owner_model
import ...Variables: qvec, bmat_gradient, gradient_strength, slew_rate, duration, variables, VariableType
import ...BuildingBlocks: GradientBlock, fixed, RFPulseBlock

struct AbstractConstantGradientBlock{T<:VariableType} <: GradientBlock
    gradient_strength :: SVector{3, T}
    duration :: T
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end

"""
    ConstantGradientBlock(gradient_strength, duration, rotate, scale)

Underlying type for any flat part in a gradient waveform.
It can overlap with an [`RFPulseBlock`], in which case it will be a `SliceSelectPulse`

The corresponding type in a fixed gradient is [`FixedConstantGradientBlock`](@ref)

Do not create this directly, it will be created by the gradient waveforms.
"""
const ConstantGradientBlock = AbstractConstantGradientBlock{VariableType}

"""
    FixedConstantGradientBlock(gradient_strength, duration, rotate, scale)

Fixed equivalent of [`ConstantGradientBlock`](@ref).
"""
const FixedConstantGradientBlock = AbstractConstantGradientBlock{Float64}

duration(cgb::AbstractConstantGradientBlock) = cgb.duration
gradient_strength(cgb::AbstractConstantGradientBlock) = cgb.gradient_strength
slew_rate(::AbstractConstantGradientBlock) = zero(SVector{3, Float64})
qvec(cgb::AbstractConstantGradientBlock) = duration(cgb) .* gradient_strength(cgb)

function bmat_gradient(cgb::AbstractConstantGradientBlock)
    grad = gradient_strength(cgb)
    return (grad .* permutedims(grad)) .* duration(cgb)^3 ./3
end

function bmat_gradient(cgb::AbstractConstantGradientBlock, qstart)
    # \int dt (qstart + t * grad)^2 = 
    #   qstart^2 * duration +
    #   qstart * grad * duration^2 +
    #   grad * grad * duration^3 / 3 +
    grad = gradient_strength(cgb)
    return (
        qstart .* permutedims(qstart) .* duration(cgb) .+
        qstart .* permutedims(grad) .* duration(cgb)^2 .+
        bmat_gradient(cgb)
    )
end

variables(::Type{<:AbstractConstantGradientBlock}) = [duration, gradient_strength]

end
