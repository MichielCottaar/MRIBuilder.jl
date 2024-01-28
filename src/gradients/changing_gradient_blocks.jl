module ChangingGradientBlocks
import StaticArrays: SVector
import ...Variables: VariableType, variables
import ...BuildingBlocks: GradientBlock
import JuMP: @constraint, @variable, Model, owner_model
import ...Variables: qvec, bmat, gradient_strength, slew_rate, duration, variables, VariableType
import ...BuildingBlocks: GradientBlock, fixed, RFPulseBlock

struct AbstractChangingGradientBlock{T<:VariableType} <: GradientBlock
    gradient_strength_start :: SVector{3, T}
    gradient_strength_end :: SVector{3, T}
    duration :: T
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end

"""
    ChangingGradientBlock(duration, gradient_strength, rotate, scale)

Underlying type for any linearly changing part in a gradient waveform.

The corresponding type in a fixed gradient is [`FixedChangingGradientBlock`](@ref)

Do not create this directly, it will be created by the gradient waveforms.
"""
const ChangingGradientBlock = AbstractChangingGradientBlock{VariableType}

"""
    FixedChangingGradientBlock(duration, gradient_strength, rotate, scale)

Fixed equivalent of [`ChangingGradientBlock`](@ref).
"""
const FixedChangingGradientBlock = AbstractChangingGradientBlock{Float64}

duration(cgb::AbstractChangingGradientBlock) = cgb.duration

slew_rate(cgb::AbstractChangingGradientBlock) = (cgb.gradient_strength_end .- cgb.gradient_strength_start) ./ duration(cgb)
qvec(cgb::AbstractChangingGradientBlock) = (cgb.gradient_strength_start .+ cgb.gradient_strength_end) .* (duration(cgb)/2)

function bmat(cgb::AbstractChangingGradientBlock, qstart)
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
        qstart .* qstart' .* duration(cgb) .+
        duration(cgb)^2 .* qstart .* (
            2 .* cgb.gradient_strength_start' .+
            cgb.gradient_strength_end') ./ 3 .+
        bmat(cgb)
    )
end

function bmat(cgb::AbstractChangingGradientBlock)
    diff = cgb.gradient_strength_end .- cgb.gradient_strength_start
    return (
        cgb.gradient_strength_start .* cgb.gradient_strength_start' ./ 3 .+
        cgb.gradient_strength_start .* diff' ./ 4 .+
        diff .* diff' ./ 10
    ) .* duration(cgb)^3
end


variables(::Type{<:AbstractChangingGradientBlock}) = [duration, qvec, bmat, slew_rate]

end
