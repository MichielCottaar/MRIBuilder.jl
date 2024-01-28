module InstantGradients
import JuMP: @constraint, @variable, Model, owner_model, AbstractJuMPScalar
import ...Variables: qvec, bmat_gradient, duration, variables, get_free_variable, VariableType
import ...BuildingBlocks: GradientBlock, fixed
import ...BuildSequences: @global_model_constructor
import ..FixedGradients: FixedInstantGradient

"""
    InstantGradientBlock(; orientation=nothing, qval=nothing, qvec=nothing, rotate=nothing, scale=nothing)

Defines an instantaneous gradient.

## Parameters
- `orientation` sets the gradient orienation (ignored if `qvec` is set). Can be set to a vector for a fixed orientation. Otherwise the orientation will be aligned with the `rotate` (if set) or fully free (if `rotate` is nothing).
- `rotate`: with which user-set parameter will this gradient be rotated (e.g., :bvec). Default is no rotation.
- `scale`: with which user-set parameter will this gradient be scaled (e.g., :bval). Default is no scaling.

## Variables
- [`qvec`](@ref): Spatial scale and direction on which spins will be dephased due to this pulsed gradient in rad/um.
- [`qval`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um.
"""
struct InstantGradientBlock <: GradientBlock
    model::Model
    qvec :: SVector{3, VariableType}
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end

@global_model_constructor InstantGradientBlock

function InstantGradientBlock(model::Model; orientation=nothing, qval=nothing, qvec=nothing, rotate=nothing, scale=nothing)
    used_qval = false
    if isnothing(qvec)
        if isnothing(orientation)
            if isnothing(rotate)
                qvec = SVector{3}(
                    get_free_variable(nothing),
                    get_free_variable(nothing),
                    get_free_variable(nothing),
                )
            else
                qvec = SVector{3}(
                    get_free_variable(qval),
                    0.,
                    0.,
                )
                used_qval = true
            end
        else
            qval = get_free_variable(qval)
            qvec = (orientation ./ norm(orientation)) * qval
            used_qval = true
        end
    end

    res = InstantGradientBlock(
        model,
        qvec,
        rotate,
        scale
    )
    if !used_qval && !isnothing(qval)
        @constraint model qval_sqr(res) == qval^2
    end
    if qval isa AbstractJuMPScalar
        @constraint model qval >= 0
    end
    return res
end


qvec(instant::InstantGradientBlock) = instant.qvec
bmat_gradient(instant::InstantGradientBlock, qstart=nothing) = zeros(3, 3)
duration(instant::InstantGradientBlock) = 0.
variables(::Type{<:InstantGradientBlock}) = [qvec, qval]

function fixed(block::InstantGradientBlock)
    return FixedInstantGradient(block.orientation, value(qval(block)))
end


end