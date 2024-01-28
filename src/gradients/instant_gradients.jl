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
- `orientation` sets the gradient orienation (ignored if `qvec` is set). Can be set to a vector for a fixed orientation. Otherwise the orientation will be aligned with the `rotate` (if set) or fully free (if `rotate` is nothing). Set to :flip to point in the inverse of the user-provided `rotate`.
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

function interpret_orientation(orientation, qval, qvec, will_rotate)
    if !isnothing(qvec)
        return (false, qvec)
    end
    if orientation == :flip
        @assert will_rotate "setting `orientation=:flip` only makes sense if the `rotate` is set as well."
        return (true, SVector{3}(
            -get_free_variable(qval),
            0.,
            0.,
        ))
    end
    if isnothing(orientation)
        if !will_rotate
            return (false, SVector{3}(
                get_free_variable(nothing),
                get_free_variable(nothing),
                get_free_variable(nothing),
            ))
        else
            return (true, SVector{3}(
                get_free_variable(qval),
                0.,
                0.,
            ))
        end
    else
        qval = get_free_variable(qval)
        return (true, (orientation ./ norm(orientation)) * qval)
    end
end

function InstantGradientBlock(model::Model; orientation=nothing, qval=nothing, qvec=nothing, rotate=nothing, scale=nothing)
    (used_qval, qvec) = interpret_orientation(orientation, qval, qvec, !isnothing(rotate))
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
    return FixedInstantGradient(qvec(block), block.rotate, block.scale)
end


end