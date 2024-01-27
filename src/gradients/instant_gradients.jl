module InstantGradients
import JuMP: @constraint, @variable, Model, owner_model
import ...Variables: qval, bval, start_time, duration, variables, get_free_variable, VariableType
import ...BuildingBlocks: GradientBlock
import ...ConcreteBlocks: to_concrete_block, AbstractConcreteBlock
import ...BuildSequences: @global_model_constructor

"""
    InstantGradientBlock(; orientation=:bvec, qval=nothing)

Defines an instantaneous gradient.

## Parameters
- `orientation` sets the gradient orienation. Can be set to a vector for a fixed orientation. Alternatively, can be set to :bvec (default) to rotate with the user-provided `bvecs` or to :neg_bvec to always be the reverse of the `bvecs`.

## Variables
- [`qval`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um.
"""
struct InstantGradientBlock <: GradientBlock
    model::Model
    orientation :: Any
    qval :: VariableType
end

@global_model_constructor InstantGradientBlock

function InstantGradientBlock(model::Model; orientation=:bvec, qval=nothing)
    res = InstantGradientBlock(
        model,
        orientation,
        get_free_variable(model, qval),
    )
    @constraint model model.qval >= 0
    return res
end


qval(instant::InstantGradientBlock) = instant.qval
bval(instant::InstantGradientBlock) = 0.
duration(instant::InstantGradientBlock) = 0.
variables(::Type{<:InstantGradientBlock}) = [qval]

"""
    ConcreteInstantGradient(orientation, qval)

Instantaneous MR gradient with no free variables.

See [`InstantGradientBlock`](@ref) for a version where [`qval`](@ref) is variable.
"""
struct ConcreteInstantGradient <: AbstractConcreteBlock
    orientation :: Any
    qval :: Number
end

function to_concrete_block(block::InstantGradientBlock)
    return ConcreteInstantGradient(block.orientation, value(qval(block)))
end


end