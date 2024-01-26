module InstantGradients
import JuMP: @constraint, @variable, VariableRef
import ...BuildingBlocks: BuildingBlock, properties, BuildingBlockPlaceholder, set_simple_constraints!, duration
import ...SequenceBuilders: SequenceBuilder, owner_model, start_time
import ...ConcreteBlocks: to_concrete_block, AbstractConcreteBlock
import ..IntegrateGradients: qval, bval

"""
    InstantGradientBlock(; orientation=:bvec, qval=nothing)

Defines an instantaneous gradient.

This is a [`BuildingBlock`](@ref) for the [`SequenceBuilder`](@ref).

## Parameters
- `orientation` sets the gradient orienation. Can be set to a vector for a fixed orientation. Alternatively, can be set to :bvec (default) to rotate with the user-provided `bvecs` or to :neg_bvec to always be the reverse of the `bvecs`.

## Variables
- [`qval`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um.
"""
struct InstantGradientBlock <: BuildingBlock
    builder::SequenceBuilder
    orientation :: Any
    qval :: VariableRef
end

InstantGradientBlock(; kwargs...) = BuildingBlockPlaceholder{InstantGradientBlock}(; kwargs...)

function InstantGradientBlock(builder::SequenceBuilder; orientation=:bvec, kwargs...)
    model = owner_model(builder)
    res = InstantGradientBlock(
        builder,
        orientation,
        @variable(model)
    )
    set_simple_constraints!(model, res, kwargs)
    return res
end


qval(instant::InstantGradientBlock) = instant.qval
bval(instant::InstantGradientBlock) = 0.
duration(instant::InstantGradientBlock) = 0.
properties(::Type{<:InstantGradientBlock}) = [qval]

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