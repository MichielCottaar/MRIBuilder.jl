module ConcreteBlocks
import JuMP: has_values
import ..BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, properties
import ..SequenceBuilders: SequenceBuilder

abstract type AbstractConcreteBlock end

"""
    ConcreteBlock(duration; pulse=nothing, gradient=nothing, rotate_bvec=false, readout_times=nothing)

A [`BuildingBlock`](@ref) that is fully defined (i.e., there are no variables to be optimised).
"""
struct ConcreteBlock <: AbstractConcreteBlock
    builder :: SequenceBuilder
    duration :: Float64
    pulse :: ConcreteShape{Float64}
    gradient :: ConcreteShape{SVector{3, Float64}}
    readout_times :: Vector{Float64}
end

ConcreteBlock(args...; kwargs...) = BuildingBlockPlaceholder{ConcreteBlock}(args...; kwargs...)

function ConcreteBlock(builder::SequenceBuilder, duration::Number; pulse=nothing, gradient=nothing, readout_times=Number[])
    ConcreteBlock(builder, duration, ConcreteShape{Float64}(pulse), ConcreteShape(SVector{3, Float64})(gradient), Float64.(readout_times))
end


"""
    ConcreteBlock(other_building_block)

Creates a [`ConcreteBlock`](@ref) from another [`BuildingBlock`](@ref).

This will raise an error if the other [`BuildingBlock`](@ref) has not been optimised yet.
If it has been optimised, then [`to_concrete_block`](@ref) will be called.
"""
function ConcreteBlock(builder::SequenceBuilder, block::BuildingBlock)
    if !has_values(block)
        error("Making `BuildingBlock` objects concrete is only possible after optimisation.")
    end
    return to_concrete_block(builder, block)
end


"""
    to_concrete_block(builder, other_building_block)

Internal function used to create [`ConcreteBlock`](@ref) from any [`BuildingBlock`](@ref).

This needs to be defined for every [`BuildingBlock`](@ref)
"""
function to_concrete_block(builder::SequenceBuilder, cb::ConcreteBlock)
    return ConcreteBlock(builder, cb.duration, cb.pulse, cb.gradient, cb.readout_times)
end

properties(::Type{<:ConcreteBlock}) = []

has_values(c::ConcreteBlock) = true

end