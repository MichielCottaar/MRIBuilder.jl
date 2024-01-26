module ConcreteBlocks
import JuMP: has_values, optimize!
import ..BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, properties, duration
import ..SequenceBuilders: SequenceBuilder, to_block, AbstractSequence

abstract type AbstractConcreteBlock <: BuildingBlock end

struct ConcreteRFPulse
    time :: Vector{Number}
    amplitude :: Vector{Number}
    phase :: Vector{Number}
end

function ConcreteRFPulse(arr::Vector)
    @assert all(length.(arr) .== 3)
    ConcreteRFPulse(
        [a[1] for a in arr],
        [a[2] for a in arr],
        [a[3] for a in arr],
    )
end

ConcreteRFPulse(::Nothing) = nothing
ConcreteRFPulse(values::Tuple{<:Vector, <:Vector, <:Vector}) = ConcreteRFPulse(values...)

struct ConcreteGradient
    time :: Vector{Number}
    Gx :: Vector{Number}
    Gy :: Vector{Number}
    Gz :: Vector{Number}
end

function ConcreteGradient(arr::Vector)
    @assert all(length.(arr) .== 4)
    ConcreteGradient(
        [a[1] for a in arr],
        [a[2] for a in arr],
        [a[3] for a in arr],
        [a[4] for a in arr],
    )
end

ConcreteGradient(::Nothing) = nothing
ConcreteGradient(values::Tuple{<:Vector, <:Vector, <:Vector, <:Vector}) = ConcreteGradient(values...)


"""
    ConcreteBlock(duration; pulse=nothing, gradient=nothing, rotate_bvec=false, readout_times=nothing)

A [`BuildingBlock`](@ref) that is fully defined (i.e., there are no variables to be optimised).
"""
struct ConcreteBlock <: AbstractConcreteBlock
    duration :: Float64
    pulse :: Union{ConcreteRFPulse, Nothing}
    gradient :: Union{ConcreteGradient, Nothing}
    readout_times :: Vector{Float64}
end

function ConcreteBlock(duration::Number; pulse=nothing, gradient=nothing, readout_times=Number[])
    ConcreteBlock(builder, duration, ConcreteRFPulse(pulse), ConcreteGradient(gradient), Float64.(readout_times))
end

to_block(::SequenceBuilder, concrete::AbstractConcreteBlock) = concrete
has_values(::AbstractConcreteBlock) = true
duration(c::AbstractConcreteBlock) = 0.
duration(c::ConcreteBlock) = c.duration


"""
    ConcreteBlock(other_building_block)

Creates a [`ConcreteBlock`](@ref) from another [`BuildingBlock`](@ref).

This will raise an error if the other [`BuildingBlock`](@ref) has not been optimised yet.
If it has been optimised, then [`to_concrete_block`](@ref) will be called.
"""
function ConcreteBlock(block::BuildingBlock)
    if !has_values(block)
        error("Making `BuildingBlock` objects concrete is only possible after optimisation.")
    end
    return to_concrete_block(block)
end


"""
    to_concrete_block(other_building_block)

Internal function used to create [`ConcreteBlock`](@ref) from any [`BuildingBlock`](@ref).

This needs to be defined for every [`BuildingBlock`](@ref)
"""
function to_concrete_block(cb::AbstractConcreteBlock)
    return cb
end

properties(::Type{<:ConcreteBlock}) = []

has_values(c::ConcreteBlock) = true


"""
    Sequence(builder::SequenceBuilder)

A fully defined sequence with no free variables.

When created from a [`SequenceBuilder`](@ref), all non-fixed variables are optimised given any constraints
and the resulting sequence is returned.
"""
struct Sequence <: AbstractSequence
    blocks :: Vector{<:AbstractConcreteBlock}
    TR :: Number
end

function Sequence(seq::SequenceBuilder)
    if !has_values(seq)
        optimize!(seq.model)
    end
    return Sequence(ConcreteBlock.(seq.blocks), value(seq.TR))
end

end