module ConcreteBlocks
import JuMP: has_values, optimize!, value
import ..BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, properties, duration
import ..SequenceBuilders: SequenceBuilder, to_block, AbstractSequence, TR, get_blocks

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
ConcreteRFPulse(values::Tuple{<:AbstractVector, <:AbstractVector, <:AbstractVector}) = ConcreteRFPulse(values...)

struct ConcreteGradient
    time :: Vector{Number}
    Gx :: Vector{Number}
    Gy :: Vector{Number}
    Gz :: Vector{Number}
end

function ConcreteGradient(arr::Vector)
    if length(arr[1]) == 4
        @assert all(length.(arr) .== 4)
        ConcreteGradient(
            [a[1] for a in arr],
            [a[2] for a in arr],
            [a[3] for a in arr],
            [a[4] for a in arr],
        )
    elseif length(arr[1]) == 2
        @assert all(length.(arr) .== 2)
        ConcreteGradient(
            [a[1] for a in arr],
            [a[2] for a in arr],
        )
    else
        error()
    end
end

function ConcreteGradient(times, gradients)
    ConcreteGradient(
        times,
        [g[1] for g in gradients],
        [g[2] for g in gradients],
        [g[3] for g in gradients],
    )
end

ConcreteGradient(::Nothing) = nothing
ConcreteGradient(values::Tuple{<:Vector, <:Vector}) = ConcreteGradient(values...)
ConcreteGradient(values::Tuple{<:Vector, <:Vector, <:Vector, <:Vector}) = ConcreteGradient(values...)


"""
    ConcreteBlock(duration; pulse=nothing, gradient=nothing, rotate_bvec=false, readout_times=nothing)

A [`BuildingBlock`](@ref) that is fully defined (i.e., there are no variables to be optimised).
"""
struct ConcreteBlock <: AbstractConcreteBlock
    builder :: AbstractSequence
    duration :: Float64
    pulse :: Union{ConcreteRFPulse, Nothing}
    gradient :: Union{ConcreteGradient, Nothing}
    readout_times :: Vector{Float64}
    rotate_gradient :: Bool
end

ConcreteBlock(args...; kwargs...) = BuildingBlockPlaceholder{ConcreteBlock}(args...; kwargs...)

function ConcreteBlock(builder::AbstractSequence, duration::Number; pulse=nothing, gradient=nothing, readout_times=Number[], rotate_gradient=false)
    ConcreteBlock(builder, duration, ConcreteRFPulse(pulse), ConcreteGradient(gradient), Float64.(readout_times), rotate_gradient)
end


has_values(c::AbstractConcreteBlock) = has_values(c.builder)
duration(c::AbstractConcreteBlock) = 0.
duration(c::ConcreteBlock) = c.duration


"""
    ConcreteBlock(sequence, other_building_block)

Creates a [`ConcreteBlock`](@ref) from another [`BuildingBlock`](@ref).

This will raise an error if the other [`BuildingBlock`](@ref) has not been optimised yet.
If it has been optimised, then [`to_concrete_block`](@ref) will be called.
"""
function ConcreteBlock(sequence::AbstractSequence, block::BuildingBlock)
    if !has_values(block)
        error("Making `BuildingBlock` objects concrete is only possible after optimisation.")
    end
    return to_concrete_block(sequence, block)
end


"""
    to_concrete_block(other_building_block)

Internal function used to create [`ConcreteBlock`](@ref) from any [`BuildingBlock`](@ref).

This needs to be defined for every [`BuildingBlock`](@ref)
"""
function to_concrete_block(sequence::AbstractSequence, cb::ConcreteBlock)
    return ConcreteBlock(sequence, cb.duration, cb.pulse, cb.gradient, cb.readout_times)
end

properties(::Type{<:AbstractConcreteBlock}) = []


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

TR(seq::Sequence) = seq.TR
get_blocks(seq::Sequence) = seq.blocks
has_values(seq::Sequence) = true

function Sequence(builder::SequenceBuilder)
    if !has_values(builder)
        optimize!(builder.model)
    end
    seq = Sequence(AbstractConcreteBlock[], value(TR(builder)))
    for block in builder.blocks
        @show block
        push!(seq.blocks, ConcreteBlock(seq, block))
    end
    return seq
end

end