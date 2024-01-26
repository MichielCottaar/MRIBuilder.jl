module ConcreteBlocks
import JuMP: has_values, optimize!, value
import ..Variables: variables, duration
import ..BuildingBlocks: BuildingBlock


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


abstract type AbstractConcreteBlock <: BuildingBlock end

"""
    ConcreteBlock(duration; pulse=nothing, gradient=nothing, rotate_bvec=false, readout_times=nothing)

A [`BuildingBlock`](@ref) that is fully defined (i.e., there are no variables to be optimised).
"""
struct ConcreteBlock <: AbstractConcreteBlock
    duration :: Float64
    pulse :: Union{ConcreteRFPulse, Nothing}
    gradient :: Union{ConcreteGradient, Nothing}
    readout_times :: Vector{Float64}
    rotate_gradient :: Bool
end

function ConcreteBlock(duration::Number; pulse=nothing, gradient=nothing, readout_times=Number[], rotate_gradient=false)
    ConcreteBlock(duration, ConcreteRFPulse(pulse), ConcreteGradient(gradient), Float64.(readout_times), rotate_gradient)
end


has_values(c::AbstractConcreteBlock) = true
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
function to_concrete_block(cb::ConcreteBlock)
    return ConcreteBlock(cb.duration, cb.pulse, cb.gradient, cb.readout_times, cb.rotate_gradient, value.(cb.start_time))
end

variables(::Type{<:AbstractConcreteBlock}) = []


end