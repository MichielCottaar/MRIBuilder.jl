"""
Defines a set of different options for MRI gradients.
"""
module PulsedGradients

import JuMP: @constraint, @variable, Model, VariableRef, owner_model, value
import StaticArrays: SVector
import ...Variables: qval, qvec, rise_time, flat_time, slew_rate, gradient_strength, variables, duration, δ, get_free_variable, VariableType
import ...BuildingBlocks: ContainerBlock, duration, set_simple_constraints!, fixed, start_time, get_children_indices
import ...BuildSequences: @global_model_constructor
import ..FixedGradients: FixedGradient
import ..ChangingGradientBlocks: ChangingGradientBlock
import ..ConstantGradientBlocks: ConstantGradientBlock


"""
    PulsedGradient(; orientation=:bvec, variables...)

Defines a trapezoidal pulsed gradient

## Parameters
- `orientation` sets the gradient orienation (ignored if `qvec` is set). Can be set to a vector for a fixed orientation. Otherwise the orientation will be aligned with the `rotate` (if set) or fully free (if `rotate` is nothing). Set to :flip to point in the inverse of the user-provided `rotate`.
- `rotate`: with which user-set parameter will this gradient be rotated (e.g., :bvec). Default is no rotation.
- `scale`: with which user-set parameter will this gradient be scaled (e.g., :bval). Default is no scaling.

## Variables
Variables can be set during construction or afterwards as an attribute.
If not set, they will be determined during the sequence optimisation.
### Timing variables
- [`rise_time`](@ref): Time of the gradient to reach from 0 to maximum in ms. If explicitly set to 0, the scanner slew rate will be ignored.
- [`flat_time`](@ref): Time that the gradient stays at maximum strength in ms.
- [`δ`](@ref): effective pulse duration (`rise_time` + `flat_time`) in ms.
- [`duration`](@ref): total pulse duration (2 * `rise_time` + `flat_time`) in ms.
### Gradient variables
- [`gradient_strength`](@ref): Maximum gradient strength achieved during the pulse in kHz/um
- [`qval`](@ref)/[`qvec`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um (given by `δ` * `gradient_strength`).

The [`bvalue`](@ref) can be constrained for multiple gradient pulses.
"""
struct PulsedGradient <: ContainerBlock
    model :: Model
    rise :: ChangingGradientBlock
    flat :: ConstantGradientBlock
    fall :: ChangingGradientBlock
    slew_rate :: SVector{3, VariableType}
    _scaling :: Union{Nothing, VariableType}
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
end

@global_model_constructor PulsedGradient

function PulsedGradient(model::Model; orientation=nothing, rise_time=nothing, flat_time=nothing, rotate=nothing, scale=nothing, kwargs...)
    if isnothing(orientation) && isnothing(rotate)
        rate_1d = nothing
        slew_rate = (
            get_free_variable(model, nothing),
            get_free_variable(model, nothing),
            get_free_variable(model, nothing),
        )
    else
        rate_1d = get_free_variable(model, nothing)
        @constraint model rate_1d >= 0
        if isnothing(orientation)
            rate_1d = get_free_variable(model, nothing)
            slew_rate = (
                rate_1d,
                0.,
                0.,
            )
        elseif orientation == :flip
            @assert !isnothing(rotate) "setting `orientation=:flip` only makes sense if the `rotate` is set as well."
            slew_rate = (
                -rate_1d,
                0.,
                0.,
            )
        else
            slew_rate = rate_1d .* (orientation ./ norm(orientation))
        end
    end
    rise_time = get_free_variable(model, rise_time)
    flat_time = get_free_variable(model, flat_time)
    grad = slew_rate .* rise_time

    res = PulsedGradient(
        model,
        ChangingGradientBlock(zeros(3), slew_rate, rise_time, rotate, scale),
        ConstantGradientBlock(grad, flat_time, rotate, scale),
        ChangingGradientBlock(grad, -slew_rate, rise_time, rotate, scale),
        slew_rate,
        rate_1d,
        rotate,
        scale
    )

    set_simple_constraints!(model, res, kwargs)
    @constraint model flat_time >= 0
    @constraint model rise_time >= 0
    return res
end

rise_time(pg::PulsedGradient) = duration(pg.rise)
flat_time(pg::PulsedGradient) = duration(pg.flat)
gradient_strength(g::PulsedGradient) = gradient_strength(g.flat)
slew_rate(g::PulsedGradient) = g.slew_rate
δ(g::PulsedGradient) = rise_time(g) + flat_time(g)
duration(g::PulsedGradient) = 2 * rise_time(g) + flat_time(g)
qvec(g::PulsedGradient) = δ(g) .* gradient_strength(g)

get_children_indices(::PulsedGradient) = (:rise, :flat, :fall)
Base.getindex(pg::PulsedGradient, symbol::Symbol) = pg[Val(symbol)]
Base.getindex(pg::PulsedGradient, ::Val{:rise}) = pg.rise
Base.getindex(pg::PulsedGradient, ::Val{:flat}) = pg.flat
Base.getindex(pg::PulsedGradient, ::Val{:fall}) = pg.fall

start_time(pg::PulsedGradient, symbol::Symbol) = start_time(pg, Val(symbol))
start_time(pg::PulsedGradient, ::Val{:rise}) = 0.
start_time(pg::PulsedGradient, ::Val{:flat}) = rise_time(pg)
start_time(pg::PulsedGradient, ::Val{:fall}) = δ(pg)


variables(::Type{<:PulsedGradient}) = [qval, δ, gradient_strength, duration, rise_time, flat_time]


function fixed(block::PulsedGradient)
    grad = value.(gradient_strength(block))
    t_rise = value(rise_time(block))
    t_d = value(δ(block))
    return FixedGradient(
        [0., t_rise, t_d, td + t_rise],
        [zeros(3), grad, grad, zeros(3)]; 
        rotate=rotate
    )
end


end