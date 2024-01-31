"""
Defines a set of different options for MRI gradients.
"""
module TrapezoidGradients

import JuMP: @constraint, @variable, Model, VariableRef, owner_model, value
import StaticArrays: SVector
import ...Variables: qvec, rise_time, flat_time, slew_rate, gradient_strength, variables, duration, δ, get_free_variable, VariableType
import ...BuildingBlocks: duration, set_simple_constraints!, fixed
import ...BuildSequences: @global_model_constructor
import ..Generic: GenericWaveform, SpecificWaveform
import ..ChangingGradientBlocks: ChangingGradientBlock
import ..ConstantGradientBlocks: ConstantGradientBlock


"""
    TrapezoidGradient(; orientation=:bvec, rotate=nothing, scale=nothing, pulse=nothing, variables...)

Defines a trapezoidal pulsed gradient

## Parameters
- `orientation` sets the gradient orienation (ignored if `qvec` is set). Can be set to a vector for a fixed orientation. Otherwise the orientation will be aligned with the `rotate` (if set) or fully free (if `rotate` is nothing). Set to :flip to point in the inverse of the user-provided `rotate`.
- `rotate`: with which user-set parameter will this gradient be rotated (e.g., :bvec). Default is no rotation.
- `scale`: with which user-set parameter will this gradient be scaled (e.g., :bval). Default is no scaling.
- `pulse`: RF pulse that will play during the flat part of the trapezoidal gradient. If you want to add dead time around the pulse, you can set it to a tuple (e.g., `pulse=(2, pulse, :min)`).

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
### Pulse variables
Any variables defined for the specific pulse added. Also:
- [`slice_thickness`](@ref): vector with the slice thickness

The [`bvalue`](@ref) can be constrained for multiple gradient pulses.
"""
struct TrapezoidGradient <: SpecificWaveform
    model :: Model
    rise_time :: VariableType
    flat_time :: VariableType
    slew_rate :: SVector{3, VariableType}
    rotate :: Union{Nothing, Symbol}
    scale :: Union{Nothing, Symbol}
    time_before_pulse :: VariableType
    pulse :: Union{Nothing, RFPulseBlock}
    time_after_pulse :: VariableType
end

@global_model_constructor TrapezoidGradient

function TrapezoidGradient(model::Model; orientation=nothing, rise_time=nothing, flat_time=nothing, rotate=nothing, scale=nothing, kwargs...)
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
    slew_rate = SVector{3}(slew_rate)
    rise_time = get_free_variable(model, rise_time)
    flat_time = get_free_variable(model, flat_time)

    res = TrapezoidGradient(
        model,
        rise_time,
        flat_time,
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

function GenericWaveform(pg::PulsedGradient)
    GenericWaveform(
        duration(pg),
        [
            ChangingGradientBlock(zeros(3), slew_rate(pg), rise_time(pg), pg.rotate, pg.scale),
            ConstantGradientBlock(gradient_strength(pg), flat_time(pg), pg.rotate, pg.scale),
            ChangingGradientBlock(gradient_strength(pg), -slew_rate(pg), rise_time(pg), pg.rotate, pg.scale),
        ],
        []
    )
end

rise_time(pg::TrapezoidGradient) = pg.rise_time
flat_time(pg::TrapezoidGradient) = pg.flat_time
gradient_strength(g::TrapezoidGradient) = slew_rate(g) .* rise_time(g)
slew_rate(g::TrapezoidGradient) = g.slew_rate
δ(g::TrapezoidGradient) = rise_time(g) + flat_time(g)
duration(g::TrapezoidGradient) = 2 * rise_time(g) + flat_time(g)
qvec(g::TrapezoidGradient) = δ(g) .* gradient_strength(g) .* 2π

get_children_indices(::TrapezoidGradient) = (:rise, :flat, :fall)
Base.getindex(pg::TrapezoidGradient, symbol::Symbol) = pg[Val(symbol)]
Base.getindex(pg::TrapezoidGradient, ::Val{:rise}) = pg.rise
Base.getindex(pg::TrapezoidGradient, ::Val{:flat}) = pg.flat
Base.getindex(pg::TrapezoidGradient, ::Val{:fall}) = pg.fall

start_time(pg::TrapezoidGradient, symbol::Symbol) = start_time(pg, Val(symbol))
start_time(pg::TrapezoidGradient, ::Val{:rise}) = 0.
start_time(pg::TrapezoidGradient, ::Val{:flat}) = rise_time(pg)
start_time(pg::TrapezoidGradient, ::Val{:fall}) = δ(pg)


variables(::Type{<:TrapezoidGradient}) = [qval, δ, gradient_strength, duration, rise_time, flat_time]


end