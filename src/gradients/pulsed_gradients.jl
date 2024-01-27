"""
Defines a set of different options for MRI gradients.
"""
module PulsedGradients

import JuMP: @constraint, @variable, Model, VariableRef, owner_model, value
import StaticArrays: SVector
import ...Variables: qval, bval, rise_time, flat_time, slew_rate, gradient_strength, variables, duration, δ, get_free_variable, VariableType
import ...BuildingBlocks: BuildingBlock, duration, set_simple_constraints!
import ...ConcreteBlocks: ConcreteBlock, to_concrete_block
import ...BuildSequences: @global_model_constructor


"""
    PulsedGradient(; orientation=:bvec, variables...)

Defines a trapezoidal pulsed gradient

This is a [`BuildingBlock`](@ref) for a [`SequenceBuilder`](@ref).

## Parameters
- `orientation` sets the gradient orienation. Can be set to a vector for a fixed orientation. Alternatively, can be set to :bvec (default) to rotate with the user-provided `bvecs` or to :neg_bvec to always be the reverse of the `bvecs`.

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
- [`qval`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um (given by `δ` * `gradient_strength`).

The [`bvalue`](@ref) can be constrained for multiple gradient pulses.
"""
mutable struct PulsedGradient <: BuildingBlock
    model :: Model
    orientation :: Any
    slew_rate :: VariableType
    rise_time :: VariableType
    flat_time :: VariableType
end

@global_model_constructor PulsedGradient

function PulsedGradient(model::Model; orientation=:bvec, slew_rate=nothing, rise_time=nothing, flat_time=nothing, kwargs...)
    res = PulsedGradient(
        model,
        orientation,
        [get_free_variable(model, value) for value in (slew_rate, rise_time, flat_time)]...
    )

    set_simple_constraints!(model, res, kwargs)
    @constraint model res.flat_time >= 0
    @constraint model res.rise_time >= 0
    @constraint model res.slew_rate >= 0
    return res
end

rise_time(pg::PulsedGradient) = pg.rise_time
flat_time(pg::PulsedGradient) = pg.flat_time
gradient_strength(g::PulsedGradient) = rise_time(g) * slew_rate(g)
slew_rate(g::PulsedGradient) = g.slew_rate
δ(g::PulsedGradient) = rise_time(g) + flat_time(g)
duration(g::PulsedGradient) = 2 * rise_time(g) + flat_time(g)
qval(g::PulsedGradient) = (g.orientation == :neg_bvec ? -1 : 1) * gradient_strength(g) * δ(g)


function bval(g::PulsedGradient, qstart=0.)
    tr = rise_time(g)
    td = δ(g)
    grad = gradient_strength(g)
    return (
        # b-value due to just the gradient
        grad * (1//60 * tr^3 - 1//12 * tr^2 * td + 1//2 * tr * td^2 + 1//3 * td^3) + 
        # b-value due to cross-term
        qstart * grad * (td * (td + tr)) +
        # b-value due to just `qstart`
        (td + tr) * qstart^2
    )
end

variables(::Type{<:PulsedGradient}) = [qval, δ, gradient_strength, duration, rise_time, flat_time, slew_rate]


function to_concrete_block(block::PulsedGradient)
    if block.orientation == :bvec
        rotate = true
        qvec = [value(qval(block)), 0., 0.]
    elseif block.orientation == :neg_bvec
        rotate = true
        qvec = [-value(qval(block)), 0., 0.]
    elseif block.orientation isa AbstractVector && size(block.orientation) == (3, )
        rotate = false
        qvec = block.orientation .* (value(qval(block)) / norm(block.orientation))
    else
        error("Gradient orientation should be :bvec, :neg_bvec or a length-3 vector, not $(block.orienation)")
    end
    t_rise = value(rise_time(block))
    t_d = value(δ(block))
    return ConcreteBlock(t_d + t_rise, gradient=[
        (0., zeros(3)),
        (t_rise, qvec),
        (t_d, qvec),
        (t_d + t_rise, zeros(3)),
    ], rotate_gradient=rotate)
end


end