"""
Defines a set of different options for MRI gradients.
"""
module PulsedGradients

import JuMP: @constraint, @variable, Model, VariableRef, owner_model, value
import StaticArrays: SVector
import ....Sequences: MRGradients
import ...BuildingBlocks: BuildingBlock, duration, properties, set_simple_constraints!, BuildingBlockPlaceholder, gradient_strength, slew_rate, to_mcmr_components
import ...SequenceBuilders: SequenceBuilder, start_time
import ..IntegrateGradients: qval, bval


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
    builder::SequenceBuilder
    orientation :: Any
    slew_rate :: VariableRef
    rise_time :: VariableRef
    flat_time :: VariableRef
end

function PulsedGradient(; kwargs...)
    return BuildingBlockPlaceholder{PulsedGradient}(; kwargs...)
end

function PulsedGradient(builder::SequenceBuilder; orientation=:bvec, kwargs...)
    model = owner_model(builder)
    res = PulsedGradient(
        builder,
        orientation,
        @variable(model),
        @variable(model),
        @variable(model)
    )

    set_simple_constraints!(model, res, kwargs)
    @constraint model flat_time(res) >= 0
    @constraint model rise_time(res) >= 0
    @constraint model slew_rate(res) >= 0
    return res
end

"""
    rise_time(pulsed_gradient)

The time from 0 till the maximum gradient strength in ms.
"""
rise_time(pg::PulsedGradient) = pg.rise_time

"""
    flat_time(pulsed_gradient)

The time spent at the maximum gradient strength in ms.
"""
flat_time(pg::PulsedGradient) = pg.flat_time

gradient_strength(g::PulsedGradient) = rise_time(g) * slew_rate(g)

slew_rate(g::PulsedGradient) = g.slew_rate

"""
    δ(pulsed_gradient)

Pulse gradient duration (`rise_time + `flat_time`).  This is the effective duration of the gradient. The real duration is longer (and given by [`duration`](@ref)).
"""
δ(g::PulsedGradient) = rise_time(g) + flat_time(g)

duration(g::PulsedGradient) = 2 * rise_time(g) + flat_time(g)

"""
    qval(gradient)

q-value at the end of the gradient (rad/um).
"""
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

properties(::Type{<:PulsedGradient}) = [qval, δ, gradient_strength, duration, rise_time, flat_time, slew_rate]


function to_mcmr_components(block::PulsedGradient)
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
    t_start = value(start_time(block))
    t_rise = value(rise_time(block))
    t_d = value(δ(block))
    return MRGradients([
        (t_start, zeros(3)),
        (t_start + t_rise, qvec),
        (t_start + t_d, qvec),
        (t_start + t_d + t_rise, zeros(3)),
    ])
end


end