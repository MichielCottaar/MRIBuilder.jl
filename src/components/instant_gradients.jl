module InstantGradients
import StaticArrays: SVector, SMatrix
import JuMP: @constraint
import ...Variables: @defvar, VariableType, variables, get_free_variable, set_simple_constraints!, make_generic, adjust_internal, adjustable, gradient_orientation
import ...BuildSequences: global_model
import ..AbstractTypes: EventComponent, GradientWaveform

"""
    InstantGradient(; orientation=nothing, group=nothing, variables...)

If the `orientation` is set an [`InstantGradient1D`](@ref) is returned, otherwise an [`InstantGradient3D`](@ref).

## Parameters
- `orientation` sets the gradient orientation as a length-3 vector. If not set, the gradient can be in any direction.
- `group`: name of the group to which this gradient belongs (used for scaling and rotating).

## Variables
- [`qval`](@ref): Spatial frequency on which spins will be dephased due to this pulsed gradient in rad/um (scalar if `orientation` is set and vector otherwise).
- [`spoiler_scale`](@ref): Length-scale on which spins will be dephased by exactly 2π in mm.
"""
abstract type InstantGradient{N} <: EventComponent end

function (::Type{InstantGradient})(; orientation=nothing, group=nothing, qval=nothing, variables...)
    if isnothing(orientation)
        res = InstantGradient3D(get_free_variable.(qval), group)
    else
        res = InstantGradient1D(get_free_variable(qval), orientation, group)
        if !(res.qval isa Number)
            @constraint global_model() res.qval >= 0
        end
    end
    set_simple_constraints!(res, variables)
    return res
end

"""
An [`InstantGradient`](@ref) with a fixed orientation.
"""
struct InstantGradient1D <: InstantGradient{1}
    qval :: VariableType
    orientation :: SVector{3, Float64}
    group :: Union{Nothing, Symbol}
end

@defvar gradient qvec(ig::InstantGradient1D) = ig.qval .* ig.orientation

"""
An [`InstantGradient`](@ref) with a variable orientation.
"""
struct InstantGradient3D <: InstantGradient{3}
    qvec :: SVector{3, VariableType}
    group :: Union{Nothing, Symbol}
end

function InstantGradient3D(; qval=[nothing, nothing, nothing], group=nothing, variables...)
    if isnothing(qval)
        qval = [nothing, nothing, nothing]
    end
    res = InstantGradient3D(get_free_variable.(qval), group)
    set_simple_constraints!(res, variables)
    return res
end

@defvar gradient qvec(ig::InstantGradient3D) = ig.qvec


variables.duration.f(::InstantGradient) = 0.
variables.effective_time.f(::InstantGradient) = 0.
@defvar gradient bmat_gradient(::InstantGradient, qstart=nothing) = zero(SMatrix{3, 3, Float64, 3})

make_generic(ig::InstantGradient) = ig

adjustable(::InstantGradient) = :gradient

gradient_orientation(ig::InstantGradient{1}) = ig.orientation

function adjust_internal(ig::InstantGradient1D; orientation=nothing, scale=1., rotation=nothing)
    if !isnothing(orientation) && !isnothing(rotation)
        error("Cannot set both the gradient orientation and rotation.")
    end
    new_orientation = isnothing(orientation) ? (isnothing(rotation) ? ig.orientation : rotation * ig.orientation) : orientation
    return InstantGradient1D(
        ig.qval * scale,
        new_orientation,
        ig.group
    )
end

function adjust_internal(ig::InstantGradient3D; scale=1., rotation=nothing)
    return InstantGradient3D(
        (
            isnothing(rotation) ? 
            (ig.qvec .* scale) : 
            (rotation * (ig.qvec .* scale))
        ),
        ig.group
    )
end

end