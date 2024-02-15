module InstantGradients
import StaticArrays: SVector, SMatrix
import ...Variables: VariableType, duration, qval, bmat_gradient, get_free_variable, set_simple_constraints!, qval, effective_time, make_generic
import ..AbstractTypes: EventComponent, GradientWaveform

"""
    InstantGradient1D(; orientation=[1, 0, 0], group=nothing, variables...)

Defines an instantaneous gradient with given fixed `orientation` (default: x-direction).

To have a variable `orientation`, see [`InstantGradient3D`](@ref).

## Parameters
- `orientation` sets the gradient orienation as a length-3 vector. If not set, the gradient can be in any direction.
- `group`: name of the group to which this gradient belongs (used for scaling and rotating).

## Variables
- [`qval`](@ref): Spatial frequency on which spins will be dephased due to this pulsed gradient in rad/um (scalar if `orientation` is set and vector otherwise).
- [`spoiler_scale`](@ref): Length-scale on which spins will be dephased by exactly 2π in mm.
"""
abstract type InstantGradient <: EventComponent end

"""
    InstantGradient1D(; orientation=[1, 0, 0], group=nothing, variables...)

Defines an instantaneous gradient with given fixed `orientation` (default: x-direction).

To have a variable `orientation`, see [`InstantGradient3D`](@ref).

## Parameters
- `orientation` sets the gradient orienation as a length-3 vector (default: x-direction).
- `group`: name of the group to which this gradient belongs (used for scaling and rotating).

## Variables
- [`qval`](@ref): Spatial frequency on which spins will be dephased due to this pulsed gradient in rad/um.
- [`spoiler_scale`](@ref): Length-scale on which spins will be dephased by exactly 2π in mm.
"""
struct InstantGradient1D <: InstantGradient
    qval :: VariableType
    orientation :: SVector{3, Number}
    group :: Union{Nothing, Symbol}
end

function InstantGradient1D(; orientation=[1, 0, 0], group=nothing, qval=nothing, variables...)
    res = InstantGradient1D(qval, orientation, group)
    set_simple_constraints!(res, variables)
    return res
end

qval(ig::InstantGradient1D) = ig.qval


"""
    InstantGradient3D(; group=nothing, variables...)

Defines an instantaneous gradient without a fixed orientation.

To have a fixed `orientation`, see [`InstantGradient1D`](@ref).

## Parameters
- `group`: name of the group to which this gradient belongs (used for scaling and rotating).

## Variables
- [`qval`](@ref): Vector of spatial frequency on which spins will be dephased due to this pulsed gradient in rad/um.
- [`spoiler_scale`](@ref): Vector with length-scale on which spins will be dephased by exactly 2π in mm.
"""
struct InstantGradient3D <: InstantGradient
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

qval(ig::InstantGradient3D) = ig.qvec


duration(::InstantGradient) = 0.
effective_time(::InstantGradient) = 0.
bmat_gradient(::InstantGradient, qstart=nothing) = zero(SMatrix{3, 3, Float64, 3})

make_generic(ig::InstantGradient) = ig

end