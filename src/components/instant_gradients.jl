module InstantGradients
import StaticArrays: SVector, SMatrix
import ...Variables: VariableType, duration, qval, bmat_gradient, get_free_variable, set_simple_constraints!, effective_time, make_generic
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
    end
    set_simple_constraints!(res, variables)
    return res
end

"""
An [`InstantGradient`](@ref) with a fixed orientation.
"""
struct InstantGradient1D <: InstantGradient{1}
    qval :: VariableType
    orientation :: SVector{3, Number}
    group :: Union{Nothing, Symbol}
end

qval(ig::InstantGradient1D) = ig.qval

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

qval(ig::InstantGradient3D) = ig.qvec


duration(::InstantGradient) = 0.
effective_time(::InstantGradient) = 0.
bmat_gradient(::InstantGradient, qstart=nothing) = zero(SMatrix{3, 3, Float64, 3})

make_generic(ig::InstantGradient) = ig

end