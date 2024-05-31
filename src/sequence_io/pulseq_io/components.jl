"""Define a list of all the [`PulseqShape`](@ref) and [`PulseqComponent`](@ref) that are used in a [`PulseqSequence`](@ref)."""
module Components

import ..Types: AnyPulseqComponent, PulseqShape, PulseqRFPulse, AnyPulseqGradient, PulseqADC, PulseqExtension, PulseqSequence

"""
    PulseqComponents(shapes, pulses, grads, adcs, extensions)

All the shapes, pulses, grads, adcs, and extensions used in a [`PulseqSequence`](@ref).

They can be provided as a dictionary from the integer ID to the object or as a vector.
"""
struct PulseqComponents
    shapes:: Dict{Int, PulseqShape}
    pulses:: Dict{Int, PulseqRFPulse}
    grads:: Dict{Int, AnyPulseqGradient}
    adc:: Dict{Int, PulseqADC}
    PulseqComponents(shapes, pulses, grads, adc) = new(
        _convert_to_dict(shapes, PulseqShape),
        _convert_to_dict(pulses, PulseqRFPulse),
        _convert_to_dict(grads, AnyPulseqGradient),
        _convert_to_dict(adc, PulseqADC),
    )
end

_convert_to_dict(d::Dict{Int, Any}, ::Type{T}) where {T} = Dict(Int, T)(d)
_convert_to_dict(d::Dict{Int, T}, ::Type{T}) where {T} = d
_convert_to_dict(vec::Vector, ::Type{T}) where {T} = Dict{Int, T}(i => v for (i, v) in enumerate(vec))


PulseqComponents() = PulseqComponents(
    PulseqShape[],
    PulseqRFPulse[],
    AnyPulseqGradient[],
    PulseqADC[],
)

"""
    add_components(comp::PulseqComponents, search_vec::Vector, component)
    add_components(comp::PulseqComponents, shape::PulseqShape)

Adds a component to the [`search_vec`](@ref), which is assumed to be the appropriate vector within the `comp`.

It will check whether the `component` is already part of the `search_vec` before adding it.
The integer ID of the position of the `component` in `search_vec` is returned.

0 is returned if `component` is `nothing`.
"""
add_components!(::PulseqComponents, search_vec::Dict{Int, <:Any}, ::Nothing) = 0
function add_components!(comp::PulseqComponents, search_vec::Dict{Int, <:T}, component::T) where {T <: AnyPulseqComponent}
    for (i, c) in search_vec
        if same_component(comp, c, component)
            return i
        end
    end
    search_vec[length(search_vec) + 1] = component
    return length(search_vec)
end

"""
    same_component(comp::PulseqComponents, a, b)

Check whether components `a` and `b` are the same when using the shapes represented in `comp`.
"""
same_component(::PulseqComponents, ::Any, ::Any) = false
function same_component(comp::PulseqComponents, a::T, b::T) where {T}
    for name in fieldnames(T)
        v1 = getfield(a, name)
        v2 = getfield(b, name)

        if v1 isa PulseqShape
            v1 = add_components!(comp, v1)
        end
        if v2 isa PulseqShape
            v2 = add_components!(comp, v2)
        end
        if v1 != v2
            return false
        end
    end
    return true
end


add_components!(comp::PulseqComponents, ::Nothing) = 0
function add_components!(comp::PulseqComponents, shape::PulseqShape)
    for (i, s) in comp.shapes
        if same_shape(s, shape)
            return i
        end
    end
    comp.shapes[length(comp.shapes) + 1] = shape
    return length(comp.shapes)
end

"""
    same_component(a, b)

Check whether shapes `a` and `b` are the same.
"""
same_shape(shape1::PulseqShape, shape2::PulseqShape) = length(shape1.samples) == length(shape2.samples) && all(shape1.samples .≈ shape2.samples)

end