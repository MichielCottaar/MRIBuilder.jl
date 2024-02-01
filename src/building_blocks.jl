module BuildingBlocks
import JuMP: has_values, value, Model, @constraint, @objective, owner_model, objective_function, optimize!, AbstractJuMPScalar
import Printf: @sprintf
import ..Variables: variables, start_time, duration, end_time, gradient_strength, slew_rate, effective_time, VariableType, qval_square

"""
Parent type for all individual components out of which a sequence can be built.

Required methods:
- [`duration`](@ref)(block, parameters): Return block duration in ms.
- [`fixed`](block): Return the equivalent fixed BuildingBlock (i.e., `FixedBlock`, `FixedPulse`, `FixedGradient`, `FixedInstantPulse`, `FixedInstantGradient`, or `InstantReadout`). These all have in common that they have no free variables and explicitly set any gradient and RF pulse profiles.
- [`variables`](@ref): A list of all functions that are used to compute variables of the building block. Any of these can be used in constraints or objective functions.
"""
abstract type BuildingBlock end

"""
Parent type for all RF pulses.

RF pulses combined with gradients, should be childrent of [`ContainerBlock`](@ref) instead.

Required methods:
- [`effective_time`](@ref)(pulse): Best approximation of time the RF pulse is applied. This is defined relative to the start of the pulse.
"""
abstract type RFPulseBlock <: BuildingBlock end

"""
Parent type for all gradient profiles.
"""
abstract type GradientBlock <: BuildingBlock end

"""
Parent type for all types combining one or more pulses/gradients.

Required methods:
- [`get_children_blocks`](@ref)(container): return all the [`BuildingBlock`](@ref) objects includes in this container with their indices.
- [`start_time`](@ref)(container, index): returns the starting time of the child corresponding to `index` relative to the start of the `container` in ms.
- `Base.getindex`(container, index): get child [`BuildingBlock`](@ref) corresponding to `index`.
"""
abstract type ContainerBlock <: BuildingBlock end


"""
    get_children_blocks(container)

Return all the [`BuildingBlock`](@ref) objects includes in this container with their indices.
"""
get_children_blocks(bb::BuildingBlock) = [(i, bb[i]) for i in get_children_indices(bb)]

"""
    get_children_indices(container)

Return the indices of all the children in a [`ContainerBlock`](@ref).

This needs to be defined for every [`ContainerBlock`](@ref).
It is not part of the external API, but is used by [`get_children_blocks`](@ref)
"""
function get_children_indices end


"""
    start_time(container, args...)

Returns the starting time of the specific [`BuildingBlock`](@ref) within the container.
The [`BuildingBlock`](@ref) is defined by one or more indices as defined below.
"""
start_time(bb::BuildingBlock) = 0.
start_time(container::ContainerBlock, index1, index2, more_indices...) = start_time(container, index1) + start_time(container[index1], index2, more_indices)

"""
    effective_time(pulse)
    effective_time(readout)
    effective_time(container, indices...)

Returns the effective time of a pulse or readout.

For a pulse, this means the timepoint at which one would place an [`InstantRFPulseBlock`](@ref) if one would want to have a similar effect.

For a reaodut, this is the time the readout passes through the zero-point in k-space (or the minimum in k-space if it does not go through zero).

The time is given with respect to the start of the pulse or readout, or to the start of a container if the pulse/readout is identified using indices.
"""
effective_time(bb::ContainerBlock, index, indices...) = start_time(bb, index) + effective_time(bb[index], indices...)

"""
    end_time(container, args...)

Returns the end time of the specific [`BuildingBlock`](@ref) within the container.
The [`BuildingBlock`](@ref) is defined by one or more indices as defined below.
"""
end_time(bb::BuildingBlock) = duration(bb::BuildingBlock)
end_time(container::ContainerBlock, index1, indices...) = start_time(container, index1) + end_time(container[index1], indices...)


"""
    to_block(object)

Function used internally to convert a wide variety of objects into [`BuildingBlock`](@ref) objects.
"""
to_block(bb::BuildingBlock) = bb


"""
    fixed(block::BuildingBlock)

Return the fixed equivalent of the `BuildingBlock`

Possible return types are `FixedSequence`, `FixedBlock`, `FixedPulse`, `FixedGradient`, `FixedInstantPulse`, `FixedInstantGradient`, or `InstantReadout`. 
These all have in common that they have no free variables and explicitly set any gradient and RF pulse profiles.
"""
function fixed end



"""
    variables(building_block)

Returns a list of function that can be called to constrain the `building_block`.
"""
variables(bb::BuildingBlock) = variables(typeof(bb))


struct BuildingBlockPrinter{T<:BuildingBlock}
    bb::T
    start_time::Union{Nothing, Number}
    spaces::Int
end



Base.show(io::IO, block::BuildingBlock) = print(io, BuildingBlockPrinter(block, nothing, 0))

function _robust_value(possible_number::VariableType)
    try
        return value(possible_number)
    catch
        return nothing
    end
end

function _robust_value(possible_vector::AbstractVector) 
    result = _robust_value.(possible_vector)
    if any(isnothing.(result))
        return nothing
    end
    return result
end

function Base.show(io::IO, printer::BuildingBlockPrinter)
    block = printer.bb
    print(io, string(typeof(block)), "(")
    variable_names = nameof.(variables(block))
    printed_duration = false
    if !isnothing(printer.start_time)
        print(io, "t=", @sprintf("%.3g", printer.start_time))

        dur = _robust_value(duration(block))
        @assert !(dur isa AbstractVector)
        if !isnothing(dur) && !iszero(dur)
            print(io, "-", @sprintf("%.3g", printer.start_time + dur))
            printed_duration = true
        end
        print(io, ", ")
    end
    for name in propertynames(block)
        ft = fieldtype(typeof(block), name)
        if (
            ft in (VariableType, Model) ||
            (ft <: AbstractVector && eltype(ft) == VariableType) ||
            string(name)[1] == '_'
        )
            continue
        end

        print(io, name, "=", repr(getproperty(block, name)), ", ")
    end

    for fn in variables(block)
        if printed_duration && fn == duration
            continue
        end
        numeric_value = _robust_value(fn(block))
        if isnothing(numeric_value)
            continue
        end
        if numeric_value isa AbstractVector
            printed_value = "[" * join(map(v -> @sprintf("%.3g", v), numeric_value), ", ") * "]"
        else
            printed_value = @sprintf("%.3g", numeric_value)
        end
        print(io, "$(nameof(fn))=$(printed_value), ")
    end
    print(io, ")")
end


"""
    set_simple_constraints!(model, block, kwargs)

Add any constraints or objective functions to the variables of a [`BuildingBlock`](@ref) in the JuMP `model`.

Each keyword argument has to match one of the functions in [`variables`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_simple_constraints!(model::Model, block::BuildingBlock, kwargs)
    to_funcs = Dict(nameof(fn) => fn for fn in variables(block))

    invert_value(value::VariableType) = 1 / value
    invert_value(value::Symbol) = invert_value(Val(value))
    invert_value(::Val{:min}) = Val(:max)
    invert_value(::Val{:max}) = Val(:min)
    invert_value(value::AbstractVector) = invert_value.(value)
    invert_value(value) = value

    for (key, value) in kwargs
        if key in keys(to_funcs)
            apply_simple_constraint!(model, to_funcs[key](block), value)
        else
            if key == :qval
                apply_simple_constraint!(model, to_funcs[:qval_square](block), value isa VariableType ? value^2 : value)
            elseif key == :slice_thickness && :inverse_slice_thickness in keys(to_funcs)
                apply_simple_constraint!(model, to_funcs[:inverse_slice_thickness](block), invert_value(value))
            elseif key == :bandwidth && :inverse_bandwidth in keys(to_funcs)
                apply_simple_constraint!(model, to_funcs[:inverse_bandwidth](block), invert_value(value))
            else
                error("Trying to set an unrecognised variable $key.")
            end
        end
    end
    nothing
end

"""
    apply_simple_constraint!(model, variable, value)

Add a single constraint or objective to the JuMP `model`.
This is an internal function used by [`set_simple_constraints`](@ref).
"""
apply_simple_constraint!(model::Model, variable, ::Nothing) = nothing
apply_simple_constraint!(model::Model, variable, value::Symbol) = apply_simple_constraint!(model, variable, Val(value))
apply_simple_constraint!(model::Model, variable, ::Val{:min}) = @objective model Min objective_function(model) + variable
apply_simple_constraint!(model::Model, variable, ::Val{:max}) = @objective model Min objective_function(model) - variable
apply_simple_constraint!(model::Model, variable, value::VariableType) = @constraint model variable == value
apply_simple_constraint!(model::Model, variable::AbstractVector, value::AbstractVector) = [apply_simple_constraint!(model, v1, v2) for (v1, v2) in zip(variable, value)]


"""
    match_blocks!(block1, block2[, property_list])

Matches the listed variables between two [`BuildingBlock`](@ref) objects.
By default all shared variables (i.e., those with the same name) are matched.
"""
function match_blocks!(block1::BuildingBlock, block2::BuildingBlock, property_list)
    model = owner_model(block1)
    @assert model == owner_model(block2)
    for fn in property_list
        @constraint model fn(block1) == fn(block2)
    end
end

function match_blocks!(block1::BuildingBlock, block2::BuildingBlock)
    property_list = intersect(variables(block1), variables(block2))
    match_blocks!(block1, block2, property_list)
end


optimize!(bb::BuildingBlock) = optimize!(owner_model(bb))
function owner_model(bb::BuildingBlock)
    if hasproperty(bb, :model)
        return bb.model
    else
        for name in propertynames(bb)
            value = getproperty(bb, name)
            if value isa AbstractJuMPScalar
                return owner_model(value)
            end
        end
    end
    error("Cannot find owner model")
end

function has_values(bb::BuildingBlock) 
    try 
        return has_values(owner_model(bb))
    catch
        # return true for building blocks without a model
        return true
    end
end

end