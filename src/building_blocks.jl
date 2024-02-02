module BuildingBlocks
import JuMP: value, Model, @constraint, @objective, objective_function, AbstractJuMPScalar
import Printf: @sprintf
import ..Variables: variables, start_time, duration, end_time, gradient_strength, slew_rate, effective_time, VariableType, qval_square
import ..BuildSequences: global_model, global_scanner, fixed
import ..Scanners: Scanner

"""
Parent type for all individual components out of which a sequence can be built.

Required methods:
- [`duration`](@ref)(block, parameters): Return block duration in ms.
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

_robust_value(possible_tuple::Tuple) = _robust_value([possible_tuple...])


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
            ft == VariableType ||
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
    set_simple_constraints!(block, kwargs)

Add any constraints or objective functions to the variables of a [`BuildingBlock`](@ref).

Each keyword argument has to match one of the functions in [`variables`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_simple_constraints!(block::BuildingBlock, kwargs)
    to_funcs = Dict(nameof(fn) => fn for fn in variables(block))

    invert_value(value::VariableType) = 1 / value
    invert_value(value::Symbol) = invert_value(Val(value))
    invert_value(::Val{:min}) = Val(:max)
    invert_value(::Val{:max}) = Val(:min)
    invert_value(value::AbstractVector) = invert_value.(value)
    invert_value(value) = value

    for (key, value) in kwargs
        if key in keys(to_funcs)
            apply_simple_constraint!(to_funcs[key](block), value)
        else
            if key == :qval
                apply_simple_constraint!(to_funcs[:qval_square](block), value isa VariableType ? value^2 : value)
            elseif key == :slice_thickness && :inverse_slice_thickness in keys(to_funcs)
                apply_simple_constraint!(to_funcs[:inverse_slice_thickness](block), invert_value(value))
            elseif key == :bandwidth && :inverse_bandwidth in keys(to_funcs)
                apply_simple_constraint!(to_funcs[:inverse_bandwidth](block), invert_value(value))
            else
                error("Trying to set an unrecognised variable $key.")
            end
        end
    end
    nothing
end

"""
    apply_simple_constraint!(variable, value)

Add a single constraint or objective to the `variable`.

`value` can be one of:
- `nothing`: do nothing
- `:min`: minimise the variable
- `:max`: maximise the variable
- `number`: fix variable to this value
- `equation`: fix variable to the result of this equation
"""
apply_simple_constraint!(variable, ::Nothing) = nothing
apply_simple_constraint!(variable, value::Symbol) = apply_simple_constraint!(variable, Val(value))
apply_simple_constraint!(variable, ::Val{:min}) = @objective global_model() Min objective_function(global_model()) + variable
apply_simple_constraint!(variable, ::Val{:max}) = @objective global_model() Min objective_function(global_model()) - variable
apply_simple_constraint!(variable, value::VariableType) = @constraint global_model() variable == value
apply_simple_constraint!(variable::AbstractVector, value::AbstractVector) = [apply_simple_constraint!(v1, v2) for (v1, v2) in zip(variable, value)]
apply_simple_constraint!(variable::Number, value::Number) = @assert variable ≈ value "Variable set to multiple incompatible values."


"""
    match_blocks!(block1, block2[, property_list])

Matches the listed variables between two [`BuildingBlock`](@ref) objects.
By default all shared variables (i.e., those with the same name) are matched.
"""
function match_blocks!(block1::BuildingBlock, block2::BuildingBlock, property_list)
    for fn in property_list
        @constraint global_model() fn(block1) == fn(block2)
    end
end

function match_blocks!(block1::BuildingBlock, block2::BuildingBlock)
    property_list = intersect(variables(block1), variables(block2))
    match_blocks!(block1, block2, property_list)
end

"""
    scanner_constraints!(building_block[, scanner])

Adds the gradient strength and slew rate constraints from a specific [`Scanner`](@ref) to a [`BuildingBlock`]{@ref}.

This is applied iteratively to each part of a `Sequence`.
"""
scanner_constraints!(building_block::BuildingBlock) = scanner_constraints!(building_block, global_scanner())

function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner)
    for func in [gradient_strength, slew_rate]
        if isfinite(func(scanner))
            scanner_constraints!(building_block, scanner, func)
        end
    end
end

function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner, func::Function)
    model = global_model()
    if func in variables(building_block)
        # apply constraint at this level
        res_bb = func(building_block)
        if res_bb isa AbstractVector
            if isnothing(building_block.rotate)
                # no rotation; apply constraint to each dimension independently
                for expr in res_bb
                    @constraint model expr <= func(scanner)
                    @constraint model expr >= -func(scanner)
                end
            else
                # with rotation: apply constraint to total squared
                total_squared = sum(map(n->n^2, res_bb))
                @constraint model total_squared <= func(scanner)^2
            end
        else
            @constraint model res_bb <= func(scanner)
            @constraint model res_bb >= -func(scanner)
        end
    elseif building_block isa ContainerBlock
        # apply constraints at lower level
        for (_, child_block) in get_children_blocks(building_block)
            scanner_constraints!(child_block, scanner, func)
        end
    end
end


function fixed(bb::BuildingBlock)
    arguments = []
    for name in propertynames(bb)
        push!(arguments, fixed(getproperty(bb, name)))
    end
    return typeof(bb)(arguments...)
end

end