module BuildingBlocks
import JuMP: has_values, value, Model, @constraint, @objective, owner_model, objective_function, optimize!, AbstractJuMPScalar
import Printf: @sprintf
import ..Scanners: Scanner
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
    scanner_constraints!([model, ]building_block, scanner)

Adds any constraints from a specific scanner to a [`BuildingBlock`]{@ref}.
"""
function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner)
    scanner_constraints!(owner_model(building_block), building_block, scanner)
end

function scanner_constraints!(model::Model, building_block::BuildingBlock, scanner::Scanner)
    for func in [gradient_strength, slew_rate]
        if isfinite(func(scanner))
            scanner_constraints!(model, building_block, scanner, func)
        end
    end
end

function scanner_constraints!(model::Model, building_block::BuildingBlock, scanner::Scanner, func::Function)
    if func in variables(building_block)
        # apply constraint at this level
        res_bb = func(building_block)
        if res_bb isa AbstractVector
            if isnothing(building_block.rotate)
                # no rotation; apply constraint to each dimension independently
                for expr in res_bb
                    @constraint model expr <= func(scanner)
                end
            else
                # with rotation: apply constraint to total squared
                total_squared = sum(map(n->n^2, res_bb))
                @constraint model total_squared <= func(scanner)^2
            end
        else
            @constraint model res_bb <= func(scanner)
        end
    elseif building_block isa ContainerBlock
        # apply constraints at lower level
        for (_, child_block) in get_children_blocks(building_block)
            scanner_constraints!(model, child_block, scanner, func)
        end
    end
end

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
            ft <: BuildingBlock ||
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

    if block isa ContainerBlock
        use_start_time = isnothing(printer.start_time) ? 0. : printer.start_time
        print(io, ":")
        for (child_index, child_block) in get_children_blocks(block)
            child_printer = BuildingBlockPrinter(
                child_block,
                _robust_value(start_time(block, child_index) + use_start_time),
                printer.spaces + 2
            )
            print(io, "\n", repeat(' ', printer.spaces + 2), "- ", child_index, ": ", child_printer)
        end
    end
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
    for (key, value) in kwargs
        if key == :qval
            apply_simple_constraint!(model, qval_square(block), value isa VariableType ? value^2 : value)
        else
            apply_simple_constraint!(model, to_funcs[key](block), value)
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