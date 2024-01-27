module BuildingBlocks
import JuMP: has_values, value, Model, @constraint, @objective, owner_model, objective_function, optimize!, AbstractJuMPScalar
import Printf: @sprintf
import ..Scanners: Scanner
import ..Variables: variables, start_time, duration, end_time, gradient_strength, slew_rate

"""
Parent type for all individual components out of which a sequence can be built.

Required methods:
- [`duration`](@ref)(block, parameters): returns block duration in ms.
- `to_concrete_block`(sequence, block): converts the block into a `ConcreteBlock`, which will be part of given `Sequence`.
- [`variables`](@ref): A list of all functions that are used to compute variables of the building block. Any of these can be used in constraints or objective functions.
"""
abstract type BuildingBlock end

"""
Parent type for all RF pulses.

RF pulses combined with gradients, should be childrent of [`ContainerBlock`](@ref) instead.
"""
abstract type RFPulseBlock <: BuildingBlock end

"""
Parent type for all gradient profiles.
"""
abstract type GradientBlock <: BuildingBlock end

"""
Parent type for all types combining one or more pulses/gradients.

Required methods:
- [`get_children_blocks`](@ref): return all the [`BuildingBlock`](@ref) objects includes in this container with their indices.
"""
abstract type ContainerBlock <: BuildingBlock end


"""
    get_children_blocks(container)

Return all the [`BuildingBlock`](@ref) objects includes in this container with their indices.
"""
function get_children_blocks end


"""
    start_time(container, args...)

Returns the starting time of the specific [`BuildingBlock`](@ref) within the container.
The [`BuildingBlock`](@ref) is defined by one or more indices as defined below.
"""
start_time(bb) = 0.

"""
    end_time(container, args...)

Returns the end time of the specific [`BuildingBlock`](@ref) within the container.
The [`BuildingBlock`](@ref) is defined by one or more indices as defined below.
"""
end_time(bb) = duration(bb)


"""
    to_block(object)

Function used internally to convert a wide variety of objects into [`BuildingBlock`](@ref) objects.
"""
to_block(bb::BuildingBlock) = bb


"""
    scanner_constraints!([model, ]building_block, scanner)

Adds any constraints from a specific scanner to a [`BuildingBlock`]{@ref}.
"""
function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner)
    scanner_constraints!(owner_model(building_block), building_block, scanner)
end

function scanner_constraints!(model::Model, building_block::BuildingBlock, scanner::Scanner)
    for func in [gradient_strength, slew_rate]
        if (func in variables(building_block)) && isfinite(func(scanner))
            @constraint model func(building_block) <= func(scanner)
        end
    end
    if building_block isa ContainerBlock
        for (_, child_block) in get_children_blocks(building_block)
            scanner_constraints!(model, child_block, scanner)
        end
    end
end

"""
    variables(building_block)

Returns a list of function that can be called to constrain the `building_block`.
"""
variables(bb::BuildingBlock) = variables(typeof(bb))


function Base.show(io::IO, block::BuildingBlock)
    print(io, string(typeof(block)), "(")
    internal_print(io, block)
    print(io, ")")
end

function internal_print(io::IO, block::BuildingBlock)
    for name in propertynames(block)
        value = getproperty(block, name)
        if value isa AbstractJuMPScalar || value isa Model || string(name)[1] == '_'
            continue
        end
        print(io, name, "=", repr(value), ", ")
    end

    if has_values(block)
        for fn in variables(block)
            if fn == duration
                continue
            end
            numeric_value = value(fn(block))
            printed_value = @sprintf("%.3g", numeric_value)
            print(io, "$(nameof(fn))=$(printed_value), ")
        end
    end
end

# The `start_time` and `end_time` functions are properly defined in "sequence_builders.jl"
function start_time end
function end_time end

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
        apply_simple_constraint!(model, to_funcs[key](block), value)
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
apply_simple_constraint!(model::Model, variable, value::Number) = @constraint model variable == value


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
    error("Cannot find owner model of $bb")
end
has_values(bb::BuildingBlock) = has_values(owner_model(bb))

end