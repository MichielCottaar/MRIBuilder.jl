module BuildingBlocks
import JuMP: value, Model, @constraint, @objective, objective_function, AbstractJuMPScalar
import ..Variables: Variables, variables, start_time, duration, end_time, gradient_strength, slew_rate, effective_time, VariableType, alternative_variables, qval_square
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
- [`get_children_blocks`](@ref)(container): return all the [`BuildingBlock`](@ref) objects includes in this container.
- [`start_time`](@ref)(container, index): returns the starting time of the child corresponding to `index` relative to the start of the `container` in ms.
- `Base.getindex`(container, index): get child [`BuildingBlock`](@ref) corresponding to `index`.
"""
abstract type ContainerBlock <: BuildingBlock end


"""
    get_children_blocks(container)

Return all the [`BuildingBlock`](@ref) objects includes in this container.
"""
get_children_blocks(bb::BuildingBlock) = [bb[i] for i in get_children_indices(bb)]

"""
    get_children_indices(container)

Return the indices of all the children in a [`ContainerBlock`](@ref).

This needs to be defined for every [`ContainerBlock`](@ref).
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
    match_blocks!(block1, block2, property_list)

Matches the listed variables between two [`BuildingBlock`](@ref) objects.
"""
function match_blocks!(block1::BuildingBlock, block2::BuildingBlock, property_list)
    for fn in property_list
        @constraint global_model() fn(block1) == fn(block2)
    end
end

"""
    scanner_constraints!(building_block[, scanner])

Adds the gradient strength and slew rate constraints from a specific [`Scanner`](@ref) to a [`BuildingBlock`]{@ref}.

This is applied iteratively to each part of a `Sequence`.
"""
function scanner_constraints!(building_block::BuildingBlock) 
    try
        scanner_constraints!(building_block, global_scanner())
    catch e
        if occursin("No valid scanner", e.msg)
            return
        end
        rethrow()
    end
end

function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner)
    for func in [gradient_strength, slew_rate]
        if isfinite(func(scanner))
            scanner_constraints!(building_block, scanner, func)
        end
    end
end

function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner, func::Function)
    model = global_model()
    try
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

    catch e
        if !(e isa VariableNotAvailable)
            rethrow()
        end
        if building_block isa ContainerBlock
            for child_block in get_children_blocks(building_block)
                scanner_constraints!(child_block, scanner, func)
            end
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


"""
    make_generic(bb::BuildingBlock)

Replaces a [`BuildingBlock`](@ref) or whole sequence with a generic version.

This replaces all functional RF pulses and gradient waveforms with their generic equivalents.
"""
function make_generic end

end