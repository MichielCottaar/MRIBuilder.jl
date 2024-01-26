module InsertedBuildingBlocks
import Printf: @sprintf
import JuMP: owner_model, has_values, Model
import ..GlobalModel: @global_model_constructor
import ..Variables: Variables, duration, start_time, variables, end_time, get_free_variable
import ..BuildingBlocks: BuildingBlock, internal_print

"""
    InsertedBuildingBlock(building_block; start_time=nothing)
    InsertedBuildingBlock(building_block, index)

Represents a specific insertion of the [`BuildingBlock`](@ref) object within a larger sequence.

If `index` is set an existing [`InsertedBuildingBlock`](@ref) is returned. Otherwise, a new one is created with the given `start_time` as constraint.
"""
struct InsertedBuildingBlock{T<:BuildingBlock}
    bb :: T
    index :: Int
    function InsertedBuildingBlock(bb, index)
        if index < 1 || index > length(bb)
            error("$index is out of range for $bb")
        end
        return new{typeof(bb)}(bb, index)
    end
end

owner_model(inserted::InsertedBuildingBlock) = owner_model(inserted.bb)
has_values(inserted::InsertedBuildingBlock) = has_values(owner_model(inserted))

@global_model_constructor InsertedBuildingBlock

function InsertedBuildingBlock(model::Model, bb::BuildingBlock; start_time=nothing)
    index = length(bb.start_time) + 1
    push!(bb.start_time, get_free_variable(model, start_time))
    InsertedBuildingBlock(bb, index)
end


function Base.show(io::IO, inserted::InsertedBuildingBlock)
    print(io, string(typeof(block)), "(")
    if has_values(inserted)
        if iszero(value(duration(block)))
            print(io, "time=", @sprintf("%.3f", value(start_time(block))), ", ")
        else
            print(
                io, "time=", 
                @sprintf("%.3f", value(start_time(block))), "-",
                @sprintf("%.3f", value(end_time(block))), ", "
            )
        end
    end
    internal_print(io, inserted.bb)
    print(io, ")")
end

for func in variables()
    if func in (start_time, end_time)
        continue
    end
    @eval Variables.$(nameof(func))(inserted::InsertedBuildingBlock, args...; kwargs...) = Variables.$(nameof(func))(inserted.bb, args...; kwargs...)
end


end