"""
Defines [`BaseSequence`](@ref) and [`Sequence`](@ref)
"""
module BaseSequences
import StaticArrays: SVector
import JuMP: @constraint
import ...Variables: get_free_variable, VariableType, variables, set_simple_constraints!, make_generic, get_gradient, get_pulse, get_gradient, @defvar, add_cost_function!
import ...BuildSequences: global_scanner
import ...Components: EventComponent, NoGradient, edge_times
import ...Scanners: Scanner, B0
import ..Abstract: ContainerBlock, start_time
import ..BuildingBlocks: Wait, BuildingBlock, BaseBuildingBlock

"""
Super-type of any sequence of non-overlapping building blocks that should be played after each other.

It contains `N` [`ContainerBlock`](@ref) objects (e.g., building blocks or other sequences).

Main interface:
- Acts as an iterable containing the blocks and sequences.
    - Indiviual blocks/sequences can be obtained using indexing.
    - If there is a finite number of repeats, the iteration will continue over all repeats.

Sub-types need to implement:
- `get_index_single_TR`: return the index assuming it is between 1 and N
"""
abstract type BaseSequence{N} <: ContainerBlock end

function Base.getindex(bs::BaseSequence{N}, index::Integer) where {N}
    if nrepeat(bs) > 0 && (index < 1 || index > length(bs))
        throw(BoundsError(bs, index))
    end
    base_index = ((index - 1) % N) + 1
    return get_index_single_TR(bs, base_index)
end

function (sequence::BaseSequence{N})(time::AbstractFloat) where {N}
    var_time = mod(time, variables.duration(sequence))
    for block in sequence
        var_time -= variables.duration(block)
        if var_time < 0
            return (var_time + variables.duration(block), block)
        end
    end
    if abs(var_time) <= 1e-6
        return (variables.duration(sequence[N]) + var_time, sequence[N])
    end
    error("Total duration of building blocks does not match duration of the sequence.")
end

Base.iterate(bs::BaseSequence) = Base.iterate(bs, 1)
Base.iterate(bs::BaseSequence{N}, index::Integer) where {N} = index > length(bs) ? nothing : (bs[index], index + 1)
Base.length(bs::BaseSequence{N}) where {N} = iszero(nrepeat(bs)) ? N : (nrepeat(bs) * N)
Base.eltype(::Type{<:BaseSequence}) = ContainerBlock
Base.keys(seq::BaseSequence{N}) where {N} = 1:N

function start_time(bs::BaseSequence{N}, index::Integer) where {N}
    nTR = div(index-1, N, RoundDown)
    if 0 < nrepeat(bs) <= nTR
        throw(BoundsError(bs, index))
    end
    base_index = ((index - 1) % N) + 1
    base_time = sum(i -> variables.duration(bs[i]), 1:base_index-1; init=0.)
    if iszero(nTR)
        return base_time
    else
        return nTR * variables.duration(bs) + base_time
    end
end

function start_time(bs::BaseSequence{N}, s::Symbol) where {N}
    index = findfirst(p -> p[1] == s, bs.blocks)
    if isnothing(index)
        error("Sequence does not contain a block with name $s")
    end
    return start_time(bs, index)
end

function gradient_orientation(seq::BaseSequence)
    return gradient_orientation(get_gradient(seq))
end

gradient_orientation(nt::NamedTuple) = map(gradient_orientation, nt)


"""
    get_index_single_TR(sequence, index)

Used internally by any [`BaseSequence`](@ref) to get a specific block.
The `index` should be between 1 and N.
It should be implemented for any sub-classes of [`BaseSequence`](@ref).
"""
get_index_single_TR(seq::BaseSequence, index::Integer) = get_index_single_TR(seq, Val(index))


"""
    nrepeat(sequence)

How often sequence should be repeated. 
"""
nrepeat(bs::BaseSequence) = 1

@defvar begin
    duration(bs::BaseSequence{0}) = 0.
    duration(bs::BaseSequence) = sum(variables.duration.(bs); init=0.)
end

function edge_times(seq::BaseSequence; tol=1e-6)
    res = Float64[]
    for (index, block) in enumerate(seq)
        append!(res, edge_times(block) .+ start_time(seq, index))
    end
    unique_res = Float64[res[1]]
    for edge in res
        if !isapprox(edge, unique_res[end], atol=tol)
            push!(unique_res, edge)
        end
    end
    return sort(unique_res)
end

for fn in (:gradient_strength, :amplitude, :phase, :frequency)
    @eval function variables.$fn.f(sequence::BaseSequence, time::AbstractFloat)
        (block_time, block) = sequence(time)
        return variables.$fn.f(block, block_time)
    end
end

for fn in (:get_gradient, :get_pulse)
    @eval function $fn(sequence::BaseSequence, time::AbstractFloat)
        (block_time, block) = sequence(time)
        return $fn(block, block_time)
    end
end
"""
    Sequence(blocks; name=:Sequence, variables...)
    Sequence(blocks...; name=:Sequence, variables...)

Defines an MRI sequence from a vector of building blocks.

## Arguments
- `blocks`: The actual building blocks that will be played in sequence. All the building blocks must be of type [`ContainerBlock`](@ref), which means that they cannot only contain actual [`BaseBuildingBlock`](@ref) objects, but also other [`BaseSequence`](@ref) objects.
    Objects of a different type are converted into a [`ContainerBlock`](@ref) internally:
    - numbers/`nothing`/`:min`/`:max` : replaced with a [`Wait`](@ref) block with the appropriate constraint/objective added to its [`variables.duration`](@ref).
    - RF pulse or readout: will be embedded within a [`BuildingBlock`](@ref) of the appropriate length

Specific named sequences might define additional variables.
"""
struct Sequence{S, N} <: BaseSequence{N}
    blocks :: SVector{N, Pair{<:Union{Symbol, Nothing}, <:ContainerBlock}}
    scanner :: Scanner
end

function Sequence(blocks::AbstractVector; name=:Sequence, scanner=nothing, vars...)
    blocks = to_block_pair.(blocks)
    res = Sequence{name, length(blocks)}(SVector{length(blocks)}(blocks), isnothing(scanner) ? global_scanner() : scanner)
    set_simple_constraints!(res, vars)
    add_cost_function!(variables.duration(res), 3)
    return res
end

Base.show(io::IO, ::Type{<:Sequence{S, N}}) where {S, N} = print(io, S, "{$N}")
B0(sequence::Sequence) = B0(sequence.scanner)

Sequence(blocks...; kwargs...) = Sequence([blocks...]; kwargs...)

get_index_single_TR(s::Sequence, i::Integer) = s.blocks[i][2]
Base.getindex(seq::Sequence, sym::Symbol) = seq[findfirst(p -> p[1] == sym, seq.blocks)]
nrepeat(::Sequence) = 0

@defvar repetition_time(seq::Sequence) = duration(seq)

to_block_pair(pair::Pair) = pair[1] => to_block(pair[2])
to_block_pair(other) = nothing => to_block(other)

"""
    to_block(block_like)

Converst object into something that can be included in the sequence:
- :min/:max/number/variable/nothing => [`Wait`](@ref).
- `building_block` or `sequence` => no change.
- RF pulse/readout => will be embedded within a [`BuildingBlock`](@ref).
"""
to_block(cb::ContainerBlock) = cb
to_block(s::Symbol) = to_block(Val(s))
to_block(s::Union{VariableType, Nothing, Val{:min}, Val{:max}}) = Wait(s)
to_block(ec::EventComponent) = BuildingBlock([NoGradient(variables.duration(ec)), (0., ec)])


function make_generic(seq::BaseSequence)
    blocks = BaseBuildingBlock[]

    function add_sequence!(to_add, sub_seq)
        for block in sub_seq
            if block isa BaseSequence
                add_sequence!(to_add, block)
            else
                push!(to_add, make_generic(block))
            end
        end
    end
    add_sequence!(blocks, seq)

    return Sequence(blocks)
end



end
