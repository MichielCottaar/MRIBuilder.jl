module Printing
import JuMP: value
import Printf: @sprintf
import ..BuildingBlocks: BuildingBlock, get_children_indices, VariableNotAvailable
import ..Alternatives: AlternativeBlocks
import ..Overlapping: AbstractOverlapping, waveform, interruptions
import ..Sequences: Sequence
import ..Variables: VariableType, duration, start_time, variables, alternative_variables, effective_time, flip_angle
import ..Pulses: GenericPulse

struct BuildingBlockPrinter{T<:BuildingBlock}
    bb::T
    start_time::Union{Nothing, Number}
    spaces::Int
end

Base.show(io::IO, block::BuildingBlock) = print(io, BuildingBlockPrinter(block, nothing, 0))
Base.show(io::IO, seq::Sequence) = print(io, BuildingBlockPrinter(seq, 0., 0))

function Base.show(io::IO, alt_printer::BuildingBlockPrinter{<:AlternativeBlocks})
    block = alt_printer.bb
    print(io, "AlternativeBlocks(", block.name, ", ", length(block), " options)")
end


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
        if (
            block isa AbstractOverlapping && (
            ft <: Union{Nothing, BuildingBlock} ||
            (ft <: AbstractVector && eltype(ft) <: BuildingBlock) ||
            name == :interruptions
        ))
            continue
        end

        print(io, name, "=", repr(getproperty(block, name)), ", ")
    end

    for fn in values(variables)
        if printed_duration && fn == duration
            continue
        end
        if fn in [fn_alt for (fn_alt, _, _, _) in values(alternative_variables)]
            continue
        end
        try
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
        catch e
            if e isa VariableNotAvailable
                continue
            end
            rethrow()
        end
    end
    print(io, ")")
    if block isa AbstractOverlapping
        print(io, ":")
        ref_start_time = isnothing(printer.start_time) ? 0. : printer.start_time
        inter = copy(interruptions(block))
        prev_time = 0.
        for (index_grad, (time, grad)) in enumerate(waveform(block))
            while length(inter) > 0 && index_grad == (inter[1].index + 1)
                to_print = popfirst!(inter)
                t_eff = ref_start_time + prev_time + to_print.time
                as_printer = BuildingBlockPrinter(to_print.object, t_eff - effective_time(to_print.object), printer.spaces + 2)
                print(io, "\n", repeat(' ', printer.spaces + 2), "- ", @sprintf("%.3g", t_eff), ": ", as_printer)
            end
            printed_grad = "[" * join(map(v -> @sprintf("%.3g", v), grad), ", ") * "]"
            print(io, "\n", repeat(' ', printer.spaces + 2), "- ", @sprintf("%.3g", ref_start_time + time), ": ", printed_grad)
            prev_time = time
        end
    end
end

function Base.show(io::IO, printer::BuildingBlockPrinter{<:Sequence})
    seq = printer.bb
    print(io, "Sequence(")
    d = _robust_value(duration(seq))
    if !isnothing(d)
        if !isnothing(printer.start_time)
            print(io, "t=", @sprintf("%.3g", printer.start_time), "-", @sprintf("%.3g", printer.start_time + d), ", ")
        else
            print(io, "duration=", @sprintf("%.3g", d), ",")
        end
    end
    TR = _robust_value(seq.TR)
    if !isnothing(TR) && isfinite(TR)
        print(io, "TR=", Int(round(TR)))
    end
    print(io, "):")

    for child_index in get_children_indices(seq)
        child_block = seq[child_index]
        child_printer = BuildingBlockPrinter(
            child_block,
            isnothing(printer.start_time) ? nothing : _robust_value(start_time(seq, child_index) + printer.start_time),
            printer.spaces + 2
        )
        print(io, "\n", repeat(' ', printer.spaces + 2), "- ", child_index, ": ", child_printer)
    end

end

function Base.show(io::IO, printer::BuildingBlockPrinter{<:GenericPulse})
    fp = printer.bb
    pp(value::Number) = @sprintf("%.3g", value)

    print(io, "GenericPulse(")
    if isnothing(printer.start_time)
        print(io, "duration=", pp(duration(fp)))
    else
        print(io, "t=", pp(printer.start_time), "-", pp(printer.start_time + duration(fp)))
    end
    print(io, ", flip_angle=", pp(flip_angle(fp)), ", effective_time=", pp(effective_time(fp)), ")")
end




end