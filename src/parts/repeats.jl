module Repeats
import ...Containers: BaseSequence, get_index_single_TR, Wait
import ...Variables: variables, VariableType, AbstractBlock, get_free_variable, set_simple_constraints!, apply_simple_constraint!, @defvar, get_pulse, get_gradient, get_readout

"""
    Repeat(block, variables...)

# Parameters
- `block`: RF pulse/gradient/ADC/sequence part that needs to be repeated.

# Variables
- `nrepeat`: number of times the block is repeated.
- `wait_time`: time between end of one block and start of the next repeat in ms.
- `interblock_delay`: time between start of one block and start of the next repeat in ms.
- `frequency`: temporal frequency of the block repeats (i.e., `1/interblock_delay`) in kHz.
- `duration`: total duration of this sequence part all repeats in ms.
"""
struct Repeat{T<:AbstractBlock} <: BaseSequence{1}
    block :: T
    wait_time :: VariableType
    nrepeat :: VariableType
end

function Repeat(block::AbstractBlock; wait_time=nothing, nrepeat=nothing, variables...)
    res = Repeat{typeof(block)}(block, get_free_variable(wait_time), get_free_variable(nrepeat, integer=true))
    set_simple_constraints!(res, variables)
    apply_simple_constraint!(res.wait_time, :>=, 0.)
    apply_simple_constraint!(res.nrepeat, :>=, 1)
    return res
end

function get_index_single_TR(r::Repeat, i::Integer)
    if isone(i % 2)
        return r.block
    else
        return Wait(r.wait_time)
    end
end

Base.length(r::Repeat) = r.nrepeat * 2 - 1

@defvar begin
    nrepeat(r::Repeat) = r.nrepeat
    wait_time(r::Repeat) = r.wait_time
    duration(r::Repeat) = r.interblock_delay * r.nrepeat - r.wait_time
    interblock_delay(r::Repeat) = r.block.duration + r.wait_time
    frequency(r::Repeat) = 1/(r.interblock_delay)
end

get_pulse(r::Repeat) = get_pulse(r.block)
get_gradient(r::Repeat) = get_gradient(r.block)
get_readout(r::Repeat) = get_readout(r.block)

end