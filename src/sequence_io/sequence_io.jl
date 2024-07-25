module SequenceIO

include("pulseq_io/pulseq_io.jl")
include("pulseq.jl")

import Serialization: serialize, deserialize
import .PulseqIO: read_pulseq, write_pulseq, PulseqSequence
import ..Containers: Sequence


function read_sequence(filename::AbstractString, args...; kwargs...)
    open(filename, "r") do io
        read_sequence(io, args...; kwargs...)
    end
end

function read_sequence(io::IO; format=nothing, B0=nothing, scanner=nothing,kwargs...)
    if isnothing(format)
        pos = position(io)
        for format in (:pulseq, :serialize)
            try
                return read_sequence(io; format=format, kwargs...)
            catch
                seek(io, pos)
            end
        end
        error("Could not read the input filename. Tried all formats (:pulseq/:serialize).")
    end
    if format == :pulseq
        return Sequence(read_pulseq(io, kwargs...), B0=B0, scanner=scanner)
    elseif format == :serialize
        return deserialize(io, kwargs...)
    else
        error("Cannot read file $filename. Extension is not recognised")
    end
end


function write_sequence(filename::AbstractString, args...; kwargs...)
    open(filename, "w") do io
        write_sequence(io, args...; kwargs...)
    end
end

function write_sequence(io::IO, sequence; format::Symbol=:pulseq)
    if format == :pulseq
        write_pulseq(io, PulseqSequence(sequence))
    elseif format == :serialize
        serialize(io, sequence)
    else
        error("Unrecognised selected format for output ($format). Only :pulseq or :serialize are supported..")
    end
end

end