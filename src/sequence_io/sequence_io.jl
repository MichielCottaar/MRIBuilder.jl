module SequenceIO

include("pulseq_io/pulseq_io.jl")

import .PulseqIO: read_pulseq
import Serialization: serialize, deserialize


function read_sequence(filename::AbstractString, args...; kwargs...)
    open(filename, "r") do io
        read_sequence(io, args...; kwargs...)
    end
end

function read_sequence(io::IO; format=nothing, kwargs...)
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
        return read_pulseq(io, kwargs...)
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

function write_sequence(io::IO, sequence; format::Symbol=:serialize)
    if format == :serialize
        serialize(io, sequence)
    else
        error("Unrecognised selected format for output ($format). Only :serialize is supported..")
    end
end

end