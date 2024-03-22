module SequenceIO

include("pulseq.jl")

import .Pulseq: read_pulseq


function read_sequence(filename::AbstractString)
    if endswith(filename, ".seq")
        return read_pulseq(filename)
    else
        error("Cannot read file $filename. Extension is not recognised")
    end
end

end