using MRIBuilder
using Test

@testset "MRIBuilder.jl" begin
    include("test_sequences.jl")
    include("test_IO.jl")
    include("test_plot.jl")
end
