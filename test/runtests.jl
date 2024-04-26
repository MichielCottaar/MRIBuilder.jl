using MRIBuilder
using Test

@testset "MRIBuilder.jl" begin
    include("test_building_blocks.jl")
    include("test_post_hoc.jl")
    include("test_components.jl")
    include("test_sequences.jl")
    include("test_IO.jl")
    include("test_plot.jl")
end
