module AllBuildingBlocks
include("base_building_blocks.jl")
include("building_blocks.jl")
include("trapezoids.jl")

import .BaseBuildingBlocks: waveform, waveform_sequence, events, BaseBuildingBlock
import .BuildingBlocks: BuildingBlock
import .Trapezoids: Trapezoid, SliceSelect, LineReadout

end