module AllBuildingBlocks
include("base_building_blocks.jl")
include("building_blocks.jl")
include("trapezoids.jl")
include("spoilt_slice_selects.jl")
include("wait_blocks.jl")

import .BaseBuildingBlocks: waveform, waveform_sequence, events, BaseBuildingBlock
import .BuildingBlocks: BuildingBlock
import .Trapezoids: Trapezoid, SliceSelect, LineReadout
import .SpoiltSliceSelects: SpoiltSliceSelect
import .WaitBlocks: Wait

end