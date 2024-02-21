module Parts
include("trapezoids.jl")
include("slice_select_rephases.jl")
include("spoilt_slice_selects.jl")
include("epi_readouts.jl")
include("helper_functions.jl")

import .Trapezoids: Trapezoid, SliceSelect, LineReadout, opposite_kspace_lines
import .SpoiltSliceSelects: SpoiltSliceSelect
import .SliceSelectRephases: SliceSelectRephase
import .EPIReadouts: EPIReadout
import .HelperFunctions: excitation_pulse, refocus_pulse, readout_event


end