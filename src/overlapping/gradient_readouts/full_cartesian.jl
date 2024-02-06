module Cartesian
import ....Variables: inverse_fov, resolution, inverse_voxel_size
import ...GradientPulses: TrapezoidGradient
import ..SingleLines: SingleLine

"""
    FullCartesianReadout(nlines_per_readout, nreadout, )
"""
struct FullCartesianReadout
    nlines_per_readout :: Int
    nreadout :: Int
    prepare :: TrapezoidGradient
    positive_readout :: SingleLine
    negative_readout :: SingleLine
    blip :: TrapezoidGradient
end

inverse_fov(fc::FullCartesianReadout) = [inverse_fov(fc.positive_readout), 1e3 * qvec(fc.blip)[2]]
resolution(fc::FullCartesianReadout) = [resolution(fc.positive_readout), fc.nlines_per_readout * fc.nreadout]
inverse_voxel_size(fc::FullCartesianReadout) = inverse_fov(fc) .* resolution(fc)



end