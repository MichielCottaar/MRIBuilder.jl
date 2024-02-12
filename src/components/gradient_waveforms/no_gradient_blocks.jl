module NoGradientBlocks
import StaticArrays: SVector, SMatrix
import ....Variables: VariableType, duration, qvec, bmat_gradient, gradient_strength, slew_rate, get_free_variable
import ...AbstractTypes: GradientWaveform
import ..ChangingGradientBlocks: split_gradient

"""
    NoGradientBlock(duration)

Part of a gradient waveform when there is no gradient active.

Usually, you do not want to create this object directly, use a `BuildingBlock` instead.
"""
struct NoGradientBlock{N} <: GradientWaveform{N}
    duration :: VariableType
    function NoGradientBlock{N}(duration)
        if !(N in (1, 3))
            error("Dimensionality of the gradient should be 1 or 3, not $N")
        end
        new(duration)
    end
end

duration(ngb::NoGradientBlock) = duration(ngb)
for func in (:qvec, :gradient_strength, :slew_rate)
    @eval $func(::NoGradientBlock{1}) = 0.
    @eval $func(::NoGradientBlock{3}) = zero(SVector{3, Float64})
end

bmat_gradient(::NoGradientBlock{1}) = 0.
bmat_gradient(::NoGradientBlock{3}) = zero(SMatrix{3, 3, Float64, 9})

bmat_gradient(::NoGradientBlock, ) = 0.
bmat_gradient(ngb::NoGradientBlock{1}, qstart::VariableType) = qstart^2 * duration(ngb)
bmat_gradient(ngb::NoGradientBlock{3}, qstart::AbstractVector{<:VariableType}) = @. qstart * permutedims(qstart) * duration(ngb)

function split_gradient(ngb::NoGradientBlock, times::VariableType...)
    durations = [times[1], [t[2] - t[1] for t in zip(times[1:end-1], times[2:end])]..., duration(ngb) - times[end]]
    @assert all(durations >= 0.)
    return [NoGradientBlock(d) for d in durations]
end

end