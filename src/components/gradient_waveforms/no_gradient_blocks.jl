module NoGradientBlocks
import StaticArrays: SVector, SMatrix
import ....Variables: VariableType, duration, qval, qvec, bmat_gradient, gradient_strength, slew_rate, get_free_variable
import ...AbstractTypes: GradientWaveform
import ..ChangingGradientBlocks: split_gradient

"""
    NoGradient{N}(duration)

Part of a gradient waveform when there is no gradient active.

`N` needs to be set to 1 if `orientation` is fixed in the gradient waveform or 3 otherwise.

Usually, you do not want to create this object directly, use a `BuildingBlock` instead.
"""
struct NoGradient{N} <: GradientWaveform{N}
    duration :: VariableType
    function NoGradient{N}(duration) where {N}
        if !(N in (1, 3))
            error("Dimensionality of the gradient should be 1 or 3, not $N")
        end
        new(duration)
    end
end

duration(ngb::NoGradient) = ngb.duration
for func in (:qval, :gradient_strength, :slew_rate)
    @eval $func(::NoGradient{1}) = 0.
    @eval $func(::NoGradient{3}) = zero(SVector{3, Float64})
end

NoGradient(duration) = NoGradient{1}(duration)

qvec(::NoGradient, index1, index2) = zero(SVector{3, Float64})
qvec(::NoGradient) = zero(SVector{3, Float64})

bmat_gradient(::NoGradient{1}) = 0.
bmat_gradient(::NoGradient{3}) = zero(SMatrix{3, 3, Float64, 9})
bmat_gradient(ngb::NoGradient{1}, qstart::VariableType) = qstart^2 * duration(ngb)
bmat_gradient(ngb::NoGradient{3}, qstart::AbstractVector{<:VariableType}) = @. qstart * permutedims(qstart) * duration(ngb)

function split_gradient(ngb::NoGradient, times::VariableType...)
    durations = [times[1], [t[2] - t[1] for t in zip(times[1:end-1], times[2:end])]..., duration(ngb) - times[end]]
    return [NoGradient(d) for d in durations]
end

end