module NoGradientBlocks
import StaticArrays: SVector, SMatrix
import ....Variables: VariableType, get_free_variable, adjustable, variables, @defvar
import ...AbstractTypes: GradientWaveform
import ..ChangingGradientBlocks: split_gradient

"""
    NoGradient(duration)

Part of a gradient waveform when there is no gradient active.

Usually, you do not want to create this object directly, use a `BuildingBlock` instead.
"""
struct NoGradient <: GradientWaveform{0}
    duration :: VariableType
end

for func in (:qval, :gradient_strength, :slew_rate)
    @eval variables.$(func).f(::NoGradient) = 0.
end

@defvar begin
    duration(ngb::NoGradient) = ngb.duration
    gradient_strength(nb::NoGradient, time::Number) = 0.
end

@defvar gradient begin
    bmat_gradient(::NoGradient) = 0.
    bmat_gradient(ngb::NoGradient, qstart::VariableType) = qstart^2 * duration(ngb)
    bmat_gradient(ngb::NoGradient, qstart::AbstractVector{<:VariableType}) = @. qstart * permutedims(qstart) * duration(ngb)
end

function split_gradient(ngb::NoGradient, times::VariableType...)
    durations = [times[1], [t[2] - t[1] for t in zip(times[1:end-1], times[2:end])]..., duration(ngb) - times[end]]
    return [NoGradient(d) for d in durations]
end

adjustable(::NoGradient) = :false

end