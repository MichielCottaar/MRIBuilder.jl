module Pathways
import LinearAlgebra: norm
import StaticArrays: SVector, SMatrix
import ..Variables: qval, bval

"""
    PathWay(sequence::Sequence, pulse_effects::Vector{:Symbol/Number}, readout_index=1)

Describes how a specific spin/isochromat might experience the sequence.

Only a single pathway through the RF pulses is considered,
so that at every point in time the spins are in one of the following four states:
- +longitudinal: initial relaxed state
- +transverse: excited state. During this time gradients will affect the [`area_under_curve`](@ref) (or [`qval`](@ref)) and [`bval`](@ref).
- -longitudinal: inverse state
- -transverse: inverse excited state. During this time all gradients will have the inverse effect compared with +transverse.

The RF pulses cause mappings between these different states as described below.

## Parameters
- `sequence`: MRI [`Sequence`](@ref) to be considered.
- `pulse_effects`: How each RF pulse affects the spins. This can be one of the following:
    - `:skip`/`:ignore`/0: This RF pulse leaves the spins unaffected.
    - `:refocus`/`:invert`/180: Flips the sign of the spin state (i.e., +longitudinal <-> -longitudinal, +transverse <-> -transverse)
    - `:excite`/90: Takes spin state one step along the following sequence +longitudinal -> +transverse -> -longitudinal -> -transverse -> +longitudinal
    - `:neg_excite`/270/-90: Inverse step compared with `:excite`.
- `readout_index`: After encountering the number of pulses as defined in `pulse_effects`, continue the `PathWay` until the readout given by `index` is reached. If set to 0 the `PathWay` is terminated immediately after the last RF pulse.

## Attributes
Over the pathway the following values are computed. Each can be accessed by calling the appropriate function:

### Timings
- [`duration_state`](@ref): The total amount of time spent in a specific state in this pathway (+longitudinal, +transverse, -longitudinal, or -transverse)
- [`duration_transverse`](@ref): The total amount of time the spins spent in the transverse plane in ms. This can be used to quantify the expected effect of T2-decay.
- [`duration_dephase`](@ref): The total amount of time the spins spent in the +transverse relative to -transverse state in ms. The absolute value of this can be used to quantify the expected effect of T2'-decay.

### Effect of gradients
- [`qvec`](@ref): Net displacement vector in k-space/q-space.
- [`qval`](@ref)/[`area_under_curve`](@ref): size of the displacement in k-space/q-space. For a spoiled pathway, this should be large compared with 1/voxel size; for unspoiled pathways it should be (close to) zero.
- [`bmat`](@ref): Net diffusion weighting due to gradients along the [`PathWay`](@ref) in matrix form.
- [`bval`](@ref): Net diffusion weighting due to gradients along the [`PathWay`](@ref) as a single number.
"""
struct PathWay
    # user provided
    sequence :: Sequence
    pulse_effects :: Vector{Union{:Symbol, :Number}}
    readout_index :: Integer

    # computed
    duration_states :: SVector{4, Float64}
    qvec :: SVector{3, Float64}
    bmat :: SMatrix{3, 3, Float64, 9}
end

"""
    duration_state(pathway::PathWay, transverse::Bool, positive::Bool)

Returns how long the [`PathWay`](@ref) spent in a specific state.

The requested state can be set using `transverse` and `positive` as follows:
- `transverse=false`, `positive=true`: +longitudinal
- `transverse=true`, `positive=true`: +transverse
- `transverse=false`, `positive=false`: -longitudinal
- `transverse=true`, `positive=false`: -transverse
"""
function duration_state(pathway::PathWay, transverse, positive)
    index = Dict([
        (false, true) => 1,
        (true, true) => 2,
        (false, false) => 3,
        (true, false) => 4,
    ])[(Bool(transverse), Bool(positive))]
    return pathway.duration_states[index]
end

"""
    duration_transverse(pathway::PathWay)

Returns the total amount of time that spins following the given [`PathWay`](@ref) spent in the transverse plane.
This determines the amount of T2-weighting as ``e^{t/T_2}``, where ``t`` is the `duration_transverse`.

Also see [`duration_dephase`](@ref) for T2'-weighting.
"""
function duration_transverse(pathway::Pathway)
    return duration_state(pathway, true, true) + duration_state(pathway, true, false)
end

"""
    duration_dephase(pathway::PathWay)

Returns the net time that spins following the given [`PathWay`](@ref) spent in the +transverse versus the -transverse state.
This determines the amount of T2'-weighting as ``e^{t/T_2'}``, where ``t`` is the `duration_dephase`.

Also see [`duration_transverse`](@ref) for T2-weighting.
"""
function duration_dephase(pathway::Pathway)
    return duration_state(pathway, true, true) - duration_state(pathway, true, false)
end


"""
    qvec(pathway::PathWay)

Return net displacement vector in k-space/q-space experienced by the spins following a specific [`PathWay`](@ref).

Only gradients active while the spins are in the transverse plane are considered.
"""
qvec(pathway::PathWay) = pathway.qvec

"""
    qval(pathway::PathWay)

Return net displacement in k-space/q-space experienced by the spins following a specific [`PathWay`](@ref).

Only gradients active while the spins are in the transverse plane are considered.
"""
qval(pathway::PathWay) = norm(qvec(pathway))

"""
    bmat(pathway::PathWay)

Return 3x3 diffusion-weighted matrix experienced by the spins following a specific [`PathWay`](@ref).

Only gradients active while the spins are in the transverse plane are considered.
"""
bmat(pathway::PathWay) = pathway.bmat

"""
    bval(pathway::PathWay)

Return size of diffusion-weighting experienced by the spins following a specific [`PathWay`](@ref).

Only gradients active while the spins are in the transverse plane are considered.
"""
bval(pathway::PathWay) = norm(bmat(pathway))

end