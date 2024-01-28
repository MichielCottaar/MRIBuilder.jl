module Pathways
import LinearAlgebra: norm
import StaticArrays: SVector, SMatrix
import ..Containers: Sequence
import ..Variables: qval, bval, area_under_curve
import ..PathwayWalkers: PathwayWalker

"""
    Pathway(sequence::Sequence, pulse_effects::Vector{:Symbol/Number}, readout_index=1)

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
Some gradients will be scaled/rotated with user-provided values (e.g., bvals/bvecs).
The area under curve, q-values, and b-values are computed separately for such gradients.
You can select which gradients to consider when accessing these values.
- [`qvec`](@ref): Net displacement vector in k-space/q-space.
- [`qval`](@ref)/[`area_under_curve`](@ref): size of the displacement in k-space/q-space. For a spoiled pathway, this should be large compared with 1/voxel size; for unspoiled pathways it should be (close to) zero.
- [`bmat`](@ref): Net diffusion weighting due to gradients along the [`Pathway`](@ref) in matrix form.
- [`bval`](@ref): Net diffusion weighting due to gradients along the [`Pathway`](@ref) as a single number.
"""
struct Pathway
    # user provided
    sequence :: Sequence
    pulse_effects :: Vector{Union{:Symbol, :Number}}
    readout_index :: Integer

    # computed
    duration_states :: Dict{Any, SVector{4, Float64}}
    qvec :: Dict{Any, SVector{3, Float64}}
    bmat :: Dict{Any, SMatrix{3, 3, Float64, 9}}
end

function Pathway(sequence::Sequence, pulse_effects::AbstractVector, readout_index::Integer)
    walker = PathwayWalker(pulse_effects, readout_index)
    walk_pathway!(walker, sequence)
    return Pathway(
        sequence,
        pulse_effects,
        readout_index,

        Dict(k => SVector{4, Float64}(v) for (k, v) in pairs(walker.duration_states)),
        Dict(k => SVector{3, Float64}(v) for (k, v) in pairs(walker.qvec)),
        Dict(k => SMatrix{3, 3, Float64, 9}(v) for (k, v) in pairs(walker.bmat)),
    )
end


"""
    duration_state(pathway::Pathway, transverse::Bool, positive::Bool)

Returns how long the [`Pathway`](@ref) spent in a specific state.

The requested state can be set using `transverse` and `positive` as follows:
- `transverse=false`, `positive=true`: +longitudinal
- `transverse=true`, `positive=true`: +transverse
- `transverse=false`, `positive=false`: -longitudinal
- `transverse=true`, `positive=false`: -transverse
"""
function duration_state(pathway::Pathway, transverse, positive)
    index = Dict([
        (false, true) => 1,
        (true, true) => 2,
        (false, false) => 3,
        (true, false) => 4,
    ])[(Bool(transverse), Bool(positive))]
    return pathway.duration_states[index]
end

"""
    duration_transverse(pathway::Pathway)

Returns the total amount of time that spins following the given [`Pathway`](@ref) spent in the transverse plane.
This determines the amount of T2-weighting as ``e^{t/T_2}``, where ``t`` is the `duration_transverse`.

Also see [`duration_dephase`](@ref) for T2'-weighting.
"""
function duration_transverse(pathway::Pathway)
    return duration_state(pathway, true, true) + duration_state(pathway, true, false)
end

"""
    duration_dephase(pathway::Pathway)

Returns the net time that spins following the given [`Pathway`](@ref) spent in the +transverse versus the -transverse state.
This determines the amount of T2'-weighting as ``e^{t/T_2'}``, where ``t`` is the `duration_dephase`.

Also see [`duration_transverse`](@ref) for T2-weighting.
"""
function duration_dephase(pathway::Pathway)
    return duration_state(pathway, true, true) - duration_state(pathway, true, false)
end


"""
    qvec(pathway::Pathway; scale=nothing, rotate=nothing)

Return net displacement vector in k-space/q-space experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by speficic `scale`/`rotate` parameters
"""
qvec(pathway::Pathway; scale=nothing, rotate=nothing) = get(pathway.qvec, (scale, rotate), zero(SVector{3, Float64}))

"""
    qval(pathway::Pathway; scale=nothing, rotate=nothing)

Return net displacement in k-space/q-space experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by speficic `scale`/`rotate` parameters
"""
qval(pathway::Pathway; scale=nothing, rotate=nothing) = norm(qvec(pathway; scale, rotate))

"""
    area_under_curve(pathway::Pathway; scale=nothing, rotate=nothing)

Return net displacement in k-space (i.e., spoiling) experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by speficic `scale`/`rotate` parameters
"""
area_under_curve(pathway::Pathway; scale=nothing, rotate=nothing) = qval(pathway; scale, rotate)


"""
    bmat(pathway::Pathway; scale=nothing, rotate=nothing)

Return 3x3 diffusion-weighted matrix experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by speficic `scale`/`rotate` parameters
"""
bmat(pathway::Pathway; scale=nothing, rotate=nothing)  = get(pathway.qvec, (scale, rotate), zero(SMatrix{3, 3, Float64, 9}))

"""
    bval(pathway::Pathway; scale=nothing, rotate=nothing)

Return size of diffusion-weighting experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by speficic `scale`/`rotate` parameters
"""
bval(pathway::Pathway; scale=nothing, rotate=nothing) = tr(bmat(pathway; scale, rotate))

end