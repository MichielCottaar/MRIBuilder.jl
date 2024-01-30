module Pathways
import LinearAlgebra: norm
import StaticArrays: SVector, SMatrix
import ..BuildingBlocks: BuildingBlock, GradientBlock, RFPulseBlock, ContainerBlock, get_children_blocks
import ..Containers: Sequence
import ..Variables: qvec, qval, bmat_gradient, VariableType, start_time, effective_time


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
    pulse_effects :: Vector{Union{Symbol, Number}}
    readout_index :: Integer

    # computed
    duration_states :: Dict{Any, SVector{4, VariableType}}
    qvec :: Dict{Any, SVector{3, VariableType}}
    bmat :: Dict{Any, SMatrix{3, 3, VariableType, 9}}
end

function Pathway(sequence::Sequence, pulse_effects::AbstractVector, readout_index::Integer=1)
    walker = PathwayWalker()
    walk_pathway!(sequence, walker, interpret_pulse_effects.(pulse_effects), Ref(readout_index))
    return Pathway(
        sequence,
        pulse_effects,
        readout_index,
        Dict(k => SVector{4, VariableType}(v) for (k, v) in pairs(walker.duration_states)),
        Dict(k => SVector{3, VariableType}(v) for (k, v) in pairs(walker.qvec)),
        Dict(k => SMatrix{3, 3, VariableType, 9}(v) for (k, v) in pairs(walker.bmat)),
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
    return pathway.duration_states[duration_state_index(transverse, positive)]
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


"""
    interpret_pulse_effects(number_or_symbol)

Interpret the various numbers and symbols that can be passed on to a `Pathway`.

The result will be one of:
- :ignore (if input is 0, :ignore, or :skip).
- :excite (if input is 90 or :excite).
- :refocus (if input is 180, :refocus, or :excite).
- :neg_excite (if input is -90, 270, or :negexcite).
"""
function interpret_pulse_effects(number::Number)
    normed = Int(number % 360)
    if normed == 0
        return :ignore
    elseif normed == 90
        return :excite
    elseif normed == 180
        return :refocus
    elseif normed == 270
        return :neg_excite
    else
        error("The pulse effect along a pathway should be divisible by 90, not $number.")
    end
end

function interpret_pulse_effects(sym::Symbol)
    mapping = Dict{Symbol, Symbol}(
        :ignore => :ignore,
        :skip => :ignore,
        :excite => :excite,
        :net_excite => :net_excite,
        :refocus => :refocus,
        :invert => :refocus,
    )
    if sym in keys(mapping)
        return mapping[sym]
    else
        all_symbols = join(keys(mapping), ", ")
        error("The pulse effect along a pathway should be one of ($all_symbols), not $sym.")
    end
end


"""
Helper structure for [`PathwayWalker`](@ref), which is itself a helper for `Pathway`.  You are deep down the rabit hole now...
    
For documentation, see that structure and [`walk_pathway`](@ref) and [`update_walker_gradient!`](@ref).
"""
mutable struct GradientTracker
    last_gradient_time :: VariableType
    qvec :: Vector{<:VariableType}
    bmat :: Matrix{<:VariableType}
end

GradientTracker() = GradientTracker(0., zeros(3), zeros(3, 3))


"""
Helper structure for `Pathway`.
    
For documentation, see that structure and [`walk_pathway`](@ref).
"""
mutable struct PathwayWalker
    last_pulse_time :: VariableType
    is_transverse :: Bool
    is_positive :: Bool
    duration_states :: Vector{<:VariableType}
    gradient_trackers :: Dict{Any, GradientTracker}
end

PathwayWalker() = PathwayWalker(
    0., false, true,
    zeros(4),
    Dict{Any, GradientTracker}()
)

"""
    walk_pathway!(bb::BuildingBlock, walker::PathwayWalker, pulse_effects::Vector, nreadout::Ref{Int}, start_time)

Computes the effect of a specific [`BuildingBlock`](@ref) (starting at `start_time`) on the [`PathwayWalker`](@ref).

For individual pulses and gradients, the following behaviour is implemented:
- If a pulse is encountered, call [`update_walker_pulse!`](@ref)`(walker, pulse_effects, pulse_effective_time)`
- If a gradient is encountered, call [`update_walker_gradient!`](@ref)(gradient, walker, gradient_start_time)

For overlapping gradients/pulses, one should first encounter the part of the gradient before the [`effective_time`](@ref) of the pulse,
then apply the pulse, and then the rest of the gradient.
This effective time can be passed on to [`update_walker_gradient!`](@ref) to allow part of the gradient waveform to be applied.

The function should return `true` if the `Pathway` has reached its end (i.e., the final readout) and `false` otherwise.
"""
function walk_pathway!(grad::GradientBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time=0.::VariableType) 
    update_walker_gradient!(grad, walker, block_start_time)
    return false
end

function walk_pathway!(pulse::RFPulseBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time=0.::VariableType) 
    update_walker_pulse!(walker, pulse_effects, block_start_time + effective_time(pulse))
    return iszero(length(pulse_effects)) && iszero(nreadout[])
end

function walk_pathway!(container::ContainerBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time=0.::VariableType)
    for (index, child) in get_children_blocks(container)
        if walk_pathway!(child, walker, pulse_effects, nreadout, block_start_time + start_time(container, index))
            return true
        end
    end
    return false
end


"""
    update_walker_pulse!(walker::PathwayWalker, pulse_effects::Vector, pulse_time)

Apply the first element of `pulse_effects` to the `walker` at the given `pulse_time`.

The following steps will be taken if the first `pulse_effect` is not `:ignore`
- if `walker.transverse` is true before the pulse, increase the `walker.bmat` by the outer product of `walker.qvec` with itself multiplied by the time since the last gradient
- update `walker.duration_states` with time since last pulse.
- update `walker.last_pulse_time`
- update `walker.is_transverse`, and `walker.is_positive` based on the first value in `pulse_effects`. 
- if `walker.is_positive` changed in the previous step than the `walker.qvec` needs to be flipped.
- remove the first element from `pulse_effects`.

If the first element is `:ignore` the only effect is that the first element is removed from `pulse_effects`.
"""
function update_walker_pulse!(walker::PathwayWalker, pulse_effects::AbstractVector{Symbol}, pulse_time::VariableType)
    if length(pulse_effects) == 0
        error("Pathway definition is invalid! Another RF pulse was encountered before the number of readouts expected from `nreadout` where detected.")
    end
    instruction = popfirst!(pulse_effects)
    if instruction == :ignore
        return
    end

    # update qvec/bmat
    if walker.is_transverse
        for gradient_tracker in values(walker.gradient_trackers)
            gradient_tracker.bmat += (
                (gradient_tracker.qvec .* gradient_tracker.qvec') .* 
                (pulse_time - gradient_tracker.last_gradient_time)
            )
            gradient_tracker.last_gradient_time = pulse_time
        end
    end
    prev_sign = walker.is_positive
    
    # update durations
    index = duration_state_index(walker.is_transverse, walker.is_positive)
    walker.duration_states[index] = walker.duration_states[index] + (pulse_time - walker.last_pulse_time)

    walker.last_pulse_time = pulse_time

    # -transverse, +longitudinal, +transverse, -longitudinal, -transverse, +longitudinal
    ordering = [(true, false), (false, true), (true, true), (false, false), (true, false)]

    if instruction == :refocus
        walker.is_positive = !walker.is_positive
    elseif instruction == :excite
        index = findfirst(isequal(walker.is_transverse, walker.is_positive), ordering)
        (walker.is_transverse, walker.is_positive) = ordering[index + 1]
    elseif instruction == :neg_excite
        index = findlast(isequal(walker.is_transverse, walker.is_positive), ordering)
        (walker.is_transverse, walker.is_positive) = ordering[index - 1]
    else
        error("Invalid pulse instruction ($instruction); This error should have been caught earlier.")
    end

    # flip qvec if needed
    if prev_sign != walker.is_positive
        for gradient_tracker in values(walker.gradient_trackers)
            gradient_tracker.qvec = -gradient_tracker.qvec
        end
    end
end

"""
    update_walker_gradient!(gradient_block::GradientBlock, walker::PathwayWalker, gradient_start_time::VariableType; overlapping_pulses=[], overlapping_readouts=[])

Update the walker's `qvec` and `bmat` based on the given `gradient_block`.

The following steps will be taken:
- Do nothing if `walker.transverse` is false
- increase the appropriate `walker.bmat` by the outer product of `walker.qvec` with itself multiplied by the time since the last gradient
- update the appropriate `walker.qvec` and `walker.bmat` based on the gradient waveform. This will require appropriate `qvec`/`bmat` functions to be defined for the gradient building block.
- update `walker.last_gradient_time` to the time at the end of the gradient.

This requires [`bmat`](@ref) and [`qvec`](@ref) to be implemented for the [`GradientBlock`](@ref).
"""
function update_walker_gradient!(gradient::GradientBlock, walker::PathwayWalker, gradient_start_time::VariableType)
    if walker.transverse
        return
    end

    # make sure the appropriate gradient tracker exists
    key = (gradient.scale, gradient.rotate)
    if !(key in keys(walker.gradient_trackers))
        walker.gradient_trackers[key] = GradientTracker()
    end
    tracker = walker.gradient_trackers[key]

    # update bmat till start of gradient
    tracker.bmat = tracker.bmat .+ (
        (tracker.qvec .* tracker.qvec') .* 
        (gradient_start_time - tracker.last_gradient_time)
    )

    tracker.last_gradient_time = gradient_start_time
    tracker.bmat = tracker.bmat .+ bmat(gradient, tracker.qvec, internal_start_time, internal_end_time)
    tracker.qvec = tracker.qvec .+ qvec(gradient, internal_start_time, internal_end_time)
end

"""
    duration_state_index(transverse, positive)

Returns the index of a specific state in the `Pathway.duration_state` vector.

This function is used internally to access that vector.
"""
function duration_state_index(transverse, positive)
    return Dict([
        (false, true) => 1,
        (true, true) => 2,
        (false, false) => 3,
        (true, false) => 4,
    ])[(Bool(transverse), Bool(positive))]
end

end