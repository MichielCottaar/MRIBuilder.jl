module Pathways
import LinearAlgebra: norm, tr
import StaticArrays: SVector, SMatrix
import ..Components: NoGradient, ConstantGradient, ChangingGradient, RFPulseComponent, ReadoutComponent, InstantGradient
import ..AllSequences: BaseSequence
import ..AllBuildingBlocks: BaseBuildingBlock, waveform, events, get_parts
import ..Variables: qvec, qval, bmat_gradient, VariableType, start_time, effective_time, duration, qval_square
import ..Alternatives: AlternativeBlocks


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
    pulse_effects :: Vector{<:Union{Symbol, Number}}
    readout_index :: Integer

    # computed
    duration_states :: SVector{4, <:VariableType}
    qvec :: Dict{Any, SVector{3, <:VariableType}}
    bmat :: Dict{Any, SMatrix{3, 3, <:VariableType, 9}}
end

function Pathway(sequence::Sequence, pulse_effects::AbstractVector, readout_index::Integer=1)
    walker = PathwayWalker()
    walk_pathway!(sequence, walker, interpret_pulse_effects.(pulse_effects), Ref(readout_index))
    return Pathway(
        sequence,
        pulse_effects,
        readout_index,
        SVector{4}(walker.duration_states),
        Dict(k => SVector{3}(v.qvec) for (k, v) in pairs(walker.gradient_trackers)),
        Dict(k => SMatrix{3, 3}(v.bmat) for (k, v) in pairs(walker.gradient_trackers)),
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

function qval_square(pathway::Pathway; kwargs...)
    vec = qvec(pathway; kwargs...)
    return vec[1]^2 + vec[2]^2 + vec[3]^2
end
qval(pathway::Pathway; kwargs...) = sqrt(qval_square(pathway; kwargs...))


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

Return 3x3 diffusion-weighted matrix experienced by the spins following a specific [`Pathway`](@ref) in rad^2 ms/um^2.

Only gradients active while the spins are in the transverse plane are considered.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by speficic `scale`/`rotate` parameters
"""
bmat(pathway::Pathway; scale=nothing, rotate=nothing)  = get(pathway.bmat, (scale, rotate), zero(SMatrix{3, 3, Float64, 9}))

"""
    bval(pathway::Pathway; scale=nothing, rotate=nothing)

Return size of diffusion-weighting experienced by the spins following a specific [`Pathway`](@ref) in rad^2 ms/um^2.

Only gradients active while the spins are in the transverse plane will contribute to the diffusion weighting.

By default gradients that are affected by user-provided `scale` or `rotate` parameters (e.g., bvals/bvecs) are ignored.
You can set `scale` and/or `rotate` to specific symbols to only consider gradients that are affected by specific `scale`/`rotate` parameters
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
    qvec :: Vector{VariableType}
    bmat :: Matrix{VariableType}
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
    duration_states :: Vector{VariableType}
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
function walk_pathway!(seq::Sequence, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}) 
    current_TR = 0
    while !(walk_pathway!(seq, walker, pulse_effects, nreadout, current_TR * seq.TR))
        current_TR += 1
    end
    return true
end

function walk_pathway!(grad::GradientBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType) 
    update_walker_gradient!(grad, walker, block_start_time)
    return false
end

function walk_pathway!(alt::AlternativeBlocks, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType) 
    walk_pathway!(alt[1], walker, pulse_effects, nreadout, block_start_time)
end

function walk_pathway!(pulse::RFPulseBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType) 
    update_walker_pulse!(walker, pulse_effects, block_start_time + effective_time(pulse))
    return iszero(length(pulse_effects)) && iszero(nreadout[])
end

function walk_pathway!(container::ContainerBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType)
    for index in get_children_indices(container)
        child = container[index]
        if walk_pathway!(child, walker, pulse_effects, nreadout, block_start_time + start_time(container, index))
            return true
        end
    end
    return false
end

walk_pathway!(wait::WaitBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType) = false

function walk_pathway!(ao::AbstractOverlapping, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType)
    current_index = nothing
    current_time = block_start_time
    for (index_inter, interruption) in enumerate(interruptions(ao))
        if interruption.object isa RFPulseBlock
            if iszero(length(pulse_effects))
                error("Pathway definition is invalid! Another RF pulse was encountered before the number of readouts expected from `nreadout` where detected.")
            end
            if pulse_effects[1] == :ignore
                popfirst!(pulse_effects)
                continue
            end
        elseif (interruption.object isa InstantReadout || interruption.object isa ADC)
            if length(pulse_effects) > 0
                continue
            end
            nreadout[] -= 1
            if nreadout[] > 0
                continue
            end
        end 
        for part in get_parts(ao, current_index, index_inter)
            update_walker_gradient!(part, walker, current_time)
            current_time = current_time + duration(part)
        end 
        if interruption.object isa RFPulseBlock
            update_walker_pulse!(walker, pulse_effects, current_time)
        end
        current_index = index_inter
        if length(pulse_effects) == 0 && nreadout[] == 0
            update_walker_till_time!(walker, current_time)
            return true
        end
    end
    for part in get_parts(ao, current_index, nothing)
        update_walker_gradient!(part, walker, current_time)
        current_time = current_time + duration(part)
    end
    return false
end

function walk_pathway!(::Union{InstantReadout, ADC}, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time=0.::VariableType)
    if length(pulse_effects) > 0
        return false
    end
    nreadout[] -= 1
    if iszero(nreadout[])
        update_walker_till_time!(walker, block_start_time)
        return true
    elseif nreadout[] < 0
        error("Pathway walker continued past the point where it should have ended. Did you start with a negative `nreadout`?")
    end
    return false
end


"""
    update_walker_till_time!(walker::PathwayWalker, new_time[, (rotate, scale)])

Updates all parts of a [`PathwayWalker`](@ref) up to the given time.

This updates the `walker.duration_states` and the `bmat` for each gradient tracker.
If `rotate` and `scale` are provided, then only the gradient tracker matching those properties will be updated.
If that gradient tracker does not exist, it will be created.

This function is used to get the `walker` up to date till the start of a gradient, pulse, or final readout.
"""
function update_walker_till_time!(walker::PathwayWalker, new_time::VariableType, gradient_key=nothing)
    # update duration state and pulse time
    index = duration_state_index(walker.is_transverse, walker.is_positive)
    walker.duration_states[index] = walker.duration_states[index] + (new_time - walker.last_pulse_time)
    walker.last_pulse_time = new_time

    if isnothing(gradient_key)
        for tracker in values(walker.gradient_trackers)
            update_gradient_tracker_till_time!(tracker, new_time)
        end
    else
        update_gradient_tracker_till_time!(walker.gradient_trackers, gradient_key, new_time)
    end
end

"""
    update_gradient_tracker_till_time!(walker::PathwayWalker, key, new_time)
    update_gradient_tracker_till_time!(tracker::GradientTracker, new_time)

Update the `bmat` for any time passed since the last update (assuming there will no gradients during that period).

The `bmat` is updated with the outer produce of `qvec` with itself multiplied by the time since the last update.

When called with the first signature the tracker will be created from scratch if a tracker with that `key` does not exist.
"""
function update_gradient_tracker_till_time!(walker::PathwayWalker, key::Tuple, new_time::VariableType)
    if !(key in keys(walker.gradient_trackers))
        walker.gradient_trackers[key] = GradientTracker()
    end
    update_gradient_tracker_till_time!(walker.gradient_trackers[key], new_time)
end

function update_gradient_tracker_till_time!(gradient_tracker::GradientTracker, new_time::VariableType)
    gradient_tracker.bmat += (
        (gradient_tracker.qvec .* permutedims(gradient_tracker.qvec)) .* 
        (new_time - gradient_tracker.last_gradient_time)
    )
    gradient_tracker.last_gradient_time = new_time
end

"""
    update_walker_pulse!(walker::PathwayWalker, pulse_effects::Vector, pulse_time)

Apply the first element of `pulse_effects` to the `walker` at the given `pulse_time`.

The following steps will be taken if the first `pulse_effect` is not `:ignore`
- if `walker.is_transverse` is true before the pulse, increase the `walker.bmat` by the outer product of `walker.qvec` with itself multiplied by the time since the last gradient
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

    update_walker_till_time!(walker, pulse_time)

    prev_sign = walker.is_positive
    # -transverse, +longitudinal, +transverse, -longitudinal, -transverse, +longitudinal
    ordering = [(true, false), (false, true), (true, true), (false, false), (true, false)]

    if instruction == :refocus
        walker.is_positive = !walker.is_positive
    elseif instruction == :excite
        index = findfirst(isequal((walker.is_transverse, walker.is_positive)), ordering)
        (walker.is_transverse, walker.is_positive) = ordering[index + 1]
    elseif instruction == :neg_excite
        index = findlast(isequal((walker.is_transverse, walker.is_positive)), ordering)
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
- Do nothing if `walker.is_transverse` is false
- increase the appropriate `walker.bmat` by the outer product of `walker.qvec` with itself multiplied by the time since the last gradient
- update the appropriate `walker.qvec` and `walker.bmat` based on the gradient waveform. This will require appropriate `qvec`/`bmat` functions to be defined for the gradient building block.
- update `walker.last_gradient_time` to the time at the end of the gradient.

This requires [`bmat_gradient`](@ref) and [`qvec`](@ref) to be implemented for the [`GradientBlock`](@ref).
"""
function update_walker_gradient!(gradient::GradientBlock, walker::PathwayWalker, gradient_start_time::VariableType)
    if !walker.is_transverse
        return
    end

    # update gradient tracker till start of gradient
    key = (gradient.scale, gradient.rotate)
    update_gradient_tracker_till_time!(walker, key, gradient_start_time)

    # update qvec/bmat during gradient
    tracker = walker.gradient_trackers[key]
    tracker.bmat = tracker.bmat .+ bmat_gradient(gradient, tracker.qvec)
    tracker.qvec = tracker.qvec .+ qvec(gradient)
    tracker.last_gradient_time = gradient_start_time + duration(gradient)
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