module Pathways
import LinearAlgebra: norm, tr
import StaticArrays: SVector, SMatrix
import ..Components: NoGradient, RFPulseComponent, ReadoutComponent, InstantGradient, GradientWaveform
import ..Containers: BaseSequence, Sequence, BaseBuildingBlock, waveform, events, waveform_sequence, start_time, AlternativeBlocks
import ..Variables: qvec, qval, qval3, bmat_gradient, VariableType, effective_time, duration, TR, bmat, bval, area_under_curve, duration_dephase, duration_transverse, VariableNotAvailable


"""
    Pathway(sequence::Sequence, pulse_effects::Vector{:Symbol/Number}, readout_index=1; group=nothing)

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
- `group`: which gradient grouping to consider for the `qvec` and `bmat`. If not set, all gradients will be considered (using their current alignment).

## Attributes
Over the pathway the following values are computed. Each can be accessed by calling the appropriate function:

### Timings
- [`duration_state`](@ref): The total amount of time spent in a specific state in this pathway (+longitudinal, +transverse, -longitudinal, or -transverse)
- [`duration_transverse`](@ref): The total amount of time the spins spent in the transverse plane in ms. This can be used to quantify the expected effect of T2-decay.
- [`duration_dephase`](@ref): The total amount of time the spins spent in the +transverse relative to -transverse state in ms. The absolute value of this can be used to quantify the expected effect of T2'-decay.

### Effect of gradients
The area under curve, q-values, and b-values are computed separately for each group of gradients (depending on the `group` keyword set during construction).
- [`qvec`](@ref): Net displacement vector in k-space/q-space.
- [`qval`](@ref)/[`area_under_curve`](@ref): size of the displacement in k-space/q-space. For a spoiled pathway, this should be large compared with 1/voxel size; for unspoiled pathways it should be (close to) zero.
- [`bmat`](@ref): Net diffusion weighting due to gradients along the `Pathway` in matrix form.
- [`bval`](@ref): Net diffusion weighting due to gradients along the `Pathway` as a single number.
"""
struct Pathway
    # user provided
    sequence :: Sequence
    pulse_effects :: Vector{<:Union{Symbol, Number}}
    readout_index :: Integer

    # computed
    duration_states :: SVector{4, <:VariableType}
    qvec :: SVector{3, <:VariableType}
    bmat :: SMatrix{3, 3, <:VariableType, 9}
end

function Pathway(sequence::Sequence, pulse_effects::AbstractVector, readout_index::Integer=1; group=nothing)
    walker = PathwayWalker()
    walk_pathway!(sequence, walker, interpret_pulse_effects.(pulse_effects), Ref(readout_index))
    tracker = walker.gradient_trackers[group]
    return Pathway(
        sequence,
        pulse_effects,
        readout_index,
        SVector{4}(walker.duration_states),
        tracker.qvec,
        tracker.bmat,
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
    qvec(pathway::Pathway)

Return net displacement vector in k-space/q-space experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

Returns a NamedTuple with the `qvec` for all gradient groups.
"""
qvec(pathway::Pathway) = pathway.qvec


"""
    area_under_curve(pathway::Pathway)

Return net displacement in k-space (i.e., spoiling) experienced by the spins following a specific [`Pathway`](@ref).

Only gradients active while the spins are in the transverse plane are considered.

Returns a NamedTuple with the `area_under_curve` for all gradient groups.
"""
area_under_curve(pathway::Pathway) = norm(qvec(pathway))


"""
    bmat(pathway::Pathway)

Return 3x3 diffusion-weighted matrix experienced by the spins following a specific [`Pathway`](@ref) in rad^2 ms/um^2.

Only gradients active while the spins are in the transverse plane are considered.

Returns a NamedTuple with the `bmat` for all gradient groups.
"""
bmat(pathway::Pathway)  = pathway.bmat

"""
    bval(pathway::Pathway)

Return size of diffusion-weighting experienced by the spins following a specific [`Pathway`](@ref) in rad^2 ms/um^2.

Only gradients active while the spins are in the transverse plane will contribute to the diffusion weighting.

Returns a NamedTuple with the `bval` for all gradient groups.
"""
bval(pathway::Pathway) = tr(bmat(pathway))


"""
    get_pathway(sequence)

Gets the main [`PathWay`](@ref) that spins are expected to experience in the sequence.

Multiple pathways might be returned as an array or (named)tuple.
"""
function get_pathway end

for fn in (:qvec, :area_under_curve, :bmat, :bval, :duration_dephase, :duration_transverse)
    @eval function $fn(seq::Sequence)
        pathway = try
            get_pathway(seq)
        catch e
            if e isa MethodError
                throw(VariableNotAvailable(typeof(seq), $fn))
            end
            rethrow()
        end
        if pathway isa Pathway
            return $fn(pathway)
        elseif pathway isa AbstractVector || pathway isa Tuple
            return $fn.(pathway)
        elseif pathway isa NumedTuple
            return NamedTuple(k => $fn(v) for (k, v) in pairs(pathway))
        end
        error("get_pathway returned unexpected type for $seq")
    end
end


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
    
For documentation, see that structure and [`walk_pathway!`](@ref) and [`update_walker_gradient!`](@ref).
"""
mutable struct GradientTracker
    last_gradient_time :: VariableType
    qvec :: Vector{VariableType}
    bmat :: Matrix{VariableType}
end

GradientTracker() = GradientTracker(0., zeros(3), zeros(3, 3))


"""
Helper structure for `Pathway`.
    
For documentation, see that structure and [`walk_pathway!`](@ref).
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
    Dict{Any, GradientTracker}(nothing => GradientTracker())
)

"""
    walk_pathway!(bb::Sequence/BuildingBlock, walker::PathwayWalker, pulse_effects::Vector, nreadout::Ref{Int}, start_time)

Computes the effect of a specific [`ContainerBlock`](@ref) (starting at `start_time`) on the [`PathwayWalker`](@ref).

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
    nwait = length(pulse_effects) + nreadout[]
    while !(walk_pathway!(seq, walker, pulse_effects, nreadout, current_TR * TR(seq)))
        new_nwait = length(pulse_effects) + nreadout[]
        if nwait == new_nwait
            not_seen = iszero(length(pulse_effects)) ? "readout" : "pulse"
            error("Pathway iterated through the whole sequence without seeing a valid $not_seen. Terminating...")
        end
        nwait = new_nwait

        current_TR += 1
    end
    return true
end

function walk_pathway!(seq::BaseSequence, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType)
    for (index, child) in enumerate(seq)
        if walk_pathway!(child, walker, pulse_effects, nreadout, block_start_time + start_time(seq, index))
            return true
        end
    end
    return false
end

function walk_pathway!(block::BaseBuildingBlock, walker::PathwayWalker, pulse_effects::Vector{Symbol}, nreadout::Ref{Int}, block_start_time::VariableType)
    current_index = nothing
    current_time = block_start_time
    for (index_inter, interruption) in events(block)

        # determine if action should be taken
        if interruption isa RFPulseComponent
            if iszero(length(pulse_effects))
                error("Pathway definition is invalid! Another RF pulse was encountered before the number of readouts expected from `nreadout` where detected.")
            end
            if pulse_effects[1] == :ignore
                popfirst!(pulse_effects)
                continue
            end
        elseif interruption isa ReadoutComponent
            if length(pulse_effects) > 0
                continue
            end
            nreadout[] -= 1
            if nreadout[] > 0
                continue
            end
        end 

        # apply gradients up till interrupt
        for (_, part) in waveform_sequence(block, current_index, index_inter)
            update_walker_gradient!(part, walker, current_time)
            current_time = current_time + duration(part)
        end 

        # apply interrupt
        if interruption isa RFPulseComponent
            update_walker_pulse!(walker, pulse_effects, current_time)
        elseif interruption isa InstantGradient
            update_walker_instant_gradient!(interruption, walker, current_time)
        end
        current_index = index_inter
        if length(pulse_effects) == 0 && nreadout[] == 0
            update_walker_till_time!(walker, current_time)
            return true
        end
    end

    # apply remaining gradients
    for (_, part) in waveform_sequence(block, current_index, nothing)
        update_walker_gradient!(part, walker, current_time)
        current_time = current_time + duration(part)
    end
    return false
end


"""
    update_walker_till_time!(walker::PathwayWalker, new_time[, gradient_group])

Updates all parts of a [`PathwayWalker`](@ref) up to the given time.

This updates the `walker.duration_states` and the `bmat` for each gradient tracker.
If `gradient_group` are provided, then only the gradient tracker matching that group will be updated.
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
function update_gradient_tracker_till_time!(walker::PathwayWalker, key::Union{Nothing, Symbol}, new_time::VariableType)
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
    update_walker_gradient!(gradient_block::GradientWaveform, walker::PathwayWalker, gradient_start_time::VariableType; overlapping_pulses=[], overlapping_readouts=[])

Update the walker's `qvec` and `bmat` based on the given `gradient_block`.

The following steps will be taken:
- Do nothing if `walker.is_transverse` is false
- increase the appropriate `walker.bmat` by the outer product of `walker.qvec` with itself multiplied by the time since the last gradient
- update the appropriate `walker.qvec` and `walker.bmat` based on the gradient waveform. This will require appropriate `qvec`/`bmat` functions to be defined for the gradient building block.
- update `walker.last_gradient_time` to the time at the end of the gradient.

This requires [`bmat_gradient`](@ref) and [`qvec`](@ref) to be implemented for the [`GradientWaveform`](@ref).
"""
update_walker_gradient!(gradient::NoGradient, walker::PathwayWalker, gradient_start_time::VariableType) = nothing

function update_walker_gradient!(gradient::GradientWaveform, walker::PathwayWalker, gradient_start_time::VariableType)
    if !walker.is_transverse
        return
    end

    # update gradient tracker till start of gradient
    for key in (isnothing(gradient.group) ? [nothing] : [nothing, gradient.group])
        update_gradient_tracker_till_time!(walker, key, gradient_start_time)

        # update qvec/bmat during gradient
        tracker = walker.gradient_trackers[key]
        tracker.bmat = tracker.bmat .+ bmat_gradient(gradient, tracker.qvec)
        tracker.qvec = tracker.qvec .+ qval3(gradient)
        tracker.last_gradient_time = gradient_start_time + duration(gradient)
    end
end

"""
    update_walker_instant_gradient!(gradient, walker)
"""
function update_walker_instant_gradient!(gradient::InstantGradient{N}, walker::PathwayWalker, gradient_start_time::VariableType) where {N}
    if N == 1
        qvec3 = qval(gradient) .* gradient.orientation
    else
        qvec3 = qval(gradient)
    end
    for key in (isnothing(gradient.group) ? [nothing] : [nothing, gradient.group])
        update_gradient_tracker_till_time!(walker, key, gradient_start_time)
        tracker = walker.gradient_trackers[key]
        tracker.qvec = tracker.qvec  .+ qvec3
    end
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