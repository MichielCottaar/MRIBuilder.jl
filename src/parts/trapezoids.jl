"""
Defines a set of different options for MRI gradients.
"""
module Trapezoids

import JuMP: @constraint
import StaticArrays: SVector
import LinearAlgebra: norm
import ...Variables: variables, @defvar, scanner_constraints!, get_free_variable, set_simple_constraints!, gradient_orientation, VariableType, get_gradient, get_pulse, get_readout, adjust, adjust_groups, adjust_internal, apply_simple_constraint!, add_cost_function!
import ...Components: ChangingGradient, ConstantGradient, RFPulseComponent, ADC
import ...Containers: BaseBuildingBlock


"""
Parent type for any `BuildingBlock` that has a trapezoidal gradient waveform.

Sub-types:
- [`Trapezoid`](@ref)
- [`SliceSelect`](@ref)
- [`LineReadout`](@ref)

The `N` indicates whether the gradient has a fixed orientation (N=1) or is free (N=3).
"""
abstract type BaseTrapezoid{N} <: BaseBuildingBlock end

"""
    Trapezoid(; orientation=nothing, group=nothing, variables...)

Defines a trapezoidal pulsed gradient

## Parameters
- `orientation` sets the gradient orientation (completely free by default). Can be set to a vector for a fixed orientation.
- `group`: assign the trapezoidal gradient to a specific group. This group will be used to scale or rotate the gradients after optimisation.

## Variables
Variables can be set during construction or afterwards as an attribute.
If not set, they will be determined during the sequence optimisation.
### Timing variables
- [`variables.rise_time`](@ref): Time of the gradient to reach from 0 to maximum in ms. If explicitly set to 0, the scanner slew rate will be ignored.
- [`variables.flat_time`](@ref): Time that the gradient stays at maximum strength in ms.
- [`variables.δ`](@ref): effective pulse duration (`rise_time` + `flat_time`) in ms.
- [`variables.duration`](@ref): total pulse duration (2 * `rise_time` + `flat_time`) in ms.
### Gradient variables
- [`variables.gradient_strength`](@ref): Maximum gradient strength achieved during the pulse in kHz/um
- [`variables.qvec`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um (given by `δ` * `gradient_strength`).

The `bvalue` can be constrained for multiple gradient pulses by creating a `Pathway`.
"""
abstract type Trapezoid{N} <: BaseTrapezoid{N} end

struct Trapezoid1D <: Trapezoid{1}
    rise_time :: VariableType
    flat_time :: VariableType
    slew_rate :: VariableType
    orientation :: SVector{3, Float64}
    group :: Union{Nothing, Symbol}
end

struct Trapezoid3D <: Trapezoid{3}
    rise_time :: VariableType
    flat_time :: VariableType
    slew_rate :: SVector{3, VariableType}
    group :: Union{Nothing, Symbol}
end

function (::Type{Trapezoid})(; orientation=nothing, rise_time=nothing, flat_time=nothing, group=nothing, slew_rate=nothing, kwargs...)
    if isnothing(orientation)
        if isnothing(slew_rate)
            slew_rate = (nothing, nothing, nothing)
        end
        res = Trapezoid3D(
            get_free_variable(rise_time),
            get_free_variable(flat_time),
            get_free_variable.(slew_rate),
            group
        )
        add_cost_function!(-sum(abs.(res.slew_rate)))
    else
        res = Trapezoid1D(
            get_free_variable(rise_time),
            get_free_variable(flat_time),
            get_free_variable(slew_rate),
            orientation,
            group
        )
        apply_simple_constraint!(res.slew_rate, :>=, 0)
        add_cost_function!(-res.slew_rate)
    end

    set_simple_constraints!(res, kwargs)

    apply_simple_constraint!(res.flat_time, :>=, 0)
    apply_simple_constraint!(res.rise_time, :>=, 0)
    scanner_constraints!(res)
    return res
end

Base.keys(::Trapezoid) = (Val(:rise), Val(:flat), Val(:fall))

Base.getindex(pg::BaseTrapezoid{N}, ::Val{:rise}) where {N} = ChangingGradient(zeros(3), variables.slew_rate(pg), variables.rise_time(pg), get_group(pg))
Base.getindex(pg::BaseTrapezoid, ::Val{:flat}) = ConstantGradient(variables.gradient_strength(pg), variables.flat_time(pg), get_group(pg))
Base.getindex(pg::BaseTrapezoid, ::Val{:fall}) = ChangingGradient(variables.gradient_strength(pg), -variables.slew_rate(pg), variables.rise_time(pg), get_group(pg))
gradient_orientation(::BaseTrapezoid{3}) = nothing
gradient_orientation(pg::BaseTrapezoid{1}) = gradient_orientation(get_gradient(pg))
gradient_orientation(pg::Trapezoid{1}) = pg.orientation

get_group(pg::Trapezoid) = pg.group
get_group(pg::BaseTrapezoid) = get_group(get_gradient(pg))

@defvar gradient begin
    rise_time(pg::Trapezoid) = pg.rise_time
    flat_time(pg::Trapezoid) = pg.flat_time
    slew_rate(g::Trapezoid1D) = g.slew_rate .* g.orientation
    slew_rate(g::Trapezoid3D) = g.slew_rate
end

"""
    rise_time(trapezoid)

Returns the rise time of a [`Trapezoid`](@ref) gradient profile in ms.
"""
variables.rise_time

"""
    flat_time(trapezoid)

Returns the flat time of a [`Trapezoid`](@ref) gradient profile in ms.
"""
variables.flat_time

@defvar gradient begin
    gradient_strength(g::Trapezoid) = variables.slew_rate(g) .* variables.rise_time(g)
    δ(g::Trapezoid) = variables.rise_time(g) + variables.flat_time(g)
end

"""
    δ(trapezoid)

Returns the effective duration of a [`Trapezoid`](@ref) gradient profile in ms.

Defined as [`variables.rise_time`](@ref) + [`variables.flat_time`](@ref).
"""
variables.δ

@defvar duration(g::BaseTrapezoid) = 2 * variables.rise_time(g) + variables.flat_time(g)

@defvar qvec(g::BaseTrapezoid, ::Nothing, ::Nothing) = variables.δ(g) .* variables.gradient_strength(g) .* 2π

adjust_groups(t::Trapezoid) = [t.group, :gradient]

function adjust_internal(trap::Trapezoid1D; orientation=nothing, scale=1., rotation=nothing)
    if !isnothing(orientation) && !isnothing(rotation)
        error("Cannot set both the gradient orientation and rotation.")
    end
    new_orientation = isnothing(orientation) ? (isnothing(rotation) ? trap.orientation : rotation * trap.orientation) : orientation
    return Trapezoid1D(
        trap.rise_time,
        trap.flat_time,
        trap.slew_rate * scale,
        new_orientation,
        trap.group
    )
end

function adjust_internal(trap::Trapezoid3D; scale=1., rotation=nothing)
    return Trapezoid3D(
        trap.rise_time,
        trap.flat_time,
        (
            isnothing(rotation) ? 
            (trap.slew_rate .* scale) : 
            (rotation * (trap.slew_rate .* scale))
        ),
        trap.group
    )
end


"""
    SliceSelect(pulse; orientation=nothing, group=nothing, variables...)

Defines a trapezoidal gradient with a pulse played out during the flat time.

Parameters and variables are identical as for [`Trapezoid`](@ref) with the addition of:

## Parameters
- `pulse`: sub-type of [`RFPulseComponent`](@ref) that describes the RF pulse.

## Variables
- `slice_thickness`: thickness of the selected slice in mm
"""
struct SliceSelect{N} <: BaseTrapezoid{N}
    trapezoid :: Trapezoid{N}
    pulse :: RFPulseComponent
    group :: Symbol
end

function SliceSelect(pulse::RFPulseComponent; orientation=nothing, rise_time=nothing, group=:FOV, slew_rate=nothing, vars...)
    res = SliceSelect(
        Trapezoid(; orientation=orientation, rise_time=rise_time, flat_time=variables.duration(pulse), slew_rate=slew_rate, group=nothing),
        pulse,
        group
    )
    set_simple_constraints!(res, vars)
    return res
end

Base.keys(::SliceSelect) = (Val(:rise), Val(:flat), Val(:pulse), Val(:fall))
Base.getindex(pg::SliceSelect, ::Val{:pulse}) = (0., pg.pulse)

@defvar pulse inverse_slice_thickness(ss::SliceSelect) = 1e3 * variables.gradient_strength_norm(ss.trapezoid) * variables.inverse_bandwidth(ss.pulse)

"""
    slice_thickness(slice_select)

Defines the slice thickness for a RF pulse with an active gradient in mm (e.g., [`SliceSelect`](@ref)).

Defines as [`variables.gradient_strength_norm`](@ref)(gradient) / [`variables.bandwidth`](@ref)(pulse)
"""
variables.slice_thickness

get_pulse(ss::SliceSelect) = ss.pulse
get_gradient(ss::SliceSelect) = ss.trapezoid
@defvar effective_time(ss::SliceSelect) = variables.effective_time(ss, :pulse)

adjust_groups(ss::SliceSelect) = [ss.group, :slice]
function adjust_internal(ss::SliceSelect; shift)
    return SliceSelect(
        ss.trapezoid,
        adjust(
            ss.pulse;
            pulse=(frequency=variables.gradient_strength_norm(ss) * shift * 1e3,)
        ),
        ss.group,
    )
end


"""
    LineReadout(adc; ramp_overlap=1., orientation=nothing, group=nothing, variables...)

Defines a trapezoidal gradient with an ADC readout overlaid.

Parameters and variables are identical as for [`Trapezoid`](@ref) with the addition of:

## Parameters
- `adc`: [`ADC`](@ref) object that describes the readout.
- `ramp_overlap`: how much the gradient ramp should overlap with the ADC. 0 for no overlap, 1 for full overlap (default: 1). Can be set to `nothing` to become a free variable.

## Variables
- [`variables.fov`](@ref): FOV of the output image along this single k-space line in mm.
- [`variables.voxel_size`](@ref): size of each voxel along this single k-space line in mm.
"""
struct LineReadout{N} <: BaseTrapezoid{N}
    trapezoid :: Trapezoid{N}
    adc :: ADC
    ramp_overlap :: VariableType
end

function LineReadout(adc::ADC; ramp_overlap=nothing, orientation=nothing, group=nothing, vars...)
    res = LineReadout(
        Trapezoid(; orientation=orientation, group=group),
        adc,
        get_free_variable(ramp_overlap)
    )
    set_simple_constraints!(res, vars)
    if !(res.ramp_overlap isa Number)
        apply_simple_constraint!(res.ramp_overlap, :>=, 0.)
        apply_simple_constraint!(res.ramp_overlap, :<=, 1.)
    end
    apply_simple_constraint!(
        res.ramp_overlap * variables.rise_time(res.trapezoid) * 2 + variables.flat_time(res.trapezoid),
        variables.duration(res.adc)
    )
    return res
end

Base.keys(::LineReadout) = (Val(:rise), Val(:adc), Val(:flat), Val(:fall))
Base.getindex(lr::LineReadout, ::Val{:adc}) = ((1 - variables.ramp_overlap(lr)) * variables.rise_time(lr), lr.adc)

@defvar begin
    ramp_overlap(lr::LineReadout) = lr.ramp_overlap
    inverse_fov(lr::LineReadout) = 1e3 * variables.dwell_time(lr.adc) * variables.gradient_strength_norm(lr.trapezoid) * lr.adc.oversample
    inverse_voxel_size(lr::LineReadout) = 1e3 * variables.duration(lr.adc) * variables.gradient_strength(lr.trapezoid)
    effective_time(lr::LineReadout) = variables.effective_time(lr, :adc)
end

"""
    fov(readout)

Defines the field of view of a readout in mm.
"""
variables.fov

"""
    voxel_size(readout)

Defines the voxel size of a readout in mm.
"""
variables.voxel_size

"""
    ramp_overlap(line_readout)

Return the fraction of the gradient ramp that overlaps with the ADC readout.

Set to 0 to ensure that the ADC is only active during the flat time of the readout.
"""
variables.ramp_overlap

get_readout(lr::LineReadout) = lr.adc
get_gradient(lr::LineReadout) = lr.trapezoid

"""
    opposite_kspace_lines(; orientation=[1, 0, 0], kwargs...)

Return a positive and negative readout of a k-space line.
"""
function opposite_kspace_lines(; orientation=[1, 0, 0], kwargs...)
    if isnothing(orientation)
        error("orientation of k-space readout should be fixed at construction.")
    end
    positive = LineReadout(ADC(); orientation=orientation, kwargs...)
    trap = positive.trapezoid
    negative = LineReadout(
        Trapezoid1D(trap.rise_time, trap.flat_time, trap.slew_rate, -orientation, trap.group),
        positive.adc,
        positive.ramp_overlap,
    )
    return (positive, negative)
end

end